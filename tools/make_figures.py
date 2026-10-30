"""Gera as figuras da aula (imagens/*.png).

Estruturas reais do PDB desenhadas em esferas (projecao ortografica, mesma
linguagem visual do curso de dinamica molecular); esquemas marcados como
ilustrativos. As poses e energias do docking vem de runs/swissdock, gerado por
tools/swissdock_cli.py, ou de resultados/ quando runs/ nao existe.

Uso: python tools/make_figures.py
"""
from __future__ import annotations

import json
import sys
import urllib.request
from pathlib import Path

import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Polygon, Rectangle
from rdkit import Chem
from rdkit.Chem import Draw

sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import (C, CACHE, CRISTAL, DESTAQUE, ERRO, GREY, INK, BAND, OUT, POSE, PROT, ROOT, SITE, WATER,  # noqa: E402
                   draw_spheres, fetch_pdb, frame_off, project, save, shade, view_basis)
from analyze_poses import poses, skeleton  # noqa: E402

A = 0.1          # angstrom -> nm (mdtraj trabalha em nm)
R_ATOM = 0.17    # raio de desenho dos atomos pesados (nm)


def resultados(nome):
    for base in (ROOT / 'runs' / 'swissdock' / nome, ROOT / 'resultados' / nome):
        if (base / 'vina_dock.pdbqt').exists():
            return base
    raise FileNotFoundError(nome)


def ligand_sdf(code, comp, chain='A'):
    path = CACHE / f'{code}_{comp}.sdf'
    if not path.exists():
        CACHE.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(
            f'https://models.rcsb.org/v1/{code}/ligand?label_comp_id={comp}&auth_asym_id={chain}&encoding=sdf', path)
    m = Chem.MolFromMolFile(str(path), sanitize=False, removeHs=False)
    return m.GetConformer().GetPositions() * A, [a.GetSymbol() for a in m.GetAtoms()]


def heavy(t, sel):
    return t.topology.select(f'({sel}) and not element H')


def fit(fig, ax, largura=12, extra=0.8):
    """Altura da figura pela proporcao dos limites do eixo (evita faixa em branco com aspect igual)."""
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    fig.set_size_inches(largura, largura * abs(y1 - y0) / abs(x1 - x0) + extra)


def legend(fig, items, ncol=None, y=-0.02, ax=None):
    """Legenda embaixo; com `ax`, ancorada ao eixo (evita vao quando o eixo tem aspecto fixo)."""
    hs = [Line2D([], [], ls='none', marker='o', ms=9, mfc=c, mec=INK, mew=0.5) for _, c in items]
    alvo = ax if ax is not None else fig
    alvo.legend(hs, [l for l, _ in items], loc='lower center', ncol=ncol or len(items), frameon=False,
                fontsize=8.5, bbox_to_anchor=(0.5, y))


def label(ax, xy, text, xytext, color):
    """Rotulo com linha de indicacao; a origem da linha fica marcada por um ponto."""
    ax.annotate(text, xy=xy, xytext=xytext, fontsize=8.5, ha='center', va='center',
                arrowprops=dict(arrowstyle='-', color=color, lw=1.1, shrinkA=2, shrinkB=0),
                bbox=dict(boxstyle='round,pad=0.2', fc='white', ec=color, lw=0.8), zorder=10)
    ax.plot(*xy, 'o', ms=6, mfc=color, mec='white', mew=1.2, zorder=11)


def front(P, idx):
    """Posicao na tela do atomo mais a frente entre os indices `idx` (o que fica visivel)."""
    idx = list(idx)
    return P[idx[int(np.argmax(P[idx, 2]))], :2]


def pocket_scene(t, lig_xyz, extra=(), cut=1.0, chain_sel='protein'):
    """Atomos pesados do receptor a menos de `cut` nm do ligante, com o ligante e camadas extras."""
    prot = heavy(t, chain_sel)
    xyz = t.xyz[0]
    d = np.min(np.linalg.norm(xyz[prot][:, None] - lig_xyz[None], axis=2), axis=1)
    near = prot[d < cut]
    return near, xyz


# ---------------------------------------------------------------- 1. modelos de ligacao
def _beads_in(poly_test, x0, x1, y0, y1, step, rng, jitter=0.25):
    xs, ys = np.meshgrid(np.arange(x0, x1 + 1e-9, step), np.arange(y0, y1 + 1e-9, step * 0.87))
    xs = xs + (np.arange(xs.shape[0])[:, None] % 2) * step / 2
    pts = np.c_[xs.ravel(), ys.ravel()]
    pts = pts + rng.normal(0, step * jitter * 0.3, pts.shape)
    return pts[poly_test(pts)]


def _receptor(ax, x, y, pocket, color=SITE, w=2.6, h=1.6, r=0.11, seed=1):
    """Receptor em corte feito de esferas, com um bolsao (poligono `pocket`, relativo ao topo)."""
    from matplotlib.path import Path as MPath
    rng = np.random.default_rng(seed)
    body = MPath(np.array([[-w / 2, -h], [w / 2, -h], [w / 2, 0], [-w / 2, 0]]))
    hole = MPath(pocket)
    test = lambda p: body.contains_points(p, radius=-0.02) & ~hole.contains_points(p, radius=0.12)
    pts = _beads_in(test, -w / 2 + r, w / 2 - r, -h + r, -r, 2 * r * 0.95, rng)
    depth = rng.uniform(0, 1, len(pts))
    P = np.c_[pts[:, 0] + x, pts[:, 1] + y, depth]
    from style import draw_spheres as ds
    lim = (ax.get_xlim(), ax.get_ylim())
    ds(ax, P, [color] * len(P), [r] * len(P), z=1)
    ax.set_xlim(*lim[0]); ax.set_ylim(*lim[1])


def _ligand(ax, x, y, shape, color=POSE, r=0.1, angle=0.0, alpha=1.0, seed=2):
    """Ligante feito de esferas dentro do poligono `shape`, girado por `angle` graus."""
    from matplotlib.path import Path as MPath
    rng = np.random.default_rng(seed)
    th = np.radians(angle)
    R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    poly = MPath(shape)
    lo, hi = shape.min(0), shape.max(0)
    pts = _beads_in(lambda p: poly.contains_points(p, radius=-0.01), lo[0] + r, hi[0] - r, lo[1] + r, hi[1] - r,
                    2 * r * 0.95, rng, jitter=0)
    pts = pts @ R.T
    for (px, py) in pts:
        ax.add_patch(Circle((x + px, y + py), r, fc=shade(color, 1 - alpha), ec=shade(INK, 0.6 + 0.4 * (1 - alpha)),
                            lw=0.3, zorder=5))


def fig_binding_models():
    fig, axs = plt.subplots(3, 1, figsize=(12, 9.6))
    sq = np.array([[-0.5, 0.2], [-0.5, -0.62], [0.5, -0.62], [0.5, 0.2]])
    vee = np.array([[-0.75, 0.2], [0, -0.75], [0.75, 0.2]])
    lig = np.array([[-0.42, 0.0], [-0.42, -0.55], [0.42, -0.55], [0.42, 0.0]])
    rows = [('Chave e fechadura (Fischer, 1894)', sq, sq, 'o bolsão já tem a forma do ligante', 'encaixe sem mudança'),
            ('Encaixe induzido (Koshland, 1958)', vee, sq, 'o bolsão tem outra forma', 'o bolsão muda ao receber o ligante'),
            ('Seleção conformacional', vee, sq, 'o receptor alterna entre A (em V) e B',
             'o ligante se liga à conformação B, que já existia')]
    for ax, (title, before, after, t1, t2) in zip(axs, rows):
        frame_off(ax); ax.set_aspect('equal'); ax.set_xlim(-1.8, 11.2); ax.set_ylim(-2.3, 1.45)
        ax.text(-1.7, 1.35, title, fontsize=10.5, fontweight='bold', va='top')
        _receptor(ax, 1.0, 0.0, before)
        if 'Seleção' in title:
            _receptor(ax, 4.2, 0.0, sq, seed=3)
            ax.add_patch(FancyArrowPatch((2.35, -0.8), (2.85, -0.8), arrowstyle='<|-|>', mutation_scale=12, color=GREY, lw=1.4))
            ax.text(1.0, -2.0, 'A', ha='center', fontsize=9); ax.text(4.2, -2.0, 'B', ha='center', fontsize=9)
            ax.text(2.6, 0.75, t1, ha='center', fontsize=8.5)
            _ligand(ax, 6.0, 0.9, lig)
        else:
            ax.text(1.0, -2.0, t1, ha='center', fontsize=8.5)
            _ligand(ax, 3.4, 0.9, lig)
            ax.text(3.4, 0.05, 'ligante', ha='center', fontsize=8.5)
        ax.add_patch(FancyArrowPatch((6.4, -0.6), (7.6, -0.6), arrowstyle='-|>', mutation_scale=14, color=INK, lw=1.2))
        _receptor(ax, 9.3, 0.0, after, seed=4)
        _ligand(ax, 9.3, 0.0, lig)
        ax.text(9.3, -2.0, t2, ha='center', fontsize=8.5)
    fig.suptitle('Três modelos de reconhecimento entre receptor e ligante (esquema)', y=0.99, fontsize=11)
    save(fig, 'fund_modelos_ligacao.png')


# ---------------------------------------------------------------- 2. vocabulario em estrutura real (8HBK)
def fig_vocabulary():
    t = md.load_pdb(str(fetch_pdb('8HBK')))
    lig = Chem.MolFromMolFile(str(ROOT / 'dados' / 'mpro' / 'ligante_cristal.sdf'), sanitize=False)
    L = lig.GetConformer().GetPositions() * A
    prot = heavy(t, 'protein')
    X = t.xyz[0]
    d = np.min(np.linalg.norm(X[prot][:, None] - L[None], axis=2), axis=1)
    site = set(prot[d < 0.45])
    box = json.loads((ROOT / 'dados' / 'mpro' / 'caixa.json').read_text())
    bc, bs = np.array(box['centro_A']) * A, np.array(box['tamanho_A']) * A
    look = X[prot].mean(0) - L.mean(0)
    B = view_basis(look, up_hint=np.array([0, 0, 1.0]))
    c0 = X[prot].mean(0)
    P = project(np.r_[X[prot], L], B, c0)
    cols = [SITE if i in site else PROT for i in prot] + [CRISTAL] * len(L)
    rad = [R_ATOM] * len(prot) + [0.2] * len(L)
    fig, ax = plt.subplots(figsize=(9, 7.4))
    draw_spheres(ax, P, cols, rad, pad=1.6)
    # caixa: projecao dos 8 vertices
    corners = np.array([[sx, sy, sz] for sx in (-.5, .5) for sy in (-.5, .5) for sz in (-.5, .5)]) * bs + bc
    Pc = project(corners, B, c0)
    for i in range(8):
        for j in range(i + 1, 8):
            if np.sum(np.abs(corners[i] - corners[j]) > 1e-6) == 1:
                ax.plot(Pc[[i, j], 0], Pc[[i, j], 1], color=INK, lw=0.9, ls=(0, (3, 2)), zorder=5)
    pl = front(P, range(len(prot), len(P)))
    xmin, xmax = P[:, 0].min(), P[:, 0].max()
    label(ax, pl, 'ligante (ensitrelvir)', (xmax + 0.9, pl[1] + 0.6), CRISTAL)
    vtx = Pc[np.argmax(Pc[:, 0] + Pc[:, 1]), :2]
    label(ax, vtx, 'caixa de busca', (xmax + 0.9, Pc[:, 1].max() + 0.4), INK)
    sidx = [k for k, i in enumerate(prot) if i in site]
    lx = P[len(prot):, 0].mean()
    sp = front(P, [k for k in sidx if P[k, 0] < lx - 0.3] or sidx)   # lado esquerdo: a linha nao cruza a do ligante
    label(ax, sp, 'sítio de ligação', (xmin - 0.9, sp[1] - 0.6), SITE)
    low = [k for k in range(len(prot)) if P[k, 1] < P[:len(prot), 1].min() + 0.6]
    label(ax, front(P, low), 'receptor (Mpro)', (xmin - 0.9, P[:, 1].min() + 0.4), GREY)
    ax.set_xlim(xmin - 2.2, xmax + 2.2)
    fit(fig, ax, largura=11, extra=0.5)
    ax.set_title('Vocabulário do docking numa estrutura real: Mpro do SARS-CoV-2 com ensitrelvir (PDB 8HBK)')
    save(fig, 'fund_vocabulario.png')


# ---------------------------------------------------------------- 3. poses, energia e RMSD (dados reais)
def _load_rows(nome):
    return json.loads((resultados(nome) / 'poses.json').read_text())


def _pose_xyz(nome, k):
    e, m = poses(resultados(nome) / 'vina_dock.pdbqt')[k - 1]
    m = skeleton(m)
    return m.GetConformer().GetPositions() * A, e


def _trp_rings(t, resids):
    out = {}
    for r in resids:
        idx = t.topology.select(f'resSeq {r} and name CD2 CE2 CE3 CZ2 CZ3 CH2')
        out[r] = t.xyz[0][idx].mean(0)
    return out


# ---------------------------------------------------------------- 4. RMSD explicado com dados reais
def fig_rmsd():
    """RMSD explicado com duas poses do DiffDock para o cap, contra o cristal 1IPC."""
    ref = skeleton(Chem.MolFromMolFile(str(ROOT / 'dados' / 'eif4e' / 'ligante_cristal.sdf'), sanitize=False, removeHs=False))
    P0 = ref.GetConformer().GetPositions()
    rows = _dd_rows('eif4e_diffdock')
    perto = min(rows, key=lambda r: r['rmsd_A'])['pose']; longe = max(rows, key=lambda r: r['rmsd_A'])['pose']
    sup = Chem.SDMolSupplier(str(ROOT / 'resultados' / 'eif4e_diffdock' / 'poses.sdf'), sanitize=False, removeHs=False)
    fig, axs = plt.subplots(1, 2, figsize=(13, 5))
    u, s, vt = np.linalg.svd(P0 - P0.mean(0)); B = np.c_[vt[0], vt[1], vt[2]]
    for ax, k in zip(axs, (perto, longe)):
        m = skeleton(sup[k - 1]); Q = m.GetConformer().GetPositions()
        match = min(ref.GetSubstructMatches(m, uniquify=False, useChirality=False, maxMatches=5000),
                    key=lambda mt: ((Q - P0[list(mt)]) ** 2).sum())
        Pm = P0[list(match)]
        a = (Pm - P0.mean(0)) @ B; b = (Q - P0.mean(0)) @ B
        for i in range(len(a)):
            ax.plot([a[i, 0], b[i, 0]], [a[i, 1], b[i, 1]], color=ERRO, lw=0.8, alpha=0.8)
        ax.scatter(a[:, 0], a[:, 1], s=40, color=CRISTAL, edgecolor=INK, lw=0.4, zorder=3)
        ax.scatter(b[:, 0], b[:, 1], s=40, color=POSE, edgecolor=INK, lw=0.4, zorder=3)
        rmsd = np.sqrt(((Q - Pm) ** 2).sum(1).mean())
        ax.set_aspect('equal'); frame_off(ax)
        ax.set_title(f'Pose {k} do DiffDock: RMSD = {rmsd:.2f} Å'.replace('.', ','), fontsize=10)
    legend(fig, [('átomo no cristal 1IPC', CRISTAL), ('mesmo átomo na pose', POSE)], y=-0.06)
    fig.text(0.5, 0.06, r'RMSD = $\sqrt{\frac{1}{N}\sum_i d_i^2}$, com $d_i$ a distância (linha vermelha) '
             'entre as duas posições do átomo $i$', ha='center', fontsize=10)
    fig.suptitle('RMSD: a distância média entre a pose e o cristal (cap do mRNA na eIF4E, átomos pesados)',
                 y=1.02, fontsize=11)
    save(fig, 'fund_rmsd.png')


# ---------------------------------------------------------------- 5. forcas (esquema)
def fig_forces():
    fig, axs = plt.subplots(2, 3, figsize=(13, 7.4))
    for ax in axs.flat:
        frame_off(ax); ax.set_aspect('equal'); ax.set_xlim(0, 6); ax.set_ylim(0, 4)
    def atom(ax, x, y, s, col, r=0.32, fc='white'):
        ax.add_patch(Circle((x, y), r, fc=col, ec=INK, lw=0.6, zorder=3))
        ax.text(x, y, s, ha='center', va='center', fontsize=9, color=fc, fontweight='bold', zorder=4)
    a = axs[0, 0]   # ligacao de hidrogenio
    atom(a, 1.4, 2, 'N', '#4c72b0'); atom(a, 2.25, 2, 'H', '#d9d9d9', 0.2, INK); atom(a, 4.3, 2, 'O', '#c44e52')
    a.plot([2.5, 3.95], [2, 2], color=INK, ls=':', lw=1.4)
    a.annotate('', xy=(4.3, 1.2), xytext=(1.4, 1.2), arrowprops=dict(arrowstyle='<->', color=GREY))
    a.text(2.85, 0.85, '2,7 a 3,2 Å entre N e O', ha='center', fontsize=8.5)
    a.set_title('Ligação de hidrogênio', fontsize=10)
    a = axs[0, 1]   # ponte salina
    atom(a, 1.5, 2, '+', '#4c72b0', 0.45); atom(a, 4.5, 2, '−', '#c44e52', 0.45)
    a.text(1.5, 1.1, 'Lys, Arg', ha='center', fontsize=8.5); a.text(4.5, 1.1, 'Asp, Glu, fosfato', ha='center', fontsize=8.5)
    a.plot([2.0, 4.0], [2, 2], color=INK, ls=':', lw=1.4)
    a.set_title('Interação iônica (ponte salina)', fontsize=10)
    a = axs[0, 2]   # hidrofobico
    for x, y in ((1.6, 2.3), (1.9, 1.7), (1.3, 1.6)):
        atom(a, x, y, '', '#9a9a9a', 0.33)
    for x, y in ((4.1, 2.3), (4.4, 1.7), (3.8, 1.6)):
        atom(a, x, y, '', '#9a9a9a', 0.33)
    for x, y in ((2.6, 3.0), (3.4, 3.2), (3.0, 3.5)):
        a.add_patch(Circle((x, y), 0.17, fc=WATER, ec=INK, lw=0.4))
    a.annotate('', xy=(3.0, 2.85), xytext=(3.0, 2.05), arrowprops=dict(arrowstyle='-|>', color=C[0]))
    a.text(3.0, 0.45, 'os grupos apolares se juntam\ne liberam as águas que os cercavam', ha='center', fontsize=8.5)
    a.set_title('Efeito hidrofóbico', fontsize=10)
    a = axs[1, 0]   # empilhamento pi
    for y in (1.4, 2.6):
        hexa = np.array([[np.cos(t), 0.35 * np.sin(t)] for t in np.linspace(0, 2 * np.pi, 7)[:-1]]) * 1.1 + [3, y]
        a.add_patch(Polygon(hexa, closed=True, fc='#e6d9a8', ec=INK, lw=0.8))
    a.annotate('', xy=(4.6, 2.6), xytext=(4.6, 1.4), arrowprops=dict(arrowstyle='<->', color=GREY))
    a.text(4.8, 2.0, '3,4 a 4 Å', fontsize=8.5, va='center')
    a.set_title('Empilhamento π (anéis aromáticos)', fontsize=10)
    a = axs[1, 1]   # cation-pi
    hexa = np.array([[np.cos(t), 0.35 * np.sin(t)] for t in np.linspace(0, 2 * np.pi, 7)[:-1]]) * 1.1 + [3, 1.4]
    a.add_patch(Polygon(hexa, closed=True, fc='#e6d9a8', ec=INK, lw=0.8))
    atom(a, 3, 2.8, '+', '#4c72b0', 0.4)
    a.plot([3, 3], [1.75, 2.4], color=INK, ls=':', lw=1.4)
    a.text(3, 0.6, 'ex.: m7G+ entre dois Trp do eIF4E', ha='center', fontsize=8.5)
    a.set_title('Cátion-π', fontsize=10)
    a = axs[1, 2]   # dessolvatacao
    atom(a, 1.6, 2.4, 'L', POSE, 0.45)
    for t in np.linspace(0, 2 * np.pi, 9)[:-1]:
        a.add_patch(Circle((1.6 + 0.8 * np.cos(t), 2.4 + 0.8 * np.sin(t)), 0.15, fc=WATER, ec=INK, lw=0.4))
    a.annotate('', xy=(3.7, 2.4), xytext=(2.7, 2.4), arrowprops=dict(arrowstyle='-|>', color=INK))
    atom(a, 4.6, 2.4, 'L', POSE, 0.45)
    a.add_patch(Rectangle((4.0, 1.65), 1.2, 0.3, fc=SITE, ec=INK, lw=0.6))
    a.text(3.0, 0.5, 'perder a água custa energia;\no ganho tem que pagar esse custo', ha='center', fontsize=8.5)
    a.set_title('Dessolvatação', fontsize=10)
    fig.suptitle('Forças que seguram o ligante no sítio (esquema)', y=0.98, fontsize=11)
    save(fig, 'fund_forcas.png')


# ---------------------------------------------------------------- 6. tipos de docking
def fig_flexibility():
    fig, axs = plt.subplots(1, 4, figsize=(18, 3.9))
    sq = np.array([[-0.5, 0.2], [-0.5, -0.62], [0.5, -0.62], [0.5, 0.2]])
    lig = np.array([[-0.42, 0.0], [-0.42, -0.55], [0.42, -0.55], [0.42, 0.0]])
    titles = ['Rígido-rígido\nnada se move', 'Rígido-flexível (semiflexível)\nligante flexível, receptor rígido',
              'Flexibilidade parcial\ncadeias laterais do sítio se movem', 'Flexível-flexível\nligante e receptor inteiros se ajustam']
    notas = ['o ligante só translada e gira', 'esferas claras: outras conformações testadas',
             'só os resíduos da borda do bolsão', 'paredes e fundo do bolsão se movem']
    for k, (ax, title) in enumerate(zip(axs, titles)):
        frame_off(ax); ax.set_aspect('equal'); ax.set_xlim(-1.8, 1.8); ax.set_ylim(-2.0, 0.55)
        _receptor(ax, 0.0, 0.0, sq, seed=5 + k)
        if k == 0:
            _ligand(ax, 0, 0, lig)
        else:
            for ang, al in ((-25, 0.35), (25, 0.35), (0, 1.0)):
                _ligand(ax, 0, 0.1 if ang else 0, lig, angle=ang, alpha=al)
        if k == 2:
            for dx in (-1, 1):
                ax.add_patch(FancyArrowPatch((dx * 0.62, 0.32), (dx * 0.62, 0.05), arrowstyle='<|-|>',
                                             mutation_scale=9, color=DESTAQUE, lw=1.3, zorder=6))
        if k == 3:
            for dx in (-1, 1):
                ax.add_patch(FancyArrowPatch((dx * 0.95, -0.35), (dx * 0.62, -0.35), arrowstyle='<|-|>',
                                             mutation_scale=10, color=DESTAQUE, lw=1.4, zorder=6))
            ax.add_patch(FancyArrowPatch((0, -1.0), (0, -0.7), arrowstyle='<|-|>', mutation_scale=10, color=DESTAQUE,
                                         lw=1.4, zorder=6))
        ax.text(0, -1.85, notas[k], ha='center', fontsize=8.5, color=DESTAQUE if k >= 2 else GREY)
        ax.set_title(title, fontsize=10)
    fig.suptitle('Tipos de docking pela flexibilidade, do mais simples ao mais completo (esquema)', y=1.04, fontsize=11)
    save(fig, 'fund_tipos_flexibilidade.png')


def fig_search():
    t = md.load_pdb(str(fetch_pdb('8HBK')))
    L = Chem.MolFromMolFile(str(ROOT / 'dados' / 'mpro' / 'ligante_cristal.sdf'), sanitize=False).GetConformer().GetPositions() * A
    prot = heavy(t, 'protein')
    X = t.xyz[0]
    B = view_basis(X[prot].mean(0) - L.mean(0)); c0 = X[prot].mean(0)
    P = project(np.r_[X[prot], L], B, c0)
    cols = [PROT] * len(prot) + [CRISTAL] * len(L)
    rad = [R_ATOM] * len(prot) + [0.18] * len(L)
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.2))
    pl = P[len(prot):, :2]
    titles = ['Baseado em molde: parte de um complexo\nconhecido e parecido', 'Local: caixa só no sítio conhecido',
              'Cego (global): a proteína inteira,\no programa procura o sítio (DiffDock, esta prática)']
    for k, (ax, title) in enumerate(zip(axs, titles)):
        draw_spheres(ax, P if k else P, cols if k != 2 else [PROT] * len(prot) + [shade(CRISTAL, 0.5)] * len(L), rad, pad=1.0)
        if k == 1:
            lo, hi = pl.min(0) - 0.8, pl.max(0) + 0.8
            ax.add_patch(Rectangle(lo, *(hi - lo), fill=False, ec=INK, lw=1.2, ls=(0, (3, 2)), zorder=6))
        if k == 2:
            lo, hi = P[:, :2].min(0) - 0.3, P[:, :2].max(0) + 0.3
            ax.add_patch(Rectangle(lo, *(hi - lo), fill=False, ec=INK, lw=1.2, ls=(0, (3, 2)), zorder=6))
        if k == 0:
            ax.text(pl[:, 0].mean(), pl[:, 1].max() + 0.9, 'ligante do molde', ha='center', fontsize=8.5,
                    bbox=dict(boxstyle='round,pad=0.2', fc='white', ec=CRISTAL), zorder=7)
        ax.set_title(title, fontsize=10)
    fig.suptitle('Tipos de docking pelo que se sabe do sítio (Mpro, PDB 8HBK)', y=0.99, fontsize=11)
    save(fig, 'fund_tipos_busca.png')


# ---------------------------------------------------------------- 7. niveis e custo
def fig_levels_cost():
    # ordem de grandeza do tempo por ligante (s); faixas indicativas, citadas no texto
    levels = [('Modelo de ML (DiffDock), GPU', 5, 60, C[4]),
              ('Docking com receptor rígido (AutoDock), 1 CPU', 30, 600, C[0]),
              ('Docking com cadeias laterais flexíveis', 600, 3600 * 2, C[6]),
              ('Ensemble docking (várias conformações da MD)', 1800, 3600 * 10, C[3]),
              ('MD do complexo + MM/GBSA, GPU', 3600 * 5, 3600 * 48, C[1]),
              ('Perturbação de energia livre (FEP), GPU', 3600 * 12, 3600 * 24 * 4, C[7])]
    fig, ax = plt.subplots(figsize=(11, 4.6))
    for k, (name, lo, hi, col) in enumerate(levels):
        y = len(levels) - 1 - k
        ax.plot([lo, hi], [y, y], color=col, lw=9, solid_capstyle='butt')
        ax.text(lo, y + 0.32, name, ha='left', va='bottom', fontsize=9)
    ax.plot(15, len(levels) - 1, 'o', color='white', mec=INK, ms=8, zorder=5)
    ax.annotate('medido nesta aula:\nDiffDock, 5 a 25 s', (15, len(levels) - 1), xytext=(4000, len(levels) - 0.75),
                ha='left', va='center', fontsize=8.5, arrowprops=dict(arrowstyle='-', color=INK, lw=0.8))
    ax.set_xscale('log'); ax.set_xlim(2, 3600 * 24 * 10); ax.set_ylim(-0.5, len(levels) + 0.1)
    ax.set_yticks([])
    for sp in ('left', 'right', 'top'):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis='x', top=False)
    ax.set_xticks([1, 60, 3600, 3600 * 24]); ax.set_xticklabels(['1 s', '1 min', '1 h', '1 dia'])
    ax.tick_params(axis='x', which='minor', bottom=False, top=False)
    ax.set_xlabel('Tempo por ligante (ordem de grandeza, escala logarítmica)')
    ax.set_title('Níveis de modelagem da ligação: mais física, mais custo')
    save(fig, 'fund_niveis_custo.png')


# ---------------------------------------------------------------- 8. SMILES
def fig_smiles():
    mols = [('Aspirina', 'CC(=O)Oc1ccccc1C(=O)O'),
            ('Ensitrelvir', (ROOT / 'dados' / 'mpro' / 'ligante.smi').read_text().split()[0]),
            ('m7GTP (cap do mRNA)', (ROOT / 'dados' / 'eif4e' / 'ligante.smi').read_text().split()[0])]
    fig, axs = plt.subplots(1, 3, figsize=(15, 5.4))
    for ax, (name, smi) in zip(axs, mols):
        m = Chem.MolFromSmiles(smi)
        img = Draw.MolToImage(m, size=(520, 420))
        ax.imshow(img); frame_off(ax)
        ax.set_title(name, fontsize=10.5)
        wrap = '\n'.join(smi[i:i + 38] for i in range(0, len(smi), 38))
        ax.text(0.5, -0.04, wrap, transform=ax.transAxes, ha='center', va='top', fontsize=8.5, family='monospace')
    fig.suptitle('A mesma molécula em desenho e em SMILES', y=1.0, fontsize=11)
    save(fig, 'fund_smiles.png')


# ---------------------------------------------------------------- 9. exemplos: resistencia, PPI
def _scene_with_ligand(ax, code, comp, chain, hl, hl_color, cut=1.0, lig_col=CRISTAL, title='', labels=()):
    t = md.load_pdb(str(fetch_pdb(code)))
    L, _ = ligand_sdf(code, comp, chain)
    prot = heavy(t, f'protein and chainid 0')
    X = t.xyz[0]
    d = np.min(np.linalg.norm(X[prot][:, None] - L[None], axis=2), axis=1)
    near = prot[d < cut]
    hli = set(t.topology.select(f'chainid 0 and resSeq {hl} and not element H'))
    B = view_basis(X[prot].mean(0) - L.mean(0)); c0 = L.mean(0)
    P = project(np.r_[X[near], L], B, c0)
    cols = [hl_color if i in hli else PROT for i in near] + [lig_col] * len(L)
    rad = [R_ATOM] * len(near) + [0.17] * len(L)
    draw_spheres(ax, P[:len(near)], cols[:len(near)], rad[:len(near)], pad=0.5)
    draw_spheres(ax, P[len(near):], cols[len(near):], rad[len(near):], pad=0.5, z=3)
    hk = [k for k, i in enumerate(near) if i in hli]
    draw_spheres(ax, P[hk], [hl_color] * len(hk), [R_ATOM] * len(hk), pad=0.5, z=4)   # destaque por cima
    ax.set_xlim(P[:, 0].min() - 0.3, P[:, 0].max() + 0.3); ax.set_ylim(P[:, 1].min() - 0.3, P[:, 1].max() + 0.3)
    ax.set_title(title, fontsize=10)
    hp = front(P, hk)
    return P, hp


def fig_resistance():
    fig, ax = plt.subplots(figsize=(10, 6.6))
    P, hp = _scene_with_ligand(ax, '1IEP', 'STI', 'A', 315, DESTAQUE, cut=0.9,
                               title='Imatinibe no sítio da quinase ABL (PDB 1IEP) e a Thr315')
    label(ax, hp, 'Thr315: na mutação T315I\nvira Ile, maior e\nsem a hidroxila', (P[:, 0].max() + 1.3, hp[1] - 0.6), DESTAQUE)
    ax.set_xlim(P[:, 0].min() - 0.3, P[:, 0].max() + 2.4)
    fit(fig, ax, largura=11, extra=0.5)
    legend(fig, [('imatinibe', CRISTAL), ('Thr315', DESTAQUE), ('ABL: átomos a menos de 9 Å do fármaco', PROT)], y=-0.08, ax=ax)
    save(fig, 'fund_resistencia_abl.png')


def _superpose(mobile, target, sel_m, sel_t):
    """Sobrepoe `mobile` em `target` pelos CA de MDM2 (Kabsch)."""
    a = mobile.xyz[0][sel_m]; b = target.xyz[0][sel_t]
    ca, cb = a.mean(0), b.mean(0)
    H = (a - ca).T @ (b - cb)
    U, _, Vt = np.linalg.svd(H)
    D = np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))])
    R = Vt.T @ D @ U.T
    return lambda x: (x - ca) @ R.T + cb


def fig_p53_mdm2():
    t1 = md.load_pdb(str(fetch_pdb('1YCR')))
    t2 = md.load_pdb(str(fetch_pdb('4HG7')))
    # residuos de MDM2 presentes nas duas estruturas
    r1 = {r.resSeq: r for r in t1.topology.chain(0).residues}
    r2 = {r.resSeq: r for r in t2.topology.chain(0).residues}
    common = sorted(set(r1) & set(r2))
    ca1 = [r1[i].atom('CA').index for i in common]; ca2 = [r2[i].atom('CA').index for i in common]
    move = _superpose(t2, t1, ca2, ca1)
    mdm2 = heavy(t1, 'chainid 0')
    pep = heavy(t1, 'chainid 1')
    hot = set(t1.topology.select('chainid 1 and resSeq 19 23 26 and not element H'))
    N, _ = ligand_sdf('4HG7', 'NUT', 'A'); N = move(N)
    X = t1.xyz[0]
    B = view_basis(X[mdm2].mean(0) - X[pep].mean(0)); c0 = X[mdm2].mean(0)
    fig, axs = plt.subplots(1, 2, figsize=(13, 6))
    P = project(np.r_[X[mdm2], X[pep]], B, c0)
    cols = [PROT] * len(mdm2) + [ERRO if i in hot else SITE for i in pep]
    draw_spheres(axs[0], P[:len(mdm2)], cols[:len(mdm2)], [R_ATOM] * len(mdm2), pad=0.5)
    draw_spheres(axs[0], P[len(mdm2):], cols[len(mdm2):], [R_ATOM] * len(pep), pad=0.5, z=3)
    hk = [len(mdm2) + j for j, i in enumerate(pep) if i in hot]
    draw_spheres(axs[0], P[hk], [ERRO] * len(hk), [R_ATOM] * len(hk), pad=0.5, z=4)
    axs[0].set_xlim(P[:, 0].min() - 0.3, P[:, 0].max() + 0.3); axs[0].set_ylim(P[:, 1].min() - 0.3, P[:, 1].max() + 0.3)
    axs[0].set_title('MDM2 com a hélice da p53 (PDB 1YCR):\nPhe19, Trp23 e Leu26 entram no sulco', fontsize=10)
    Q = project(np.r_[X[mdm2], N], B, c0)
    draw_spheres(axs[1], Q[:len(mdm2)], [PROT] * len(mdm2), [R_ATOM] * len(mdm2), pad=0.5)
    draw_spheres(axs[1], Q[len(mdm2):], [CRISTAL] * len(N), [0.17] * len(N), pad=0.5, z=3)
    axs[1].set_xlim(P[:, 0].min() - 0.3, P[:, 0].max() + 0.3); axs[1].set_ylim(P[:, 1].min() - 0.3, P[:, 1].max() + 0.3)
    axs[1].set_title('MDM2 com nutlina-3a (PDB 4HG7, sobreposta):\no fármaco ocupa o mesmo sulco', fontsize=10)
    legend(fig, [('MDM2', PROT), ('hélice da p53', SITE), ('Phe19, Trp23, Leu26', ERRO), ('nutlina-3a', CRISTAL)])
    fig.suptitle('Proteína-proteína e proteína-fármaco no mesmo sítio: MDM2 e p53', y=1.0, fontsize=11)
    save(fig, 'fund_p53_mdm2.png')


def fig_pd1_pdl1():
    """Docking proteina-proteina com um exemplo de imunoterapia: PD-1 e PD-L1 (PDB 4ZQK)."""
    t = md.load_pdb(str(fetch_pdb('4ZQK')))
    a = heavy(t, 'chainid 0 and protein'); b = heavy(t, 'chainid 1 and protein')
    X = t.xyz[0]
    d = np.min(np.linalg.norm(X[a][:, None] - X[b][None], axis=2), axis=1)
    ia = set(a[d < 0.45])
    d2 = np.min(np.linalg.norm(X[b][:, None] - X[a][None], axis=2), axis=1)
    ib = set(b[d2 < 0.45])
    axis = X[b].mean(0) - X[a].mean(0)
    B = view_basis(np.cross(axis, [0, 0, 1.0]), up_hint=axis / np.linalg.norm(axis)); c0 = X[np.r_[a, b]].mean(0)
    P = project(X[np.r_[a, b]], B, c0)
    cols = [C[0] if i in ia else '#b8c9e0' for i in a] + [C[5] if i in ib else '#ecc7a8' for i in b]
    fig, ax = plt.subplots(figsize=(7.6, 7.2))
    draw_spheres(ax, P, cols, [R_ATOM] * len(P), pad=0.6)
    ax.set_title('Docking proteína-proteína: PD-1 e PD-L1 (PDB 4ZQK)', fontsize=10.5)
    legend(fig, [('PD-L1 (proteína A)', '#b8c9e0'), ('PD-1 (proteína B)', '#ecc7a8'), ('interface do PD-L1', C[0]),
                 ('interface do PD-1', C[5])], ncol=2, y=0.0)
    save(fig, 'fund_pd1_pdl1.png')


# ---------------------------------------------------------------- 10. docking precisa de MD (esquema)
def _blob(ax, cx, cy, width, depth=0.55, scale=1.0, color=SITE):
    th = np.linspace(-np.pi, np.pi, 400)
    r = 1 + 0.07 * np.sin(5 * th + 0.6) + 0.04 * np.cos(9 * th) - depth * np.exp(-(th / width) ** 2)
    ax.fill(cx + scale * r * np.cos(th), cy + scale * r * np.sin(th), color=color, alpha=0.35, lw=0)
    ax.plot(cx + scale * r * np.cos(th), cy + scale * r * np.sin(th), color=C[0], lw=1.4)


def fig_md_docking():
    fig, ax = plt.subplots(figsize=(14, 4.6))
    frame_off(ax); ax.set_aspect('equal'); ax.set_xlim(-1.4, 18.6); ax.set_ylim(-2.6, 1.8)
    widths = [0.05, 0.18, 0.4, 0.12]
    for k, w in enumerate(widths):
        cx = 0.6 + 2.6 * k
        _blob(ax, cx, 0, w)
        ax.text(cx, 1.35, f'quadro {k + 1} da MD', ha='center', fontsize=8.5)
        if k == 2:
            ax.add_patch(Circle((cx + 0.72, 0), 0.17, fc=POSE, ec=INK, lw=0.6, zorder=3))
    ax.text(4.5, -1.5, '1. a MD gera várias conformações do receptor;\n2. o docking roda em cada uma (ensemble docking)',
            ha='center', fontsize=9)
    ax.add_patch(FancyArrowPatch((10.4, 0), (11.6, 0), arrowstyle='-|>', mutation_scale=14, color=INK))
    x = np.linspace(0, 1, 200)
    ox, oy, w, h = 12.2, -1.0, 4.4, 2.2
    ax.add_patch(Rectangle((ox, oy), w, h, fill=False, ec=INK, lw=0.8))
    ax.plot(ox + w * x, oy + 0.25 + 0.12 * h * (1 - np.exp(-x * 8)) + 0.03 * np.sin(40 * x), color=C[3], lw=1.6)
    ax.plot(ox + w * x, oy + 0.25 + 0.75 * h * np.clip((x - 0.35) * 2.5, 0, 1) + 0.04 * np.sin(33 * x), color=ERRO, lw=1.6,
            ls=(0, (6, 3)))
    ax.text(ox + w / 2, oy - 0.35, 'tempo de MD do complexo', ha='center', fontsize=8.5)
    ax.text(ox + 0.1, oy + h - 0.1, 'RMSD do ligante', ha='left', va='top', fontsize=8.5)
    ax.text(ox + w + 0.15, oy + 0.45, 'pose estável', ha='left', va='center', fontsize=8.5, color=C[3])
    ax.text(ox + w + 0.15, oy + h - 0.25, 'pose instável:\no ligante sai', ha='left', va='center', fontsize=8.5, color=ERRO)
    ax.text(ox + w / 2, oy + h + 0.25, '3. a MD testa a pose escolhida', ha='center', fontsize=9)
    fit(fig, ax, largura=14, extra=0.6)
    fig.suptitle('Como a dinâmica molecular entra antes e depois do docking (esquema)', y=1.0, fontsize=11)
    save(fig, 'fund_md_docking.png')


# ---------------------------------------------------------------- 11. metricas com os dados da pratica
RT = 0.593   # kcal/mol a 298 K


def _interactions(nome):
    return json.loads((resultados(nome) / 'interacoes.json').read_text())


def _pairs(rows):
    return {(r['tipo'], r['residuo']) for r in rows}


def fig_metrics():
    """Tres numeros para julgar uma pose, com as poses do DiffDock nos dois sistemas."""
    ra, rb = _dd_rows('eif4e_diffdock'), _dd_rows('mpro_diffdock')
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.6))
    a = axs[0]
    a.axvspan(-4, -1.5, color=shade(ERRO, 0.8), lw=0); a.axvspan(-1.5, 0, color=shade(C[1], 0.8), lw=0)
    a.axvspan(0, 1.5, color=shade(C[3], 0.75), lw=0)
    for x, t in ((-2.75, 'baixa'), (-0.75, 'moderada'), (0.75, 'alta')):
        a.text(x, 0.94, t, ha='center', va='top', fontsize=9, transform=a.get_xaxis_transform())
    for v, lab, cor, y in ((ra[0]['confianca'], 'eIF4E, pose 1', C[7], 0.62), (rb[0]['confianca'], 'Mpro, pose 1', C[0], 0.38)):
        a.plot(v, y, 'o', color=cor, mec=INK, ms=9, transform=a.get_xaxis_transform())
        a.text(v - 0.12, y, f'{lab}: {v:.2f}'.replace('.', ',').replace('-', '−'), ha='right', va='center', fontsize=8.5,
               transform=a.get_xaxis_transform())
    a.set_xlim(-4, 1.5); a.set_yticks([]); a.set_xlabel('Confiança do DiffDock (maior = melhor)')
    a.set_title('Confiança: a nota do próprio modelo', fontsize=10)
    b = axs[1]
    b.axvspan(0, 2, color=shade(C[3], 0.75), lw=0); b.axvspan(2, 3, color=shade(C[1], 0.8), lw=0)
    b.axvspan(3, 10, color=shade(ERRO, 0.8), lw=0)
    for x, t in ((1, 'correta'), (2.5, 'quase'), (6.5, 'errada')):
        b.text(x, 0.94, t, ha='center', va='top', fontsize=9, transform=b.get_xaxis_transform())
    errada = max(rb[:10], key=lambda r: r['rmsd_nucleo_A'])
    for r, lab, cor, y in ((ra[0]['rmsd_A'], 'eIF4E\npose 1', C[7], 0.62), (rb[0]['rmsd_nucleo_A'], 'Mpro\npose 1', C[0], 0.3),
                           (errada['rmsd_nucleo_A'], f"Mpro\npose {errada['pose']}", C[0], 0.62)):
        b.plot(r, y, 'o', color=cor, mec=INK, ms=9, transform=b.get_xaxis_transform())
        b.text(r + 0.2, y, f'{lab}: {r:.1f} Å'.replace('.', ','), fontsize=8.5, va='center', transform=b.get_xaxis_transform())
    b.set_xlim(0, 10); b.set_yticks([]); b.set_xlabel('RMSD em relação ao cristal (Å)')
    b.set_title('RMSD: só existe quando há um cristal', fontsize=10)
    c = axs[2]
    vals = [(ra[0]['contatos_do_cristal'] / ra[0]['contatos_no_cristal'], 'eIF4E pose 1', C[7]),
            (rb[0]['contatos_do_cristal'] / rb[0]['contatos_no_cristal'], 'Mpro pose 1', C[0]),
            (errada['contatos_do_cristal'] / errada['contatos_no_cristal'], f"Mpro pose {errada['pose']}", C[0])]
    y = np.arange(len(vals))[::-1]
    c.barh(y, [100 * v for v, _, _ in vals], color=[k for _, _, k in vals], edgecolor=INK, lw=0.5, height=0.55)
    for yy, (v, _, _) in zip(y, vals):
        c.text(100 * v + 2, yy, f'{100 * v:.0f}%', va='center', fontsize=9)
    c.set_yticks(y); c.set_yticklabels([n for _, n, _ in vals], fontsize=9); c.tick_params(axis='y', right=False)
    c.set_xlim(0, 115); c.set_xlabel('Contatos do cristal que a pose repete (%)')
    c.set_title('Interações: quais contatos se repetem', fontsize=10)
    fig.suptitle('Três números para julgar uma pose do DiffDock', y=1.02, fontsize=11)
    fig.tight_layout(w_pad=2.5)
    save(fig, 'fund_metricas.png')


# ---------------------------------------------------------------- 12. cavidades e tipos de sitio
def _prankweb_raw(code):
    """Resposta do PrankWeb com os residuos de cada cavidade, guardada em runs/pdb."""
    path = CACHE / f'{code}_pockets.json'
    if not path.exists():
        CACHE.mkdir(parents=True, exist_ok=True)
        url = f'https://prankweb.cz/api/v2/prediction/v4-conservation-hmm/{code}/public/prediction.json'
        path.write_bytes(urllib.request.urlopen(url).read())
    return json.loads(path.read_text())


def _pockets(code):
    return json.loads((ROOT / 'resultados' / 'cavidades' / f'{code}.json').read_text())['cavidades']


def _pocket_atoms(t, pocket, chain=0):
    """Atomos pesados dos residuos da cavidade (indices 'A_145' do PrankWeb)."""
    nums = sorted({int(r.split('_')[1]) for r in pocket['residues']}) if 'residues' in pocket else []
    sel = ' '.join(str(n) for n in nums)
    return t.topology.select(f'chainid {chain} and resSeq {sel} and not element H') if nums else np.array([], int)


def _scene(t, chain=0, lig=None):
    prot = heavy(t, f'chainid {chain} and protein')
    X = t.xyz[0]
    centro = X[prot].mean(0)
    look = centro - (lig if lig is not None else X[prot][np.argmin(X[prot][:, 2])])
    return prot, X, view_basis(look), centro


def fig_cavidades():
    """Mpro: centro de cada cavidade prevista pelo P2Rank; so a primeira e o sitio real."""
    code = '8HBK'
    t = md.load_pdb(str(fetch_pdb(code)))
    L = Chem.MolFromMolFile(str(ROOT / 'dados' / 'mpro' / 'ligante_cristal.sdf'), sanitize=False).GetConformer().GetPositions() * A
    pk = _pockets(code)
    prot = heavy(t, 'chainid 0 and protein')
    X = t.xyz[0]
    B = view_basis(X[prot].mean(0) - L.mean(0)); c0 = X[prot].mean(0)
    P = project(np.r_[X[prot], L], B, c0)
    fig, ax = plt.subplots(figsize=(11, 7))
    draw_spheres(ax, P[:len(prot)], [PROT] * len(prot), [R_ATOM] * len(prot), pad=1.0)
    draw_spheres(ax, P[len(prot):], [POSE] * len(L), [0.17] * len(L), pad=1.0, z=3)
    xmin, xmax, ymin, ymax = P[:, 0].min(), P[:, 0].max(), P[:, 1].min(), P[:, 1].max()
    esq, dir_ = [], []
    for meta in pk:
        c = project(np.array([meta['centro']]) * A, B, c0)[0]
        col = CRISTAL if meta['distancia_A'] <= 5 else ERRO
        # raio proporcional ao numero de residuos da cavidade
        ax.add_patch(Circle(c[:2], 0.1 * np.sqrt(meta['n_residuos']), fc=col, alpha=0.5, ec=col, lw=1.4, zorder=4))
        (esq if c[0] < 0 else dir_).append((meta, c))
    for grupo, x_lab in ((esq, xmin - 2.2), (dir_, xmax + 2.2)):
        grupo.sort(key=lambda g: -g[1][1])
        passo = (ymax - ymin) / max(len(grupo), 1)
        for i, (meta, c) in enumerate(grupo):
            col = CRISTAL if meta['distancia_A'] <= 5 else ERRO
            texto = (f"cavidade {meta['rank']}: sítio ativo\nescore {meta['score']:.1f}" if meta['distancia_A'] <= 5
                     else f"cavidade {meta['rank']}: escore {meta['score']:.1f}\na {meta['distancia_A']:.0f} Å do ligante")
            label(ax, c[:2], texto, (x_lab, ymax - 0.5 - i * passo), col)
    ax.set_xlim(xmin - 5.0, xmax + 5.0); ax.set_ylim(ymin - 0.6, ymax + 0.6)
    fit(fig, ax, largura=13, extra=0.6)      # altura pela proporcao dos dados: sem faixa em branco
    legend(fig, [('cavidade que é o sítio real', CRISTAL), ('cavidade sem ligante conhecido', ERRO),
                 ('ensitrelvir no cristal', POSE), ('Mpro', PROT)], y=-0.07, ax=ax)
    ax.set_title('O P2Rank encontra 8 cavidades na Mpro; só a primeira é o sítio ativo (PDB 8HBK).\n'
                 'O tamanho do círculo acompanha o número de resíduos da cavidade', fontsize=10.5)
    save(fig, 'fund_cavidades.png')


def fig_tipos_sitio():
    """Quatro tipos de sitio em estruturas reais."""
    fig, axs = plt.subplots(1, 4, figsize=(17, 5.4))
    # (a) ortosterico e (b) alosterico: ABL com nilotinibe (sitio do ATP) e asciminibe (bolso miristoil)
    t = md.load_pdb(str(fetch_pdb('5MO4')))
    nil, _ = ligand_sdf('5MO4', 'NIL', 'A')
    ay7, _ = ligand_sdf('5MO4', 'AY7', 'A')
    prot = heavy(t, 'chainid 0 and protein')
    X = t.xyz[0]
    B = view_basis(np.cross(ay7.mean(0) - nil.mean(0), [0, 0, 1.0]), up_hint=(ay7.mean(0) - nil.mean(0)))
    c0 = X[prot].mean(0)
    P = project(np.r_[X[prot], nil, ay7], B, c0)
    n = len(prot)
    for ax, sel, col, title in ((axs[0], slice(n, n + len(nil)), CRISTAL,
                                'Ortostérico: nilotinibe no sítio do ATP\n(ABL, PDB 5MO4)'),
                                (axs[1], slice(n + len(nil), None), DESTAQUE,
                                 'Alostérico: asciminibe no bolso do miristoil,\nna mesma proteína')):
        draw_spheres(ax, P[:n], [PROT] * n, [R_ATOM] * n, pad=0.5)
        draw_spheres(ax, P[sel], [col] * len(P[sel]), [0.18] * len(P[sel]), pad=0.5, z=3)
        ax.set_xlim(P[:, 0].min() - 0.3, P[:, 0].max() + 0.3); ax.set_ylim(P[:, 1].min() - 0.3, P[:, 1].max() + 0.3)
        ax.set_title(title, fontsize=9.5)
    # (c) criptico: KRAS com sotorasibe
    t2 = md.load_pdb(str(fetch_pdb('6OIM')))
    mov, _ = ligand_sdf('6OIM', 'MOV', 'A')
    gdp, _ = ligand_sdf('6OIM', 'GDP', 'A')
    p2 = heavy(t2, 'chainid 0 and protein')
    X2 = t2.xyz[0]
    B2 = view_basis(X2[p2].mean(0) - mov.mean(0)); c2 = X2[p2].mean(0)
    Q = project(np.r_[X2[p2], mov, gdp], B2, c2)
    m = len(p2)
    draw_spheres(axs[2], Q[:m], [PROT] * m, [R_ATOM] * m, pad=0.5)
    draw_spheres(axs[2], Q[m:m + len(mov)], [C[6]] * len(mov), [0.18] * len(mov), pad=0.5, z=3)
    draw_spheres(axs[2], Q[m + len(mov):], [GREY] * len(gdp), [0.16] * len(gdp), pad=0.5, z=3)
    axs[2].set_xlim(Q[:, 0].min() - 0.3, Q[:, 0].max() + 0.3); axs[2].set_ylim(Q[:, 1].min() - 0.3, Q[:, 1].max() + 0.3)
    axs[2].set_title('Críptico: o bolso switch II do KRAS G12C\nse abre com o sotorasibe (PDB 6OIM)', fontsize=9.5)
    # (d) interface proteina-proteina: MDM2 com a helice da p53
    t3 = md.load_pdb(str(fetch_pdb('1YCR')))
    mdm2 = heavy(t3, 'chainid 0'); pep = heavy(t3, 'chainid 1')
    X3 = t3.xyz[0]
    B3 = view_basis(X3[mdm2].mean(0) - X3[pep].mean(0)); c3 = X3[mdm2].mean(0)
    R = project(np.r_[X3[mdm2], X3[pep]], B3, c3)
    draw_spheres(axs[3], R[:len(mdm2)], [PROT] * len(mdm2), [R_ATOM] * len(mdm2), pad=0.5)
    draw_spheres(axs[3], R[len(mdm2):], [SITE] * len(pep), [R_ATOM] * len(pep), pad=0.5, z=3)
    axs[3].set_xlim(R[:, 0].min() - 0.3, R[:, 0].max() + 0.3); axs[3].set_ylim(R[:, 1].min() - 0.3, R[:, 1].max() + 0.3)
    axs[3].set_title('Interface proteína-proteína: o sulco do MDM2\nque recebe a hélice da p53 (PDB 1YCR)', fontsize=9.5)
    legend(fig, [('ligante ortostérico', CRISTAL), ('ligante alostérico', DESTAQUE), ('inibidor em bolso críptico', C[6]),
                 ('GDP (sítio do nucleotídeo)', GREY), ('parceiro proteico', SITE), ('receptor', PROT)], ncol=3, y=-0.04)
    fig.suptitle('Tipos de sítio: onde um ligante pode se encaixar', y=1.0, fontsize=11)
    save(fig, 'fund_tipos_sitio.png')


# ---------------------------------------------------------------- 13. fluxo da pratica
def fig_fluxo():
    """Os nove passos da pratica com o DiffDock, com a ferramenta de cada um."""
    passos = [('1', 'Inspecionar\nos cristais', 'RCSB PDB', SITE),
              ('2', 'Abrir o\nDiffDock', 'NVIDIA', C[7]),
              ('3', 'eIF4E: enviar\nproteína e cap', 'DiffDock', C[7]),
              ('4', 'Ler poses\ne confiança', 'arquivo SDF', PROT),
              ('5', 'Comparar com\no cristal 1IPC', 'Mol*', CRISTAL),
              ('6', 'Mpro: rodar\no exemplo', 'DiffDock', C[0]),
              ('7', 'Sobrepor ao\ncristal 8HBK', 'Mol*', CRISTAL),
              ('8', 'Confiança\ncontra RMSD', 'tabela', CRISTAL),
              ('9', 'Ver as\ninterações', 'PLIP', DESTAQUE)]
    fig, ax = plt.subplots(figsize=(16, 2.6))
    frame_off(ax)
    ax.set_xlim(-0.6, len(passos) - 0.4); ax.set_ylim(-1.5, 0.6)
    for i, (num, texto, ferramenta, cor) in enumerate(passos):
        ax.add_patch(FancyBboxPatch((i - 0.42, -0.45), 0.84, 0.9, boxstyle='round,pad=0,rounding_size=0.08',
                                    fc='white', ec=cor, lw=1.6, zorder=2))
        ax.add_patch(Circle((i - 0.3, 0.33), 0.1, fc=cor, ec='none', zorder=3))
        ax.text(i - 0.3, 0.33, num, ha='center', va='center', fontsize=8, color='white', fontweight='bold', zorder=4)
        ax.text(i, -0.07, texto, ha='center', va='center', fontsize=9, zorder=4)
        ax.text(i, -0.75, ferramenta, ha='center', va='center', fontsize=8, color=cor, zorder=4)
        if i < len(passos) - 1:
            ax.add_patch(FancyArrowPatch((i + 0.46, 0), (i + 0.54, 0), arrowstyle='-|>', mutation_scale=11,
                                         color=GREY, lw=1.2, zorder=1))
    ax.annotate('', xy=(2.0, -1.15), xytext=(4.0, -1.15), arrowprops=dict(arrowstyle='-', color=C[7], lw=1.2))
    ax.text(3.0, -1.3, 'sistema A: eIF4E + cap', ha='center', va='top', fontsize=8.5, color=C[7])
    ax.annotate('', xy=(5.0, -1.15), xytext=(6.0, -1.15), arrowprops=dict(arrowstyle='-', color=C[0], lw=1.2))
    ax.text(5.5, -1.3, 'sistema B: Mpro + análogo', ha='center', va='top', fontsize=8.5, color=C[0])
    fig.suptitle('Os nove passos da prática e a ferramenta de cada um', y=1.02, fontsize=11)
    save(fig, 'pratica_fluxo.png')


# ---------------------------------------------------------------- energia de ligacao: dG e ddG
RT = 0.593        # kcal/mol a 25 C
ERRO_VINA = 2.85  # erro padrao da funcao do Vina no conjunto de treino (Trott e Olson, 2010)


def _poco(ax, x0, largura, profundidade, cor, n=200):
    """Perfil de energia livre com um poco gaussiano; devolve o ponto mais fundo."""
    x = np.linspace(x0, x0 + largura, n)
    c, w = x0 + largura / 2, largura / 7
    y = -profundidade * np.exp(-((x - c) ** 2) / (2 * w ** 2))
    ax.plot(x, y, color=INK, lw=1.4, zorder=2)
    ax.fill_between(x, y, 0, color=cor, alpha=0.18, lw=0, zorder=1)
    return c, -profundidade


def _ligante(ax, x, y, cor=POSE):
    ax.scatter([x - 0.12, x, x + 0.12], [y + 0.32, y + 0.38, y + 0.32], s=60, color=cor, edgecolor=INK,
               lw=0.5, zorder=5)


def fig_dg_poco():
    fig, ax = plt.subplots(figsize=(12, 5.2))
    casos = ((0.3, 'ligação forte', 9.0, C[0]), (7.6, 'ligação fraca', 4.5, C[7]))
    ax.axhline(0, color=GREY, lw=1, ls=(0, (4, 3)), zorder=0)
    ax.text(14.4, 0.25, 'receptor e ligante separados: nível de referência', ha='right', fontsize=9.5, color=GREY)
    for x0, nome, prof, cor in casos:
        c, fundo = _poco(ax, x0, 5.2, prof, cor)
        _ligante(ax, c, fundo)
        ax.text(c, 1.1, nome, ha='center', fontsize=11)
        ax.plot([c, c + 2.9], [fundo, fundo], color=cor, lw=0.8, ls=(0, (2, 2)))
        ax.add_patch(FancyArrowPatch((c + 2.85, 0), (c + 2.85, fundo), arrowstyle='-|>', mutation_scale=13,
                                     lw=1.3, color=cor, zorder=4))
        ax.text(c + 3.0, fundo / 2, 'ΔG mais\nnegativo' if prof > 6 else 'ΔG menos\nnegativo', va='center',
                fontsize=10, color=cor)
    ax.text(7.2, -12.2, 'Poço mais fundo: o ligante precisa de mais energia para sair, e passa mais tempo ligado.',
            ha='center', fontsize=10)
    ax.set_xlim(0, 14.5); ax.set_ylim(-12.8, 1.9)
    ax.set_ylabel('Energia livre')
    ax.tick_params(axis='x', bottom=False, top=False, labelbottom=False)
    ax.tick_params(axis='y', left=False, right=False, labelleft=False)
    ax.set_title('ΔG de ligação: a profundidade do poço (esquema)')
    save(fig, 'fund_dg_poco.png')


def fig_ddg():
    fig, a = plt.subplots(figsize=(10, 5.0))
    a.axhline(0, color=GREY, lw=1, ls=(0, (4, 3)), zorder=0)
    a.text(10.9, 0.15, 'separados', ha='right', va='bottom', fontsize=9.5, color=GREY)
    c1, f1 = _poco(a, 0.3, 4.6, 9.0, CRISTAL)
    c2, f2 = _poco(a, 5.6, 4.6, 5.5, ERRO)
    _ligante(a, c1, f1); _ligante(a, c2, f2)
    a.text(c1, 0.7, 'proteína selvagem', ha='center', fontsize=10.5)
    a.text(c2, 0.7, 'proteína mutante', ha='center', fontsize=10.5)
    a.plot([c1, c2 + 1.3], [f1, f1], color=CRISTAL, lw=0.9, ls=(0, (2, 2)))
    a.add_patch(FancyArrowPatch((c2 + 1.2, f1), (c2 + 1.2, f2), arrowstyle='-|>', mutation_scale=13, lw=1.3,
                                color=ERRO, zorder=4))
    a.text(c2 + 1.4, (f1 + f2) / 2, 'ΔΔG > 0', va='center', fontsize=10.5, color=ERRO)
    a.text(5.5, -10.6, 'ΔΔG = ΔG(mutante) − ΔG(selvagem)', ha='center', fontsize=10.5)
    a.text(5.5, -11.5, 'positivo: a mutação deixa o poço mais raso e enfraquece a ligação', ha='center', fontsize=9.5)
    a.set_xlim(0, 11); a.set_ylim(-12, 1.4)
    a.set_ylabel('Energia livre')
    a.tick_params(axis='x', bottom=False, top=False, labelbottom=False)
    a.tick_params(axis='y', left=False, right=False, labelleft=False)
    a.set_title('ΔΔG: o mesmo ligante em duas proteínas (esquema)')
    save(fig, 'fund_ddg.png')


# ---------------------------------------------------------------- conceitos basicos: receptor, ligante, pose
def _heavy_lig(path):
    m = Chem.MolFromMolFile(str(path), sanitize=False, removeHs=False)
    keep = [a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum() > 1]
    return m.GetConformer().GetPositions()[keep] * A


def _separado(ax, code, lig_path, lig_col, gap):
    """Receptor inteiro e o ligante afastado do sitio, na direcao que sai do bolsao."""
    t = md.load_pdb(str(fetch_pdb(code)))
    prot = heavy(t, 'chainid 0 and protein')
    X = t.xyz[0][prot]
    L = _heavy_lig(lig_path)
    pc, lc = X.mean(0), L.mean(0)
    out = (lc - pc) / np.linalg.norm(lc - pc)
    look = np.cross(out, [0, 0, 1.0])
    if np.linalg.norm(look) < 0.1:
        look = np.cross(out, [1.0, 0, 0])
    B = view_basis(look, up_hint=out)
    L2 = L + out * gap
    P = project(np.r_[X, L2], B, pc)
    draw_spheres(ax, P[:len(X)], [PROT] * len(X), [R_ATOM] * len(X), pad=0.6)
    draw_spheres(ax, P[len(X):], [lig_col] * len(L2), [0.17] * len(L2), pad=0.6, z=3)
    sitio = project(L, B, pc)
    ax.add_patch(FancyArrowPatch((P[len(X):, 0].mean(), P[len(X):, 1].min() - 0.15),
                                 (sitio[:, 0].mean(), sitio[:, 1].max() + 0.25), arrowstyle='-|>',
                                 mutation_scale=16, color=INK, lw=1.4, zorder=6))
    ax.set_xlim(P[:, 0].min() - 0.6, P[:, 0].max() + 0.6); ax.set_ylim(P[:, 1].min() - 0.4, P[:, 1].max() + 0.5)
    return P, len(X), len(L)


def fig_receptor_ligante():
    fig, axs = plt.subplots(1, 2, figsize=(12, 7.2))
    casos = [(axs[0], '1IPC', ROOT / 'dados' / 'eif4e' / 'ligante_cristal.sdf', CRISTAL, 'eIF4E', 'cap do mRNA (m7GTP)',
              'Ligante natural: a eIF4E reconhece o cap do mRNA\n(PDB 1IPC, célula saudável)'),
             (axs[1], '8HBK', ROOT / 'dados' / 'mpro' / 'ligante_cristal.sdf', POSE, 'Mpro do SARS-CoV-2', 'ensitrelvir',
              'Fármaco: o ensitrelvir bloqueia a protease do vírus\n(PDB 8HBK, infecção)')]
    for ax, code, lig, col, rec, nome, titulo in casos:
        P, nx, nl = _separado(ax, code, lig, col, gap=2.4)
        pl = P[nx:, :2]
        label(ax, (pl[:, 0].max(), pl[:, 1].mean()), f'ligante: {nome}\n{nl} átomos pesados',
              (pl[:, 0].max() + 1.6, pl[:, 1].mean()), col)
        low = [k for k in range(nx) if P[k, 1] < P[:nx, 1].min() + 0.6]
        label(ax, front(P, low), f'receptor: {rec}\n{nx} átomos pesados', (P[:, 0].max() + 0.4, P[:nx, 1].min() + 0.6), GREY)
        ax.set_xlim(P[:, 0].min() - 0.4, max(P[:, 0].max(), pl[:, 0].max() + 1.6) + 1.4)
        ax.set_title(titulo, fontsize=10.5)
    fig.suptitle('Receptor: a macromolécula que tem o sítio. Ligante: a molécula que se encaixa nele.', y=1.0, fontsize=11.5)
    fig.tight_layout(w_pad=0.5)
    save(fig, 'fund_receptor_ligante.png')


def _mol_2d(ax, Q, bonds, cor, alpha=1.0, z=3, s=46):
    for i, j in bonds:
        ax.plot(Q[[i, j], 0], Q[[i, j], 1], color=shade(INK, 1 - alpha), lw=1.6, zorder=z, solid_capstyle='round')
    ax.scatter(Q[:, 0], Q[:, 1], s=s, color=shade(cor, 1 - alpha), edgecolor=shade(INK, 1 - alpha), lw=0.5, zorder=z + 1)


def fig_pose():
    """Pose = posicao + orientacao + conformacao, com o ensitrelvir do cristal."""
    m = Chem.MolFromMolFile(str(ROOT / 'dados' / 'mpro' / 'ligante_cristal.sdf'))
    m = Chem.RemoveHs(m)
    from rdkit.Chem import rdMolTransforms as T
    X0 = m.GetConformer().GetPositions()
    u, s_, vt = np.linalg.svd(X0 - X0.mean(0)); B = vt.T
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in m.GetBonds()]
    # ligacao rotavel cujo lado movel tem ~10 atomos
    patt = Chem.MolFromSmarts('[!$([NH])&!D1&!$(*#*)]-&!@[!D1&!$(*#*)]')
    melhor = None
    for i, j in m.GetSubstructMatches(patt):
        em = Chem.RWMol(m); em.RemoveBond(i, j)
        frags = Chem.GetMolFrags(em.GetMol(), asMols=False)
        lado = [f for f in frags if j in f][0]
        if melhor is None or abs(len(lado) - 9) < abs(melhor[2] - 9):
            melhor = (i, j, len(lado))
    i, j, _ = melhor
    a_ = [n.GetIdx() for n in m.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx() != j][0]
    d_ = [n.GetIdx() for n in m.GetAtomWithIdx(j).GetNeighbors() if n.GetIdx() != i][0]
    mt = Chem.Mol(m); conf = mt.GetConformer()
    T.SetDihedralDeg(conf, a_, i, j, d_, T.GetDihedralDeg(conf, a_, i, j, d_) + 140)
    Xt = conf.GetPositions()
    c = X0.mean(0)
    th = np.radians(60); R = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
    proj = lambda X: (X - c) @ B
    ref = proj(X0)
    casos = [('Pose de referência (cristal)', ref, CRISTAL),
             ('Outra posição: translação\n(3 números: x, y, z)', ref + np.array([4.0, 2.0, 0]), POSE),
             ('Outra orientação: rotação\n(3 números: três ângulos)', ref @ R.T, POSE),
             ('Outra conformação: torção\n(1 número por ligação rotável)', proj(Xt), POSE)]
    fig, axs = plt.subplots(1, 4, figsize=(17, 4.4))
    for k, (ax, (titulo, Q, cor)) in enumerate(zip(axs, casos)):
        if k:
            _mol_2d(ax, ref, bonds, CRISTAL, alpha=0.3, z=1, s=30)
        _mol_2d(ax, Q, bonds, cor)
        if k == 3:
            ax.plot(*Q[[i, j], :2].T, color=DESTAQUE, lw=4, zorder=5, solid_capstyle='round')
        ax.set_aspect('equal'); frame_off(ax)
        ax.set_xlim(-10, 13); ax.set_ylim(-7, 8.5)
        ax.set_title(titulo, fontsize=10.5)
    legend(fig, [('ensitrelvir no cristal (8HBK)', CRISTAL), ('nova pose', POSE), ('ligação girada', DESTAQUE)], y=-0.04)
    fig.suptitle('Pose = posição + orientação + conformação do ligante dentro do receptor', y=1.02, fontsize=11.5)
    save(fig, 'fund_pose.png')


TORS = {'ensitrelvir': 6, 'm7GTP': 11}   # TORSDOF do Vina em resultados/*/vina_dock.pdbqt


def fig_espaco_conformacional():
    fig, (a, b) = plt.subplots(1, 2, figsize=(13.5, 4.6), gridspec_kw={'width_ratios': [1.15, 1]})
    nomes = list(TORS)
    y = np.arange(len(nomes))[::-1]
    for yy, n in zip(y, nomes):
        partes = [(3, C[0], 'translação'), (3, C[6], 'rotação'), (TORS[n], C[5], 'torções')]
        x0 = 0
        for v, cor, lab in partes:
            a.barh(yy, v, left=x0, color=cor, edgecolor=INK, lw=0.5, height=0.5, label=lab if yy == y[0] else None)
            a.text(x0 + v / 2, yy, str(v), ha='center', va='center', fontsize=10, color='white')
            x0 += v
        a.text(x0 + 0.3, yy, f'{x0} números por pose', va='center', fontsize=10)
    a.set_yticks(y); a.set_yticklabels(nomes, fontsize=10.5); a.tick_params(axis='y', right=False)
    a.set_xlim(0, 25); a.set_xlabel('Graus de liberdade do ligante')
    a.legend(loc='upper right', ncol=3, fontsize=9, bbox_to_anchor=(1.0, 1.0))
    a.set_ylim(-0.6, 1.9)
    a.set_title('Cada pose é descrita por 6 + n números')
    comb = [3 ** TORS[n] for n in nomes]
    b.bar(range(len(nomes)), comb, color=[C[0], C[7]], edgecolor=INK, lw=0.5, width=0.5)
    for k, v in enumerate(comb):
        b.text(k, v * 1.4, f'$3^{{{TORS[nomes[k]]}}}$ = {v:,}'.replace(',', ' '), ha='center', va='bottom', fontsize=10)
    b.set_yscale('log'); b.set_ylim(1, 1e7)
    b.set_xticks(range(len(nomes))); b.set_xticklabels(nomes, fontsize=10.5); b.tick_params(axis='x', top=False)
    b.set_ylabel('Conformações possíveis')
    b.set_title('Testando só 3 ângulos por ligação rotável (contagem ilustrativa)')
    fig.suptitle('Espaço conformacional: cada ligação rotável multiplica o que o programa precisa percorrer', y=1.02,
                 fontsize=11.5)
    fig.tight_layout(w_pad=3)
    save(fig, 'fund_espaco_conformacional.png')


def fig_busca_pontuacao():
    x = np.linspace(0, 10, 600)
    real = -6 * np.exp(-((x - 6.6) ** 2) / 0.35) - 3.2 * np.exp(-((x - 2.4) ** 2) / 0.6) + 0.15 * np.sin(3 * x)
    func = -4.4 * np.exp(-((x - 6.6) ** 2) / 0.35) - 5.0 * np.exp(-((x - 2.4) ** 2) / 0.6) + 0.15 * np.sin(3 * x)
    fig, (a, b) = plt.subplots(1, 2, figsize=(14, 4.8), sharey=True)
    rng = np.random.default_rng(4)
    vis = np.sort(np.r_[rng.uniform(0.3, 4.5, 9), rng.uniform(8.2, 9.8, 3)])
    a.plot(x, func, color=INK, lw=1.6)
    a.scatter(vis, np.interp(vis, x, func), s=60, color=POSE, edgecolor=INK, lw=0.5, zorder=4)
    k = np.argmin(np.interp(vis, x, func))
    a.annotate('melhor pose encontrada', (vis[k], np.interp(vis[k], x, func)), xytext=(3.4, -6.6), fontsize=9.5,
               ha='left', arrowprops=dict(arrowstyle='-', color=INK, lw=0.8))
    a.annotate('poço que a busca\nnunca visitou', (6.75, np.interp(6.75, x, func)), xytext=(7.6, -2.0), fontsize=9.5, ha='left',
               arrowprops=dict(arrowstyle='-', color=ERRO, lw=0.8), color=ERRO)
    a.set_title('Busca: o programa só conhece as poses que visitou')
    b.plot(x, real, color=CRISTAL, lw=1.6, ls=(0, (5, 3)), label='energia real (o que a natureza faz)')
    b.plot(x, func, color=INK, lw=1.6, label='função de pontuação (o que o programa calcula)')
    b.annotate('a função prefere\neste poço', (2.4, func[np.argmin(np.abs(x - 2.4))]), xytext=(0.2, -7.6), fontsize=9.5,
               arrowprops=dict(arrowstyle='-', color=INK, lw=0.8))
    b.annotate('a pose certa\nestá aqui', (6.6, real.min()), xytext=(7.6, -6.4), fontsize=9.5, color=CRISTAL,
               arrowprops=dict(arrowstyle='-', color=CRISTAL, lw=0.8))
    b.legend(loc='upper right', fontsize=9)
    b.set_title('Pontuação: a nota da função não é a energia real')
    for ax in (a, b):
        ax.set_xlim(0, 10); ax.set_ylim(-8.4, 2.4)
        ax.set_xlabel('Poses possíveis (espaço conformacional)')
        ax.tick_params(axis='x', bottom=False, top=False, labelbottom=False)
    a.set_ylabel('Pontuação (mais baixo = melhor)')
    fig.suptitle('Todo programa de docking tem duas partes: um algoritmo de busca e uma função de pontuação (esquema)',
                 y=1.02, fontsize=11.5)
    fig.tight_layout(w_pad=2)
    save(fig, 'fund_busca_pontuacao.png')


# ---------------------------------------------------------------- tipos de docking: graus de liberdade reais (Mpro)
def _graus_mpro():
    t = md.load_pdb(str(fetch_pdb('8HBK')))
    a = t.atom_slice(t.topology.select('chainid 0 and protein'))
    sitio = {41, 49, 142, 145, 163, 165, 166, 189}
    chis = [getattr(md, f'compute_chi{i}')(a)[0] for i in range(1, 5)]
    chi_sitio = sum(1 for idx in chis for q in idx if a.topology.atom(q[1]).residue.resSeq in sitio)
    bb = len(md.compute_phi(a)[0]) + len(md.compute_psi(a)[0])
    n_lig = 6 + TORS['ensitrelvir']
    return [('Rígido-rígido', 6), ('Ligante flexível,\nreceptor rígido', n_lig),
            ('+ cadeias laterais\ndo sítio flexíveis', n_lig + chi_sitio),
            ('Flexível-flexível\n(proteína inteira)', n_lig + bb + sum(len(c) for c in chis))], chi_sitio, len(sitio)


def fig_graus_tipos():
    tipos, chi_sitio, n_res = _graus_mpro()
    ganho = ['só a forma: o ligante e o bolsão já se encaixam',
             'a molécula dobra e gira suas ligações',
             f'{n_res} resíduos do bolsão ajustam seus {chi_sitio} ângulos χ',
             'alças e domínios se movem: encaixe induzido real']
    custo = ['segundos', 'minutos', 'dezenas de minutos a horas', 'horas a dias (dinâmica)']
    fig, ax = plt.subplots(figsize=(13, 4.8))
    y = np.arange(len(tipos))[::-1]
    cores = [GREY, C[0], C[6], C[7]]
    ax.barh(y, [v for _, v in tipos], color=cores, edgecolor=INK, lw=0.5, height=0.55)
    for yy, (_, v), g, c_ in zip(y, tipos, ganho, custo):
        ax.text(v * 1.15, yy + 0.12, f'{v:,}'.replace(',', ' ') + ' graus de liberdade', va='center', fontsize=10)
        ax.text(v * 1.15, yy - 0.2, f'{g}; custo: {c_}', va='center', fontsize=9, color=GREY)
    ax.set_xscale('log'); ax.set_xlim(3, 4e5)
    ax.set_yticks(y); ax.set_yticklabels([n for n, _ in tipos], fontsize=10); ax.tick_params(axis='y', right=False)
    ax.set_xlabel('Números que o programa ajusta por pose (escala logarítmica)')
    ax.set_title('Mpro + ensitrelvir: quanto mais flexível o modelo, mais fiel e mais caro')
    save(fig, 'fund_graus_tipos.png')


# ---------------------------------------------------------------- cada programa tem a sua escala
# ---------------------------------------------------------------- pratica: contexto da eIF4E e o resultado do AC
def fig_cap_mrna():
    """Esquema: mRNA com cap 5', eIF4E no cap, eIF4F recrutando o ribossomo; 4E-BP como freio."""
    fig, (a, b) = plt.subplots(1, 2, figsize=(15, 4.0), gridspec_kw={'width_ratios': [1, 1]})
    for ax in (a, b):
        frame_off(ax); ax.set_xlim(0, 10); ax.set_ylim(0.3, 5.0); ax.set_aspect('equal')
    # mRNA
    a.plot([1.2, 9.6], [1.6, 1.6], color=INK, lw=3, solid_capstyle='round', zorder=2)
    a.add_patch(Rectangle((4.2, 1.35), 3.2, 0.5, fc=shade(C[0], 0.4), ec=INK, lw=0.8, zorder=3))
    a.text(5.8, 1.6, 'região codificadora', ha='center', va='center', fontsize=9, zorder=4)
    a.text(9.6, 1.1, 'AAAA…', ha='right', fontsize=9.5)
    a.add_patch(Circle((1.0, 1.6), 0.28, fc=CRISTAL, ec=INK, lw=0.8, zorder=4))
    a.text(1.0, 0.95, "cap 5'\n(m7G-ppp)", ha='center', va='top', fontsize=9)
    # eIF4F
    a.add_patch(Circle((1.0, 2.35), 0.45, fc=DESTAQUE, ec=INK, lw=0.8, zorder=4))
    a.text(1.0, 2.35, 'eIF4E', ha='center', va='center', fontsize=8.5, color='white', zorder=5)
    a.add_patch(FancyBboxPatch((1.5, 2.05), 1.9, 0.6, boxstyle='round,pad=0.05', fc=shade(C[1], 0.3), ec=INK, lw=0.8, zorder=3))
    a.text(2.45, 2.35, 'eIF4G', ha='center', va='center', fontsize=9, zorder=4)
    a.add_patch(Circle((2.2, 3.0), 0.35, fc=shade(C[6], 0.2), ec=INK, lw=0.8, zorder=4))
    a.text(2.2, 3.0, 'eIF4A', ha='center', va='center', fontsize=8, zorder=5)
    a.add_patch(FancyBboxPatch((3.2, 2.7), 1.6, 1.3, boxstyle='round,pad=0.1', fc=shade(GREY, 0.4), ec=INK, lw=0.8, zorder=3))
    a.text(4.0, 3.35, 'ribossomo\n(40S)', ha='center', va='center', fontsize=9, zorder=4)
    a.add_patch(FancyArrowPatch((4.9, 2.2), (6.5, 2.2), arrowstyle='-|>', mutation_scale=14, color=INK, lw=1.2))
    a.text(5.7, 2.45, 'varre até o AUG', ha='center', fontsize=9)
    a.text(0.4, 4.7, 'Célula saudável: a eIF4E segura o cap e o complexo eIF4F traz o ribossomo', fontsize=10)
    # freio e cancer
    b.add_patch(Circle((2.0, 3.2), 0.5, fc=DESTAQUE, ec=INK, lw=0.8))
    b.text(2.0, 3.2, 'eIF4E', ha='center', va='center', fontsize=8.5, color='white')
    b.add_patch(FancyBboxPatch((2.55, 2.95), 1.2, 0.5, boxstyle='round,pad=0.05', fc=shade(ERRO, 0.4), ec=INK, lw=0.8))
    b.text(3.15, 3.2, '4E-BP', ha='center', va='center', fontsize=9)
    b.text(2.6, 2.2, 'freio: a 4E-BP ocupa o lugar\nda eIF4G; a mTOR solta o freio', ha='center', fontsize=9)
    for k, dx in enumerate((6.0, 6.8, 7.6, 8.4, 6.4, 7.2, 8.0)):
        yy = 3.4 if k < 4 else 2.8
        b.add_patch(Circle((dx, yy), 0.32, fc=DESTAQUE, ec=INK, lw=0.6))
    b.text(7.2, 2.0, 'câncer: eIF4E em excesso traduz\nmais mRNAs de proliferação', ha='center', fontsize=9)
    b.text(0.3, 4.7, 'Regulação (freio) e doença (excesso)', fontsize=10)
    fig.suptitle("O cap 5' do mRNA e a eIF4E: o primeiro passo da tradução (esquema)", y=1.0, fontsize=11.5)
    fig.tight_layout(w_pad=2)
    save(fig, 'pratica_cap_mrna.png')


def fig_eif4e_sitio():
    fig, ax = plt.subplots(figsize=(10, 6.6))
    P, hp = _scene_with_ligand(ax, '1IPC', 'MGP', 'A', '56 102', DESTAQUE, cut=0.9,
                               title='O cap no sítio da eIF4E (PDB 1IPC): a guanina metilada entre Trp56 e Trp102')
    label(ax, hp, 'Trp56 e Trp102:\no anel do cap fica\nentre os dois (cátion-π)', (P[:, 0].max() + 1.3, hp[1] + 0.4),
          DESTAQUE)
    ax.set_xlim(P[:, 0].min() - 0.3, P[:, 0].max() + 2.4)
    fit(fig, ax, largura=11, extra=0.5)
    legend(fig, [('cap (m7GTP)', CRISTAL), ('Trp56 e Trp102', DESTAQUE), ('eIF4E: átomos a menos de 9 Å do cap', PROT)],
           y=-0.08, ax=ax)
    save(fig, 'pratica_eif4e_sitio.png')


def _dd_rows(nome):
    d = json.loads((ROOT / 'resultados' / nome / 'comparacao_cristal.json').read_text())
    return d['poses'] if isinstance(d, dict) else d


def _sdf_xyz(path, k):
    m = Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False)[k - 1]
    return skeleton(m).GetConformer().GetPositions() * A


def _cena_sitio(ax, t, L, Q, near, X, B, c0, destaque, titulo):
    P = project(np.r_[X[near], L, Q], B, c0)
    cols = [DESTAQUE if i in destaque else PROT for i in near] + [CRISTAL] * len(L) + [POSE] * len(Q)
    rad = [R_ATOM] * len(near) + [0.16] * (len(L) + len(Q))
    draw_spheres(ax, P[:len(near)], cols[:len(near)], rad[:len(near)], pad=0.5)
    draw_spheres(ax, P[len(near):], cols[len(near):], rad[len(near):], pad=0.5, z=3)
    ax.set_xlim(P[:, 0].min() - 0.3, P[:, 0].max() + 0.3); ax.set_ylim(P[:, 1].min() - 0.3, P[:, 1].max() + 0.3)
    ax.set_title(titulo.replace('.', ',').replace('-', '−'), fontsize=10)


def fig_diffdock_eif4e():
    """eIF4E: cristal, pose 1 e pose 20 do DiffDock."""
    t = md.load_pdb(str(fetch_pdb('1IPC')))
    L = _heavy_lig(ROOT / 'dados' / 'eif4e' / 'ligante_cristal.sdf')
    near, X = pocket_scene(t, L, cut=0.9)
    trp = set(t.topology.select('resSeq 56 102 and not element H'))
    rows = _dd_rows('eif4e_diffdock')
    B = view_basis(X[heavy(t, 'protein')].mean(0) - L.mean(0)); c0 = L.mean(0)
    sdf = ROOT / 'resultados' / 'eif4e_diffdock' / 'poses.sdf'
    fig, axs = plt.subplots(1, 3, figsize=(16, 5.8))
    for ax, k in ((axs[0], 0), (axs[1], 1), (axs[2], 20)):
        Q = _sdf_xyz(sdf, k) if k else np.zeros((0, 3))
        tit = ('Cristal (1IPC): a guanina\nentre Trp56 e Trp102' if not k else
               f"Pose {k}: confiança {rows[k - 1]['confianca']:.2f}\nRMSD {rows[k - 1]['rmsd_A']:.1f} Å")
        _cena_sitio(ax, t, L, Q, near, X, B, c0, trp, tit)
    legend(fig, [('cap no cristal', CRISTAL), ('pose do DiffDock', POSE), ('Trp56 e Trp102', DESTAQUE),
                 ('eIF4E: átomos a menos de 9 Å do cap', PROT)])
    fig.suptitle('eIF4E + cap: a pose de maior e a de menor confiança do DiffDock contra o cristal', y=1.0, fontsize=11)
    save(fig, 'pratica_diffdock_eif4e.png')


def fig_diffdock_mpro():
    """Mpro: ensitrelvir no cristal 8HBK e poses do analogo (receptor OpenFold sobreposto a 8HBK)."""
    t = md.load_pdb(str(fetch_pdb('8HBK')))
    L = _heavy_lig(ROOT / 'dados' / 'mpro' / 'ligante_cristal.sdf')
    near, X = pocket_scene(t, L, cut=0.9)
    rows = _dd_rows('mpro_diffdock')
    B = view_basis(X[heavy(t, 'protein')].mean(0) - L.mean(0)); c0 = L.mean(0)
    sdf = ROOT / 'resultados' / 'mpro_diffdock' / 'poses_no_8hbk.sdf'
    errada = max(rows[:10], key=lambda r: r['rmsd_nucleo_A'])['pose']
    fig, axs = plt.subplots(1, 3, figsize=(16, 5.8))
    for ax, k in ((axs[0], 0), (axs[1], 1), (axs[2], errada)):
        Q = _sdf_xyz(sdf, k) if k else np.zeros((0, 3))
        tit = ('Cristal (8HBK): o ensitrelvir\nno sítio ativo' if not k else
               f"Análogo, pose {k}: confiança {rows[k - 1]['confianca']:.2f}\nRMSD do núcleo {rows[k - 1]['rmsd_nucleo_A']:.1f} Å")
        _cena_sitio(ax, t, L, Q, near, X, B, c0, set(), tit)
    legend(fig, [('ensitrelvir no cristal', CRISTAL), ('análogo, pose do DiffDock', POSE),
                 ('Mpro: átomos a menos de 9 Å do ligante', PROT)])
    fig.suptitle('Mpro: o análogo docado na proteína predita (OpenFold), sobreposta ao cristal 8HBK', y=1.0, fontsize=11)
    save(fig, 'pratica_diffdock_mpro.png')


def fig_analogo():
    """Ensitrelvir e o analogo do exemplo da NVIDIA, com o nucleo comum destacado."""
    from rdkit.Chem import rdFMCS, AllChem
    ens = Chem.MolFromSmiles((ROOT / 'dados' / 'mpro' / 'ligante.smi').read_text().split()[0])
    ana = Chem.RemoveHs(Chem.MolFromMolFile(str(ROOT / 'dados' / 'mpro_analogo' / 'ligante_analogo_3d.sdf')))
    ana = Chem.MolFromSmiles(Chem.MolToSmiles(ana))
    mcs = Chem.MolFromSmarts(rdFMCS.FindMCS([ens, ana], ringMatchesRingOnly=True, completeRingsOnly=True, timeout=20).smartsString)
    AllChem.Compute2DCoords(mcs)
    fig, axs = plt.subplots(1, 2, figsize=(14, 5.4))
    for ax, m, nome in ((axs[0], ens, 'Ensitrelvir (no cristal 8HBK)'), (axs[1], ana, 'Análogo do exemplo da NVIDIA')):
        AllChem.GenerateDepictionMatching2DStructure(m, mcs)
        hit = m.GetSubstructMatch(mcs)
        diff = [a.GetIdx() for a in m.GetAtoms() if a.GetIdx() not in hit]
        img = Draw.MolToImage(m, size=(640, 480), highlightAtoms=diff,
                              highlightColor=tuple(v for v in plt.matplotlib.colors.to_rgb(shade(ERRO, 0.4))))
        ax.imshow(img); frame_off(ax)
        ax.set_title(f'{nome}: {m.GetNumHeavyAtoms()} átomos pesados', fontsize=10.5)
    fig.text(0.5, 0.02, f'Núcleo comum: {mcs.GetNumAtoms()} átomos. Em vermelho, o que muda: o metiltriazol vira um grupo '
             'hidrazona com bromo, e o cloro do indazol vira flúor.', ha='center', fontsize=10)
    fig.suptitle('O ligante do exemplo é um análogo: o RMSD contra o cristal só vale no núcleo comum', y=1.0, fontsize=11.5)
    save(fig, 'pratica_analogo.png')


def fig_diffdock_confianca():
    """Confianca do DiffDock contra o RMSD, nos dois sistemas."""
    fig, axs = plt.subplots(1, 2, figsize=(14, 4.9), sharey=True)
    for ax, nome, chave, titulo, cor in ((axs[0], 'eif4e_diffdock', 'rmsd_A', 'eIF4E + cap (receptor do cristal 1IPC)', C[7]),
                                         (axs[1], 'mpro_diffdock', 'rmsd_nucleo_A',
                                          'Mpro + análogo (receptor predito, OpenFold)', C[0])):
        rows = _dd_rows(nome)
        c = np.array([r['confianca'] for r in rows]); r = np.array([r[chave] for r in rows])
        ax.axvspan(0, 2, color=BAND, lw=0, zorder=0)
        ax.axhline(-1.5, color=GREY, ls=(0, (4, 3)), lw=1)
        ax.scatter(r, c, s=46, color=cor, edgecolor=INK, lw=0.4, zorder=3)
        ax.scatter(r[0], c[0], s=150, facecolor='none', edgecolor=INK, lw=1.4, zorder=4)
        ax.annotate('pose 1', (r[0], c[0]), xytext=(12, 4), textcoords='offset points', fontsize=9)
        ax.text(13.8, -1.45, 'acima de −1,5: confiança moderada', ha='right', va='bottom', fontsize=8.5, color=GREY)
        ax.set_xlim(0, 14); ax.set_xlabel('RMSD do cap (Å)' if chave == 'rmsd_A' else 'RMSD do núcleo comum (Å)')
        ax.set_title(f'{titulo}: {(r < 2).sum()} de {len(r)} abaixo de 2 Å', fontsize=10)
    axs[0].set_ylabel('Confiança do DiffDock (maior = melhor)')
    axs[0].text(1.0, 0.03, 'RMSD < 2 Å', transform=axs[0].get_xaxis_transform(), ha='center', fontsize=8.5)
    fig.suptitle('DiffDock contra o cristal: a confiança acompanha o acerto?', y=1.02, fontsize=11.5)
    fig.tight_layout(w_pad=2)
    save(fig, 'pratica_diffdock_confianca.png')


def fig_interacoes_diffdock():
    """Contatos do cristal e da pose 1 do DiffDock, nos dois sistemas."""
    tipo_cor = {'ligação de hidrogênio': C[0], 'ponte salina': ERRO, 'hidrofóbica': C[1], 'empilhamento π': C[4],
                'cátion-π': C[6], 'halogênio': C[3]}
    fig, axs = plt.subplots(1, 2, figsize=(15, 3.6), gridspec_kw={'width_ratios': [1, 1.3]})
    for ax, nome, kc, titulo in ((axs[0], 'eif4e_diffdock', 'cristal', 'eIF4E + cap'),
                                 (axs[1], 'mpro_diffdock', 'cristal (ensitrelvir, 8HBK)', 'Mpro: ensitrelvir e análogo')):
        it = json.loads((ROOT / 'resultados' / nome / 'interacoes.json').read_text())
        keys = [kc, 'pose 1']
        res = sorted({r['residuo'] for k in keys for r in it[k]}, key=lambda s: int(''.join(ch for ch in s if ch.isdigit())))
        for j, k in enumerate(keys):
            for i, rn in enumerate(res):
                tipos = sorted({r['tipo'] for r in it[k] if r['residuo'] == rn})
                ax.add_patch(Rectangle((i - 0.45, j - 0.45), 0.9, 0.9, fc='white', ec=BAND, lw=0.8))
                for ti, tp in enumerate(tipos):
                    w = 0.9 / len(tipos)
                    ax.add_patch(Rectangle((i - 0.45 + ti * w, j - 0.45), w, 0.9, fc=tipo_cor.get(tp, GREY), ec='white', lw=0.5))
        ax.set_xlim(-0.6, len(res) - 0.4); ax.set_ylim(1.6, -0.6)
        ax.set_xticks(range(len(res))); ax.set_xticklabels(res, rotation=60, ha='right', fontsize=9)
        ax.set_yticks([0, 1]); ax.set_yticklabels(['cristal', 'DiffDock, pose 1'], fontsize=9.5)
        ax.tick_params(top=False, right=False, length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_title(titulo, fontsize=10.5)
    hs = [Rectangle((0, 0), 1, 1, fc=c_) for c_ in tipo_cor.values()]
    fig.legend(hs, list(tipo_cor), loc='lower center', ncol=6, frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.2))
    fig.suptitle('Contatos (PLIP) do cristal e da pose 1 do DiffDock, uma coluna por resíduo', y=1.04, fontsize=11)
    save(fig, 'pratica_interacoes_diffdock.png')


def main():
    fig_binding_models(); fig_vocabulary(); fig_forces(); fig_flexibility(); fig_search(); fig_levels_cost()
    fig_dg_poco(); fig_ddg()
    fig_smiles(); fig_resistance(); fig_p53_mdm2(); fig_pd1_pdl1(); fig_md_docking()
    fig_cavidades(); fig_tipos_sitio(); fig_fluxo()
    fig_rmsd(); fig_metrics()
    fig_receptor_ligante(); fig_pose(); fig_espaco_conformacional(); fig_busca_pontuacao(); fig_graus_tipos()
    fig_cap_mrna(); fig_eif4e_sitio()
    fig_diffdock_eif4e(); fig_diffdock_mpro(); fig_analogo(); fig_diffdock_confianca(); fig_interacoes_diffdock()

if __name__ == '__main__':
    if len(sys.argv) > 1:
        for nome in sys.argv[1:]:
            globals()[nome]()
    else:
        main()
