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


def fig_poses_scatter():
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    for nome, col, mk, lab in (('mpro', C[0], 'o', 'Mpro + ensitrelvir'), ('eif4e', C[7], 's', 'eIF4E + m7GTP')):
        rows = _load_rows(nome)
        e = [r['energia_kcal_mol'] for r in rows]; r = [r['rmsd_A'] for r in rows]
        ax.scatter(r, e, s=46, color=col, marker=mk, edgecolor=INK, lw=0.4, label=lab, zorder=3)
        ax.annotate('pose 1', (r[0], e[0]), xytext=(10, 4), textcoords='offset points', fontsize=8.5, color=col)
    ax.axvspan(0, 2, color=BAND, lw=0, zorder=0)
    ax.text(1.0, ax.get_ylim()[1] - 0.15, 'RMSD < 2 Å:\npose correta', ha='center', va='top', fontsize=8.5)
    ax.set_xlim(0, None)
    ax.set_xlabel('RMSD em relação à pose do cristal (Å)')
    ax.set_ylabel('Energia do Vina (kcal/mol)')
    ax.legend(loc='lower right')
    ax.set_title('Todas as poses da prática: a melhor energia nem sempre é a pose certa')
    save(fig, 'pratica_energia_rmsd.png')


def _pose_xyz(nome, k):
    e, m = poses(resultados(nome) / 'vina_dock.pdbqt')[k - 1]
    m = skeleton(m)
    return m.GetConformer().GetPositions() * A, e


def fig_mpro_result():
    t = md.load_pdb(str(fetch_pdb('8HBK')))
    L = Chem.MolFromMolFile(str(ROOT / 'dados' / 'mpro' / 'ligante_cristal.sdf'), sanitize=False).GetConformer().GetPositions() * A
    near, X = pocket_scene(t, L, cut=0.9)
    rows = _load_rows('mpro')
    wrong = max(rows[1:6], key=lambda r: r['rmsd_A'])['pose']
    fig, axs = plt.subplots(1, 2, figsize=(13, 5.8))
    look = X[heavy(t, 'protein')].mean(0) - L.mean(0)
    B = view_basis(look); c0 = L.mean(0)
    for ax, k, title in ((axs[0], 1, f"Pose 1: {rows[0]['energia_kcal_mol']:.1f} kcal/mol, RMSD {rows[0]['rmsd_A']:.1f} Å"),
                         (axs[1], wrong, f"Pose {wrong}: {rows[wrong - 1]['energia_kcal_mol']:.1f} kcal/mol, "
                                         f"RMSD {rows[wrong - 1]['rmsd_A']:.1f} Å")):
        Q, _ = _pose_xyz('mpro', k)
        P = project(np.r_[X[near], L, Q], B, c0)
        cols = [PROT] * len(near) + [CRISTAL] * len(L) + [POSE] * len(Q)
        rad = [R_ATOM] * len(near) + [0.16] * len(L) + [0.16] * len(Q)
        # ligantes por cima do bolsao, para nao sumirem atras dos atomos do receptor
        draw_spheres(ax, P[:len(near)], cols[:len(near)], rad[:len(near)], pad=0.5)
        draw_spheres(ax, P[len(near):], cols[len(near):], rad[len(near):], pad=0.5, z=3)
        ax.set_xlim(P[:, 0].min() - 0.3, P[:, 0].max() + 0.3); ax.set_ylim(P[:, 1].min() - 0.3, P[:, 1].max() + 0.3)
        ax.set_title(title, fontsize=10)
    legend(fig, [('ensitrelvir no cristal (8HBK)', CRISTAL), ('pose do docking (SwissDock, Vina)', POSE),
                 ('Mpro: átomos a menos de 9 Å do ligante', PROT)])
    fig.suptitle('Prática A: o redocking do ensitrelvir na Mpro', y=1.0, fontsize=11)
    save(fig, 'pratica_mpro_poses.png')


def _trp_rings(t, resids):
    out = {}
    for r in resids:
        idx = t.topology.select(f'resSeq {r} and name CD2 CE2 CE3 CZ2 CZ3 CH2')
        out[r] = t.xyz[0][idx].mean(0)
    return out


def fig_eif4e_result():
    t = md.load_pdb(str(fetch_pdb('1IPC')))
    L = Chem.MolFromMolFile(str(ROOT / 'dados' / 'eif4e' / 'ligante_cristal.sdf'), sanitize=False).GetConformer().GetPositions() * A
    near, X = pocket_scene(t, L, cut=0.9)
    trp = set(t.topology.select('resSeq 56 102 and not element H'))
    rows = _load_rows('eif4e')
    best = min(rows, key=lambda r: r['rmsd_A'])['pose']
    look = X[heavy(t, 'protein')].mean(0) - L.mean(0)
    B = view_basis(look); c0 = L.mean(0)
    fig, axs = plt.subplots(1, 3, figsize=(16, 5.8))
    panels = [(axs[0], None, 'Cristal (1IPC): a guanina metilada\nentre Trp56 e Trp102'),
              (axs[1], 1, f"Pose 1: {rows[0]['energia_kcal_mol']:.1f} kcal/mol, RMSD {rows[0]['rmsd_A']:.1f} Å"),
              (axs[2], best, f"Pose {best} (a mais próxima): {rows[best - 1]['energia_kcal_mol']:.1f} kcal/mol,\n"
                             f"RMSD {rows[best - 1]['rmsd_A']:.1f} Å")]
    for ax, k, title in panels:
        Q = _pose_xyz('eif4e', k)[0] if k else np.zeros((0, 3))
        P = project(np.r_[X[near], L, Q], B, c0)
        cols = [DESTAQUE if i in trp else PROT for i in near] + [CRISTAL] * len(L) + [POSE] * len(Q)
        rad = [R_ATOM] * len(near) + [0.16] * (len(L) + len(Q))
        draw_spheres(ax, P[:len(near)], cols[:len(near)], rad[:len(near)], pad=0.5)
        draw_spheres(ax, P[len(near):], cols[len(near):], rad[len(near):], pad=0.5, z=3)
        ax.set_xlim(P[:, 0].min() - 0.3, P[:, 0].max() + 0.3); ax.set_ylim(P[:, 1].min() - 0.3, P[:, 1].max() + 0.3)
        ax.set_title(title, fontsize=10)
    legend(fig, [('m7GTP no cristal', CRISTAL), ('pose do docking (Vina)', POSE), ('Trp56 e Trp102', DESTAQUE),
                 ('eIF4E: átomos a menos de 9 Å do ligante', PROT)])
    fig.suptitle('Prática B: o Vina não reproduz o cap entre os dois triptofanos', y=1.0, fontsize=11)
    save(fig, 'pratica_eif4e_poses.png')


# ---------------------------------------------------------------- 4. RMSD explicado com dados reais
def fig_rmsd():
    ref = skeleton(Chem.MolFromMolFile(str(ROOT / 'dados' / 'mpro' / 'ligante_cristal.sdf'), sanitize=False, removeHs=False))
    P0 = ref.GetConformer().GetPositions()
    rows = _load_rows('mpro')
    fig, axs = plt.subplots(1, 2, figsize=(13, 5))
    B = None
    for ax, k in zip(axs, (1, max(rows[1:6], key=lambda r: r['rmsd_A'])['pose'])):
        e, m = poses(resultados('mpro') / 'vina_dock.pdbqt')[k - 1]
        m = skeleton(m); Q = m.GetConformer().GetPositions()
        match = min(ref.GetSubstructMatches(m, uniquify=False, useChirality=False),
                    key=lambda mt: ((Q - P0[list(mt)]) ** 2).sum())
        Pm = P0[list(match)]
        if B is None:
            u, s, vt = np.linalg.svd(P0 - P0.mean(0)); B = np.c_[vt[0], vt[1], vt[2]]
        a = (Pm - P0.mean(0)) @ B; b = (Q - P0.mean(0)) @ B
        for i in range(len(a)):
            ax.plot([a[i, 0], b[i, 0]], [a[i, 1], b[i, 1]], color=ERRO, lw=0.7, alpha=0.7)
        ax.scatter(a[:, 0], a[:, 1], s=40, color=CRISTAL, edgecolor=INK, lw=0.4, zorder=3)
        ax.scatter(b[:, 0], b[:, 1], s=40, color=POSE, edgecolor=INK, lw=0.4, zorder=3)
        rmsd = np.sqrt(((Q - Pm) ** 2).sum(1).mean())
        ax.set_aspect('equal'); frame_off(ax)
        ax.set_title(f'Pose {k}: RMSD = {rmsd:.2f} Å', fontsize=10)
    legend(fig, [('átomo no cristal', CRISTAL), ('mesmo átomo na pose', POSE)], y=-0.06)
    fig.text(0.5, 0.06, r'RMSD = $\sqrt{\frac{1}{N}\sum_i d_i^2}$, com $d_i$ a distância (linha vermelha) '
             'entre as duas posições do átomo $i$', ha='center', fontsize=10)
    fig.suptitle('RMSD: a distância média entre a pose do docking e a do cristal (ensitrelvir, átomos pesados)',
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
    fig, axs = plt.subplots(1, 3, figsize=(15, 3.8))
    sq = np.array([[-0.5, 0.2], [-0.5, -0.62], [0.5, -0.62], [0.5, 0.2]])
    lig = np.array([[-0.42, 0.0], [-0.42, -0.55], [0.42, -0.55], [0.42, 0.0]])
    titles = ['Rígido-rígido: nada se move', 'Ligante flexível, receptor rígido\n(Vina, SwissDock, DockThor)',
              'Flexível-flexível: os dois se ajustam\n(Rosetta, HADDOCK, encaixe induzido)']
    for k, (ax, title) in enumerate(zip(axs, titles)):
        frame_off(ax); ax.set_aspect('equal'); ax.set_xlim(-1.8, 1.8); ax.set_ylim(-2.0, 0.9)
        _receptor(ax, 0.0, 0.0, sq, seed=5 + k)
        if k == 0:
            _ligand(ax, 0, 0, lig)
        else:
            for ang, al in ((-25, 0.35), (25, 0.35), (0, 1.0)):
                _ligand(ax, 0, 0.1 if ang else 0, lig, angle=ang, alpha=al)
        if k == 2:
            for dx in (-1, 1):
                ax.add_patch(FancyArrowPatch((dx * 0.95, -0.35), (dx * 0.62, -0.35), arrowstyle='<|-|>',
                                             mutation_scale=10, color=DESTAQUE, lw=1.4, zorder=6))
            ax.text(0, -1.85, 'as paredes do bolsão também se movem', ha='center', fontsize=8.5, color=DESTAQUE)
        if k == 1:
            ax.text(0, -1.85, 'esferas claras: outras orientações testadas', ha='center', fontsize=8.5, color=GREY)
        ax.set_title(title, fontsize=10)
    fig.suptitle('Tipos de docking pela flexibilidade (esquema)', y=1.02, fontsize=11)
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
    titles = ['Baseado em molde: parte de um complexo\nconhecido e parecido', 'Local: caixa só no sítio conhecido\n(esta prática)',
              'Cego (global): caixa na proteína inteira,\no programa procura o sítio']
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
              ('Docking com receptor rígido (Vina), 1 CPU', 30, 600, C[0]),
              ('Docking com cadeias laterais flexíveis', 600, 3600 * 2, C[6]),
              ('Ensemble docking (várias conformações da MD)', 1800, 3600 * 10, C[3]),
              ('MD do complexo + MM/GBSA, GPU', 3600 * 5, 3600 * 48, C[1]),
              ('Perturbação de energia livre (FEP), GPU', 3600 * 12, 3600 * 24 * 4, C[7])]
    fig, ax = plt.subplots(figsize=(11, 4.6))
    for k, (name, lo, hi, col) in enumerate(levels):
        y = len(levels) - 1 - k
        ax.plot([lo, hi], [y, y], color=col, lw=9, solid_capstyle='butt')
        ax.text(lo, y + 0.32, name, ha='left', va='bottom', fontsize=9)
    ax.plot(122, len(levels) - 2, 'o', color='white', mec=INK, ms=8, zorder=5)
    ax.annotate('medido nesta aula:\nMpro no SwissDock, 2 min', (122, len(levels) - 2), xytext=(4000, len(levels) - 1.55),
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


def fig_smiles_reading():
    """Aspirina: cada atomo do desenho numerado na ordem em que aparece na SMILES."""
    smi = 'CC(=O)Oc1ccccc1C(=O)O'
    m = Chem.MolFromSmiles(smi)
    for a in m.GetAtoms():
        a.SetProp('atomNote', str(a.GetIdx() + 1))
    img = Draw.MolToImage(m, size=(560, 440))
    fig, (a, b) = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw={'width_ratios': [1, 1.25]})
    a.imshow(img); frame_off(a)
    a.set_title('Aspirina: átomos numerados na ordem da SMILES', fontsize=10)
    frame_off(b); b.set_xlim(0, 10); b.set_ylim(0, 10)
    b.text(0.2, 9.0, smi, fontsize=15, family='monospace')
    rules = [('C, O, N', 'átomo alifático'), ('c, n, o', 'átomo aromático (minúscula)'),
             ('( )', 'ramificação'), ('=', 'ligação dupla'),
             ('1 ... 1', 'abre e fecha um anel'), ('[n+], [O-]', 'carga entre colchetes'),
             ('@, @@', 'quiralidade')]
    for k, (sym, mean) in enumerate(rules):
        b.text(0.2, 7.4 - 1.0 * k, sym, fontsize=11, family='monospace', color=C[0])
        b.text(3.0, 7.4 - 1.0 * k, mean, fontsize=10)
    b.set_title('Como ler', fontsize=10)
    fig.suptitle('SMILES: a molécula escrita como uma linha de texto (Weininger, 1988)', y=1.0, fontsize=11)
    save(fig, 'fund_smiles_leitura.png')


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


def fig_barnase_barstar():
    t = md.load_pdb(str(fetch_pdb('1BRS')))
    a = heavy(t, 'chainid 0'); b = heavy(t, 'chainid 3')
    X = t.xyz[0]
    d = np.min(np.linalg.norm(X[a][:, None] - X[b][None], axis=2), axis=1)
    ia = set(a[d < 0.45])
    d2 = np.min(np.linalg.norm(X[b][:, None] - X[a][None], axis=2), axis=1)
    ib = set(b[d2 < 0.45])
    axis = X[b].mean(0) - X[a].mean(0)
    B = view_basis(np.cross(axis, [0, 0, 1.0]), up_hint=axis / np.linalg.norm(axis)); c0 = X[np.r_[a, b]].mean(0)
    P = project(X[np.r_[a, b]], B, c0)
    cols = [shade(C[0], 0.0) if i in ia else '#b8c9e0' for i in a] + [C[5] if i in ib else '#ecc7a8' for i in b]
    fig, ax = plt.subplots(figsize=(7.6, 7.2))
    draw_spheres(ax, P, cols, [R_ATOM] * len(P), pad=0.6)
    ax.set_title('Docking proteína-proteína: barnase e barstar (PDB 1BRS)', fontsize=10.5)
    legend(fig, [('barnase', '#b8c9e0'), ('barstar', '#ecc7a8'), ('interface da barnase', C[0]), ('interface da barstar', C[5])],
           ncol=2, y=0.0)
    save(fig, 'fund_barnase_barstar.png')


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
    rows_m, rows_e = _load_rows('mpro'), _load_rows('eif4e')
    best_e = min(rows_e, key=lambda r: r['rmsd_A'])
    fig, axs = plt.subplots(2, 2, figsize=(13.5, 8.6))
    # (a) energia e Kd estimado
    a = axs[0, 0]
    g = np.linspace(-13, -3, 300)
    a.axvspan(-13, -9, color=shade(C[3], 0.75), lw=0); a.axvspan(-9, -6, color=shade(C[1], 0.8), lw=0)
    a.axvspan(-6, -3, color=shade(ERRO, 0.8), lw=0)
    for x, t in ((-11, 'forte'), (-7.5, 'moderada'), (-4.5, 'fraca')):
        a.text(x, 0.92, t, ha='center', va='top', fontsize=9, transform=a.get_xaxis_transform())
    for e, lab, col, dy in ((rows_m[0]['energia_kcal_mol'], 'Mpro, pose 1', C[0], 0.62),
                            (rows_e[0]['energia_kcal_mol'], 'eIF4E, pose 1', C[7], 0.38)):
        kd = np.exp(e / RT)
        txt = f'{kd * 1e9:.0f} nM' if kd < 1e-6 else f'{kd * 1e6:.0f} µM'
        a.axvline(e, color=col, lw=1.6)
        a.text(e + 0.12, dy, f'{lab}\n{e:.1f} kcal/mol\nKd ≈ {txt}', fontsize=8.5, va='center',
               transform=a.get_xaxis_transform(), color=col)
    a.set_xlim(-13, -3); a.set_yticks([])
    a.set_xlabel('Energia do docking (kcal/mol); mais negativa = ligação mais forte')
    a.set_title('Energia: Kd ≈ exp(ΔG/RT), uma estimativa grosseira', fontsize=10)
    # (b) RMSD
    b = axs[0, 1]
    b.axvspan(0, 2, color=shade(C[3], 0.75), lw=0); b.axvspan(2, 3, color=shade(C[1], 0.8), lw=0)
    b.axvspan(3, 10, color=shade(ERRO, 0.8), lw=0)
    for x, t in ((1, 'correta'), (2.5, 'quase'), (6.5, 'errada')):
        b.text(x, 0.92, t, ha='center', va='top', fontsize=9, transform=b.get_xaxis_transform())
    marks = [(rows_m[0]['rmsd_A'], 'Mpro\npose 1', C[0], 0.6), (rows_m[1]['rmsd_A'], 'Mpro\npose 2', C[0], 0.25),
             (rows_e[0]['rmsd_A'], 'eIF4E\npose 1', C[7], 0.6),
             (best_e['rmsd_A'], f"eIF4E\npose {best_e['pose']}", C[7], 0.25)]
    for r, lab, col, y in marks:
        b.plot(r, y, 'o', color=col, mec=INK, ms=9, transform=b.get_xaxis_transform())
        b.text(r + 0.15, y, f'{lab}: {r:.1f} Å', fontsize=8.5, va='center', transform=b.get_xaxis_transform())
    b.set_xlim(0, 10); b.set_yticks([])
    b.set_xlabel('RMSD em relação ao cristal (Å); menor = mais parecida')
    b.set_title('RMSD: só existe quando há uma pose de referência', fontsize=10)
    # (c) recuperacao das interacoes do cristal
    c = axs[1, 0]
    labels, vals, cols = [], [], []
    for nome, col, poses_ in (('mpro', C[0], (1, 2, 3)), ('eif4e', C[7], (1, 2, 3, best_e['pose']))):
        it = _interactions(nome); ref = _pairs(it['cristal'])
        for k in poses_:
            key = f'pose {k}'
            if key in it:
                labels.append(f"{'Mpro' if nome == 'mpro' else 'eIF4E'} pose {k}")
                vals.append(100 * len(ref & _pairs(it[key])) / len(ref)); cols.append(col)
    y = np.arange(len(labels))[::-1]
    c.barh(y, vals, color=cols, edgecolor=INK, lw=0.5, height=0.6)
    for yy, v in zip(y, vals):
        c.text(v + 1.5, yy, f'{v:.0f}%', va='center', fontsize=8.5)
    c.set_yticks(y); c.set_yticklabels(labels, fontsize=9)
    c.set_xlim(0, 110); c.set_xlabel('Interações do cristal recuperadas pela pose (%)')
    c.tick_params(axis='y', right=False)
    c.set_title('Interações: a pose repete os contatos do cristal?', fontsize=10)
    # (d) eficiencia do ligante
    d = axs[1, 1]
    n_m = Chem.MolFromMolFile(str(ROOT / 'dados' / 'mpro' / 'ligante_cristal.sdf')).GetNumHeavyAtoms()
    n_e = Chem.MolFromMolFile(str(ROOT / 'dados' / 'eif4e' / 'ligante_cristal.sdf'), sanitize=False).GetNumHeavyAtoms()
    le = [(-rows_m[0]['energia_kcal_mol'] / n_m, f'ensitrelvir\n{n_m} átomos pesados', C[0]),
          (-rows_e[0]['energia_kcal_mol'] / n_e, f'm7GTP\n{n_e} átomos pesados', C[7])]
    d.bar([0, 1], [v for v, _, _ in le], color=[c_ for _, _, c_ in le], edgecolor=INK, lw=0.5, width=0.55)
    for k, (v, _, _) in enumerate(le):
        d.text(k, v + 0.008, f'{v:.2f}', ha='center', fontsize=9)
    d.axhline(0.3, color=INK, ls=(0, (4, 3)), lw=1)
    d.text(1.45, 0.305, 'referência usual: 0,3', ha='right', va='bottom', fontsize=8.5)
    d.set_xticks([0, 1]); d.set_xticklabels([t for _, t, _ in le], fontsize=9)
    d.set_xlim(-0.6, 1.6); d.set_ylim(0, 0.4)
    d.set_ylabel('kcal/mol por átomo pesado')
    d.tick_params(axis='x', top=False)
    d.set_title('Eficiência do ligante: energia dividida pelo tamanho', fontsize=10)
    fig.suptitle('Métricas do docking com os resultados da prática', y=0.995, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save(fig, 'fund_metricas.png')


def fig_fingerprint():
    """Impressao digital de interacoes: residuo x (cristal, poses), cor pelo tipo."""
    tipo_cor = {'ligação de hidrogênio': C[0], 'ponte salina': ERRO, 'hidrofóbica': C[1], 'empilhamento π': C[4],
                'cátion-π': C[6], 'halogênio': C[3], 'ponte de água': WATER}
    fig, axs = plt.subplots(1, 2, figsize=(14, 4.8), gridspec_kw={'width_ratios': [1, 1]})
    best_e = min(_load_rows('eif4e'), key=lambda r: r['rmsd_A'])['pose']
    for ax, nome, title, keys in ((axs[0], 'mpro', 'Mpro + ensitrelvir', ['cristal', 'pose 1', 'pose 2', 'pose 3']),
                                  (axs[1], 'eif4e', 'eIF4E + m7GTP', ['cristal', 'pose 1', 'pose 2', f'pose {best_e}'])):
        it = _interactions(nome)
        res = sorted({r['residuo'] for k in keys for r in it.get(k, [])}, key=lambda s: int(''.join(ch for ch in s if ch.isdigit())))
        for j, k in enumerate(keys):
            for i, rname in enumerate(res):
                tipos = sorted({r['tipo'] for r in it.get(k, []) if r['residuo'] == rname})
                ax.add_patch(Rectangle((i - 0.45, j - 0.45), 0.9, 0.9, fc='white', ec=BAND, lw=0.8))
                for t_i, t in enumerate(tipos):
                    w = 0.9 / len(tipos)
                    ax.add_patch(Rectangle((i - 0.45 + t_i * w, j - 0.45), w, 0.9, fc=tipo_cor[t], ec='white', lw=0.5))
        ax.set_xlim(-0.6, len(res) - 0.4); ax.set_ylim(len(keys) - 0.4, -0.6)
        ax.set_xticks(range(len(res))); ax.set_xticklabels(res, rotation=60, ha='right', fontsize=8.5)
        ax.set_yticks(range(len(keys))); ax.set_yticklabels(keys, fontsize=9)
        ax.tick_params(top=False, right=False, length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_title(title, fontsize=10)
    hs = [Rectangle((0, 0), 1, 1, fc=c_) for c_ in tipo_cor.values()]
    fig.legend(hs, list(tipo_cor), loc='lower center', ncol=7, frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, -0.08))
    fig.suptitle('Impressão digital das interações (PLIP): cada coluna é um resíduo do receptor', y=1.02, fontsize=11)
    save(fig, 'pratica_interacoes.png')


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


def fig_pocket_ranking():
    """Escore das cavidades previstas em quatro estruturas; marcadas as que sao sitios reais."""
    casos = [('8HBK', 'Mpro + ensitrelvir', {1: 'sítio ativo'}),
             ('5MO4', 'ABL + nilotinibe e asciminibe', {1: 'sítio do ATP', 3: 'bolso do miristoil'}),
             ('6OIM', 'KRAS G12C + sotorasibe', {1: 'sítio do GDP', 2: 'bolso switch II'}),
             ('4OBE', 'KRAS sem inibidor', {1: 'sítio do GDP (cópia B)', 2: 'sítio do GDP', 3: 'switch II ainda raso'})]
    fig, axs = plt.subplots(1, 4, figsize=(16, 4.2), sharey=True)
    for ax, (code, titulo, reais) in zip(axs, casos):
        pk = _pockets(code)
        x = [r['rank'] for r in pk]; y = [r['score'] for r in pk]
        cols = [CRISTAL if r['rank'] in reais else ERRO for r in pk]
        ax.bar(x, y, color=cols, edgecolor=INK, lw=0.5, width=0.7)
        for r in pk:
            if r['rank'] in reais:
                ax.text(r['rank'], r['score'] + 0.6, reais[r['rank']], ha='center', fontsize=8, rotation=90, va='bottom')
        ax.set_xticks(x); ax.tick_params(top=False, right=False)
        ax.set_xlabel('posto da cavidade')
        ax.set_title(f'{titulo}\n({code})', fontsize=9.5)
        ax.set_ylim(0, 44)
    axs[0].set_ylabel('escore do P2Rank')
    legend(fig, [('cavidade com ligante conhecido', CRISTAL), ('cavidade sem ligante no cristal', ERRO)], y=-0.06)
    fig.suptitle('Nem toda cavidade é um sítio: cavidades previstas pelo P2Rank (PrankWeb)', y=1.03, fontsize=11)
    save(fig, 'fund_cavidades_ranking.png')


def main():
    fig_binding_models(); fig_vocabulary(); fig_forces(); fig_flexibility(); fig_search(); fig_levels_cost()
    fig_smiles(); fig_smiles_reading(); fig_resistance(); fig_p53_mdm2(); fig_barnase_barstar(); fig_md_docking()
    fig_cavidades(); fig_tipos_sitio(); fig_pocket_ranking()
    fig_poses_scatter(); fig_mpro_result(); fig_eif4e_result(); fig_rmsd(); fig_metrics(); fig_fingerprint()


if __name__ == '__main__':
    main()
