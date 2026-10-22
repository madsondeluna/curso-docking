"""Estilo das figuras do curso, o mesmo do curso de dinamica molecular.

Moldura completa com ticks para dentro, titulo centrado, fundo branco; esferas
ordenadas por profundidade para estruturas (projecao ortografica).
"""
from __future__ import annotations

import urllib.request
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.colors import to_rgb
from matplotlib.patches import Circle

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'imagens'
CACHE = ROOT / 'runs' / 'pdb'
DPI = 300

INK = '#0d1321'
BAND = '#dde1e9'
GREY = '#9a9a9a'
C = ['#3973b1', '#9f8322', '#9e527f', '#4c985f', '#745ba5', '#ba6f3e', '#1990ad', '#ac5551']
LINESTYLES = ['-', (0, (6, 3)), (0, (1, 2))]

# cores das estruturas
PROT = '#c9c2b4'        # receptor
SITE = '#8fa6c4'        # residuos do sitio
CRISTAL = '#4c985f'     # ligante do cristal
POSE = '#dd8452'        # pose do docking
DESTAQUE = '#8172b2'    # residuo em destaque
ERRO = '#c44e52'
WATER = '#a6cfe3'

STYLE = {
    'axes.spines.top': True, 'axes.spines.right': True,
    'axes.edgecolor': INK, 'axes.linewidth': 0.8, 'axes.grid': False,
    'axes.titlelocation': 'center', 'axes.titlesize': 11, 'axes.titleweight': 'normal',
    'axes.labelcolor': INK, 'axes.labelsize': 10,
    'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
    'xtick.color': INK, 'ytick.color': INK, 'text.color': INK,
    'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white',
    'legend.frameon': False, 'legend.fontsize': 9,
    'font.family': 'sans-serif', 'font.sans-serif': ['Public Sans', 'Helvetica Neue', 'Arial', 'DejaVu Sans'],
    'font.size': 10,
}
plt.rcParams.update(STYLE)


def save(fig, name, margem=28):
    """Grava a figura e recorta a borda em branco que sobra quando o eixo tem aspecto fixo."""
    OUT.mkdir(exist_ok=True)
    caminho = OUT / name
    fig.savefig(caminho, dpi=DPI, bbox_inches='tight')
    plt.close(fig)
    try:
        from PIL import Image
        im = Image.open(caminho).convert('RGB')
        a = np.asarray(im.convert('L'))
        tinta = np.argwhere(a < 245)
        if len(tinta):
            y0, x0 = tinta.min(0); y1, x1 = tinta.max(0)
            caixa = (max(int(x0) - margem, 0), max(int(y0) - margem, 0),
                     min(int(x1) + margem + 1, im.width), min(int(y1) + margem + 1, im.height))
            im.crop(caixa).save(caminho)
    except ImportError:
        pass
    print('  ', name)


def frame_off(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def fetch_pdb(code):
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f'{code}.pdb'
    if not path.exists():
        urllib.request.urlretrieve(f'https://files.rcsb.org/download/{code}.pdb', path)
    return path


def shade(color, back):
    """Clareia o que esta ao fundo (back = 0 frente, 1 fundo)."""
    c = np.array(to_rgb(color))
    return tuple(c + (1 - c) * 0.55 * back)


def draw_spheres(ax, P, colors, radii, pad=0.4, z=2):
    """Projecao ortografica: P (n, 3) com colunas (tela x, tela y, profundidade); fundo primeiro."""
    d = P[:, 2]
    dn = (d - d.min()) / max(np.ptp(d), 1e-6)
    order = np.argsort(dn)
    patches = [Circle((P[k, 0], P[k, 1]), radii[k]) for k in order]
    faces = [shade(colors[k], 1 - dn[k]) for k in order]
    ax.add_collection(PatchCollection(patches, facecolors=faces, edgecolors=[shade(INK, 0.6)] * len(order),
                                      linewidths=0.25, match_original=False, zorder=z))
    ax.set_aspect('equal'); frame_off(ax)
    ax.set_xlim(P[:, 0].min() - pad, P[:, 0].max() + pad); ax.set_ylim(P[:, 1].min() - pad, P[:, 1].max() + pad)


def view_basis(look, up_hint=None):
    """Base de tela (x, y) e profundidade; quem esta mais perto do observador tem profundidade maior.
    `look` aponta do observador para a cena."""
    w = look / np.linalg.norm(look)
    ref = up_hint if up_hint is not None else (np.array([0, 0, 1.0]) if abs(w[2]) < 0.9 else np.array([1.0, 0, 0]))
    u = np.cross(ref, w); u /= np.linalg.norm(u)
    v = np.cross(w, u)
    return np.c_[u, v, -w]


def project(xyz, basis, center):
    return (xyz - center) @ basis
