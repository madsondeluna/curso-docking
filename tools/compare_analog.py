"""Poses do DiffDock para o analogo do ensitrelvir na Mpro predita (OpenFold) contra o cristal 8HBK.

O receptor do exemplo e um modelo do OpenFold, em outro sistema de coordenadas. O
script sobrepoe o modelo a cadeia A da 8HBK pelos carbonos alfa, leva as poses junto,
mede o RMSD no nucleo que o analogo divide com o ensitrelvir (maior subestrutura comum)
e compara os contatos (PLIP) com os do ensitrelvir no cristal.

Grava em resultados/mpro_diffdock/: poses_no_8hbk.sdf, receptor_openfold_no_8hbk.pdb,
interacoes.json e comparacao_cristal.json.

Uso: python tools/compare_analog.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import rdFMCS

from interactions import complex_pdb, ligand_lines, plip

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / 'resultados' / 'mpro_diffdock'
ENTRADA = ROOT / 'dados' / 'mpro_analogo'


def ca(path):
    """Carbonos alfa por numero de residuo."""
    out = {}
    for l in Path(path).read_text().splitlines():
        if l.startswith('ATOM') and l[12:16].strip() == 'CA':
            out[int(l[22:26])] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return out


def kabsch(P, Q):
    """Rotacao R e translacao t que levam P sobre Q (P @ R.T + t)."""
    pc, qc = P.mean(0), Q.mean(0)
    u, _, vt = np.linalg.svd((P - pc).T @ (Q - qc))
    dd = np.sign(np.linalg.det(vt.T @ u.T))
    R = vt.T @ np.diag([1, 1, dd]) @ u.T
    return R, qc - pc @ R.T


def move_pdb(path, R, t):
    out = []
    for l in Path(path).read_text().splitlines():
        if l.startswith(('ATOM', 'HETATM')):
            x = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]) @ R.T + t
            l = f'{l[:30]}{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}{l[54:]}'
        out.append(l)
    return '\n'.join(out) + '\n'


def main():
    cristal_pdb = ROOT / 'dados' / 'mpro' / 'receptor.pdb'          # cadeia A da 8HBK
    a, b = ca(ENTRADA / 'receptor_openfold.pdb'), ca(cristal_pdb)
    comum = sorted(set(a) & set(b))
    R, t = kabsch(np.array([a[i] for i in comum]), np.array([b[i] for i in comum]))
    rmsd_ca = float(np.sqrt((((np.array([a[i] for i in comum]) @ R.T + t) - np.array([b[i] for i in comum])) ** 2)
                            .sum(1).mean()))
    rec = RES / 'receptor_openfold_no_8hbk.pdb'
    rec.write_text(move_pdb(ENTRADA / 'receptor_openfold.pdb', R, t))

    # poses levadas para o referencial da 8HBK
    w = Chem.SDWriter(str(RES / 'poses_no_8hbk.sdf'))
    lista = []
    for m in Chem.SDMolSupplier(str(RES / 'poses.sdf'), removeHs=False):
        c = m.GetConformer()
        X = c.GetPositions() @ R.T + t
        for i in range(m.GetNumAtoms()):
            c.SetAtomPosition(i, X[i].tolist())
        w.write(m)
        lista.append(m)
    w.close()

    # nucleo comum com o ensitrelvir do cristal
    ens = Chem.MolFromMolFile(str(ROOT / 'dados' / 'mpro' / 'ligante_cristal.sdf'))
    analogo = Chem.RemoveHs(lista[0])
    mcs = Chem.MolFromSmarts(rdFMCS.FindMCS([ens, analogo], timeout=20, ringMatchesRingOnly=True,
                                            completeRingsOnly=True).smartsString)
    E = ens.GetConformer().GetPositions()
    ie = ens.GetSubstructMatch(mcs)

    # contatos: ensitrelvir no cristal contra cada pose no receptor predito
    cr = plip(complex_pdb(cristal_pdb, ligand_lines(next(pybel.readfile('sdf', str(ROOT / 'dados' / 'mpro' /
                                                                                     'ligante_cristal.sdf'))))))
    inter = {'cristal (ensitrelvir, 8HBK)': cr}
    pb = list(pybel.readfile('sdf', str(RES / 'poses_no_8hbk.sdf')))
    ref = {(r['tipo'], r['residuo']) for r in cr}
    out = []
    for k, (m, mp) in enumerate(zip(lista, pb), 1):
        mh = Chem.RemoveHs(m)
        Q = mh.GetConformer().GetPositions()
        d = min(np.sqrt(((Q[list(ia)] - E[list(ie)]) ** 2).sum(1).mean())
                for ia in mh.GetSubstructMatches(mcs, uniquify=False))
        rows = plip(complex_pdb(rec, ligand_lines(mp)))
        inter[f'pose {k}'] = rows
        out.append({'pose': k, 'confianca': float(m.GetProp('_Name').split()[-1]), 'rmsd_nucleo_A': round(float(d), 2),
                    'contatos_do_cristal': len(ref & {(r['tipo'], r['residuo']) for r in rows}),
                    'contatos_no_cristal': len(ref)})
    (RES / 'interacoes.json').write_text(json.dumps(inter, indent=1, ensure_ascii=False) + '\n')
    resumo = {'rmsd_ca_openfold_8hbk_A': round(rmsd_ca, 2), 'residuos_alinhados': len(comum),
              'atomos_nucleo_comum': mcs.GetNumAtoms(), 'atomos_ensitrelvir': ens.GetNumHeavyAtoms(),
              'atomos_analogo': analogo.GetNumHeavyAtoms(), 'poses': out}
    (RES / 'comparacao_cristal.json').write_text(json.dumps(resumo, indent=1, ensure_ascii=False) + '\n')
    print({k: v for k, v in resumo.items() if k != 'poses'})
    for r in out:
        print(r)


if __name__ == '__main__':
    main()
