"""Cada pose do Vina contra o cristal: RMSD do ligante inteiro, RMSD so da base e contatos recuperados.

A base e a guanina metilada (aneis com nitrogenio e o metil do N7), a parte que o
sitio reconhece. Os contatos vem de interacoes.json (tools/interactions.py).
Grava resultados/<sistema>/comparacao_cristal.json.

Uso: python tools/compare_crystal.py eif4e
     python tools/compare_crystal.py eif4e_diffdock --receptor eif4e
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from rdkit import Chem

from analyze_poses import poses, skeleton

ROOT = Path(__file__).resolve().parents[1]


def base_atoms(ref):
    """Atomos dos aneis com nitrogenio e os metis presos a eles."""
    base = set()
    for ring in ref.GetRingInfo().AtomRings():
        if any(ref.GetAtomWithIdx(i).GetSymbol() == 'N' for i in ring):
            base |= set(ring)
    for i in list(base):
        for n in ref.GetAtomWithIdx(i).GetNeighbors():
            if n.GetSymbol() == 'C' and n.GetDegree() == 1:
                base.add(n.GetIdx())
    return sorted(base)


def main():
    ap = argparse.ArgumentParser(description='Poses do Vina contra o cristal')
    ap.add_argument('sistema')
    ap.add_argument('--receptor', help='sistema de dados/ quando difere do nome da pasta de resultados')
    a = ap.parse_args()
    res = ROOT / 'resultados' / a.sistema
    ref = skeleton(Chem.MolFromMolFile(str(ROOT / 'dados' / (a.receptor or a.sistema) / 'ligante_cristal.sdf'),
                                       sanitize=False, removeHs=False))
    P0 = ref.GetConformer().GetPositions()
    base = base_atoms(ref)
    it = json.loads((res / 'interacoes.json').read_text())
    cristal = {(r['tipo'], r['residuo']) for r in it['cristal']}
    out = []
    if (res / 'vina_dock.pdbqt').exists():
        lista = poses(res / 'vina_dock.pdbqt')
    else:   # DiffDock: SDF, a nota fica no titulo de cada bloco ('pose_1 confianca -0.528')
        lista = [(float(m.GetProp('_Name').split()[-1]), m)
                 for m in Chem.SDMolSupplier(str(res / 'poses.sdf'), sanitize=False, removeHs=False)]
    nota = 'energia_kcal_mol' if (res / 'vina_dock.pdbqt').exists() else 'confianca'
    for k, (e, m) in enumerate(lista, 1):
        m = skeleton(m)
        Q = m.GetConformer().GetPositions()
        # mapeamento atomo a atomo que minimiza o desvio (simetria do grafo)
        d, mt = min(((((Q - P0[list(mt)]) ** 2).sum(1), mt)
                     for mt in ref.GetSubstructMatches(m, uniquify=False, useChirality=False, maxMatches=5000)),
                    key=lambda x: x[0].sum())
        pos = {r: q for q, r in enumerate(mt)}
        contatos = {(r['tipo'], r['residuo']) for r in it.get(f'pose {k}', [])}
        out.append({'pose': k, nota: e, 'rmsd_A': round(float(np.sqrt(d.mean())), 2),
                    'rmsd_base_A': round(float(np.sqrt(np.mean([d[pos[i]] for i in base]))), 2),
                    'contatos_do_cristal': len(cristal & contatos), 'contatos_no_cristal': len(cristal),
                    'sanduiche_trp56_trp102': {('empilhamento π', 'TRP56'), ('empilhamento π', 'TRP102')} <= contatos})
    (res / 'comparacao_cristal.json').write_text(json.dumps(out, indent=1, ensure_ascii=False) + '\n')
    for r in out:
        print(r)


if __name__ == '__main__':
    main()
