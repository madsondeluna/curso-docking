"""Energia e RMSD de cada pose do Vina contra a pose do cristal.

O RMSD usa so os atomos pesados, sem realinhar (receptor e cristal estao no
mesmo referencial) e considera a simetria do ligante (RDKit CalcRMS). O
mapeamento entre os atomos usa so a conectividade dos atomos pesados.

Uso: python tools/analyze_poses.py eif4e [mpro eif4e_ac] [--fonte runs/swissdock]

Aceita as poses do Vina (vina_dock.pdbqt) e do Attracting Cavities (result.dock4).
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from openbabel import openbabel as ob

ob.obErrorLog.SetOutputLevel(0)

ROOT = Path(__file__).resolve().parents[1]


def skeleton(m):
    """Grafo dos atomos pesados: ligacoes simples, sem carga nem aromaticidade.
    PDBQT e SDF do cristal descrevem ordens de ligacao de jeitos diferentes; a
    conectividade basta para o mapeamento atomo a atomo do RMSD."""
    rw = Chem.RWMol(m)
    for i in sorted((a.GetIdx() for a in rw.GetAtoms() if a.GetAtomicNum() == 1), reverse=True):
        rw.RemoveAtom(i)
    for a in rw.GetAtoms():
        a.SetFormalCharge(0); a.SetIsAromatic(False); a.SetNoImplicit(True); a.SetNumExplicitHs(0)
    for b in rw.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE); b.SetIsAromatic(False)
    m = rw.GetMol()
    m.UpdatePropertyCache(strict=False)
    return m


def poses(pdbqt):
    """Cada MODEL do PDBQT como molecula RDKit sem hidrogenios, com a energia do Vina."""
    text = Path(pdbqt).read_text()
    energias = [float(x) for x in re.findall(r'REMARK VINA RESULT:\s+(-?\d+\.\d+)', text)]
    mols = [Chem.MolFromMolBlock(m.write('mol'), sanitize=False) for m in pybel.readfile('pdbqt', str(pdbqt))]
    return list(zip(energias, mols))


def poses_ac(dock4):
    """Poses do Attracting Cavities (result.dock4): blocos REMARK + ATOM, energia SwissParam (SP-dG)."""
    out = []
    for block in Path(dock4).read_text().split('REMARK FILE_NAME')[1:]:
        e = float(re.search(r'REMARK SP-dG\s*:\s*(-?\d+\.?\d*)', block).group(1))
        atoms = [l for l in block.splitlines() if l.startswith(('ATOM', 'HETATM'))]
        m = pybel.readstring('pdb', '\n'.join(atoms) + '\nEND\n')
        out.append((e, Chem.MolFromMolBlock(m.write('mol'), sanitize=False)))
    return out


def main():
    ap = argparse.ArgumentParser(description='Energia e RMSD das poses contra o cristal')
    ap.add_argument('sistemas', nargs='+')
    ap.add_argument('--fonte', default='runs/swissdock')
    a = ap.parse_args()
    for nome in a.sistemas:
        d = ROOT / 'dados' / nome.split('_')[0]
        ref = skeleton(Chem.MolFromMolFile(str(d / 'ligante_cristal.sdf'), sanitize=False, removeHs=False))
        rows = []
        base = ROOT / a.fonte / nome
        lista = poses(base / 'vina_dock.pdbqt') if (base / 'vina_dock.pdbqt').exists() else poses_ac(base / 'result.dock4')
        for k, (e, m) in enumerate(lista, 1):
            m = skeleton(m)
            rows.append({'pose': k, 'energia_kcal_mol': e, 'rmsd_A': round(rdMolAlign.CalcRMS(m, ref), 2)})
        out = ROOT / a.fonte / nome / 'poses.json'
        out.write_text(json.dumps(rows, indent=1) + '\n')
        best = min(rows, key=lambda r: r['rmsd_A'])
        print(nome, 'pose 1:', rows[0], '| menor RMSD:', best)
        for r in rows:
            print(f"  {r['pose']:>2} {r['energia_kcal_mol']:>7.2f} {r['rmsd_A']:>6.2f}")


if __name__ == '__main__':
    main()
