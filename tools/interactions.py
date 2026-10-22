"""Interacoes proteina-ligante (PLIP) no cristal e nas poses do docking.

Monta o complexo receptor + ligante (cristal ou pose) num PDB e roda o PLIP,
o mesmo programa do servidor web usado na pratica. Grava
resultados/<sistema>/interacoes.json com as interacoes por residuo.

Uso: python tools/interactions.py mpro eif4e [--poses 1 2 3]
"""
from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path

from openbabel import pybel
from openbabel import openbabel as ob
from plip.structure.preparation import PDBComplex

ob.obErrorLog.SetOutputLevel(0)
ROOT = Path(__file__).resolve().parents[1]
TIPOS = {'hbonds_ldon': 'ligação de hidrogênio', 'hbonds_pdon': 'ligação de hidrogênio',
         'saltbridge_lneg': 'ponte salina', 'saltbridge_pneg': 'ponte salina',
         'hydrophobic_contacts': 'hidrofóbica', 'pistacking': 'empilhamento π',
         'pication_laro': 'cátion-π', 'pication_paro': 'cátion-π', 'halogen_bonds': 'halogênio',
         'water_bridges': 'ponte de água', 'metal_complexes': 'metal'}


def complex_pdb(receptor, ligand_pdb_lines):
    rec = [l for l in Path(receptor).read_text().splitlines() if l.startswith('ATOM')]
    return '\n'.join(rec + ['TER'] + ligand_pdb_lines + ['END']) + '\n'


def ligand_lines(mol):
    """Ligante (pybel) como HETATM com resname LIG, cadeia Z."""
    out = []
    for line in mol.write('pdb').splitlines():
        if line.startswith(('ATOM', 'HETATM')):
            line = 'HETATM' + line[6:17] + 'LIG Z 900' + line[26:]
            out.append(line)
        elif line.startswith('CONECT'):
            out.append(line)
    return out


def plip(pdb_text):
    with tempfile.NamedTemporaryFile('w', suffix='.pdb', delete=False) as f:
        f.write(pdb_text)
    mol = PDBComplex()
    mol.load_pdb(f.name)
    mol.analyze()
    site = next(v for k, v in mol.interaction_sets.items() if k.startswith('LIG'))
    rows = []
    for attr, nome in TIPOS.items():
        for it in getattr(site, attr, []):
            rows.append({'tipo': nome, 'residuo': f'{it.restype}{it.resnr}'})
    return rows


def main():
    ap = argparse.ArgumentParser(description='Interacoes PLIP do cristal e das poses')
    ap.add_argument('sistemas', nargs='+')
    ap.add_argument('--poses', nargs='+', type=int, default=[1, 2, 3])
    a = ap.parse_args()
    for nome in a.sistemas:
        d = ROOT / 'dados' / nome
        res = ROOT / 'runs' / 'swissdock' / nome
        res = res if res.exists() else ROOT / 'resultados' / nome
        out = {}
        cristal = next(pybel.readfile('sdf', str(d / 'ligante_cristal.sdf')))
        out['cristal'] = plip(complex_pdb(d / 'receptor.pdb', ligand_lines(cristal)))
        poses = list(pybel.readfile('pdbqt', str(res / 'vina_dock.pdbqt')))
        for k in a.poses:
            out[f'pose {k}'] = plip(complex_pdb(d / 'receptor.pdb', ligand_lines(poses[k - 1])))
        (res / 'interacoes.json').write_text(json.dumps(out, indent=1, ensure_ascii=False) + '\n')
        for key, rows in out.items():
            resumo = {}
            for r in rows:
                resumo.setdefault(r['tipo'], set()).add(r['residuo'])
            print(nome, key, {k: sorted(v) for k, v in resumo.items()})


if __name__ == '__main__':
    main()
