"""Prepara as entradas da pratica a partir do PDB (dados/<sistema>/).

Para cada sistema: receptor limpo (so a cadeia de proteina, sem agua nem
ligantes), pose cristalografica do ligante (referencia para o RMSD), SMILES
em pH 7,4, ligante 3D com hidrogenios (SDF e MOL2, ponto de partida para o
docking) e a caixa de busca centrada no ligante do cristal.

Uso: python tools/prepare_inputs.py
"""
from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import numpy as np
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / 'runs' / 'pdb'
DADOS = ROOT / 'dados'

# sistema: (codigo PDB, cadeia, codigo do ligante, folga da caixa em A)
SISTEMAS = {
    'eif4e': ('1IPC', 'A', 'MGP', 8.0),
    'mpro': ('8HBK', 'A', '7YY', 8.0),
}


def fetch(url, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        urllib.request.urlretrieve(url, path)
    return path


def receptor(pdb_path, chain, out):
    """So ATOM da cadeia pedida; na alternativa de conformacao, fica a primeira (A ou vazia)."""
    keep = []
    for line in pdb_path.read_text().splitlines():
        if line.startswith('ATOM') and line[21] == chain and line[16] in ' A':
            keep.append(line[:16] + ' ' + line[17:])
    out.write_text('\n'.join(keep + ['TER', 'END']) + '\n')
    return len(keep)


def smiles_ph74(code):
    """SMILES do dicionario do PDB, protonado em pH 7,4 pelo Open Babel."""
    d = json.loads(urllib.request.urlopen(f'https://data.rcsb.org/rest/v1/core/chemcomp/{code}').read())
    smi = next(x['descriptor'] for x in d['pdbx_chem_comp_descriptor']
               if x['program'] == 'OpenEye OEToolkits' and x['type'] == 'SMILES_CANONICAL')
    mol = pybel.readstring('smi', smi)
    mol.OBMol.AddHydrogens(False, True, 7.4)
    can = Chem.MolToSmiles(Chem.MolFromSmiles(mol.write('can').split()[0]))
    return smi, can


def ligand_3d(smi, sdf_out, mol2_out, seed=7):
    mol = Chem.AddHs(Chem.MolFromSmiles(smi))
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    AllChem.MMFFOptimizeMolecule(mol)
    Chem.MolToMolFile(mol, str(sdf_out))
    pybel.readfile('sdf', str(sdf_out)).__next__().write('mol2', str(mol2_out), overwrite=True)
    return mol.GetNumAtoms(), Chem.rdMolDescriptors.CalcNumRotatableBonds(mol), Chem.GetFormalCharge(mol)


def main():
    resumo = {}
    for nome, (code, chain, lig, folga) in SISTEMAS.items():
        d = DADOS / nome
        d.mkdir(parents=True, exist_ok=True)
        pdb = fetch(f'https://files.rcsb.org/download/{code}.pdb', CACHE / f'{code}.pdb')
        n_at = receptor(pdb, chain, d / 'receptor.pdb')
        ref = fetch(f'https://models.rcsb.org/v1/{code}/ligand?label_comp_id={lig}&auth_asym_id={chain}&encoding=sdf',
                    d / 'ligante_cristal.sdf')
        xyz = Chem.MolFromMolFile(str(ref), removeHs=False, sanitize=False).GetConformer().GetPositions()
        centro = xyz.mean(0)
        tamanho = np.ceil(np.ptp(xyz, axis=0) + 2 * folga)
        smi_pdb, smi = smiles_ph74(lig)
        (d / 'ligante.smi').write_text(f'{smi} {lig}\n')
        n_lig, rot, carga = ligand_3d(smi, d / 'ligante_3d.sdf', d / 'ligante_3d.mol2')
        caixa = {'centro_A': [round(float(v), 2) for v in centro], 'tamanho_A': [int(v) for v in tamanho]}
        (d / 'caixa.json').write_text(json.dumps(caixa, indent=2) + '\n')
        resumo[nome] = {'pdb': code, 'cadeia': chain, 'ligante': lig, 'atomos_receptor': n_at,
                        'atomos_ligante_com_H': n_lig, 'ligacoes_rotaveis': rot, 'carga_pH74': carga,
                        'smiles_pH74': smi, **caixa}
        print(nome, json.dumps(resumo[nome], ensure_ascii=False))
    (DADOS / 'resumo.json').write_text(json.dumps(resumo, indent=2, ensure_ascii=False) + '\n')


if __name__ == '__main__':
    main()
