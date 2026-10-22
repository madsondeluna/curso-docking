"""Copia os resultados do SwissDock para resultados/ (versionado) e empacota dados/.

resultados/<sistema>/: poses do Vina (PDBQT e SDF das 3 melhores), energias e
RMSD (poses.json), interacoes (interacoes.json) e os parametros da sessao.
dados/<sistema>.zip: receptor, ligantes, SMILES e caixa, para baixar com um clique.

Uso: python tools/export_results.py
"""
from __future__ import annotations

import shutil
import zipfile
from pathlib import Path

from openbabel import pybel
from openbabel import openbabel as ob

ob.obErrorLog.SetOutputLevel(0)
ROOT = Path(__file__).resolve().parents[1]
KEEP = ('vina_dock.pdbqt', 'vina_dock.pdb', 'poses.json', 'interacoes.json', 'parameters')


def complexo(run, out, pose=1):
    """Receptor + a pose pedida num unico PDB, com o ligante como HETATM (formato que o PLIP espera)."""
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from interactions import complex_pdb, ligand_lines
    mol = list(pybel.readfile('pdbqt', str(run / 'vina_dock.pdbqt')))[pose - 1]
    (out / f'complexo_pose{pose}.pdb').write_text(complex_pdb(run / 'receptor.pdb', ligand_lines(mol)))


def main():
    for run in sorted((ROOT / 'runs' / 'swissdock').glob('[!._]*')):
        if not (run / 'vina_dock.pdbqt').exists():
            continue
        out = ROOT / 'resultados' / run.name
        out.mkdir(parents=True, exist_ok=True)
        for f in KEEP:
            if (run / f).exists():
                shutil.copy(run / f, out / f)
        complexo(run, out)
        sdf = pybel.Outputfile('sdf', str(out / 'poses_top3.sdf'), overwrite=True)
        for k, m in enumerate(pybel.readfile('pdbqt', str(run / 'vina_dock.pdbqt'))):
            if k == 3:
                break
            m.title = f'{run.name} pose {k + 1}'
            sdf.write(m)
        sdf.close()
        print('resultados/', run.name)
    for d in sorted((ROOT / 'dados').glob('[!._]*')):
        if d.is_dir():
            with zipfile.ZipFile(ROOT / 'dados' / f'{d.name}.zip', 'w', zipfile.ZIP_DEFLATED) as z:
                for f in sorted(d.glob('[!._]*')):
                    z.write(f, f'{d.name}/{f.name}')
            print('dados/', f'{d.name}.zip')


if __name__ == '__main__':
    main()
