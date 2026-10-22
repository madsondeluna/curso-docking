"""Roda a pratica no SwissDock (AutoDock Vina) pela API de linha de comando.

Mesmas entradas e caixa que o aluno usa na interface web: SMILES em pH 7,4,
receptor de dados/<sistema>/receptor.pdb e caixa de dados/<sistema>/caixa.json.
Guarda o log e os resultados em runs/swissdock/<sistema>/.

Uso: python tools/swissdock_cli.py eif4e mpro [--metodo vina|ac] [--exhaust N] [--smiles S --tag T]
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import time
import zipfile
from pathlib import Path
from urllib.parse import quote

ROOT = Path(__file__).resolve().parents[1]
API = 'https://swissdock.ch:8443'


def curl(*args):
    return subprocess.run(['curl', '-s', '-g', *args], capture_output=True, text=True, timeout=600).stdout


def run(nome, exhaust, metodo='vina', smiles=None, tag=None):
    d = ROOT / 'dados' / nome
    out = ROOT / 'runs' / 'swissdock' / (tag or nome)
    out.mkdir(parents=True, exist_ok=True)
    log = []
    smi = smiles or (d / 'ligante.smi').read_text().split()[0]
    vina = '&Vina' if metodo == 'vina' else ''
    r = curl(f'{API}/preplig?mySMILES={quote(smi, safe="")}&pH=7.4{vina}'); log.append(r)
    sess = re.search(r'sessionNumber=(\d+)', r).group(1)
    r = curl('-F', f'myTarget=@{d / "receptor.pdb"}', f'{API}/preptarget?sessionNumber={sess}'); log.append(r)
    box = json.loads((d / 'caixa.json').read_text())
    c = '_'.join(str(v) for v in box['centro_A']); s = '_'.join(str(v) for v in box['tamanho_A'])
    extra = '' if metodo == 'vina' else '&cavity=70&ric=2'      # padroes do Attracting Cavities
    r = curl(f'{API}/setparameters?sessionNumber={sess}&exhaust={exhaust}{extra}&boxCenter={c}&boxSize={s}'); log.append(r)
    r = curl(f'{API}/startdock?sessionNumber={sess}'); log.append(r)
    t0 = time.time()
    while True:
        time.sleep(30)
        r = curl(f'{API}/checkstatus?sessionNumber={sess}')
        if 'calculation is finished' in r.lower():
            break
        if 'error' in r.lower():
            raise RuntimeError(r)
        if time.time() - t0 > 3600:
            raise TimeoutError(r)
    log.append(r)
    zf = out / 'results.zip'
    curl(f'{API}/retrievesession?sessionNumber={sess}', '-o', str(zf))
    with zipfile.ZipFile(zf) as z:
        z.extractall(out)
    (out / 'log.txt').write_text(f'sessao {sess}, metodo {metodo}, exhaust {exhaust}, {time.time() - t0:.0f} s\n\n' + '\n\n'.join(log))
    print(nome, 'sessao', sess, f'{time.time() - t0:.0f} s')


def main():
    ap = argparse.ArgumentParser(description='Docking da pratica no SwissDock (Vina)')
    ap.add_argument('sistemas', nargs='+')
    ap.add_argument('--exhaust', type=int, default=None, help='Vina: 8; Attracting Cavities: 90')
    ap.add_argument('--metodo', choices=('vina', 'ac'), default='vina')
    ap.add_argument('--smiles', help='outro ligante no lugar de dados/<sistema>/ligante.smi')
    ap.add_argument('--tag', help='subpasta de saida em runs/swissdock')
    a = ap.parse_args()
    exhaust = a.exhaust or (8 if a.metodo == 'vina' else 90)
    for nome in a.sistemas:
        run(nome, exhaust, a.metodo, a.smiles, a.tag)


if __name__ == '__main__':
    main()
