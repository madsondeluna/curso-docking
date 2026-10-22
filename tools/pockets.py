"""Previsao de cavidades no PrankWeb (P2Rank) e comparacao com os sitios reais.

Para cada estrutura: submete o codigo PDB ao PrankWeb, baixa as cavidades
previstas e mede a distancia do centro de cada uma ao atomo mais proximo de
cada ligante do cristal. Assim da para dizer quais cavidades sao sitios de ligacao
conhecidos e quais sao so buracos na superficie.

Uso: python tools/pockets.py 8HBK 1IPC 5MO4 6OIM 4OBE
"""
from __future__ import annotations

import argparse
import json
import time
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
API = 'https://prankweb.cz/api/v2/prediction/v4-conservation-hmm'
CACHE = ROOT / 'runs' / 'pdb'
# ligantes ignorados: cristalizacao, tampao, ions e cofatores comuns
IGNORAR = {'HOH', 'GOL', 'EDO', 'PEG', 'PG4', 'SO4', 'PO4', 'CL', 'NA', 'K', 'MG', 'ZN', 'CA', 'MES', 'TRS',
           'ACT', 'DMS', 'IOD', 'NO3', 'FMT', 'ACE', 'NH2'}


def get(url):
    with urllib.request.urlopen(url) as r:
        return r.read()


def predict(code):
    """Submete (se preciso) e espera a previsao ficar pronta."""
    for _ in range(120):
        info = json.loads(get(f'{API}/{code}'))
        if info['status'] == 'successful':
            return json.loads(get(f'{API}/{code}/public/prediction.json'))
        if info['status'] == 'failed':
            raise RuntimeError(f'{code}: {info}')
        time.sleep(10)
    raise TimeoutError(code)


def ligands(code):
    """Atomos pesados de cada ligante (HETATM, sem agua nem aditivos) do PDB, por cadeia."""
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f'{code}.pdb'
    if not path.exists():
        urllib.request.urlretrieve(f'https://files.rcsb.org/download/{code}.pdb', path)
    out = {}
    for line in path.read_text().splitlines():
        if line.startswith('HETATM') and line[17:20].strip() not in IGNORAR and line[76:78].strip() != 'H':
            key = (line[17:20].strip(), line[21], line[22:26].strip())
            out.setdefault(key, []).append([float(line[30 + 8 * i:38 + 8 * i]) for i in range(3)])
    return {k: np.array(v) for k, v in out.items()}


def main():
    ap = argparse.ArgumentParser(description='Cavidades previstas (PrankWeb) x ligantes do cristal')
    ap.add_argument('codigos', nargs='+')
    a = ap.parse_args()
    out_dir = ROOT / 'resultados' / 'cavidades'
    out_dir.mkdir(parents=True, exist_ok=True)
    for code in a.codigos:
        pred = predict(code)
        ligs = ligands(code)
        rows = []
        for p in pred['pockets']:
            c = np.array([float(x) for x in p['center']])
            dist = {f'{k[0]} {k[1]}{k[2]}': round(float(np.linalg.norm(v - c, axis=1).min()), 1) for k, v in ligs.items()}
            perto = min(dist, key=dist.get) if dist else None
            rows.append({'rank': int(p['rank']), 'score': float(p['score']), 'prob': float(p['probability']),
                         'centro': [round(float(x), 2) for x in c], 'n_residuos': len(p['residues']),
                         'ligante_mais_proximo': perto, 'distancia_A': dist.get(perto)})
        (out_dir / f'{code}.json').write_text(json.dumps({'pdb': code, 'cavidades': rows}, indent=1, ensure_ascii=False) + '\n')
        print(f'{code}: {len(rows)} cavidades, ligantes {sorted(f"{k[0]}" for k in ligs)}')
        for r in rows[:6]:
            print(f"  rank {r['rank']:>2} score {r['score']:>6.2f}  {r['ligante_mais_proximo']} a {r['distancia_A']} Å")


if __name__ == '__main__':
    main()
