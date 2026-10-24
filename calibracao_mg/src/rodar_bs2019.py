"""Roda a fatorial de Battisti & Sentelhas (2019): 462 pontos × 3 solos × 6 semeaduras × (potencial, atingível).

Cada ponto gera 36 rodadas de 33 safras, salvas em resultados/bs2019/sim/<lon>_<lat>.csv (checkpoint:
pontos já salvos são pulados). Só roda pontos cujo clima já foi baixado.

Uso: .venv/bin/python src/rodar_bs2019.py [n_processos]
"""
import os
import shutil
import sys
import time
import warnings
from multiprocessing import current_process, get_context
from pathlib import Path

import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))
import bs2019 as bs  # noqa: E402
from DSSATTools.run import DSSAT  # noqa: E402

RAIZ = Path(__file__).resolve().parents[1]
DIR_SIM = RAIZ / "resultados" / "bs2019" / "sim"
DIR_SIM_V1 = RAIZ / "resultados" / "bs2019" / "sim_v1_fao56_arg150"
DIR_RUN = RAIZ / ".run"                     # caminho curto: o DSSAT (Fortran) não aceita caminhos longos
_DSSAT = None


def _iniciar():
    global _DSSAT
    pasta = DIR_RUN / f"p{current_process()._identity[0]:02d}"
    shutil.rmtree(pasta, ignore_errors=True)
    _DSSAT = DSSAT(str(pasta))


def arq_clima(lon, lat):
    return RAIZ / "dados" / "brutos" / "xavier" / f"xavier_no_{lon:+.2f}_{lat:+.2f}_1980_2013.csv.gz"


def rodar_ponto(p):
    p = pd.Series(p)
    saida = DIR_SIM / f"{p.lon:+.1f}_{p.lat:+.1f}.csv"
    if saida.exists():
        return "pulado"
    clima = pd.read_csv(arq_clima(p.lon, p.lat), index_col=0, parse_dates=True)
    clima = clima.interpolate(limit=5).ffill().bfill()
    ws = bs.estacao(clima, p.lat, p.lon, elev=p.elev)
    # a potencial não depende do balanço hídrico nem do solo: reaproveita a da v1 quando existir
    v1 = DIR_SIM_V1 / saida.name
    pot = pd.read_csv(v1).query("cenario == 'potencial'") if v1.exists() else None
    cenarios = (True,) if pot is not None else (True, False)
    partes = [bs.simular(_DSSAT, ws, solo, p.lat, p.lon, sem, agua=agua)
              for solo in bs.SOLOS_S1 for sem in bs.SEMEADURAS for agua in cenarios]
    pd.concat(partes + ([pot] if pot is not None else [])).to_csv(saida, index=False)
    return "ok"


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else os.cpu_count()
    DIR_SIM.mkdir(parents=True, exist_ok=True)
    DIR_RUN.mkdir(exist_ok=True)
    pontos = pd.read_csv(RAIZ / "dados" / "processados" / "pontos_bs2019_digitalizados.csv")
    pontos = pontos[[arq_clima(r.lon, r.lat).exists() for r in pontos.itertuples()]]
    print(f"{len(pontos)} pontos com clima; {n} processos")
    t0 = time.time()
    with get_context("fork").Pool(n, initializer=_iniciar) as pool:
        for k, st in enumerate(pool.imap_unordered(rodar_ponto, pontos.to_dict("records"), chunksize=1), 1):
            if k % 10 == 0 or k == len(pontos):
                print(f"{k}/{len(pontos)}  {time.time() - t0:.0f} s", flush=True)
