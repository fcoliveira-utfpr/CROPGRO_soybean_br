"""Testa as suposições da replicação numa amostra de pontos: ET (FAO-56 × Priestley-Taylor) e água do argiloso (150 × 120 mm).

Uso: .venv/bin/python src/sensibilidade_bs2019.py
"""
import shutil
import sys
import warnings
from multiprocessing import current_process, get_context
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))
import bs2019 as bs  # noqa: E402
import rodar_bs2019 as rb  # noqa: E402
from DSSATTools.run import DSSAT  # noqa: E402

SAIDA = rb.RAIZ / "resultados" / "bs2019" / "sensibilidade.csv"
# argiloso com 120 mm em 0–80 cm: SDUL do horizonte 30–200 cm reduzido para que (DUL−LL)·80 cm = 120 mm
DUL_ARG_120 = round(.284 + (120 - (0.162 * 15 + 0.155 * 15) * 10) / 500, 3)     # = 0,429
CENARIOS = {"base": dict(evapo="F", arg120=False), "ET_PT": dict(evapo="R", arg120=False),
            "arg120": dict(evapo="F", arg120=True), "ET_PT+arg120": dict(evapo="R", arg120=True)}
_D = None


def _ini():
    global _D
    pasta = rb.DIR_RUN / f"s{current_process()._identity[0]:02d}"
    shutil.rmtree(pasta, ignore_errors=True)
    _D = DSSAT(str(pasta))


def rodar(args):
    p, nome = pd.Series(args[0]), args[1]
    cfg = CENARIOS[nome]
    bs.EVAPO = cfg["evapo"]
    bs.SOLOS_S1["argiloso"][2] = (81, 12, 7, .622, DUL_ARG_120 if cfg["arg120"] else .488, .284, 1.30, 53)
    clima = pd.read_csv(rb.arq_clima(p.lon, p.lat), index_col=0, parse_dates=True).interpolate(limit=5).ffill().bfill()
    ws = bs.estacao(clima, p.lat, p.lon, elev=p.elev)
    solos = ["argiloso"] if cfg["arg120"] else list(bs.SOLOS_S1)
    out = [bs.simular(_D, ws, s, p.lat, p.lon, sem, agua=True) for s in solos for sem in bs.SEMEADURAS]
    return pd.concat(out).assign(cenario_sens=nome)


if __name__ == "__main__":
    print("DUL 30–200 cm do argiloso para 120 mm:", DUL_ARG_120)
    pts = pd.read_csv(rb.RAIZ / "dados" / "processados" / "pontos_bs2019_digitalizados.csv")
    amostra = pts[pts.zona_bs.isin([1, 3, 5, 6, 7, 8])].groupby("zona_bs", group_keys=False).apply(
        lambda g: g.sample(min(4, len(g)), random_state=1))
    tarefas = [(r, c) for r in amostra.to_dict("records") for c in CENARIOS]
    with get_context("fork").Pool(4, initializer=_ini) as pool:
        res = pd.concat(pool.map(rodar, tarefas))
    res = res.merge(pts[["lon", "lat", "zona_bs"]], on=["lon", "lat"])
    res.to_csv(SAIDA, index=False)
    t = res.groupby(["zona_bs", "solo", "cenario_sens"]).HWAM.mean().unstack().round(0)
    print(t)
