"""Análise de sensibilidade de Morris do CROPGRO-Soybean sobre a configuração validada (Battisti & Sentelhas 2019, v2).

Parâmetros: 13 coeficientes de cultivar/ecótipo (faixas: envoltória da literatura brasileira e dos MG 5–10 do
DSSAT, resultados/01_faixas_coeficientes.csv) + fator da água disponível + profundidade de raiz.
Ambientes: 2 pontos por zona do artigo (maior área de soja), solo dominante do MapBiomas, semeaduras de
15/out, 15/nov e 15/dez, 33 safras. Saídas por rodada: produtividade atingível média e CV, R1 e R7 (DAS).

Uso: .venv/bin/python src/morris.py [trajetorias] [n_processos]
"""
import os
import shutil
import sys
import time
import warnings
from datetime import date
from multiprocessing import current_process, get_context
from pathlib import Path

import numpy as np
import pandas as pd
from SALib.sample import morris as amostra_morris

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))
import bs2019 as bs  # noqa: E402
import rodar_bs2019 as rb  # noqa: E402
from DSSATTools.filex import (Field, Planting, SCGeneral, SCManagement, SCMethods, SCOptions,  # noqa: E402
                              SCOutputs, SimulationControls)
from DSSATTools.run import DSSAT  # noqa: E402

RAIZ = rb.RAIZ
DIR = RAIZ / "resultados" / "morris"
SEMEADURAS = ["15-out", "15-nov", "15-dez"]

# (nome, mínimo, máximo, onde): "cul" e "eco" vão para a cultivar; "solo" para o perfil
PARAMETROS = [
    ("csdl", 11.5, 13.0, "cul"), ("ppsen", 0.20, 0.385, "cul"), ("em-fl", 16.5, 30.0, "cul"),
    ("fl-sd", 9.0, 18.5, "cul"), ("sd-pm", 25.0, 38.0, "cul"), ("lfmax", 0.93, 1.43, "cul"),
    ("slavr", 300, 480, "cul"), ("sizlf", 137, 340, "cul"), ("wtpsd", 0.115, 0.21, "cul"),
    ("sfdur", 12.0, 26.0, "cul"), ("sdpdv", 1.7, 2.4, "cul"),
    ("fl-vs", 9.0, 44.0, "eco"), ("r1ppo", 0.30, 0.60, "eco"),
    ("fator_ad", 0.6, 1.4, "solo"), ("prof_raiz", 40, 120, "solo"),
]
PROBLEMA = {"num_vars": len(PARAMETROS), "names": [p[0] for p in PARAMETROS],
            "bounds": [[p[1], p[2]] for p in PARAMETROS]}
PROF_RAIZ_BASE = 70     # profundidade em que o SRGF da Tabela S2 chega a 0,01
_D, _WS = None, {}
_CAMADAS_ORIG, _CULTIVAR_ORIG = bs.camadas_solo, bs.cultivar   # originais, antes das substituições em rodar()


def ambientes():
    """2 pontos por zona do artigo, os de maior área de soja, com o solo dominante do MapBiomas."""
    pts = pd.read_csv(RAIZ / "dados" / "processados" / "pontos_bs2019_digitalizados.csv")
    fr = pd.read_csv(RAIZ / "dados" / "processados" / "fracoes_solo_mapbiomas_soja2016.csv")
    d = pts.merge(fr, on=["lon", "lat"]).dropna(subset=["zona_bs"])
    d["solo"] = d[["f_francoarenoso", "f_argiloarenoso", "f_argiloso"]].idxmax(axis=1).str[2:]
    return d.sort_values("soja_ha", ascending=False).groupby("zona_bs").head(2).sort_values(["zona_bs", "lat"])


def camadas(nome_solo, fator_ad, prof_raiz):
    """Perfil do artigo com a água disponível escalada e o SRGF esticado/encolhido para a profundidade de raiz."""
    c = _CAMADAS_ORIG(nome_solo)
    c["sdul"] = np.minimum(c.slll + fator_ad * (c.sdul - c.slll), c.ssat - 0.01).round(3)
    meio = (np.r_[0, c.slb[:-1]] + c.slb) / 2
    base_meio = (np.r_[0, bs.BASE_CAMADAS[:-1]] + np.array(bs.BASE_CAMADAS)) / 2
    z = meio * PROF_RAIZ_BASE / prof_raiz                      # profundidade equivalente no perfil-base
    c["srgf"] = np.interp(z, base_meio, bs.SRGF_S2, right=0.0).round(3)
    c.loc[meio <= 15, "srgf"] = 1.0
    return c


def _cultivar(mg, valores):
    c = _CULTIVAR_ORIG(mg)
    for (nome, _, _, onde), v in zip(PARAMETROS, valores):
        if onde == "cul":
            c[nome] = float(v)
        elif onde == "eco":
            c["eco#"][nome] = float(v)
    return c


def _ini():
    global _D
    pasta = rb.DIR_RUN / f"mo{current_process()._identity[0]:02d}"
    shutil.rmtree(pasta, ignore_errors=True)
    _D = DSSAT(str(pasta))


def _estacao(p):
    k = (p["lon"], p["lat"])
    if k not in _WS:
        clima = pd.read_csv(rb.arq_clima(*k), index_col=0, parse_dates=True).interpolate(limit=5).ffill().bfill()
        _WS[k] = bs.estacao(clima, p["lat"], p["lon"], elev=p["elev"])
    return _WS[k]


def rodar(tarefa):
    i, valores, p = tarefa
    v = dict(zip(PROBLEMA["names"], valores))
    ws = _estacao(p)
    mg = bs.mg_por_latitude(p["lat"])
    cam = camadas(p["solo"], v["fator_ad"], v["prof_raiz"])
    # o bs2019 monta solo, condições iniciais e cultivar por camadas_solo() e cultivar(): substitui pelas perturbadas
    bs.camadas_solo = lambda _n, cam=cam: cam
    bs.cultivar = lambda _mg, valores=valores: _cultivar(_mg, valores)
    saidas = []
    for sem in SEMEADURAS:
        try:
            s = bs.simular(_D, ws, p["solo"], p["lat"], p["lon"], sem, agua=True, mg=mg)
            hw = s.HWAM
            r1 = (bs_data(s.ADAT) - bs_data(s.PDAT)).dt.days
            r7 = (bs_data(s.MDAT) - bs_data(s.PDAT)).dt.days
            saidas.append(dict(amostra=i, lon=p["lon"], lat=p["lat"], zona=p["zona_bs"], semeadura=sem,
                               yw=hw.mean(), cv=100 * hw.std() / hw.mean(), r1=r1.median(), r7=r7.median()))
        except Exception as erro:
            saidas.append(dict(amostra=i, lon=p["lon"], lat=p["lat"], zona=p["zona_bs"], semeadura=sem, erro=str(erro)[:80]))
    return saidas


def bs_data(col):
    return pd.to_datetime(col.astype("Int64").astype(str), format="%Y%j", errors="coerce")



if __name__ == "__main__":
    r = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    n = int(sys.argv[2]) if len(sys.argv) > 2 else os.cpu_count()
    DIR.mkdir(parents=True, exist_ok=True)
    X = amostra_morris.sample(PROBLEMA, N=r, num_levels=4, seed=42)
    np.save(DIR / "amostras.npy", X)
    amb = ambientes()
    amb.to_csv(DIR / "ambientes.csv", index=False)
    tarefas = [(i, x, p) for p in amb.to_dict("records") for i, x in enumerate(X)]
    print(f"{len(X)} conjuntos × {len(amb)} pontos × {len(SEMEADURAS)} semeaduras = {len(tarefas) * 3} rodadas; {n} processos", flush=True)
    t0, res = time.time(), []
    with get_context("fork").Pool(n, initializer=_ini) as pool:
        for k, out in enumerate(pool.imap_unordered(rodar, tarefas, chunksize=4), 1):
            res += out
            if k % 200 == 0 or k == len(tarefas):
                print(f"{k}/{len(tarefas)}  {time.time() - t0:.0f} s", flush=True)
                pd.DataFrame(res).to_csv(DIR / "saidas.csv", index=False)
    pd.DataFrame(res).to_csv(DIR / "saidas.csv", index=False)
