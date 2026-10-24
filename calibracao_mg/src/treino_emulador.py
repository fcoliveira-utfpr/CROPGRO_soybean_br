"""Gera o conjunto de treino do emulador do CROPGRO-Soybean (configuração v2 de Battisti & Sentelhas, 2019).

- 10 parâmetros selecionados pelo Morris, em hipercubo latino (faixas do Morris); demais fixos no artigo.
- 40 pontos: os 4 de maior área de soja em cada zona do artigo.
- Em cada rodada, solo (um dos 3 do artigo) e semeadura (15/out, 15/nov, 15/dez) sorteados.
- Saída por safra (33 por rodada): produtividade atingível, dias até R1 e R7.

Checkpoint: resultados/emulador/treino/parte_<n>.csv; ao reiniciar, as tarefas já salvas são puladas.

Uso: .venv/bin/python src/treino_emulador.py [n_amostras] [n_processos]
"""
import os
import shutil
import sys
import time
import warnings
from multiprocessing import current_process, get_context
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import qmc

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))
import bs2019 as bs  # noqa: E402
import morris as mo  # noqa: E402
import rodar_bs2019 as rb  # noqa: E402
from DSSATTools.run import DSSAT  # noqa: E402

RAIZ = rb.RAIZ
DIR = RAIZ / "resultados" / "emulador"
DIR_TREINO = DIR / "treino"
SELECIONADOS = ["csdl", "ppsen", "em-fl", "fl-sd", "sd-pm", "lfmax", "sfdur", "sizlf", "fator_ad", "prof_raiz"]
PARAMETROS = [p for p in mo.PARAMETROS if p[0] in SELECIONADOS]
SOLOS = list(bs.SOLOS_S1)
SEMEADURAS = ["15-out", "15-nov", "15-dez"]
SEMENTE = 2026
_D, _WS = None, {}


def pontos_treino(n_por_zona=4):
    pts = pd.read_csv(RAIZ / "dados" / "processados" / "pontos_bs2019_digitalizados.csv")
    fr = pd.read_csv(RAIZ / "dados" / "processados" / "fracoes_solo_mapbiomas_soja2016.csv")
    d = pts.merge(fr[["lon", "lat", "soja_ha"]], on=["lon", "lat"]).dropna(subset=["zona_bs"])
    return d.sort_values("soja_ha", ascending=False).groupby("zona_bs").head(n_por_zona).reset_index(drop=True)


def amostras(n):
    lhs = qmc.LatinHypercube(d=len(PARAMETROS), seed=SEMENTE).random(n)
    lo, hi = np.array([p[1] for p in PARAMETROS]), np.array([p[2] for p in PARAMETROS])
    return pd.DataFrame(qmc.scale(lhs, lo, hi), columns=[p[0] for p in PARAMETROS])


def tarefas(n_amostras):
    X, pts = amostras(n_amostras), pontos_treino()
    rng = np.random.default_rng(SEMENTE)
    lista = []
    for j, p in pts.iterrows():
        for i in range(len(X)):
            lista.append(dict(id=j * len(X) + i, amostra=i, ponto=j, solo=SOLOS[rng.integers(3)],
                              semeadura=SEMEADURAS[rng.integers(3)], **p[["lon", "lat", "elev", "zona_bs"]].to_dict(),
                              **X.iloc[i].to_dict()))
    return X, pts, lista


def _ini():
    global _D
    pasta = rb.DIR_RUN / f"em{current_process()._identity[0]:02d}"
    shutil.rmtree(pasta, ignore_errors=True)
    _D = DSSAT(str(pasta))


def _estacao(t):
    k = (t["lon"], t["lat"])
    if k not in _WS:
        clima = pd.read_csv(rb.arq_clima(*k), index_col=0, parse_dates=True).interpolate(limit=5).ffill().bfill()
        ws = bs.estacao(clima, t["lat"], t["lon"], elev=t["elev"])
        texto = ws._write_wth()                  # o arquivo de clima é igual em todas as rodadas do ponto:
        ws._write_wth = lambda texto=texto: texto  # gera uma vez e reaproveita (economiza ~1,3 s por rodada)
        _WS[k] = ws
    return _WS[k]


def rodar(t):
    mg = bs.mg_por_latitude(t["lat"])
    cam = mo.camadas(t["solo"], t["fator_ad"], t["prof_raiz"])
    valores = {p[0]: t[p[0]] for p in PARAMETROS}

    def cultivar(_mg):
        c = mo._CULTIVAR_ORIG(_mg)
        for nome, _, _, onde in PARAMETROS:
            if onde == "cul":
                c[nome] = float(valores[nome])
        return c

    bs.camadas_solo = lambda _n: cam
    bs.cultivar = cultivar
    try:
        s = bs.simular(_D, _estacao(t), t["solo"], t["lat"], t["lon"], t["semeadura"], agua=True, mg=mg)
    except Exception as erro:
        return [dict(id=t["id"], erro=str(erro)[:80])]
    pdat = mo.bs_data(s.PDAT)
    return pd.DataFrame({"id": t["id"], "safra": s.safra, "hwam": s.HWAM, "cwam": s.CWAM,
                         "r1": (mo.bs_data(s.ADAT) - pdat).dt.days, "r7": (mo.bs_data(s.MDAT) - pdat).dt.days,
                         "prcp": s.PRCP, "epcp": s.EPCP}).to_dict("records")


def ja_feitas():
    feitos = set()
    for f in DIR_TREINO.glob("parte_*.csv"):
        feitos |= set(pd.read_csv(f, usecols=["id"]).id.unique())
    return feitos


if __name__ == "__main__":
    n_am = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    n = int(sys.argv[2]) if len(sys.argv) > 2 else os.cpu_count()
    DIR_TREINO.mkdir(parents=True, exist_ok=True)
    X, pts, lista = tarefas(n_am)
    X.to_csv(DIR / "amostras_lhs.csv", index_label="amostra")
    pts.to_csv(DIR / "pontos_treino.csv", index_label="ponto")
    pd.DataFrame(lista).drop(columns=[p[0] for p in PARAMETROS]).to_csv(DIR / "tarefas.csv", index=False)
    feitos = ja_feitas()
    pend = [t for t in lista if t["id"] not in feitos]
    # intercala os pontos: um resultado parcial já cobre os 40 pontos (o clima de cada ponto fica em cache)
    pend.sort(key=lambda t: (t["amostra"], t["ponto"]))
    parte = len(list(DIR_TREINO.glob("parte_*.csv")))
    print(f"{len(lista)} tarefas | {len(feitos)} já feitas | {len(pend)} pendentes | {n} processos", flush=True)
    t0, buf = time.time(), []
    with get_context("fork").Pool(n, initializer=_ini) as pool:
        for k, out in enumerate(pool.imap_unordered(rodar, pend, chunksize=8), 1):
            buf += out
            if k % 1000 == 0 or k == len(pend):
                pd.DataFrame(buf).to_csv(DIR_TREINO / f"parte_{parte:04d}.csv", index=False)
                parte, buf = parte + 1, []
                taxa = k / (time.time() - t0)
                print(f"{k}/{len(pend)}  {time.time() - t0:.0f} s  faltam ~{(len(pend) - k) / taxa / 3600:.1f} h", flush=True)
