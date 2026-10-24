"""Emulador de fenologia do CROPGRO-Soybean: (coeficientes, GMR, local, semeadura, clima) → dias até E, R1 e R7.

Entradas de ambiente (não dependem dos coeficientes): latitude, altitude, dia da semeadura, fotoperíodo na semeadura
e aos 30 e 60 dias, e temperatura média e mínima em 0–30, 30–60, 60–90 e 90–130 dias após a semeadura.
Validação: locais inteiros fora do treino (GroupKFold).
"""
from datetime import date, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import GroupKFold

import treino_fenologia as tf

DIR = tf.DIR
COEFS = ["csdl", "ppsen", "em-fl", "fl-sd", "sd-pm", "r1ppo"]
SAIDAS = ["e", "r1", "r5", "r7"]      # r5 (início do enchimento) só existe nos lotes ≥ 4 do treino; cada saída é treinada nas linhas em que existe
JANELAS = [(0, 30), (30, 60), (60, 90), (90, 130)]


def fotoperiodo(lat, doy, altitude_sol=-6.0):
    """Comprimento do dia (h) com crepúsculo civil, como no CROPGRO."""
    dec = np.radians(23.45) * np.sin(2 * np.pi * (284 + doy) / 365)
    phi = np.radians(lat)
    cosw = (np.sin(np.radians(altitude_sol)) - np.sin(phi) * np.sin(dec)) / (np.cos(phi) * np.cos(dec))
    return 24 / np.pi * np.arccos(np.clip(cosw, -1, 1))


def data_semeadura(safra, doy_sem):
    d = date(2001, 9, 10) + timedelta(days=int(doy_sem))
    return date(int(safra) + (d.month < 7), d.month, d.day)


def variaveis_ambiente(linhas, locais):
    """linhas: DataFrame com idx_local, safra, doy_sem. Devolve as variáveis de ambiente (mesmo índice), vetorizado
    com somas acumuladas por local."""
    partes = []
    base = pd.Timestamp("2001-09-10")
    for i, g in linhas.groupby("idx_local"):
        loc = locais.loc[i]
        c = pd.read_csv(tf.arq_clima(loc), index_col=0, parse_dates=True).interpolate(limit=5).ffill().bfill()
        cs_med = np.r_[0, np.cumsum(((c.Tmax + c.Tmin) / 2).values)]
        cs_min = np.r_[0, np.cumsum(c.Tmin.values)]
        md = base + pd.to_timedelta(g.doy_sem.astype(int), unit="D")
        ano = g.safra.astype(int).values + (md.dt.month.values < 7)
        datas = pd.to_datetime({"year": ano, "month": md.dt.month.values, "day": md.dt.day.values})
        k0 = c.index.get_indexer(datas)
        v = pd.DataFrame(index=g.index)
        v["lat"], v["elev"] = loc.lat, loc.elev
        doy = datas.dt.dayofyear.values
        for a in (0, 30, 60):
            v[f"fp_{a}"] = fotoperiodo(loc.lat, (doy + a - 1) % 365 + 1)
        n = len(c)
        for a, b in JANELAS:
            i0, i1 = np.clip(k0 + a, 0, n), np.clip(k0 + b, 1, n)
            dur = np.maximum(i1 - i0, 1)
            v[f"tmed_{a}"] = (cs_med[i1] - cs_med[i0]) / dur
            v[f"tmin_{a}"] = (cs_min[i1] - cs_min[i0]) / dur
        partes.append(v)
    return pd.concat(partes).loc[linhas.index]


def montar_tabela():
    locais = pd.read_csv(DIR / "locais.csv", index_col="idx_local")
    tar = pd.read_csv(DIR / "tarefas.csv")
    s = pd.concat([pd.read_csv(f) for f in sorted((DIR / "treino").glob("parte_*.csv"))], ignore_index=True)
    if "erro" in s:
        s = s[s.erro.isna()].drop(columns="erro")
    if "r1ppo" not in tar:
        tar["r1ppo"] = np.nan
    tar["r1ppo"] = tar.r1ppo.fillna(pd.Series(tf.r1ppo_padrao(tar.gmr.values), index=tar.index))   # lotes 1–5: ecótipo
    d = s.merge(tar[["id", "idx_local", "doy_sem", "gmr"] + COEFS], on="id")
    amb = variaveis_ambiente(d[["idx_local", "safra", "doy_sem"]], locais)
    return pd.concat([d, amb], axis=1)


def colunas(d=None):
    amb = ["lat", "elev", "doy_sem", "fp_0", "fp_30", "fp_60"] + [f"{v}_{a}" for a, _ in JANELAS for v in ("tmed", "tmin")]
    return COEFS + ["gmr", "mg_int"] + amb


def preparar(d):
    d = d.copy()
    d["mg_int"] = np.clip(np.round(d.gmr), 5, 9)
    return d


def modelo():
    return HistGradientBoostingRegressor(max_iter=700, learning_rate=0.07, max_leaf_nodes=63, l2_regularization=1.0, random_state=0)


def validar(d, n_folds=5):
    d = preparar(d)
    X, g = d[colunas()], d.idx_local
    res, preds = [], {}
    for s in SAIDAS:
        y = d[s]
        ok = y.notna()
        p = pd.Series(np.nan, index=d.index)
        for tr, tt in GroupKFold(n_folds).split(X[ok], y[ok], g[ok]):
            itr, itt = X[ok].index[tr], X[ok].index[tt]
            p.loc[itt] = modelo().fit(X.loc[itr], y.loc[itr]).predict(X.loc[itt])
        e = (p - y)[ok]
        res.append(dict(saida=s, rmse=float(np.sqrt((e ** 2).mean())), mae=float(e.abs().mean()), vies=float(e.mean()),
                        r2=float(1 - (e ** 2).sum() / ((y[ok] - y[ok].mean()) ** 2).sum())))
        preds[s] = p
    return pd.DataFrame(res), pd.DataFrame(preds)


def treinar(d):
    d = preparar(d)
    return {s: modelo().fit(d[colunas()][d[s].notna()], d[s].dropna()) for s in SAIDAS}
