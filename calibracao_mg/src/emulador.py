"""Emulador do CROPGRO-Soybean: aprende (parâmetros, ambiente, clima da safra) → (produtividade, R1, R7).

Entradas: os 10 parâmetros selecionados pelo Morris, MG, solo, dia da semeadura, latitude e resumos do clima
do Xavier em 3 janelas a partir da semeadura (0–30, 30–60 e 60–110 dias): chuva, ETo, Tmax, Tmin e Rs.
Modelo: gradient boosting por histogramas (scikit-learn), um por saída.
Validação: pontos inteiros fora do treino (GroupKFold por ponto).
"""
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import GroupKFold

import bs2019 as bs
import rodar_bs2019 as rb
import treino_emulador as te

RAIZ = rb.RAIZ
JANELAS = [(0, 30), (30, 60), (60, 110)]
SAIDAS = ["hwam", "r1", "r7"]
CODIGO_SOLO = {"francoarenoso": 0, "argiloarenoso": 1, "argiloso": 2}


def clima_janelas(lon, lat, safras, semeadura):
    """Resumos climáticos por safra nas janelas após a semeadura."""
    c = pd.read_csv(rb.arq_clima(lon, lat), index_col=0, parse_dates=True).interpolate(limit=5).ffill().bfill()
    mes, dia = bs.SEMEADURAS[semeadura]
    linhas = []
    for s in safras:
        ini = pd.Timestamp(date(int(s) + (mes < 9), mes, dia))
        r = {"safra": s}
        for a, b in JANELAS:
            w = c.loc[ini + pd.Timedelta(days=a): ini + pd.Timedelta(days=b - 1)]
            r.update({f"chuva_{a}": w.pr.sum(), f"eto_{a}": w.ET.sum(), f"tmax_{a}": w.Tmax.mean(),
                      f"tmin_{a}": w.Tmin.mean(), f"rs_{a}": w.Rs.mean()})
        r["deficit_total"] = sum(r[f"chuva_{a}"] - r[f"eto_{a}"] for a, _ in JANELAS)
        linhas.append(r)
    return pd.DataFrame(linhas)


def montar_tabela(dir_treino=te.DIR_TREINO):
    """Junta saídas por safra, parâmetros, ambiente e clima numa tabela de treino."""
    saidas = pd.concat([pd.read_csv(f) for f in sorted(Path(dir_treino).glob("parte_*.csv"))], ignore_index=True)
    if "erro" in saidas:
        saidas = saidas[saidas.erro.isna()].drop(columns="erro")
    tar = pd.read_csv(te.DIR / "tarefas.csv")
    X = pd.read_csv(te.DIR / "amostras_lhs.csv")
    d = saidas.merge(tar, on="id").merge(X, on="amostra")
    clima = []
    for (lon, lat, sem), g in d.groupby(["lon", "lat", "semeadura"]):
        clima.append(clima_janelas(lon, lat, sorted(g.safra.unique()), sem).assign(lon=lon, lat=lat, semeadura=sem))
    d = d.merge(pd.concat(clima), on=["lon", "lat", "semeadura", "safra"])
    d["mg"] = d.lat.map(bs.mg_por_latitude)
    d["solo_cod"] = d.solo.map(CODIGO_SOLO)
    d["doy_semeadura"] = d.semeadura.map(lambda s: date(2001, *bs.SEMEADURAS[s]).timetuple().tm_yday)
    return d


def colunas_entrada(d):
    clima = [c for c in d.columns if any(c.startswith(p) for p in ("chuva_", "eto_", "tmax_", "tmin_", "rs_"))]
    return te.SELECIONADOS + ["mg", "solo_cod", "doy_semeadura", "lat", "deficit_total"] + clima


def treinar(d, saida, **kw):
    X, y = d[colunas_entrada(d)], d[saida]
    ok = y.notna()
    m = HistGradientBoostingRegressor(max_iter=kw.get("max_iter", 600), learning_rate=kw.get("lr", 0.06),
                                      max_leaf_nodes=kw.get("folhas", 63), l2_regularization=1.0, random_state=0)
    return m.fit(X[ok], y[ok])


def validar(d, saida, n_folds=5, **kw):
    """Erro em pontos fora do treino: RMSE, viés e R² por fold e no total."""
    X, y, g = d[colunas_entrada(d)], d[saida], d["ponto"]
    ok = y.notna()
    X, y, g = X[ok], y[ok], g[ok]
    pred = pd.Series(np.nan, index=y.index)
    for tr, te_ in GroupKFold(n_folds).split(X, y, g):
        m = HistGradientBoostingRegressor(max_iter=kw.get("max_iter", 600), learning_rate=kw.get("lr", 0.06),
                                          max_leaf_nodes=kw.get("folhas", 63), l2_regularization=1.0, random_state=0)
        m.fit(X.iloc[tr], y.iloc[tr])
        pred.iloc[te_] = m.predict(X.iloc[te_])
    erro = pred - y
    return pred, dict(saida=saida, n=len(y), rmse=float(np.sqrt((erro ** 2).mean())), vies=float(erro.mean()),
                      r2=float(1 - (erro ** 2).sum() / ((y - y.mean()) ** 2).sum()))


# ── Segunda etapa: clima e balanço hídrico nas fases da planta ───────────────────────────────
DIAS_ANTES = 60          # o balde começa 60 dias antes da semeadura (pousio, Kc = 0,3)
DIAS_DEPOIS = 200        # cobre o ciclo mais longo
FRACAO_R5 = 0.4          # R5 ≈ R1 + 0,4 (R7 − R1)


def agua_raiz(solo, fator_ad, prof_raiz):
    """Água disponível (mm) nas camadas alcançadas pela raiz (SRGF > 0,005)."""
    import morris as mo
    c = mo.camadas(solo, fator_ad, prof_raiz)
    esp = np.diff(np.r_[0, c.slb]) * 10
    return float(((c.sdul - c.slll) * esp * (c.srgf > 0.005)).sum())


def _kc(dia, r1, r5, r7):
    """Coeficiente de cultura: 0,3 antes da semeadura, 0,5 → 1,05 até R1, 1,05 até R5, cai a 0,6 em R7."""
    kc = np.where(dia < 0, 0.3, 0.5 + 0.55 * np.clip(dia / np.maximum(r1, 1), 0, 1))
    kc = np.where(dia > r5, 1.05 - 0.45 * np.clip((dia - r5) / np.maximum(r7 - r5, 1), 0, 1), kc)
    return np.where(dia > r7, 0.0, kc)


def variaveis_fenologicas(d, col_r1="r1_pred", col_r7="r7_pred"):
    """Para cada safra: chuva, ETo, temperatura, radiação e razão ETa/ETc (balde) em S–R1, R1–R5 e R5–R7."""
    if "w_raiz" not in d:
        combos = d[["solo", "fator_ad", "prof_raiz"]].drop_duplicates()
        combos["w_raiz"] = [agua_raiz(*r) for r in combos.itertuples(index=False)]
        d = d.merge(combos, on=["solo", "fator_ad", "prof_raiz"], how="left")
    saida = []
    dias = np.arange(-DIAS_ANTES, DIAS_DEPOIS)
    for (lon, lat, sem), g in d.groupby(["lon", "lat", "semeadura"]):
        c = pd.read_csv(rb.arq_clima(lon, lat), index_col=0, parse_dates=True).interpolate(limit=5).ffill().bfill()
        mes, dia = bs.SEMEADURAS[sem]
        idx0 = {s: c.index.get_indexer([pd.Timestamp(date(int(s) + (mes < 9), mes, dia))])[0] for s in g.safra.unique()}
        pos = np.array([idx0[s] for s in g.safra])[:, None] + dias[None, :]
        pos = np.clip(pos, 0, len(c) - 1)
        P, E = c.pr.values[pos], c.ET.values[pos]
        T, Tn, R = c.Tmax.values[pos], c.Tmin.values[pos], c.Rs.values[pos]
        r1, r7 = g[col_r1].values[:, None], g[col_r7].values[:, None]
        r5 = r1 + FRACAO_R5 * (r7 - r1)
        etc = _kc(dias[None, :], r1, r5, r7) * E
        W = g.w_raiz.values
        S = 0.5 * W.copy()
        eta = np.zeros_like(etc)
        for k in range(len(dias)):
            S = np.minimum(W, S + P[:, k])
            a = etc[:, k] * np.clip(S / (0.5 * np.maximum(W, 1)), 0, 1)   # ETa cai abaixo de 50% da água disponível
            a = np.minimum(a, S)
            S -= a
            eta[:, k] = a
        f = {}
        for nome, (a, b) in {"veg": (0, r1), "flor": (r1, r5), "ench": (r5, r7)}.items():
            m = (dias[None, :] >= a) & (dias[None, :] < b)
            n = np.maximum(m.sum(1), 1)
            f[f"chuva_{nome}"] = (P * m).sum(1)
            f[f"eto_{nome}"] = (E * m).sum(1)
            f[f"tmax_{nome}"] = (T * m).sum(1) / n
            f[f"tmin_{nome}"] = (Tn * m).sum(1) / n
            f[f"rs_{nome}"] = (R * m).sum(1) / n
            f[f"razao_et_{nome}"] = (eta * m).sum(1) / np.maximum((etc * m).sum(1), 0.1)
        saida.append(pd.DataFrame(f, index=g.index))
    return d.join(pd.concat(saida))


def colunas_entrada_fen(d):
    fases = [c for c in d.columns if c.split("_")[-1] in ("veg", "flor", "ench")]
    return te.SELECIONADOS + ["mg", "solo_cod", "doy_semeadura", "lat", "w_raiz", "r1_pred", "r7_pred"] + fases
