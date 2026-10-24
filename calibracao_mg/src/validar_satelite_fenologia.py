"""Validação independente da fenologia calibrada com a fenologia MODIS (GMR desconhecido no pixel).

Para cada ponto × safra com satélite: a semeadura é estimada pelo início do verde (SOS, 20% da amplitude) menos uma
defasagem L_sem; o emulador, com os coeficientes θ(GMR) e um GMR regional, prevê a duração emergência→R7, que é
comparada com a duração observada SOS→queda de 50% do NDVI (EOS50) mais uma diferença constante (estimada por
mínimos quadrados, porque as definições não coincidem). Métricas: RMSE, correlação total, espacial (médias por
ponto) e interanual (anomalias dentro de cada ponto).

GMR regional: Battisti & Sentelhas (2019) por faixa de latitude (5.8 / 6.8 / 7.8), ou o centro do ciclo médio de
Kaster & Farias (2011) por macrorregião.
"""
from pathlib import Path

import numpy as np
import pandas as pd

import bs2019 as bs
import emulador_fenologia as ef

RAIZ = Path(__file__).resolve().parents[1]
L_SEMEADURA = 15          # dias entre a semeadura e o SOS (emergência ~5 d + verde visível ~10 d), SUPOSIÇÃO


def gmr_kaster(lat, lon):
    """Centro da faixa de GMR do ciclo médio por macrorregião sojícola (aproximação geográfica)."""
    if lat <= -27:
        return 6.9            # MR 1 Sul: 6.4–7.4
    if lat <= -22:
        return 7.2            # MR 2 Centro-Sul: 6.8–7.6
    if lat <= -15 and lon > -50:
        return 7.9            # MR 3 Sudeste: 7.6–8.2
    if lat <= -10:
        return 8.2            # MR 4 Centro-Oeste: 7.9–8.5
    return 9.0                # MR 5 Nordeste/Norte: 8.7–9.3


def tabela_satelite(locais_sat):
    """Satélite nos locais de treino tipo 'satelite' (o emulador conhece o clima deles)."""
    sat = pd.read_csv(RAIZ / "resultados" / "fenologia_satelite" / "fenologia_modis_pontos.csv").dropna(subset=["sos_p50", "eos50_p50"])
    sat = sat.merge(locais_sat.reset_index()[["idx_local", "lon", "lat", "elev"]], on=["lon", "lat"])
    ref = pd.to_datetime((sat.ano_colheita - 1).astype(str) + "-07-01")
    sos = ref + pd.to_timedelta(sat.sos_p50, unit="D")
    sem = sos - pd.Timedelta(days=L_SEMEADURA)
    sat["safra"] = sat.ano_colheita - 1
    sat["doy_sem"] = (sem - pd.to_datetime(sat.safra.astype(str) + "-09-10")).dt.days.clip(0, 171)
    sat["dur_obs"] = sat.eos50_p50 - sat.sos_p50
    return sat


def prever_duracao(sat, locais, coef_fun, regra_gmr="bs"):
    X = ef.variaveis_ambiente(sat[["idx_local", "safra", "doy_sem"]], locais)
    X["doy_sem"] = sat.doy_sem.values
    gmr = np.array([bs.mg_por_latitude(la) if regra_gmr == "bs" else gmr_kaster(la, lo) for la, lo in zip(sat.lat, sat.lon)])
    X["gmr"], X["mg_int"] = gmr, np.clip(np.round(gmr), 5, 9)
    for k, v in coef_fun(gmr).items():
        X[k] = v
    X = X[ef.colunas()]
    import joblib
    m = {s: joblib.load(ef.DIR / f"modelo_{s}.joblib") for s in ef.SAIDAS}
    return m["r7"].predict(X) - m["e"].predict(X), gmr


def metricas(sat, dur_pred):
    d = sat.assign(pred=dur_pred)
    offset = (d.dur_obs - d.pred).mean()
    d["pred_aj"] = d.pred + offset
    e = d.pred_aj - d.dur_obs
    esp = d.groupby("idx_local")[["dur_obs", "pred_aj"]].mean()
    an = d.assign(ao=d.dur_obs - d.groupby("idx_local").dur_obs.transform("mean"),
                  ap=d.pred_aj - d.groupby("idx_local").pred_aj.transform("mean"))
    return dict(n=len(d), offset=offset, rmse=float(np.sqrt((e ** 2).mean())), r_total=d.dur_obs.corr(d.pred_aj),
                r_espacial=esp.dur_obs.corr(esp.pred_aj), r_interanual=an.ao.corr(an.ap),
                tendencia_obs=10 * np.polyfit(d.safra, d.dur_obs, 1)[0], tendencia_pred=10 * np.polyfit(d.safra, d.pred_aj, 1)[0]), d
