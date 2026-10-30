"""Etapa 3: indicadores por safra a partir das simulações (resultados/sim/<id>.parquet) e do clima diário da célula.

Por safra × célula × solo × GMR × decêndio (rodadas com água):
  ry          Yw/Yp da mesma safra (Yp: rodada sem estresse hídrico, solo principal);
  seca_veg/flor/ench  estresse hídrico no crescimento por fase (0–1, OVERVIEW);
  calor_flor  dias com Tmax > 34 °C entre R1 e R5 (OVERVIEW); calor36: dias com Tmax ≥ 36 °C entre R1 e R5 (clima);
  geada       Tmin ≤ 2 °C (ar, Xavier) entre a emergência e o R7 corrigido;
  chuva_colh  chuva (mm) nos 10 dias após o R8 corrigido;
  r8, r8c     R8 simulado (R7 + 11 d) e corrigido pela data de semeadura (fenologia_correcao, k central);
  safrinha    colheita (R8 corrigido) + 3 d até a data-limite de semeadura do milho safrinha da célula.
Data-limite do milho safrinha (por célula, do clima 1980–2025): última data de semeadura entre 1/jan e 15/mar com
P(Tmin ≤ 2 °C entre 50 e 120 dias após a semeadura) ≤ 20% (geada no florescimento/enchimento do milho).

Saída: resultados/indicadores_safra.parquet e resultados/safrinha_limite.csv.
Uso: ../calibracao_mg/.venv/bin/python src/etapa3_indicadores.py
"""
import numpy as np
import pandas as pd

import config as cfg
import fenologia_correcao as fc
from etapa2_simular import ARQ_CEL, DIR_SIM, arq_clima

DELTA8 = 11
TMIN_GEADA = 2.0
TMAX_CALOR = 36.0
CHUVA_DIAS = 10
P_GEADA_SAFRINHA = 0.2
JANELA_GEADA_MILHO = (50, 120)
ULTIMA_SAFRINHA = (3, 15)


def datas(yyddd):
    return pd.to_datetime(yyddd.astype("Int64").astype(str), format="%Y%j", errors="coerce")


def limite_safrinha(clima):
    """Dia do ano (contado a partir de 1/jan) da última semeadura de milho safrinha com risco de geada ≤ 20%."""
    tmin = clima.Tmin
    anos = range(1981, 2025)
    ultimo = pd.Timestamp(2001, *ULTIMA_SAFRINHA).dayofyear
    lim = 1
    for doy in range(1, ultimo + 1):
        geada = []
        for a in anos:
            s = pd.Timestamp(a, 1, 1) + pd.Timedelta(days=doy - 1)
            w = tmin[s + pd.Timedelta(days=JANELA_GEADA_MILHO[0]):s + pd.Timedelta(days=JANELA_GEADA_MILHO[1])]
            geada.append(w.min() <= TMIN_GEADA)
        if np.mean(geada) <= P_GEADA_SAFRINHA:
            lim = doy
    return lim


def janela(serie, ini, fim, func):
    """Aplica func à série diária entre ini e fim (inclusive) para cada par de datas."""
    v = serie.to_numpy()
    idx0 = serie.index[0]
    i = ((ini - idx0).dt.days).to_numpy()
    j = ((fim - idx0).dt.days).to_numpy()
    out = np.full(len(i), np.nan)
    ok = ~(np.isnan(i.astype(float)) | np.isnan(j.astype(float)))
    for k in np.flatnonzero(ok):
        a, b = int(i[k]), int(j[k])
        if 0 <= a <= b < len(v):
            out[k] = func(v[a:b + 1])
    return out


def indicadores_celula(cel, k):
    s = pd.read_parquet(DIR_SIM / f"{cel.id_celula}.parquet")
    s["ano"] = s.PDAT // 1000
    yp = s[~s.agua][["gmr", "semeadura", "ano", "HWAM", "MDAT"]].rename(columns={"HWAM": "yp", "MDAT": "mdat_p"})
    w = s[s.agua].merge(yp, on=["gmr", "semeadura", "ano"], how="left")
    w["ry"] = np.where(w.yp > 0, w.HWAM / w.yp, np.nan)
    clima = pd.read_csv(arq_clima(cel.lon, cel.lat), index_col=0, parse_dates=True).interpolate(limit=5).ffill().bfill()
    pdat, edat, r7 = datas(w.PDAT), datas(w.EDAT), datas(w.MDAT)
    r7c = fc.corrigir(r7, pdat, k)
    w["r1"] = datas(w.ADAT)
    r5 = w.r1 + pd.to_timedelta(w.dias_flor, unit="D")
    w["calor36"] = janela(clima.Tmax, w.r1, r5, lambda x: float(np.sum(x >= TMAX_CALOR)))
    w["r8"] = r7 + pd.Timedelta(days=DELTA8)
    w["r8c"] = r7c + pd.Timedelta(days=DELTA8)
    w["geada"] = janela(clima.Tmin, edat, r7c, lambda x: float(x.min() <= TMIN_GEADA))
    w["tmin_min"] = janela(clima.Tmin, edat, r7c, np.min)
    w["chuva_colh"] = janela(clima.pr, w.r8c, w.r8c + pd.Timedelta(days=CHUVA_DIAS - 1), np.sum)
    w["chuva_colh_sc"] = janela(clima.pr, w.r8, w.r8 + pd.Timedelta(days=CHUVA_DIAS - 1), np.sum)
    lim = limite_safrinha(clima)
    ano_col = w.r8c.dt.year
    dlim = pd.to_datetime(ano_col.astype("Int64").astype(str) + "-01-01", errors="coerce") + pd.to_timedelta(lim - 1, unit="D")
    w["safrinha"] = ((w.r8c + pd.Timedelta(days=3)) <= dlim).astype(float)
    w["safrinha_sc"] = ((w.r8 + pd.Timedelta(days=3)) <= dlim).astype(float)
    w["ciclo_r8c"] = (w.r8c - pdat).dt.days
    w = w.rename(columns={"agua_cr_veg": "seca_veg", "agua_cr_flor": "seca_flor", "agua_cr_ench": "seca_ench",
                          "n_tmax34_flor": "calor_flor", "HWAM": "yw"})
    cols = ["id_celula", "solo", "gmr", "semeadura", "ano", "yw", "yp", "ry", "seca_veg", "seca_flor", "seca_ench",
            "calor_flor", "calor36", "geada", "tmin_min", "chuva_colh", "chuva_colh_sc", "r1", "r8", "r8c", "ciclo_r8c",
            "safrinha", "safrinha_sc", "PRCP", "ETCP"]
    return w[cols], dict(id_celula=cel.id_celula, safrinha_limite_doy=lim,
                         safrinha_limite=(pd.Timestamp(2001, 1, 1) + pd.Timedelta(days=lim - 1)).strftime("%d/%m"))


if __name__ == "__main__":
    k = fc.carregar()["k"]
    c = pd.read_csv(ARQ_CEL)
    c = c[[(DIR_SIM / f"{i}.parquet").exists() for i in c.id_celula]]
    partes, lims = [], []
    for n, cel in enumerate(c.itertuples(), 1):
        d, l = indicadores_celula(cel, k)
        partes.append(d)
        lims.append(l)
        if n % 20 == 0:
            print(f"{n}/{len(c)}", flush=True)
    pd.concat(partes, ignore_index=True).to_parquet(cfg.DIR_RES / "indicadores_safra.parquet", index=False)
    pd.DataFrame(lims).to_csv(cfg.DIR_RES / "safrinha_limite.csv", index=False)
    print(f"{len(c)} células; {sum(len(p) for p in partes)} linhas")
