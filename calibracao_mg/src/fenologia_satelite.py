"""Fenologia da soja por satélite nos pontos de clima: MODIS NDVI (MOD13Q1 + MYD13Q1, 250 m) sobre pixels de soja
do MapBiomas.

Em cada ponto e safra (julho a junho), sorteia até N pixels de soja (MapBiomas Cobertura col. 10, classe 39, ano da
colheita) na célula de 1° do ponto e baixa a série de NDVI com qualidade e dia real de aquisição. Por pixel: suavização,
amplitude e datas em que o NDVI cruza 20% da amplitude na subida (início do ciclo), o pico e 50% e 20% na descida.
Por ponto e safra: mediana e quartis entre pixels.
"""
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import ee
import numpy as np
import pandas as pd

import gee

COLECOES = ["MODIS/061/MOD13Q1", "MODIS/061/MYD13Q1"]
DIR = gee.RAIZ / "resultados" / "fenologia_satelite"
AMPLITUDE_MIN = 0.25           # descarta pixels sem ciclo nítido de soja


def amostrar_pixels(lon, lat, ano_colheita, n=150, lado=1.0, semente=1):
    soja = ee.Image(gee.LULC).select(f"classification_{ano_colheita}").eq(gee.CLASSE_SOJA).selfMask()
    reg = ee.Geometry.Rectangle([lon - lado / 2, lat - lado / 2, lon + lado / 2, lat + lado / 2])
    fc = soja.sample(region=reg, scale=250, numPixels=n * 4, seed=semente, geometries=True).limit(n)
    return [f["geometry"]["coordinates"] for f in fc.getInfo()["features"]]


def series_ndvi(coords, ano_colheita, tentativas=4):
    ini, fim = f"{ano_colheita - 1}-07-01", f"{ano_colheita}-07-15"
    col = ee.ImageCollection(COLECOES[0]).merge(ee.ImageCollection(COLECOES[1])).filterDate(ini, fim) \
        .select(["NDVI", "SummaryQA", "DayOfYear"])
    for k in range(tentativas):
        try:
            linhas = col.getRegion(ee.Geometry.MultiPoint(coords), 250).getInfo()
            break
        except Exception:
            if k == tentativas - 1:
                raise
            time.sleep(15 * (k + 1))
    d = pd.DataFrame(linhas[1:], columns=linhas[0]).dropna(subset=["NDVI"])
    d["pixel"] = d.longitude.round(5).astype(str) + "_" + d.latitude.round(5).astype(str)
    comp = pd.to_datetime(d.time, unit="ms")
    # dia real de aquisição: DayOfYear do pixel, no ano do composto (ajusta virada de ano)
    ano = comp.dt.year - ((d.DayOfYear < 20) & (comp.dt.dayofyear > 340)).astype(int) * -1
    d["data"] = pd.to_datetime(ano.astype(str), format="%Y") + pd.to_timedelta(d.DayOfYear - 1, unit="D")
    d["ndvi"] = d.NDVI * 1e-4
    return d[d.SummaryQA <= 1][["pixel", "data", "ndvi"]].sort_values(["pixel", "data"])


def metricas_pixel(s, ano_colheita):
    """Datas de cruzamento de limiares numa série diária interpolada e suavizada (mediana móvel de 15 dias).

    O pico é procurado só na janela da soja (1/nov a 31/mar); o mínimo anterior, nos 120 dias antes do pico; o
    mínimo posterior, nos 100 dias depois. Assim culturas de inverno e a safrinha não contaminam as datas.
    """
    s = s.groupby("data").ndvi.max()
    if len(s) < 15:
        return None
    dia = s.resample("D").mean().interpolate(limit=40).rolling(15, center=True, min_periods=5).median().dropna()
    janela = dia[f"{ano_colheita - 1}-11-01":f"{ano_colheita}-03-31"]
    if janela.empty:
        return None
    pico = janela.idxmax()
    vmax = janela.max()
    antes = dia[pico - pd.Timedelta(days=120):pico]
    depois = dia[pico:pico + pd.Timedelta(days=100)]
    if len(antes) < 30 or len(depois) < 30:
        return None
    t_min_ini, b_ini = antes.idxmin(), antes.min()
    b_fim = depois.min()
    amp = vmax - max(b_ini, b_fim)
    if amp < AMPLITUDE_MIN:
        return None
    subida = dia[t_min_ini:pico]
    sos = subida[subida >= b_ini + 0.2 * (vmax - b_ini)].index.min()
    eos50 = depois[depois <= b_fim + 0.5 * (vmax - b_fim)].index.min()
    eos = depois[depois <= b_fim + 0.2 * (vmax - b_fim)].index.min()
    return dict(sos=sos, pico=pico, eos50=eos50, eos=eos, ndvi_max=vmax, amplitude=amp)


def fenologia_ponto_ano(lon, lat, ano_colheita, n=150):
    coords = amostrar_pixels(lon, lat, ano_colheita, n=n)
    if len(coords) < 10:
        return dict(lon=lon, lat=lat, ano_colheita=ano_colheita, n_pixels=0)   # registra para não repetir
    d = series_ndvi(coords, ano_colheita)
    m = [r for r in (metricas_pixel(g, ano_colheita) for _, g in d.groupby("pixel")) if r]
    if len(m) < 10:
        return dict(lon=lon, lat=lat, ano_colheita=ano_colheita, n_pixels=len(m))
    m = pd.DataFrame(m)
    ref = pd.Timestamp(f"{ano_colheita - 1}-07-01")
    out = dict(lon=lon, lat=lat, ano_colheita=ano_colheita, n_pixels=len(m))
    for c in ["sos", "pico", "eos50", "eos"]:
        dias = (m[c] - ref).dt.days
        out.update({f"{c}_p25": dias.quantile(.25), f"{c}_p50": dias.median(), f"{c}_p75": dias.quantile(.75)})
    out["ciclo_sos_eos50_p50"] = ((m.eos50 - m.sos).dt.days).median()
    return out


def rodar(pontos, anos, n=150, threads=16):
    DIR.mkdir(parents=True, exist_ok=True)
    arq = DIR / "fenologia_modis_pontos.csv"
    feitos = set()
    if arq.exists() and arq.stat().st_size > 0:
        f = pd.read_csv(arq)
        feitos = set(zip(f.lon, f.lat, f.ano_colheita))
    tarefas = [(r.lon, r.lat, a) for r in pontos.itertuples() for a in anos if (r.lon, r.lat, a) not in feitos]
    print(f"{len(tarefas)} tarefas ponto × safra ({len(feitos)} já feitas)", flush=True)
    res, t0 = [], time.time()
    with ThreadPoolExecutor(threads) as ex:
        fut = {ex.submit(fenologia_ponto_ano, *t, n): t for t in tarefas}
        for k, f in enumerate(as_completed(fut), 1):
            try:
                r = f.result()
                if r:
                    res.append(r)
            except Exception as erro:
                print("falha", fut[f], str(erro)[:100], flush=True)
            if (k % 100 == 0 or k == len(tarefas)) and res:
                novo = pd.DataFrame(res)
                if arq.exists() and arq.stat().st_size > 0:
                    novo = pd.concat([pd.read_csv(arq), novo], ignore_index=True)
                novo.to_csv(arq, index=False)
                res = []
                print(f"{k}/{len(tarefas)}  {time.time() - t0:.0f} s", flush=True)
