"""Etapa 1: grade de 0,25° do Paraná com soja, clima, perfis de solo e AWC da ANA.

1. Células de 0,25° que tocam o PR (limites GAUL); área de soja do MapBiomas (col. 10, ANO_SOJA) dentro do estado;
   ficam as células com mais de MIN_SOJA_HA → dados/processados/celulas.csv (com altitude média SRTM).
2. Perfis de solo MapBiomas Solo + SoilGrids por classe textural (1–3) sob soja, recortados no PR, com a área de cada
   classe → dados/processados/perfis_solo.csv (formato de calibracao_mg/src/solo_mapbiomas.perfis_celulas).
3. AWC da ANA/UFPR (2021; m³/m³) média sob a soja de cada classe textural de cada célula → coluna awc_ana em
   perfis_solo.csv (para a chave opcional de reescalonamento do DUL e para a classe AD do ZARC).
4. Limpeza da textura (passo limpar): em algumas células o MapBiomas não tem dado em parte dos pixels de uma camada
   e a soma areia + silte + argila fica < 100% (argiloarenosos: só 90–100 cm; 17 perfis francoarenosos pequenos: quase
   todo o perfil). Camadas com soma ≥ 50% são normalizadas para 100%; abaixo de 40 cm, camada sem dado copia a de cima;
   perfil com camada sem dado em 0–40 cm é descartado (coluna valido). Acrescenta a AD do ZARC (equação da IN SPA/MAPA
   nº 1/2022, textura de 0–40 cm) e a classe AD1–AD6.
5. Clima do Xavier (BR-DWGD) 1980–2025 no pixel de 0,1° que contém o centro da célula (cache compartilhado em
   calibracao_mg/dados/brutos/xavier/xavier_px_*.csv.gz).

Uso: ../calibracao_mg/.venv/bin/python src/etapa1_grade.py [passo ...]   (padrão: todos; passos: celulas solo awc limpar clima)
"""
import json
import sys

import ee
import numpy as np
import pandas as pd

import config as cfg
import gee
import solo_mapbiomas as sm

LADO = 0.25
ANO_SOJA = 2023
MIN_SOJA_HA = 1000
AWC_ANA = "projects/fcoliveira/assets/AWC_br"
DIR = cfg.DIR_DADOS / "processados"
ARQ_CEL, ARQ_SOLO = DIR / "celulas.csv", DIR / "perfis_solo.csv"


def geometria_pr():
    g = json.load(open(cfg.CALIB / "dados" / "processados" / "brasil_ufs_gaul.geojson"))
    f = next(f for f in g["features"] if f["properties"]["ADM1_NAME"] in ("Paraná", "Parana"))
    return ee.Geometry(f["geometry"])


def _feicoes(df):
    return [ee.Feature(ee.Geometry.Rectangle([x - LADO / 2, y - LADO / 2, x + LADO / 2, y + LADO / 2]),
                       {"lon": float(x), "lat": float(y)}) for x, y in zip(df.lon, df.lat)]


def celulas():
    pr = geometria_pr()
    x0, y0, x1, y1 = -54.75, -26.75, -48.0, -22.5
    cand = pd.DataFrame([(x + LADO / 2, y + LADO / 2) for x in np.arange(x0, x1, LADO) for y in np.arange(y0, y1, LADO)],
                        columns=["lon", "lat"])
    soja = ee.Image(gee.LULC).select(f"classification_{ANO_SOJA}").eq(gee.CLASSE_SOJA) \
        .multiply(ee.Image.pixelArea()).divide(1e4).rename("soja_ha").clip(pr)
    fc = ee.FeatureCollection(_feicoes(cand)).filterBounds(pr)
    r = soja.reduceRegions(fc, ee.Reducer.sum(), 120, tileScale=4).getInfo()["features"]
    c = pd.DataFrame([dict(lon=f["properties"]["lon"], lat=f["properties"]["lat"], soja_ha=f["properties"].get("sum") or 0)
                      for f in r])
    c = c[c.soja_ha > MIN_SOJA_HA].reset_index(drop=True)
    c = c.merge(gee.elevacao(c, lado=LADO), on=["lon", "lat"])
    c.insert(0, "id_celula", [f"C{i:03d}" for i in range(len(c))])
    c.round(4).to_csv(ARQ_CEL, index=False)
    print(f"{len(c)} células; soja {c.soja_ha.sum() / 1e6:.2f} Mha; altitude {c.elev.min():.0f}–{c.elev.max():.0f} m")


def solo():
    c = pd.read_csv(ARQ_CEL)
    pr = geometria_pr()
    # perfis_celulas usa a soja de um ano; recorta no PR trocando temporariamente a imagem de cobertura
    lulc = gee.LULC
    gee.LULC = ee.Image(lulc).clip(pr)
    try:
        p = sm.perfis_celulas(c, ano=ANO_SOJA, lado=LADO, escala=250)
    finally:
        gee.LULC = lulc
    p = p.merge(c[["id_celula", "lon", "lat"]], on=["lon", "lat"])
    p.to_csv(ARQ_SOLO, index=False)
    print(f"{len(p)} perfis em {p.id_celula.nunique()} células;",
          p.groupby("solo").area_ha.sum().div(p.area_ha.sum()).round(3).to_dict())


def awc():
    p = pd.read_csv(ARQ_SOLO)
    c = pd.read_csv(ARQ_CEL)
    pr = geometria_pr()
    soja = ee.Image(gee.LULC).select(f"classification_{ANO_SOJA}").eq(gee.CLASSE_SOJA).clip(pr)
    awc = ee.FeatureCollection(AWC_ANA).filter(ee.Filter.notNull(["AWC"])).reduceToImage(["AWC"], ee.Reducer.first())
    img = awc.rename("awc").addBands(sm.classe_solo()).updateMask(soja)
    red = ee.Reducer.mean().group(groupField=1, groupName="classe")
    linhas = []
    for i in range(0, len(c), 40):
        fc = img.reduceRegions(ee.FeatureCollection(_feicoes(c.iloc[i:i + 40])), red, 250, tileScale=8).getInfo()
        for f in fc["features"]:
            for g in f["properties"].get("groups", []):
                linhas.append(dict(lon=f["properties"]["lon"], lat=f["properties"]["lat"], classe=int(g["classe"]),
                                   awc_ana=g["mean"]))
    a = pd.DataFrame(linhas)
    p = p.drop(columns=["awc_ana"], errors="ignore").merge(a, on=["lon", "lat", "classe"], how="left")
    p.to_csv(ARQ_SOLO, index=False)
    print(f"AWC da ANA em {p.awc_ana.notna().sum()}/{len(p)} perfis;",
          p.groupby("solo").awc_ana.describe()[["mean", "min", "max"]].round(3).to_string())


def limpar():
    import zarc_proxy as zp
    p = pd.read_csv(ARQ_SOLO)
    valido = []
    for i, r in p.iterrows():
        ok = True
        for k, c in enumerate(sm.CAMADAS_MB):
            soma = r[f"sand_{c}"] + r[f"silt_{c}"] + r[f"clay_{c}"]
            if soma >= 50:
                for v in ("sand", "silt", "clay"):
                    p.at[i, f"{v}_{c}"] = r[f"{v}_{c}"] * 100 / soma
            elif k >= 4:
                for v in ("sand", "silt", "clay"):
                    p.at[i, f"{v}_{c}"] = p.at[i, f"{v}_{sm.CAMADAS_MB[k - 1]}"]
            else:
                ok = False
        valido.append(ok)
    p["valido"] = valido
    p["ad_zarc"] = [zp.ad_perfil(r) if ok else np.nan for (_, r), ok in zip(p.iterrows(), valido)]
    p["classe_ad"] = [zp.classe_ad(a)[0] if ok else "" for a, ok in zip(p.ad_zarc, valido)]
    p.to_csv(ARQ_SOLO, index=False)
    print(f"{(~p.valido).sum()} perfis descartados ({p.area_ha[~p.valido].sum() / p.area_ha.sum():.2%} da área);",
          p[p.valido].groupby("classe_ad").area_ha.sum().div(p.area_ha[p.valido].sum()).round(3).to_dict())


def clima():
    c = pd.read_csv(ARQ_CEL)
    falhas = gee.baixar_series(c, 1980, 2025, no=False, n_threads=8)
    print(f"clima: {len(c) - len(falhas)}/{len(c)} séries; falhas: {falhas[:5]}")


if __name__ == "__main__":
    gee.iniciar()
    DIR.mkdir(parents=True, exist_ok=True)
    passos = sys.argv[1:] or ["celulas", "solo", "awc", "limpar", "clima"]
    for p in passos:
        print(f"── {p}", flush=True)
        globals()[p]()
