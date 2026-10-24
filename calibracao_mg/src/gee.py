"""Acesso ao Google Earth Engine: clima do Xavier (BR-DWGD), soja e solo do MapBiomas."""
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import ee
import numpy as np
import pandas as pd

PROJETO_GEE = "mapbiomas-brazil"
XAVIER = "projects/ee-alexandrexavier/assets/BR-DWGD"
XAVIER_BANDAS = ["Tmax", "Tmin", "Rs", "pr", "RH", "u2", "ET"]
XAVIER_PASSO = 0.1                                 # grau
XAVIER_ORIGEM = (-73.9, 5.4)                       # canto superior esquerdo da grade (lon, lat)

LULC = "projects/mapbiomas-public/assets/brazil/lulc/collection10/mapbiomas_brazil_collection10_coverage_v2"
CLASSE_SOJA = 39

SOLO_C3 = "projects/mapbiomas-public/assets/brazil/soil/collection3/mapbiomas_brazil_collection3_soil_"

RAIZ = Path(__file__).resolve().parents[1]
DIR_XAVIER = RAIZ / "dados" / "brutos" / "xavier"


def iniciar():
    try:
        ee.Initialize(project=PROJETO_GEE)
    except Exception:
        ee.Authenticate()
        ee.Initialize(project=PROJETO_GEE)


def centro_celula_xavier(lon, lat):
    """Centro da célula de 0,1° do Xavier que contém (lon, lat)."""
    x0, y0 = XAVIER_ORIGEM
    i = np.floor((np.asarray(lon) - x0) / XAVIER_PASSO)
    j = np.floor((y0 - np.asarray(lat)) / XAVIER_PASSO)
    return np.round(x0 + (i + 0.5) * XAVIER_PASSO, 2), np.round(y0 - (j + 0.5) * XAVIER_PASSO, 2)


def _escalar(img):
    """Converte os valores inteiros do asset para unidades físicas (BAND_x_MULT, BAND_x_ADD)."""
    bandas = [img.select(b).multiply(ee.Number(img.get(f"BAND_{b}_MULT"))).add(ee.Number(img.get(f"BAND_{b}_ADD")))
              for b in XAVIER_BANDAS]
    return ee.Image.cat(bandas).rename(XAVIER_BANDAS).copyProperties(img, ["system:time_start"])


def _preencher(img, raio=3):
    """Preenche pixels mascarados com a média focal dos vizinhos válidos (na projeção nativa)."""
    return img.unmask(img.focalMean(raio, "square", "pixels")).copyProperties(img, img.propertyNames())


def _serie_trecho(coords, ini, fim, preencher=False, tentativas=4):
    """Séries diárias em uma lista de centros de pixel; devolve a média entre eles."""
    col = ee.ImageCollection(XAVIER).filterDate(ini, fim)
    if preencher:
        col = col.map(_preencher)
    col = col.map(_escalar)
    for k in range(tentativas):
        try:
            linhas = col.getRegion(ee.Geometry.MultiPoint(coords), 10_000).getInfo()
            break
        except Exception:
            if k == tentativas - 1:
                raise
            time.sleep(15 * (k + 1))
    df = pd.DataFrame(linhas[1:], columns=linhas[0])
    df["data"] = pd.to_datetime(df["time"], unit="ms").dt.normalize()
    df[XAVIER_BANDAS] = df[XAVIER_BANDAS].apply(pd.to_numeric)      # pixels mascarados (fora do Brasil, lagos) → NaN
    return df.groupby("data")[XAVIER_BANDAS].mean()


def pixels_no(lon, lat, n=2):
    """Centros do bloco n × n de pixels do Xavier em volta de um nó da grade (ex.: -54,5; -24,5)."""
    d = XAVIER_PASSO * (np.arange(n) - (n - 1) / 2)
    return [[round(lon + dx, 2), round(lat + dy, 2)] for dx in d for dy in d]


def serie_xavier(lon, lat, ano_ini, ano_fim, no=False, anos_por_trecho=5, cache=True):
    """Série diária do Xavier com cache em CSV.

    no=False: pixel de 0,1° que contém (lon, lat).
    no=True: média dos 4 pixels em volta do nó (lon, lat) — usado nos pontos x,5° de Battisti & Sentelhas (2019).
    """
    if no:
        coords, tag = pixels_no(lon, lat), "no"
    else:
        lon, lat = centro_celula_xavier(lon, lat)
        coords, tag = [[float(lon), float(lat)]], "px"
    arq = DIR_XAVIER / f"xavier_{tag}_{lon:+.2f}_{lat:+.2f}_{ano_ini}_{ano_fim}.csv.gz"
    if cache and arq.exists():
        return pd.read_csv(arq, index_col=0, parse_dates=True)
    # nó/pixel sem dado (fronteira, lagos): média focal dos vizinhos válidos num raio de 3 pixels
    preencher = bool(_serie_trecho(coords, f"{ano_ini}-01-01", f"{ano_ini}-01-02").isna().all(axis=None))
    partes = [_serie_trecho(coords, f"{a}-01-01", f"{min(a + anos_por_trecho, ano_fim + 1)}-01-01", preencher)
              for a in range(ano_ini, ano_fim + 1, anos_por_trecho)]
    df = pd.concat(partes).sort_index()
    df = df[~df.index.duplicated()]
    arq.parent.mkdir(parents=True, exist_ok=True)
    df.round(3).to_csv(arq)
    return df


def baixar_series(pontos, ano_ini, ano_fim, no=False, n_threads=8):
    """Baixa em paralelo as séries de vários pontos (DataFrame com colunas lon, lat)."""
    falhas = []
    with ThreadPoolExecutor(n_threads) as ex:
        fut = {ex.submit(serie_xavier, r.lon, r.lat, ano_ini, ano_fim, no): r.Index for r in pontos.itertuples()}
        for n, f in enumerate(as_completed(fut), 1):
            try:
                f.result()
            except Exception as erro:
                falhas.append((fut[f], str(erro)[:120]))
            if n % 25 == 0 or n == len(fut):
                print(f"{n}/{len(fut)} pontos ({len(falhas)} falhas)")
    return falhas


def fracao_soja_celulas(lon_centros, lat_centros, lado, ano, escala=300):
    """Fração da célula (lado em graus) ocupada por soja no MapBiomas, para uma lista de centros."""
    soja = ee.Image(LULC).select(f"classification_{ano}").eq(CLASSE_SOJA)
    feicoes = [ee.Feature(ee.Geometry.Rectangle([x - lado / 2, y - lado / 2, x + lado / 2, y + lado / 2]),
                          {"lon": float(x), "lat": float(y)}) for x, y in zip(lon_centros, lat_centros)]
    saida = []
    for i in range(0, len(feicoes), 200):
        fc = soja.reduceRegions(ee.FeatureCollection(feicoes[i:i + 200]), ee.Reducer.mean(), escala, tileScale=4)
        saida += [f["properties"] for f in fc.getInfo()["features"]]
    return pd.DataFrame(saida).rename(columns={"mean": "frac_soja"})


def elevacao(pontos, lado=0.2):
    """Altitude média (SRTM 90 m) no quadrado de lado `lado`° em volta de cada ponto."""
    dem = ee.Image("CGIAR/SRTM90_V4")
    fc = ee.FeatureCollection([ee.Feature(ee.Geometry.Rectangle([x - lado / 2, y - lado / 2, x + lado / 2, y + lado / 2]),
                                          {"lon": float(x), "lat": float(y)}) for x, y in zip(pontos.lon, pontos.lat)])
    r = dem.reduceRegions(fc, ee.Reducer.mean(), 900, tileScale=4).getInfo()["features"]
    return pd.DataFrame([f["properties"] for f in r]).rename(columns={"mean": "elev"})
