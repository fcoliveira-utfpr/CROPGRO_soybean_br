"""Solo do MapBiomas (coleção 3) nas áreas de soja: frações dos 3 solos de Battisti & Sentelhas (2019) por célula.

Cada pixel de soja (MapBiomas Cobertura, classe 39) recebe a argila média de 0–80 cm (8 camadas de 10 cm) e é
classificado no solo do artigo de argila mais próxima: franco-arenoso (19%), argilo-arenoso (41%) ou
argiloso (75–81%). Limites nos pontos médios: < 30% → franco-arenoso; 30–59% → argilo-arenoso; ≥ 59% → argiloso.
"""
import ee
import pandas as pd

import gee

LIMITES_ARGILA = (30, 59)
CLASSES = {1: "francoarenoso", 2: "argiloarenoso", 3: "argiloso"}


def argila_0_80():
    col = ee.ImageCollection(gee.SOLO_C3 + "clay_fraction_v1").sort("system:index")
    return col.limit(8).mean().rename("argila")          # camadas 000_010 … 070_080


def classe_solo():
    a = argila_0_80()
    return a.lt(LIMITES_ARGILA[0]).multiply(1).add(a.gte(LIMITES_ARGILA[0]).And(a.lt(LIMITES_ARGILA[1])).multiply(2)) \
        .add(a.gte(LIMITES_ARGILA[1]).multiply(3)).rename("classe")


def fracoes_celulas(pontos, ano=2016, lado=1.0, escala=120):
    """Área de soja (ha) em cada classe de solo, por célula centrada em cada ponto."""
    soja = ee.Image(gee.LULC).select(f"classification_{ano}").eq(gee.CLASSE_SOJA)
    img = ee.Image.pixelArea().divide(1e4).addBands(classe_solo()).updateMask(soja)
    feicoes = [ee.Feature(ee.Geometry.Rectangle([x - lado / 2, y - lado / 2, x + lado / 2, y + lado / 2]),
                          {"lon": float(x), "lat": float(y)}) for x, y in zip(pontos.lon, pontos.lat)]
    red = ee.Reducer.sum().group(groupField=1, groupName="classe")
    linhas = []
    for i in range(0, len(feicoes), 50):
        fc = img.reduceRegions(ee.FeatureCollection(feicoes[i:i + 50]), red, escala, tileScale=8).getInfo()
        for f in fc["features"]:
            pr = f["properties"]
            for g in pr.get("groups", []):
                linhas.append(dict(lon=pr["lon"], lat=pr["lat"], solo=CLASSES[int(g["classe"])], area_ha=g["sum"]))
    df = pd.DataFrame(linhas, columns=["lon", "lat", "solo", "area_ha"])
    tab = df.pivot_table(index=["lon", "lat"], columns="solo", values="area_ha", aggfunc="sum").fillna(0)
    tab = tab.reindex(columns=list(CLASSES.values()), fill_value=0)
    tab["soja_ha"] = tab.sum(1)
    for c in CLASSES.values():
        tab[f"f_{c}"] = tab[c] / tab.soja_ha
    return tab.reset_index()


SOILGRIDS = "projects/soilgrids-isric/"
CAMADAS_MB = ["000_010", "010_020", "020_030", "030_040", "040_050", "050_060", "060_070", "070_080", "080_090", "090_100"]
CAMADAS_SG = ["0-5cm", "5-15cm", "15-30cm", "30-60cm", "60-100cm"]
ANO_CARBONO = 2016


def _imagem_perfil():
    """Bandas: argila/areia/silte (%) nas 10 camadas do MapBiomas, densidade (g/cm³) e carbono orgânico (%)
    do SoilGrids nas 5 camadas até 1 m, profundidade até 50% de pedregosidade (cm) e estoque de carbono
    orgânico do MapBiomas em 0–30 cm (t/ha), usado para escalar o carbono do SoilGrids."""
    bandas = []
    for v in ["clay", "sand", "silt"]:
        col = ee.ImageCollection(gee.SOLO_C3 + f"{v}_fraction_v1").sort("system:index").toBands()
        bandas.append(col.rename([f"{v}_{c}" for c in CAMADAS_MB]))
    bd = ee.Image(SOILGRIDS + "bdod_mean").select([f"bdod_{c}_mean" for c in CAMADAS_SG]).divide(100)   # cg/cm³ → g/cm³
    oc = ee.Image(SOILGRIDS + "soc_mean").select([f"soc_{c}_mean" for c in CAMADAS_SG]).divide(100)     # dg/kg → %
    bandas += [bd.rename([f"bd_{c}" for c in CAMADAS_SG]), oc.rename([f"oc_{c}" for c in CAMADAS_SG])]
    ped = ee.ImageCollection(gee.SOLO_C3 + "stoniness_v2").filter(ee.Filter.eq("system:index", "soil_depth_stoniness_50_v2"))
    bandas.append(ped.first().rename("prof_pedreg50"))
    bandas.append(ee.Image(gee.SOLO_C3 + "carbon_v1").select(f"carbon_{ANO_CARBONO}").rename("cos_mb_0_30"))   # t/ha
    return ee.Image.cat(bandas)


def perfis_celulas(pontos, ano=2016, lado=1.0, escala=250):
    """Perfil médio dos pixels de soja de cada célula, por classe de argila (1–3), com a área de cada classe."""
    soja = ee.Image(gee.LULC).select(f"classification_{ano}").eq(gee.CLASSE_SOJA)
    perfil = _imagem_perfil()
    nomes = perfil.bandNames().getInfo()
    img = perfil.addBands(ee.Image.pixelArea().divide(1e4).rename("area_ha")).addBands(classe_solo()).updateMask(soja)
    # média nas bandas do perfil e soma na banda de área (entradas separadas)
    red = ee.Reducer.mean().forEach(nomes).combine(ee.Reducer.sum().forEach(["area_ha"]), sharedInputs=False)
    linhas = []
    for c in (1, 2, 3):
        im = img.updateMask(img.select("classe").eq(c)).select(nomes + ["area_ha"])
        for i in range(0, len(pontos), 40):
            sub = pontos.iloc[i:i + 40]
            fc = ee.FeatureCollection([ee.Feature(ee.Geometry.Rectangle([x - lado / 2, y - lado / 2, x + lado / 2, y + lado / 2]),
                                                  {"lon": float(x), "lat": float(y)}) for x, y in zip(sub.lon, sub.lat)])
            r = im.reduceRegions(fc, red, escala, tileScale=8).getInfo()
            for f in r["features"]:
                pr = f["properties"]
                if pr.get("area_ha"):
                    linhas.append({**pr, "classe": c, "solo": CLASSES[c]})
    return pd.DataFrame(linhas)
