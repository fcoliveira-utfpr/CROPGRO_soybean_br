"""Fase 3 do plano Evandro: crescimento do CROPGRO-Soybean nos experimentos de Silva (2021, tese ESALQ, cap. 3).

Simula um tratamento da Tab. 3.4 com:
  - fenologia: θ(GMR) da mediana da posterior da v13 no GMR de calibração do experimento (+ R1PPO(GMR)), FL-SH e
    FL-LF da Tab. S2 de B&S (2019) no MG mais próximo — a mesma cultivar da cadeia de produtividade;
  - crescimento: os 7 coeficientes calibráveis (LFMAX, SLAVR, SIZLF, XFRT, WTPSD, SFDUR, SDPDV) passados em `cres`;
    os demais da BRS 284 (bs.COEF_S2_COMUNS);
  - solo: perfil do apêndice F da tese (SBUILD; LL, DUL, SAT, densidade, argila, silte, SRGF);
  - clima: Xavier na célula do experimento (validar_evandro.CELULA); semeadura e densidade reais (Tab. 3.1),
    espaçamento 0,5 m; início 30 d antes da semeadura com 50% da água disponível (como na tese);
  - irrigado (I) = água e N desligados (≈ Yp, como a cadeia); sequeiro (R) = água ligada, N desligado (≈ Yw).
Saídas diárias (PlantGro) para comparar com as séries digitalizadas (crescimento_tese_cap3.csv).

Uso (teste): DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1 .venv/bin/python src/crescimento_silva.py
"""
import shutil
import sys
import warnings
from datetime import timedelta
from multiprocessing import current_process
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parent
sys.path.insert(0, str(AQUI))
sys.path.insert(0, str(RAIZ / "tse"))
import bs2019 as bs  # noqa: E402
import comparar_tse as ct  # noqa: E402
import rodar_bs2019 as rb  # noqa: E402
import treino_fenologia as tf  # noqa: E402
import validar_evandro as ve  # noqa: E402
from DSSATTools.filex import (Field, InitialConditions, InitialConditionsLayer, Planting, SCGeneral,  # noqa: E402
                              SCManagement, SCMethods, SCOptions, SCOutputs, SimulationControls)
from DSSATTools.run import DSSAT  # noqa: E402
from DSSATTools.soil import SoilLayer, SoilProfile  # noqa: E402

EV = RAIZ / "dados" / "literatura" / "evandro_silva"
OUT = RAIZ / "resultados" / "evandro" / "crescimento"

# Coeficientes calibrados e faixas (Morris/literatura; incluem os valores do apêndice G e da BRS 284)
CRES = ["lfmax", "slavr", "sizlf", "xfrt", "wtpsd", "sfdur", "sdpdv"]
FAIXAS = {"lfmax": (0.93, 1.45), "slavr": (300, 500), "sizlf": (140, 260), "xfrt": (0.70, 1.00),
          "wtpsd": (0.14, 0.22), "sfdur": (15, 28), "sdpdv": (1.8, 2.5)}
BRS284 = {k: bs.COEF_S2_COMUNS[k] for k in CRES}

# Densidade (pl/m²): Tab. 3.1 diz "plant m-1", mas TE-1 (22) = TE-100/TE-50 da Tab. 2.1 em pl/m² → lida como pl/m².
DENSIDADE = {"Piracicaba-SP": 28, "Sorriso-MT": 30, "Cruz Alta-RS": 35, "Tupanciretã-RS": 35, "Teresina-PI": 22}
ESPACAMENTO = 50

# Apêndice F: base da camada (cm), densidade (g/cm³), argila, silte (%), SRGF, LL, DUL, SAT (cm³/cm³)
SOLO_F = {
    "Piracicaba-SP": [(10, 1.31, 46, 21, 1.0, .235, .363, .508), (20, 1.31, 46, 21, 1.0, .235, .363, .508),
                      (30, 1.31, 46, 21, .8, .235, .342, .508), (40, 1.32, 57, 18, .8, .264, .432, .534),
                      (50, 1.32, 57, 18, .6, .264, .432, .534), (60, 1.35, 56, 19, .5, .304, .432, .520),
                      (120, 1.35, 56, 19, .3, .304, .432, .520), (150, 1.37, 53, 20, .08, .283, .422, .526)],
    "Sorriso-MT": [(10, 1.55, 50, 38, 1.0, .290, .445, .501), (20, 1.55, 50, 38, 1.0, .290, .445, .501),
                   (30, 1.55, 50, 38, .8, .290, .445, .501), (40, 1.43, 44, 11, .8, .251, .413, .534),
                   (50, 1.43, 44, 11, .6, .251, .413, .534), (60, 1.39, 38, 10, .5, .212, .362, .540),
                   (120, 1.39, 38, 10, .3, .212, .362, .540), (150, 1.39, 38, 10, .08, .212, .362, .540)],
    "Cruz Alta-RS": [(10, 1.59, 17, 8, 1.0, .117, .202, .377), (20, 1.59, 17, 8, 1.0, .117, .202, .377),
                     (30, 1.60, 19, 6, .9, .126, .208, .373), (40, 1.60, 19, 6, .8, .126, .208, .373),
                     (50, 1.60, 19, 6, .7, .127, .209, .373), (60, 1.60, 19, 6, .6, .127, .209, .373),
                     (120, 1.57, 23, 11, .2, .144, .233, .384), (150, 1.57, 23, 11, .067, .144, .233, .384)],
    "Tupanciretã-RS": [(10, 1.56, 32, 9, 1.0, .181, .374, .487), (20, 1.56, 32, 9, 1.0, .181, .374, .487),
                       (30, 1.56, 32, 9, .9, .181, .374, .487), (40, 1.56, 36, 8, .8, .199, .393, .487),
                       (50, 1.56, 36, 8, .7, .199, .393, .487), (60, 1.56, 36, 8, .6, .199, .393, .487),
                       (120, 1.56, 38, 7, .2, .209, .302, .487), (150, 1.56, 38, 7, .067, .209, .302, .487)],
    "Teresina-PI": [(10, 1.51, 11, 20, 1.0, .081, .173, .355), (20, 1.51, 11, 20, 1.0, .081, .173, .368),
                    (30, 1.51, 11, 20, .9, .081, .183, .368), (40, 1.49, 26, 26, .8, .090, .183, .370),
                    (50, 1.49, 26, 26, .7, .090, .193, .370), (60, 1.49, 26, 26, .6, .080, .193, .370),
                    (120, 1.49, 26, 26, .2, .080, .193, .370), (150, 1.50, 25, 24, .08, .070, .203, .370),
                    (180, 1.50, 25, 24, .08, .070, .203, .370)],
}
CODIGO = {"Piracicaba-SP": "SVPI000001", "Sorriso-MT": "SVSO000001", "Cruz Alta-RS": "SVCA000001",
          "Tupanciretã-RS": "SVTU000001", "Teresina-PI": "SVTE000001"}
# SUPOSIÇÃO (a tese não informa): carbono orgânico, pH e Ksat por textura; CN, albedo, SLU1, SLDR como em bs2019
FRACAO_AGUA_INICIAL = 0.5      # tese, cap. 3: 50% da água disponível, 30 d antes da semeadura
VAR_SERIE = {"biomassa_kgha": "CWAD", "folha_kgha": "LWAD", "haste_kgha": "SWAD", "graos_kgha": "GWAD", "iaf": "LAID"}

_D = None
TH = None          # mediana da posterior da v13 (ct.posterior)
_WS = {}


def perfil(local):
    lat, lon = ve.CELULA[local][1], ve.CELULA[local][0]
    camadas = []
    for base, bd, arg, sil, srgf, ll, dul, sat in SOLO_F[local]:
        oc = (1.2 if base <= 20 else 0.6 if base <= 60 else 0.3) * (0.6 if arg < 25 else 1.0)
        ks = 3.0 if arg < 25 else 1.2
        camadas.append(SoilLayer(slb=base, slll=ll, sdul=dul, ssat=sat, sbdm=bd, sloc=round(oc, 2), slcl=arg, slsi=sil,
                                 slhw=5.3, srgf=srgf, ssks=ks, slni=round(oc / 12, 3)))
    textura = "SL" if SOLO_F[local][0][2] < 25 else "C"
    return SoilProfile(table=camadas, name=CODIGO[local], soil_data_source="Silva 2021 ApF", soil_clasification=textura,
                       soil_series_name=local[:10], site="Brasil", country="Brazil", lat=lat, long=lon, scom="R",
                       salb=bs.ALBEDO, slu1=bs.SLU1, sldr=bs.SLDR, slro=bs.CN, slnf=1.0, slpf=1.0)


def condicoes_iniciais(local, dia):
    return InitialConditions(pcr="SB", icdat=dia, icres=0, icren=0, icrip=0, icrid=0, icrn=1, icre=1,
                             table=[InitialConditionsLayer(icbl=c[0], sh2o=round(c[5] + FRACAO_AGUA_INICIAL * (c[6] - c[5]), 3),
                                                           snh4=0.5, sno3=2.0) for c in SOLO_F[local]])


def estacao(local, ano):
    """Clima do Xavier na célula do experimento, de 1/jan do ano anterior à semeadura a 31/dez do seguinte."""
    k = (local, ano)
    if k not in _WS:
        lon, lat, elev = ve.CELULA[local]
        c = pd.read_csv(tf.arq_clima(dict(lon=lon, lat=lat, elev=elev, arq="px")), index_col=0, parse_dates=True) \
            .interpolate(limit=5).ffill().bfill()
        c = c[f"{ano - 1}-01-01":f"{ano + 1}-12-31"]
        ws = bs.estacao(c, lat, lon, elev=elev)
        texto = ws._write_wth()
        ws._write_wth = lambda texto=texto: texto
        _WS[k] = ws
    return _WS[k]


def theta():
    global TH
    if TH is None:
        TH = ct.posterior()
    return TH


def cultivar(gmr, cres):
    coefs = {ct.CHAVE[K]: v for K, v in ct.coef_linha(theta(), gmr).items()}
    c = tf.cultivar(dict(gmr=gmr, **coefs))
    mg = min(bs.COEF_S2, key=lambda m: abs(m - gmr))
    c["fl-sh"], c["fl-lf"] = bs.COEF_S2[mg][3], bs.COEF_S2[mg][6]
    for k in CRES:
        c[k] = float(cres[k])
    return c


def cultivar_apg(nome):
    """Cultivar completa do apêndice G (fenologia + crescimento), ecótipo do MG informado na tese (como Silva)."""
    from DSSATTools.crop import Soybean
    g = pd.read_csv(EV / "coeficientes_tese_apG.csv").set_index("coef")[nome]
    c = Soybean(f"99000{int(round(g['gmr_informado']))}")
    for k in g.index.drop("gmr_informado"):
        c[k.lower()] = float(g[k])
    return c


def tratamentos():
    """Tab. 3.4 + local, GMR de calibração e semeadura (experimentos_tese.csv)."""
    p = pd.read_csv(EV / "produtividade_tese_tab3_4.csv")
    e = pd.read_csv(EV / "experimentos_tese.csv")[["experimento", "local", "gmr_calibracao", "data_semeadura",
                                                   "data_colheita", "cultivar"]]
    return p.merge(e, on="experimento")


def series_obs():
    return pd.read_csv(EV / "crescimento_tese_cap3.csv")


def _ini():
    global _D
    try:
        k = current_process()._identity[0]
    except IndexError:
        k = 0
    pasta = rb.DIR_RUN / f"cs{k:02d}"
    shutil.rmtree(pasta, ignore_errors=True)
    _D = DSSAT(str(pasta))


def simular(t, cres, diario=False, cul=None):
    """Um tratamento (linha de tratamentos()) com os coeficientes de crescimento `cres`.
    Devolve (resumo: dict, série diária: DataFrame ou None)."""
    if _D is None:
        _ini()
    pdate = pd.Timestamp(t["data_semeadura"]).date()
    sdate = pdate - timedelta(days=30)
    agua = t["agua"] == "R"
    ctrl = SimulationControls(
        general=SCGeneral(sdate=sdate, nyers=1),
        options=SCOptions(water="Y" if agua else "N", nitro="N", symbi="N", co2="M"),
        methods=SCMethods(evapo=bs.EVAPO), management=SCManagement(irrig="N", ferti="N", harvs="M"),
        outputs=SCOutputs(grout="Y" if diario else "N", waout="N", caout="N"))
    campo = Field(id_field="SVXV0001", wsta=estacao(t["local"], pdate.year + (pdate.month >= 7)), id_soil=perfil(t["local"]),
                  flob=0, fldt="DR000", fldd=0, flds=0, sltx="SC", sldp=180)
    dens = DENSIDADE[t["local"]]
    plant = Planting(pdate=pdate, ppop=dens, ppoe=dens, plrs=ESPACAMENTO, pldp=bs.PROF_SEMEADURA)
    res = dict(experimento=t["experimento"], tratamento=t["tratamento"])
    try:
        _D.run_treatment(field=campo, cultivar=cul if cul is not None else cultivar(float(t["gmr_calibracao"]), cres), planting=plant,
                         initial_conditions=condicoes_iniciais(t["local"], sdate), simulation_controls=ctrl, verbose=False)
        s = bs._ler_summary(_D.output_files["Summary"]).iloc[0]
        res.update(hwam=float(s.get("HWAM", np.nan)), cwam=float(s.get("CWAM", np.nan)),
                   laix=float(s.get("LAIX", np.nan)), r7=(tf.mo.bs_data(pd.Series([s.MDAT])) -
                                                          tf.mo.bs_data(pd.Series([s.PDAT]))).dt.days.iloc[0])
        serie = None
        if diario:
            import curvas_paragominas as cp
            pg = cp.ler_tabela(_D.output_files["PlantGro"])
            serie = pd.DataFrame({"das": pd.to_numeric(pg.DAP), **{v: pd.to_numeric(pg[c], errors="coerce")
                                                                       for v, c in VAR_SERIE.items()}})
            serie["iaf"] = serie.pop("iaf")
        return res, serie
    except Exception as e:      # noqa: BLE001
        res.update(hwam=np.nan, erro=str(e)[:80])
        return res, None


if __name__ == "__main__":
    import time
    T = tratamentos()
    print(T.head(3).to_string())
    for _, t in T[T.experimento.isin(["PI-2", "TE-1"])].iterrows():
        t0 = time.time()
        r, s = simular(t, BRS284, diario=True)
        print(r, f"{time.time() - t0:.2f} s", None if s is None else s.tail(2).to_string())
