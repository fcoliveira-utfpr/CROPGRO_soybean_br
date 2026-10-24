"""Configuração de Battisti & Sentelhas (2019, Field Crops Res. 240:95–105) para o CROPGRO-Soybean.

Fontes no artigo: solos da Tabela S1, coeficientes e SRGF da Tabela S2, manejo da seção 2.2.
O que o artigo não informa está marcado como SUPOSIÇÃO.
"""
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
from DSSATTools.crop import Soybean
from DSSATTools.filex import (Field, InitialConditions, InitialConditionsLayer, Planting,
                              SCGeneral, SCManagement, SCMethods, SCOptions, SCOutputs,
                              SimulationControls)
from DSSATTools.run import DSSAT
from DSSATTools.soil import SoilLayer, SoilProfile
from DSSATTools.weather import WeatherRecord, WeatherStation

# ── Manejo (seção 2.2) ────────────────────────────────────────────────────────────────────────
SEMEADURAS = {"15-set": (9, 15), "15-out": (10, 15), "15-nov": (11, 15),
              "15-dez": (12, 15), "15-jan": (1, 15), "15-fev": (2, 15)}
ANO_INI, NYERS = 1980, 33                 # safras 1980/81 a 2012/13
INICIO = (3, 15)                          # seis meses antes da primeira semeadura (15/set)
POPULACAO, ESPACAMENTO = 28, 45           # pl/m², cm
PROF_SEMEADURA = 4                        # SUPOSIÇÃO (cm)
FRACAO_AGUA_INICIAL = 0.75                # da água disponível máxima
EVAPO = "R"                               # SUPOSIÇÃO: Priestley-Taylor (padrão do DSSAT); FAO-56 deu estresse hídrico fraco demais

# ── Solos (Tabela S1) ─────────────────────────────────────────────────────────────────────────
# O texto do artigo diz 120 mm para o argiloso; a Tabela S1 dá 150 mm (DUL 0,488 em 30–200 cm).
# Usamos 120 mm, que reproduz melhor o estresse hídrico do artigo (resultados/bs2019/sensibilidade.csv).
# camadas: 0–15, 15–30, 30–200 cm → argila, silte, areia (%), SAT, DUL, LL (cm³/cm³), BD (Mg/m³), Ksat (mm/h)
SOLOS_S1 = {
    "argiloso": [(75, 16, 9, .578, .458, .296, 1.15, 34), (77, 15, 8, .582, .453, .298, 1.15, 11),
                 (81, 12, 7, .622, .429, .284, 1.30, 53)],   # DUL 0,488 na Tabela S1 → 0,429 para 120 mm (texto)
    "argiloarenoso": [(41, 11, 48, .438, .26, .16, 1.35, 67), (41, 11, 48, .438, .26, .16, 1.35, 30),
                      (41, 11, 48, .389, .26, .16, 1.50, 67)],
    "francoarenoso": [(19, 2, 79, .387, .13, .07, 1.55, 150), (19, 2, 79, .387, .13, .07, 1.55, 94),
                      (19, 2, 79, .343, .13, .07, 1.70, 150)],
}
CODIGO_SOLO = {"argiloso": "BSCL000001", "argiloarenoso": "BSSC000001", "francoarenoso": "BSSL000001"}
TEXTURA_DSSAT = {"argiloso": "C", "argiloarenoso": "SC", "francoarenoso": "SL"}
# SUPOSIÇÃO: carbono orgânico (%) por horizonte, a partir dos solos de Battisti et al. (2017, Tabela 2)
SLOC = {"argiloso": (1.16, 0.7, 0.4), "argiloarenoso": (1.0, 0.6, 0.3), "francoarenoso": (0.62, 0.4, 0.2)}
# Camadas DSSAT = camadas do SRGF da Tabela S2 (raiz efetiva até 70 cm) + 150–200 cm
BASE_CAMADAS = [5, 15, 30, 40, 50, 60, 70, 100, 125, 150, 200]
SRGF_S2 = [1.00, 1.00, 0.18, 0.09, 0.05, 0.02, 0.01, 0.0, 0.0, 0.0, 0.0]
CN, ALBEDO, SLU1, SLDR = 73, 0.13, 6.0, 0.60   # CN de Battisti et al. (2017); demais SUPOSIÇÃO (drenagem livre)

# ── Cultivares (Tabela S2): MG por faixa de latitude ──────────────────────────────────────────
COEF_S2 = {  # CSDL PPSEN EM-FL FL-SH FL-SD SD-PM FL-LF
    5.8: (12.68, .306, 17.5, 6.0, 16.0, 26.0, 24.0),
    6.8: (12.43, .315, 19.5, 6.5, 17.5, 29.0, 26.0),
    7.8: (12.17, .325, 21.0, 6.5, 18.0, 30.0, 26.0),
}
COEF_S2_COMUNS = {"lfmax": 1.12, "slavr": 340, "sizlf": 185, "xfrt": 1.0, "wtpsd": 0.21, "sfdur": 23,
                  "sdpdv": 2.3, "podur": 12, "thrsh": 74, "sdpro": 0.4, "sdlip": 0.2}
ECO_S2 = {"fl-vs": 26.0, "rhght": 1.0}
BASE_CULTIVAR = {5.8: "990006", 6.8: "990007", 7.8: "990008"}   # SUPOSIÇÃO: ecótipo do MG inteiro mais próximo


def mg_por_latitude(lat):
    """23–32°S → 5.8; 15–22,9°S → 6.8; 4–14,9°S → 7.8 (ao norte de 4°S também 7.8)."""
    return 5.8 if lat <= -23 else 6.8 if lat <= -15 else 7.8


def cultivar(mg):
    c = Soybean(BASE_CULTIVAR[mg])
    for k, v in zip(["csdl", "ppsen", "em-fl", "fl-sh", "fl-sd", "sd-pm", "fl-lf"], COEF_S2[mg]):
        c[k] = v
    for k, v in COEF_S2_COMUNS.items():
        c[k] = v
    for k, v in ECO_S2.items():
        c["eco#"][k] = v
    return c


def camadas_solo(nome):
    linhas = []
    for base, srgf in zip(BASE_CAMADAS, SRGF_S2):
        h = 0 if base <= 15 else 1 if base <= 30 else 2
        argila, silte, _, sat, dul, ll, bd, ksat = SOLOS_S1[nome][h]
        linhas.append(dict(slb=base, slll=ll, sdul=dul, ssat=sat, sbdm=bd, sloc=SLOC[nome][h], slcl=argila,
                           slsi=silte, slhw=5.3, srgf=srgf, ssks=round(ksat / 10, 2), slni=round(SLOC[nome][h] / 12, 3)))
    return pd.DataFrame(linhas)


def agua_disponivel(nome, prof=80):
    c = camadas_solo(nome)
    esp = np.diff(np.r_[0, c.slb]).clip(max=None)
    topo = np.r_[0, c.slb[:-1]]
    esp_util = np.clip(prof - topo, 0, esp)
    return float(((c.sdul - c.slll) * esp_util * 10).sum())


def perfil_solo(nome, lat, lon):
    c = camadas_solo(nome)
    return SoilProfile(table=[SoilLayer(**r) for r in c.to_dict("records")], name=CODIGO_SOLO[nome],
                       soil_data_source="BS2019 S1", soil_clasification=TEXTURA_DSSAT[nome],
                       soil_series_name=nome, site="Brasil", country="Brazil", lat=lat, long=lon,
                       scom="R", salb=ALBEDO, slu1=SLU1, sldr=SLDR, slro=CN, slnf=1.0, slpf=1.0)


def condicoes_iniciais(nome, dia):
    c = camadas_solo(nome)
    sh2o = c.slll + FRACAO_AGUA_INICIAL * (c.sdul - c.slll)
    return InitialConditions(pcr="SB", icdat=dia, icres=0, icren=0, icrip=0, icrid=0, icrn=1, icre=1,
                             table=[InitialConditionsLayer(icbl=b, sh2o=round(w, 3), snh4=0.5, sno3=2.0)
                                    for b, w in zip(c.slb, sh2o)])


# ── Clima ─────────────────────────────────────────────────────────────────────────────────────
def _es(t):
    return 0.6108 * np.exp(17.27 * t / (t + 237.3))


def estacao(clima, lat, lon, elev, insi="BSXV"):
    """WTH a partir do Xavier: Rs (MJ/m²/d), Tmax, Tmin, chuva, ponto de orvalho (da UR média) e vento (km/d)."""
    d = clima.copy()
    ea = d.RH.clip(5, 100) / 100 * (_es(d.Tmax) + _es(d.Tmin)) / 2          # FAO-56, eq. 19 com UR média
    x = np.log(ea / 0.6108)
    d["dewp"] = np.minimum(237.3 * x / (17.27 - x), d.Tmin + 2)
    d["wind"] = d.u2 * 86.4
    troca = d.Tmax < d.Tmin
    d.loc[troca, ["Tmax", "Tmin"]] = d.loc[troca, ["Tmin", "Tmax"]].values
    regs = [WeatherRecord(date=i.date(), srad=round(r.Rs, 1), tmax=round(r.Tmax, 1), tmin=round(r.Tmin, 1),
                          rain=round(max(r.pr, 0), 1), dewp=round(r.dewp, 1), wind=round(r.wind, 0))
            for i, r in d.iterrows()]
    return WeatherStation(table=regs, lat=lat, long=lon, insi=insi, elev=elev, refht=2, wndht=2)


# ── Simulação ─────────────────────────────────────────────────────────────────────────────────
SAIDAS = ["PDAT", "EDAT", "ADAT", "MDAT", "HWAM", "CWAM", "PRCP", "ETCP", "EPCP", "IRCM"]


def _ler_summary(txt):
    linhas = txt.splitlines()
    i = next(k for k, l in enumerate(linhas) if l.startswith("@"))
    cab = linhas[i].replace("@", " ").split()
    dados = [l.split() for l in linhas[i + 1:] if l.strip()]
    df = pd.DataFrame([d for d in dados if len(d) == len(cab)], columns=cab)
    return df[[c for c in SAIDAS if c in df]].apply(pd.to_numeric, errors="coerce").replace(-99, np.nan)


def simular(dssat, ws, nome_solo, lat, lon, semeadura, agua=True, mg=None):
    """Roda as 33 safras de uma combinação ponto × solo × semeadura; agua=False dá a produtividade potencial."""
    mg = mg or mg_por_latitude(lat)
    mes, dia = SEMEADURAS[semeadura]
    pdate = date(ANO_INI + (mes < 9), mes, dia)
    sdate = date(ANO_INI, *INICIO)
    ctrl = SimulationControls(
        general=SCGeneral(sdate=sdate, nyers=NYERS),
        options=SCOptions(water="Y" if agua else "N", nitro="N", symbi="N", co2="M"),
        methods=SCMethods(evapo=EVAPO),
        management=SCManagement(irrig="N", ferti="N", harvs="M"),
        outputs=SCOutputs(grout="N", waout="N", caout="N"))
    campo = Field(id_field="BSXV0001", wsta=ws, id_soil=perfil_solo(nome_solo, lat, lon), flob=0,
                  fldt="DR000", fldd=0, flds=0, sltx=TEXTURA_DSSAT[nome_solo], sldp=200)
    plant = Planting(pdate=pdate, ppop=POPULACAO, ppoe=POPULACAO, plrs=ESPACAMENTO, pldp=PROF_SEMEADURA)
    dssat.run_treatment(field=campo, cultivar=cultivar(mg), planting=plant,
                        initial_conditions=condicoes_iniciais(nome_solo, sdate),
                        simulation_controls=ctrl, verbose=False)
    s = _ler_summary(dssat.output_files["Summary"])
    s.insert(0, "safra", ANO_INI + np.arange(len(s)))
    return s.assign(lon=lon, lat=lat, mg=mg, solo=nome_solo, semeadura=semeadura,
                    cenario="atingivel" if agua else "potencial")
