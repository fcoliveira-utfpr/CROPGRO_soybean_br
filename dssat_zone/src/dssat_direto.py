"""Execução direta do DSSAT 4.8: um FileX com vários tratamentos (solo × cultivar × semeadura × água) e 45 safras por
tratamento, sem o DSSATTools no laço. O DSSATTools só fornece o executável, os caminhos e a escrita do solo e do clima.

Os parâmetros de manejo, solo (superfície), condição inicial e ET são os da configuração validada na replicação de
Battisti & Sentelhas (calibracao_mg/src/bs2019.py).
"""
import re
import shutil
import subprocess
from datetime import date

import numpy as np
import pandas as pd

import config as cfg
import bs2019 as bs
from DSSATTools import run as dr
from DSSATTools.soil import SoilLayer, SoilProfile

INSI = "ZNPR"
COLS_SOLO = ["slb", "slll", "sdul", "ssat", "sbdm", "sloc", "slcl", "slsi", "slhw", "srgf", "ssks", "slni"]


def yyddd(d):
    return f"{d.year % 100:02d}{d.timetuple().tm_yday:03d}"


def preparar_pasta(pasta):
    """Pasta de rodada com DSSATPRO, cultivares brasileiras e subpasta de clima."""
    shutil.rmtree(pasta, ignore_errors=True)
    (pasta / "Weather").mkdir(parents=True)
    shutil.copy(cfg.CUL_BR, pasta / "SBGRO048.CUL")
    shutil.copy(cfg.ECO_BR, pasta / "SBGRO048.ECO")
    (pasta / dr.CONFILE).write_text(f"WED    {pasta / 'Weather'}\nMSB    {pasta} dscsm048 CRGRO{dr.VERSION}\n"
                                    f"CRD    {dr.CRD_PATH}\nPSD    {dr.DSSAT_HOME}Pest\n"
                                    f"SLD    {dr.SLD_PATH}\nSTD    {dr.STD_PATH}\n")


def escrever_clima(pasta, clima, lat, lon, elev):
    """Arquivo .WTH com toda a série; devolve o código WSTA (INSI + ano inicial + número de anos)."""
    ws = bs.estacao(clima, lat, lon, elev=elev, insi=INSI)
    a0, a1 = clima.index[0].year, clima.index[-1].year
    wsta = f"{INSI}{a0 % 100:02d}{a1 - a0 + 1:02d}"
    (pasta / "Weather" / f"{wsta}.WTH").write_text(ws._write_wth())
    return wsta


def escrever_solos(pasta, solos, lat, lon):
    """SOIL.SOL com vários perfis. solos: {id_solo (10 car.): (camadas DataFrame, textura DSSAT)}."""
    blocos = []
    for k, (sid, (cam, tex)) in enumerate(solos.items()):
        sp = SoilProfile(table=[SoilLayer(**r) for r in cam[COLS_SOLO].to_dict("records")], name=sid,
                         soil_data_source="MapBiomas", soil_clasification=tex, soil_series_name=sid, site="Parana",
                         country="Brazil", lat=lat, long=lon, scom="R", salb=bs.ALBEDO, slu1=bs.SLU1, sldr=bs.SLDR,
                         slro=bs.CN, slnf=1.0, slpf=1.0)
        txt = sp._write_sol()
        blocos.append(txt if k == 0 else txt.split("\n", 1)[1])      # cabeçalho *SOILS só uma vez
    (pasta / "SOIL.SOL").write_text("\n".join(blocos))


def escrever_filex(pasta, wsta, solos, tratamentos, nome="ZNPR0001.SBX"):
    """tratamentos: lista de dicts com solo (id), gmr, semeadura (mês, dia) e agua (bool). Devolve o DataFrame dos
    tratamentos com TRNO. O início da simulação (e a condição inicial) é o mesmo para todos (cfg.INICIO)."""
    t = pd.DataFrame(tratamentos)
    lv = {c: {v: i + 1 for i, v in enumerate(dict.fromkeys(t[c]))} for c in ("solo", "gmr", "semeadura", "agua")}
    t["TRNO"] = np.arange(1, len(t) + 1)
    sdate = date(cfg.ANO_INI, *cfg.INICIO)
    L = ["*EXP.DETAILS: ZNPR0001SB ZONEAMENTO SOJA PR", "", "*TREATMENTS                        -------------FACTOR LEVELS------------",
         "@N R O C TNAME.................... CU FL SA IC MP MI MF MR MC MT ME MH SM"]
    for r in t.itertuples():
        nomet = f"{r.solo[-4:]} G{r.gmr:.1f} {r.semeadura[1]:02d}/{r.semeadura[0]:02d} {'W' if r.agua else 'P'}"
        cu, fl, mp, sm = lv["gmr"][r.gmr], lv["solo"][r.solo], lv["semeadura"][r.semeadura], lv["agua"][r.agua]
        L.append(f"{r.TRNO:2d} 1 0 0 {nomet:<25.25s}{cu:3d}{fl:3d}  0{fl:3d}{mp:3d}  0  0  0  0  0  0  0{sm:3d}")
    L += ["", "*CULTIVARS", "@C CR INGENO CNAME"]
    L += [f"{i:2d} SB {cfg.cultivar_id(g)} BR RMG {g:.1f}" for g, i in lv["gmr"].items()]
    L += ["", "*FIELDS", "@L ID_FIELD WSTA....  FLSA  FLOB  FLDT  FLDD  FLDS  FLST SLTX   SLDP ID_SOIL    FLNAME"]
    for s, i in lv["solo"].items():
        L.append(f"{i:2d} ZNPR{i:04d} {wsta}   -99     0 DR000     0     0   -99 {solos[s][1]:<5s}   200 {s:<10s} -99")
    L += ["", "*INITIAL CONDITIONS",]
    for s, i in lv["solo"].items():
        cam = solos[s][0]
        L += ["@C   PCR ICDAT  ICRT  ICND  ICRN  ICRE  ICWD ICRES ICREN ICREP ICRIP ICRID ICNAME",
              f"{i:2d}    SB {yyddd(sdate)}   -99   -99     1     1   -99     0  0.00   -99     0     0 -99",
              "@C  ICBL  SH2O  SNH4  SNO3"]
        for b, ll, du in zip(cam.slb, cam.slll, cam.sdul):
            L.append(f"{i:2d} {int(b):5d} {ll + bs.FRACAO_AGUA_INICIAL * (du - ll):5.3f}   0.5   2.0")
    L += ["", "*PLANTING DETAILS",
          "@P PDATE EDATE  PPOP  PPOE  PLME  PLDS  PLRS  PLRD  PLDP  PLWT  PAGE  PENV  PLPH  SPRL                        PLNAME"]
    for (m, d), i in lv["semeadura"].items():
        p = date(cfg.ANO_INI + (m < cfg.INICIO[0]), m, d)       # jan–fev: ano seguinte ao início
        L.append(f"{i:2d} {yyddd(p)}   -99 {bs.POPULACAO:5.1f} {bs.POPULACAO:5.1f}     S     R {bs.ESPACAMENTO:5d}     0 "
                 f"{bs.PROF_SEMEADURA:5.1f}   -99   -99   -99   -99     0                           -99")
    L += ["", "*SIMULATION CONTROLS"]
    for agua, i in lv["agua"].items():
        L += ["@N GENERAL     NYERS NREPS START SDATE RSEED SNAME.................... SMODEL",
              f"{i:2d} GE          {cfg.NYERS:5d}     1     S {yyddd(sdate)}  2150 {'AGUA' if agua else 'POTENCIAL':<25s} CRGRO",
              "@N OPTIONS     WATER NITRO SYMBI PHOSP POTAS DISES  CHEM  TILL   CO2",
              f"{i:2d} OP              {'Y' if agua else 'N'}     N     N     N     N     N     N     N     M",
              "@N METHODS     WTHER INCON LIGHT EVAPO INFIL PHOTO HYDRO NSWIT MESOM MESEV MESOL",
              f"{i:2d} ME              M     M     E     {bs.EVAPO}     R     C     R     1     P     R     1",
              "@N MANAGEMENT  PLANT IRRIG FERTI RESID HARVS",
              f"{i:2d} MA              R     N     N     D     M",
              "@N OUTPUTS     FNAME OVVEW SUMRY FROPT GROUT CAOUT WAOUT NIOUT MIOUT DIOUT VBOSE CHOUT OPOUT FMOPT",
              f"{i:2d} OU              N     Y     Y     1     N     N     N     N     N     N     Y     N     N     A",
              "", "@  AUTOMATIC MANAGEMENT",
              "@N PLANTING    PFRST PLAST PH2OL PH2OU PH2OD PSTMX PSTMN",
              f"{i:2d} PL          {yyddd(sdate)} {yyddd(sdate)}    40   100    30    30    10",
              "@N IRRIGATION  IMDEP ITHRL ITHRU IROFF IMETH IRAMT IREFF",
              f"{i:2d} IR             30    50   100 IB001 IB001    10  1.00",
              "@N NITROGEN    NMDEP NMTHR NAMNT NCODE NAOFF",
              f"{i:2d} NI             30    50    25 IB001 IB001",
              "@N RESIDUES    RIPCN RTIME RIDEP",
              f"{i:2d} RE            100     1    20",
              "@N HARVEST     HFRST HLAST HPCNP HPCNR",
              f"{i:2d} HA              0 {yyddd(sdate)}   100     0", ""]
    (pasta / nome).write_text("\n".join(L) + "\n")
    return t


def ler_summary(pasta):
    txt = (pasta / "Summary.OUT").read_text(encoding="latin-1").splitlines()
    i = next(k for k, l in enumerate(txt) if l.startswith("@"))
    h = txt[i].replace("@", " ")
    cab = h.split()
    # colunas de texto (TNAM, FNAM...) têm espaços: lê pelas posições do cabeçalho (fim de cada nome)
    pos = [m.end() for m in re.finditer(r"\S+", h)]
    ini = [0] + pos[:-1]
    linhas = [[l[a:b].strip() for a, b in zip(ini, pos)] for l in txt[i + 1:] if l.strip()]
    df = pd.DataFrame(linhas, columns=cab)
    num = ["RUNNO", "TRNO", "SDAT", "PDAT", "EDAT", "ADAT", "MDAT", "HDAT", "HWAM", "CWAM", "PRCP", "ETCP", "EPCP"]
    return df[[c for c in num if c in df]].apply(pd.to_numeric, errors="coerce").replace(-99, np.nan)


# Tabela "ENVIRONMENTAL AND STRESS FACTORS" do OVERVIEW.OUT: fases do CROPGRO e colunas (na ordem do arquivo)
FASES_OV = {"Emergence -First Flower": "veg", "First Flower-First Seed": "flor", "First Seed - Phys. Mat.": "ench"}
COLS_OV = ["dias", "tmax", "tmin", "tmed", "srad", "fotoper", "co2", "chuva", "et", "etp", "n_tmin0", "n_tmin2",
           "n_tmax30", "n_tmax32", "n_tmax34", "n_chuva", "agua_fs", "agua_cr", "n_fs", "n_cr", "p_fs", "p_cr"]
USAR_OV = ["dias", "tmax", "tmin", "chuva", "n_tmin0", "n_tmin2", "n_tmax32", "n_tmax34", "agua_fs", "agua_cr"]


def ler_overview(pasta):
    """Uma linha por rodada (RUNNO): variáveis de USAR_OV por fase (veg, flor, ench), com sufixo _fase."""
    linhas, run, dentro = [], 0, False                 # o número do *RUN reinicia após 999: conta as rodadas
    with open(pasta / "OVERVIEW.OUT", encoding="latin-1") as f:
        for l in f:
            if l.startswith("*RUN"):
                run, dentro = run + 1, False
            elif l.startswith("*ENVIRONMENTAL AND STRESS"):
                dentro = True
            elif dentro and l.startswith("*"):
                dentro = False
            elif dentro and l[1:24].strip() in FASES_OV:
                v = l[26:].split()
                if len(v) == len(COLS_OV):
                    linhas.append([run, FASES_OV[l[1:24].strip()]] + [float(x) for x in v])
    d = pd.DataFrame(linhas, columns=["RUNNO", "fase"] + COLS_OV)
    d = d.pivot(index="RUNNO", columns="fase", values=USAR_OV)
    d.columns = [f"{v}_{fa}" for v, fa in d.columns]
    return d.reset_index()


def rodar(pasta, nome="ZNPR0001.SBX", modo="A"):
    for f in pasta.glob("*.OUT"):
        f.unlink()
    r = subprocess.run([dr.BIN_PATH, modo, nome], cwd=pasta, capture_output=True, text=True,
                       env={"DSSAT_HOME": dr.DSSAT_HOME})
    if not (pasta / "Summary.OUT").exists():
        erro = (pasta / "ERROR.OUT").read_text()[:800] if (pasta / "ERROR.OUT").exists() else r.stdout[-800:]
        raise RuntimeError(f"DSSAT falhou: {erro}")
    s = ler_summary(pasta)
    if (pasta / "OVERVIEW.OUT").exists():
        s = s.merge(ler_overview(pasta), on="RUNNO", how="left")
        (pasta / "OVERVIEW.OUT").unlink()
    return s.drop(columns="RUNNO")


def rodar_tratamentos(pasta, wsta, solos, tratamentos, max_trt=99, tentativas=3):
    """Divide os tratamentos em FileX de até 99 (o TRNO tem duas colunas) e devolve o Summary unido aos fatores.
    Um Summary incompleto (o DSSAT interrompido, p. ex. por falta de memória com vários processos) é refeito."""
    partes = []
    for k in range(0, len(tratamentos), max_trt):
        t = escrever_filex(pasta, wsta, solos, tratamentos[k:k + max_trt])
        for n in range(tentativas):
            s = rodar(pasta)
            if len(s) == len(t) * cfg.NYERS:
                break
        else:
            raise RuntimeError(f"Summary com {len(s)} linhas; esperadas {len(t) * cfg.NYERS} ({tentativas} tentativas)")
        partes.append(s.merge(t, on="TRNO").drop(columns="TRNO"))
    return pd.concat(partes, ignore_index=True)
