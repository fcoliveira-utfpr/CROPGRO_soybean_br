"""Prepara as observações de campo da calibração no formato que o TSE lê (FileX + FileA + FileT do DSSAT).

O TSE calibra uma cultivar por vez, com experimentos (FileX) e datas observadas dos estádios (FileA). Por isso:
- só entram os experimentos de campo (453 registros, 209 tratamentos); as médias multilocais e os catálogos não são
  experimentos e ficam para a avaliação fora da amostra;
- cada tratamento vai para a classe de GMR mais próxima (5,0 a 9,5, de 0,5 em 0,5), que é a cultivar BRGxxx do
  SBGRO048_BR (valores de partida da v13; crescimento da BRS 284) — uma instalação simulada do DSSAT por classe;
- configuração das simulações igual à do treino do emulador (treino_fenologia.py): clima Xavier do local, solo
  argiloarenoso de B&S, sem balanço hídrico e de N, início 30 dias antes da semeadura, uma safra por tratamento;
- FileA: ADAT (R1), PDFT (R5) e MDAT (R7). R8 observado sem R7 vira MDAT = R8 − 11 d (δ8 fixo da v13), para o TSE
  ver a mesma relação R7–R8 que a calibração bayesiana; observações contadas da emergência viram dias após a
  semeadura somando a emergência simulada (que não depende dos coeficientes calibrados);
- FileT: o TSE exige um; aqui só uma linha por tratamento sem dado (-99), porque o modo é o fenológico.

Saída: tse/raiz_<classe>/DSSAT48 (fora do git) e tse/tratamentos_tse.csv (tratamento → arquivo, TRNO, observações).
Uso: DELTA8_FIXO=11 .venv/bin/python tse/preparar_tse.py
"""
import shutil
import subprocess
import sys
from datetime import date, timedelta
from pathlib import Path

import numpy as np
import pandas as pd

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parent / "src"))
import bs2019 as bs  # noqa: E402
import treino_fenologia as tf  # noqa: E402
from DSSATTools import run as dr  # noqa: E402

sys.path.insert(0, str(AQUI.parent.parent / "dssat_zone" / "src"))
import dssat_direto as dd  # noqa: E402  (leitor do Summary.OUT por posição do cabeçalho)

DELTA8 = 11
CLASSES = np.arange(5.0, 9.51, 0.5)
BASE = Path(dr.STATIC_PATH)
BR = AQUI.parent / "produtos" / "dados"


def yyyyddd(d):
    return f"{d.year}{d.timetuple().tm_yday:03d}"


def yyddd(d):
    return f"{d.year % 100:02d}{d.timetuple().tm_yday:03d}"


def instalar(raiz, curto):
    """Instalação simulada C:\\DSSAT48 (arquivos .CDE e DSSATPRO do Windows para o TSE; DSSATPRO.L48 para o DSSAT)."""
    r = raiz / "DSSAT48"
    shutil.rmtree(raiz, ignore_errors=True)
    for d in ("Tools/TSE", "Soybean", "Weather", "SOIL", "GENOTYPE"):
        (r / d).mkdir(parents=True)
    for f in BASE.glob("*.CDE"):
        shutil.copy(f, r)
    for f in ("DSSATPRO.v48", "DSCSM048.CTR"):
        shutil.copy(BASE / f, r)
    for d in ("StandardData", "Pest", "Default"):
        shutil.copytree(BASE / d, r / d)
    for f in (BASE / "Genotype").iterdir():
        shutil.copy(f, r / "GENOTYPE")
    shutil.copy(BR / "SBGRO048_BR.CUL", r / "GENOTYPE" / "SBGRO048.CUL")
    shutil.copy(BR / "SBGRO048_BR.ECO", r / "GENOTYPE" / "SBGRO048.ECO")
    for f in (AQUI / "raiz" / "DSSAT48" / "Tools" / "TSE").iterdir():      # eventos fenológicos padrão e ícones
        shutil.copy(f, r / "Tools" / "TSE")
    (r / "DSCSM048.EXE").write_text("#!/bin/sh\n# substituto: o tse_mac.py chama o executável do macOS\n")
    (r / "DSSATPRO.L48").write_text(f"WED    {curto}/Weather\nMSB    {curto} dscsm048 CRGRO048\nCRD    {curto}/GENOTYPE\n"
                                     f"PSD    {curto}/Pest\nSLD    {curto}/SOIL\nSTD    {curto}/StandardData\n")
    return r


def tratamentos():
    o = pd.read_csv(AQUI.parent / "resultados" / "emulador_fenologia" / "observacoes_fenologia.csv")
    e = o[o.tipo == "experimento"]
    chave = ["fonte", "local", "idx_local", "safra", "data_semeadura", "gmr", "base", "lon", "lat", "elev"]
    t = e.pivot_table(index=chave, columns="estadio", values="dias_ve", aggfunc="first").reset_index()
    t["classe"] = np.clip(np.round(t.gmr * 2) / 2, 5.0, 9.5)
    t["semeadura"] = pd.to_datetime(t.safra.astype(str) + "-" + t.data_semeadura, errors="coerce")
    jan_jul = t.semeadura.dt.month < 7                 # safra = ano de início; semeaduras de jan–jun no ano seguinte
    t.loc[jan_jul, "semeadura"] = t.loc[jan_jul, "semeadura"] + pd.DateOffset(years=1)
    return t


def escrever_filex(nome, insi, wsta, solo_id, trat, cultivar):
    L = [f"*EXP.DETAILS: {nome[:8]}SB TSE {trat.local.iloc[0][:20]} {trat.safra.iloc[0]}", "",
         "*TREATMENTS                        -------------FACTOR LEVELS------------",
         "@N R O C TNAME.................... CU FL SA IC MP MI MF MR MC MT ME MH SM"]
    for k, r in enumerate(trat.itertuples(), 1):
        nomet = f"G{r.gmr:.1f} {r.semeadura:%d/%m/%Y}"
        L.append(f"{k:2d} 1 0 0 {nomet:<25.25s}  1  1  0  0{k:3d}  0  0  0  0  0  0  0{k:3d}")
    L += ["", "*CULTIVARS", "@C CR INGENO CNAME", f" 1 SB {cultivar} BR RMG",
          "", "*FIELDS", "@L ID_FIELD WSTA....  FLSA  FLOB  FLDT  FLDD  FLDS  FLST SLTX   SLDP ID_SOIL    FLNAME",
          f" 1 {insi}0001 {wsta}   -99     0 DR000     0     0   -99 SC      200 {solo_id:<10s} -99",
          "", "*PLANTING DETAILS",
          "@P PDATE EDATE  PPOP  PPOE  PLME  PLDS  PLRS  PLRD  PLDP  PLWT  PAGE  PENV  PLPH  SPRL                        PLNAME"]
    for k, r in enumerate(trat.itertuples(), 1):
        L.append(f"{k:2d} {yyddd(r.semeadura)}   -99 {bs.POPULACAO:5.1f} {bs.POPULACAO:5.1f}     S     R {bs.ESPACAMENTO:5d}     0 "
                 f"{bs.PROF_SEMEADURA:5.1f}   -99   -99   -99   -99     0                           -99")
    L += ["", "*SIMULATION CONTROLS"]
    for k, r in enumerate(trat.itertuples(), 1):
        sdate = r.semeadura - timedelta(days=30)
        L += ["@N GENERAL     NYERS NREPS START SDATE RSEED SNAME.................... SMODEL",
              f"{k:2d} GE              1     1     S {yyddd(sdate)}  2150 TSE FENOLOGIA              CRGRO",
              "@N OPTIONS     WATER NITRO SYMBI PHOSP POTAS DISES  CHEM  TILL   CO2",
              f"{k:2d} OP              N     N     N     N     N     N     N     N     M",
              "@N METHODS     WTHER INCON LIGHT EVAPO INFIL PHOTO HYDRO NSWIT MESOM MESEV MESOL",
              f"{k:2d} ME              M     M     E     {bs.EVAPO}     R     C     R     1     P     R     1",
              "@N MANAGEMENT  PLANT IRRIG FERTI RESID HARVS",
              f"{k:2d} MA              R     N     N     N     M",
              "@N OUTPUTS     FNAME OVVEW SUMRY FROPT GROUT CAOUT WAOUT NIOUT MIOUT DIOUT VBOSE CHOUT OPOUT FMOPT",
              f"{k:2d} OU              N     Y     Y     1     Y     N     N     N     N     N     N     N     N     A",
              "", "@  AUTOMATIC MANAGEMENT",                # exigido em cada nível (senão erro IPSIM 5010)
              "@N PLANTING    PFRST PLAST PH2OL PH2OU PH2OD PSTMX PSTMN",
              f"{k:2d} PL          {yyddd(sdate)} {yyddd(sdate)}    40   100    30    30    10",
              "@N IRRIGATION  IMDEP ITHRL ITHRU IROFF IMETH IRAMT IREFF",
              f"{k:2d} IR             30    50   100 IB001 IB001    10  1.00",
              "@N NITROGEN    NMDEP NMTHR NAMNT NCODE NAOFF",
              f"{k:2d} NI             30    50    25 IB001 IB001",
              "@N RESIDUES    RIPCN RTIME RIDEP",
              f"{k:2d} RE            100     1    20",
              "@N HARVEST     HFRST HLAST HPCNP HPCNR",
              f"{k:2d} HA              0 {yyddd(sdate)}   100     0", ""]
    return "\n".join(L) + "\n"


def emergencia(pasta, nome):
    """Dias da semeadura à emergência simulada, por TRNO (Summary.OUT), com a cultivar de partida."""
    for f in pasta.glob("*.OUT"):
        f.unlink()
    curto = (pasta.parent / "DSSATPRO.L48").read_text().split()[1].rsplit("/", 1)[0]
    subprocess.run([dr.BIN_PATH, "A", nome], cwd=pasta, env={"DSSAT_HOME": curto + "/"}, capture_output=True)
    s = dd.ler_summary(pasta)
    return dict(zip(s.TRNO.astype(int), (bs_data(s.EDAT) - bs_data(s.PDAT)).dt.days))


def bs_data(x):
    return pd.to_datetime(x.astype("Int64").astype(str), format="%Y%j", errors="coerce")


if __name__ == "__main__":
    t = tratamentos()
    L = pd.read_csv(AQUI.parent / "resultados" / "emulador_fenologia" / "locais.csv", index_col="idx_local")
    linhas = []
    for c in CLASSES:
        tc = t[t.classe == c]
        if tc.empty:
            continue
        cod = f"{int(round(c * 10)):03d}"
        curto = f"/tmp/tse{cod}"
        raiz = AQUI / f"raiz_{cod}"
        r = instalar(raiz, curto)
        Path(curto).unlink(missing_ok=True)
        Path(curto).symlink_to(r)
        solos = []
        for (idx, safra), g in tc.groupby(["idx_local", "safra"]):
            loc = L.loc[idx]
            insi = f"T{idx:03d}"
            clima = pd.read_csv(tf.arq_clima(loc), index_col=0, parse_dates=True).interpolate(limit=5).ffill().bfill()
            wth = r / "Weather" / f"{insi}8046.WTH"
            if not wth.exists():
                ws = bs.estacao(clima, loc.lat, loc.lon, elev=loc.elev, insi=insi)
                wth.write_text(ws._write_wth())
            wsta = f"{insi}8046"
            perfil = bs.perfil_solo("argiloarenoso", loc.lat, loc.lon)
            solo_id = bs.CODIGO_SOLO["argiloarenoso"]
            if solo_id not in [s[0] for s in solos]:
                solos.append((solo_id, perfil._write_sol()))
            nome = f"B{idx:03d}{safra % 100:02d}01"
            g = g.sort_values("semeadura").reset_index(drop=True)
            (r / "Soybean" / f"{nome}.SBX").write_text(escrever_filex(nome, insi, wsta, solo_id, g, f"BRG{cod}"))
            g["arquivo"], g["TRNO"] = f"{nome}.SBX", np.arange(1, len(g) + 1)
            linhas.append(g.assign(dir=str(r / "Soybean")))
        (r / "SOIL" / "SOIL.SOL").write_text("\n".join([solos[0][1]] + [s[1].split("\n", 1)[1] for s in solos[1:]]))
    tt = pd.concat(linhas, ignore_index=True)
    # emergência simulada (para as observações contadas da emergência) e FileA/FileT
    for (d, arq), g in tt.groupby(["dir", "arquivo"]):
        pasta = Path(d)
        shutil.copy(pasta.parent / "DSSATPRO.L48", pasta)
        em = emergencia(pasta, arq)
        tt.loc[g.index, "e_sim"] = g.TRNO.map(em).values
        A, T = ["@TRNO  ADAT  PDFT  MDAT"], ["@TRNO   DATE  LAID"]      # FileA: campos de 6 colunas, datas AADDD
        for r_ in tt.loc[g.index].itertuples():
            desl = r_.e_sim if r_.base == "emergencia" else 0
            def dt(v):
                return yyddd(r_.semeadura + timedelta(days=int(round(v + desl)))) if pd.notna(v) else "  -99"
            r7 = r_.r7 if pd.notna(getattr(r_, "r7", np.nan)) else (r_.r8 - DELTA8 if pd.notna(getattr(r_, "r8", np.nan)) else np.nan)
            A.append(f"{r_.TRNO:5d} {dt(getattr(r_, 'r1', np.nan)):>5s} {dt(getattr(r_, 'r5', np.nan)):>5s} {dt(r7):>5s}")
            T.append(f"{r_.TRNO:6d} {yyddd(r_.semeadura + timedelta(days=30))}   -99")
        base = arq[:8]
        (pasta / f"{base}.SBA").write_text(f"*EXP. DATA (A): {base}SB TSE\n\n" + "\n".join(A) + "\n")
        (pasta / f"{base}.SBT").write_text(f"*EXP. DATA (T): {base}SB TSE\n\n" + "\n".join(T) + "\n")
        for f in pasta.glob("*.OUT"):
            f.unlink()
        (pasta / "DSSATPRO.L48").unlink()
    cols = ["classe", "dir", "arquivo", "TRNO", "fonte", "local", "safra", "semeadura", "gmr", "base", "e_sim",
            "r1", "r5", "r7", "r8"]
    tt[[c for c in cols if c in tt]].to_csv(AQUI / "tratamentos_tse.csv", index=False)
    print(tt.groupby("classe").agg(tratamentos=("TRNO", "size"), arquivos=("arquivo", "nunique")).to_string())
    print("emergência simulada (d):", tt.e_sim.describe().round(1).to_dict())
