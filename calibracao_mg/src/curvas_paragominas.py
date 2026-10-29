"""Curvas observadas × simuladas em Paragominas-PA (Sampaio et al., 2020): cultivar padrão × calibrada (29/10/2026).

Roda o próprio experimento PARA1801.SBX (14 tratamentos: GMR 7,7 e 8,8, três semeaduras e quatro densidades em 2018)
no DSSAT com três versões de cada cultivar, trocando só a linha *CULTIVARS:
  padrao  — cultivares genéricas do DSSAT do MG inteiro mais próximo (990008 e 990009);
  sampaio — coeficientes calibrados por Sampaio et al. (RB0002 e RB0003; ecótipos SB0701 e SB0801, como no artigo);
  br      — fenologia calibrada como função do GMR (v13) + crescimento da BRS 284 (as cultivares de produtos/).
Clima: BR-DWGD (Xavier) no pixel do experimento (o arquivo original não acompanha o experimento).

Saídas em resultados/curvas/: paragominas_series.csv (PlantGro), paragominas_avaliacao.csv (Evaluate: observado ×
simulado por tratamento) e paragominas_obs_series.csv (séries observadas do .SBT).
Uso: DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1 .venv/bin/python src/curvas_paragominas.py
"""
import re
import shutil
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))
import bs2019 as bs  # noqa: E402
import gerar_produtos as gp  # noqa: E402
import treino_fenologia as tf  # noqa: E402
from DSSATTools import run as dr  # noqa: E402

RAIZ = tf.RAIZ
EXP = RAIZ / "dados" / "brutos" / "dssat_exemplos_br"
OUT = RAIZ / "resultados" / "curvas"
RUN = RAIZ / ".run" / "cp"
LON, LAT, ELEV = -46.95, -3.05, 122
GEN = gp.CUL_ORIG


def preparar_pasta():
    """Cria a pasta de rodada com uma execução qualquer do DSSATTools (DSSATPRO e diretórios) e copia os arquivos."""
    shutil.rmtree(RUN, ignore_errors=True)
    dr.DSSAT(str(RUN))
    (RUN / "Weather").mkdir(exist_ok=True)
    for f in ["PARA1801.SBX", "PARA1801.SBT", "PARA1801.SBA", "BR.SOL"]:
        shutil.copy(EXP / f, RUN / f)
    # configuração que o DSSATTools grava antes de cada rodada (diretórios de clima, genótipo, pragas, solo, padrões)
    (RUN / dr.CONFILE).write_text(f"WED    {RUN / 'Weather'}\nMSB    {RUN} dscsm048 CRGRO{dr.VERSION}\n"
                                  f"CRD    {dr.CRD_PATH}\nPSD    {Path(dr.DSSAT_HOME) / 'Pest'}\n"
                                  f"SLD    {dr.SLD_PATH}\nSTD    {dr.STD_PATH}\n")
    clima = pd.read_csv(tf.arq_clima(dict(lon=LON, lat=LAT, elev=ELEV, arq="px")), index_col=0, parse_dates=True) \
        .interpolate(limit=5).ffill().bfill()
    for ano in (2018, 2019):
        c = clima[f"{ano}-01-01":f"{ano}-12-31"]
        ws = bs.estacao(c, LAT, LON, elev=ELEV, insi="PARA")
        (RUN / "Weather" / f"PARA{str(ano)[2:]}01.WTH").write_text(ws._write_wth())


def arquivos_genotipo():
    """SBGRO048.CUL/.ECO completos do DSSAT + linhas das cultivares de comparação."""
    cul = (GEN / "SBGRO048.CUL").read_text(encoding="latin-1").rstrip("\n").split("\n")
    eco = (GEN / "SBGRO048.ECO").read_text(encoding="latin-1").rstrip("\n").split("\n")
    # Sampaio et al.: ecótipos SB0777/SB0888 não existem no DSSAT 4.8 → SB0701/SB0801, como no artigo
    for l in list(cul):
        if l.startswith("RB0002") or l.startswith("RB0003"):
            novo = "SP" + l[2:6] + l[6:30] + ("SB0701" if l.startswith("RB0002") else "SB0801") + l[36:]
            cul.append(novo)
    th = np.median(gp.posterior(), axis=0)
    for g, var, ecoid in [(7.7, "BRP077", "BREP77"), (8.8, "BRP088", "BREP88")]:
        c = gp.cultivar_completa(th, g)
        cul.append(gp.linha_cul(var, f"BR RMG {g:.1f}", ecoid, c))
        mg = int(np.clip(round(g), 5, 9))
        eco.append(gp.linha_eco(ecoid, f"BR RMG {g:.1f}", gp.ler_eco(f"SB0{mg}01"), c["r1ppo"], bs.ECO_S2["fl-vs"]))
    (RUN / "SBGRO048.CUL").write_text("\n".join(cul) + "\n", encoding="latin-1")
    (RUN / "SBGRO048.ECO").write_text("\n".join(eco) + "\n", encoding="latin-1")


VERSOES = {"padrao": ("990008", "990009"), "sampaio": ("SP0002", "SP0003"), "br": ("BRP077", "BRP088")}


def ler_tabela(txt, marcador="@"):
    """Lê as tabelas de um arquivo de saída do DSSAT (uma por tratamento), devolvendo um DataFrame com TRNO."""
    partes, trno, cab = [], None, None
    for l in txt.splitlines():
        m = re.match(r"\*RUN\s+(\d+)", l) or re.match(r"\s*TREATMENT\s+(\d+)", l)
        if m:
            trno = int(m.group(1))
        if l.startswith(marcador):
            cab = l[1:].split()
            continue
        if cab and l.strip() and not l.startswith(("*", "!", "$")) and l[:1] == " ":
            v = l.split()
            if len(v) == len(cab):
                partes.append([trno] + v)
    return pd.DataFrame(partes, columns=["TRNO"] + cab) if partes else pd.DataFrame()


def rodar(versao):
    sbx = (EXP / "PARA1801.SBX").read_text(encoding="latin-1")
    c77, c88 = VERSOES[versao]
    sbx = re.sub(r"( 1 SB )RB0002", rf"\g<1>{c77}", sbx)
    sbx = re.sub(r"( 2 SB )RB0003", rf"\g<1>{c88}", sbx)
    (RUN / "PARA1801.SBX").write_text(sbx, encoding="latin-1")
    pgs, evs = [], []
    for trno in range(1, 15):              # modo C: um tratamento por vez (o modo A não acha o DSSATPRO desta instalação)
        for f in RUN.glob("*.OUT"):
            f.unlink()
        r = subprocess.run([dr.BIN_PATH, "C", "PARA1801.SBX", str(trno)], cwd=RUN, capture_output=True, text=True,
                           env={"DSSAT_HOME": dr.DSSAT_HOME})
        if not (RUN / "PlantGro.OUT").exists():
            raise RuntimeError(f"DSSAT falhou ({versao}, trt {trno}): {(RUN / 'ERROR.OUT').read_text()[:600] if (RUN / 'ERROR.OUT').exists() else r.stdout[-500:]}")
        pg = ler_tabela((RUN / "PlantGro.OUT").read_text(encoding="latin-1")).assign(TRNO=trno)
        ev = ler_tabela((RUN / "Evaluate.OUT").read_text(encoding="latin-1")).assign(TRNO=trno)
        pgs.append(pg)
        evs.append(ev)
    return pd.concat(pgs).assign(versao=versao), pd.concat(evs).assign(versao=versao)


def obs_series():
    t = ler_tabela((EXP / "PARA1801.SBT").read_text(encoding="latin-1"))
    t["TRNO"] = pd.to_numeric(t["TRNO"]) if "TRNO" in t else t.iloc[:, 1]
    return t


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    preparar_pasta()
    arquivos_genotipo()
    pgs, evs = [], []
    for v in VERSOES:
        pg, ev = rodar(v)
        pgs.append(pg)
        evs.append(ev)
        print(v, len(pg), "linhas PlantGro;", len(ev), "tratamentos no Evaluate", flush=True)
    pd.concat(pgs).to_csv(OUT / "paragominas_series.csv", index=False)
    pd.concat(evs).to_csv(OUT / "paragominas_avaliacao.csv", index=False)
    # séries observadas: o .SBT tem uma tabela única com TRNO na primeira coluna
    linhas = [l for l in (EXP / "PARA1801.SBT").read_text(encoding="latin-1").splitlines() if l.strip() and not l.startswith(("*", "!"))]
    cab = linhas[0][1:].split()
    obs = pd.DataFrame([l.split() for l in linhas[1:] if len(l.split()) == len(cab)], columns=cab)
    obs.to_csv(OUT / "paragominas_obs_series.csv", index=False)
    ev = pd.concat(evs)
    cols = [c for c in ["ADATS", "ADATM", "MDATS", "MDATM", "HWAMS", "HWAMM", "LAIXS", "LAIXM"] if c in ev.columns]
    print(ev[["versao", "TRNO"] + cols].to_string(index=False))
