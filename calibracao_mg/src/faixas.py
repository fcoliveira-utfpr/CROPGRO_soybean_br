"""Faixas de variação dos coeficientes para a análise de sensibilidade.

Regra: para cada coeficiente, a faixa é a envoltória (mínimo e máximo) de três fontes:
  1. linhas MINIMA/MAXIMA do SBGRO048.CUL;
  2. cultivares e grupos genéricos do .CUL/.ECO com MG dentro da faixa brasileira;
  3. coeficientes publicados para cultivares brasileiras (dados/literatura).
"""
from pathlib import Path

import numpy as np
import pandas as pd

from genotipo import CLASSE_CUL, COEF_CUL, ler_eco, limites_dssat, tabela_genotipos

RAIZ = Path(__file__).resolve().parents[1]
ARQ_LITERATURA = RAIZ / "dados" / "literatura" / "coeficientes_literatura.csv"

MG_BRASIL = (5.0, 10.0)                 # do RS (MG ~5) ao MATOPIBA/Norte (MG ~9–10)

# Coeficientes do .ECO que mudam com o hábito de crescimento e com o MG
COEF_ECO_CALIB = {"FL-VS": "P", "R1PPO": "P"}

# Coeficientes que entram na análise: P e G do .CUL + os dois do .ECO
COEF_ANALISE = [c for c in COEF_CUL if CLASSE_CUL[c] in "PG"] + list(COEF_ECO_CALIB)


def ler_literatura(arq=ARQ_LITERATURA):
    return pd.read_csv(arq)


def calcular_faixas(mg=MG_BRASIL):
    gen = tabela_genotipos()
    gen = gen[gen["MG"].between(*mg)]
    lit = ler_literatura()
    lim = limites_dssat()

    linhas = []
    for c in COEF_ANALISE:
        fontes = {"cul_mg": gen[c].dropna(), "literatura": lit[c].dropna() if c in lit else pd.Series(dtype=float)}
        if c in lim:
            fontes["dssat_minmax"] = lim[c]
        if c == "FL-VS":                # hábito: determinado (9) a indeterminado (26)
            fontes["habito"] = ler_eco().loc[["DETERM", "INDETE"], c]
        todos = pd.concat(fontes.values())
        linhas.append({
            "coef": c,
            "classe": CLASSE_CUL.get(c, COEF_ECO_CALIB.get(c)),
            "arquivo": "ECO" if c in COEF_ECO_CALIB else "CUL",
            "min": todos.min(),
            "max": todos.max(),
            "dssat_min": lim[c].get("MINIMA", np.nan) if c in lim else np.nan,
            "dssat_max": lim[c].get("MAXIMA", np.nan) if c in lim else np.nan,
            "cul_mg_min": fontes["cul_mg"].min(),
            "cul_mg_max": fontes["cul_mg"].max(),
            "lit_min": fontes["literatura"].min(),
            "lit_max": fontes["literatura"].max(),
            "n_literatura": len(fontes["literatura"]),
        })
    df = pd.DataFrame(linhas).set_index("coef")
    df["amplitude_rel"] = (df["max"] - df["min"]) / ((df["max"] + df["min"]) / 2)
    return df
