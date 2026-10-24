"""Leitura dos arquivos de genótipo do CROPGRO-Soybean (SBGRO048.CUL e .ECO)."""
from pathlib import Path
import re

import pandas as pd

import DSSATTools

DIR_GENOTIPO = Path(DSSATTools.__file__).parent / "dssat-csm-os" / "Data" / "Genotype"
ARQ_CUL = DIR_GENOTIPO / "SBGRO048.CUL"
ARQ_ECO = DIR_GENOTIPO / "SBGRO048.ECO"

COEF_CUL = ["CSDL", "PPSEN", "EM-FL", "FL-SH", "FL-SD", "SD-PM", "FL-LF", "LFMAX", "SLAVR",
            "SIZLF", "XFRT", "WTPSD", "SFDUR", "SDPDV", "PODUR", "THRSH", "SDPRO", "SDLIP"]
COEF_ECO = ["THVAR", "PL-EM", "EM-V1", "V1-JU", "JU-R0", "PM06", "PM09", "LNGSH", "R7-R8",
            "FL-VS", "TRIFL", "RWDTH", "RHGHT", "R1PPO", "OPTBI", "SLOBI"]

# Classificação do próprio .CUL: P = fenologia, G = crescimento, N = não usado na calibração
CLASSE_CUL = dict(zip(COEF_CUL, "P P P N P P N G G G N G G G N N N N".split()))

# MG numérico para os grupos abaixo de 1 e para cultivares com MG fracionário no nome
MG_TEXTO = {"000": -2.0, "00": -1.0, "0": 0.0}


def _mg_do_nome(nome):
    m = re.search(r"MG\s*(\d)\.?(\d)", nome)            # "MG 3.6", "MG77PA"
    if m:
        return float(f"{m.group(1)}.{m.group(2)}")
    m = re.search(r"\((\d{1,2}|0{1,3})\)", nome)         # "BRAGG (7)", "ALTONA (00)"
    if m:
        return MG_TEXTO.get(m.group(1), float(m.group(1)))
    m = re.match(r"M GROUP\s+(\S+)$", nome)              # "M GROUP   5", "M GROUP 000"
    if m:
        txt = m.group(1)
        return MG_TEXTO.get(txt, float(txt) if txt.isdigit() else None)
    m = re.search(r"\s(0{1,3})$", nome)                  # "MAJ. 000", "LABR.  00"
    if m:
        return MG_TEXTO[m.group(1)]
    return None


def ler_eco(arq=ARQ_ECO):
    linhas = []
    for l in Path(arq).read_text(encoding="latin-1").splitlines():
        if not l.strip() or l[0] in "*!@":
            continue
        resto = l[24:].split()
        linhas.append([l[:6].strip(), l[7:24].strip(), *resto[:2 + len(COEF_ECO)]])
    df = pd.DataFrame(linhas, columns=["ECO#", "ECONAME", "MG_ECO", "TM", *COEF_ECO])
    df[["MG_ECO", "TM", *COEF_ECO]] = df[["MG_ECO", "TM", *COEF_ECO]].astype(float)
    return df.set_index("ECO#")


def ler_cul(arq=ARQ_CUL):
    linhas = []
    for l in Path(arq).read_text(encoding="latin-1").splitlines():
        if not l.strip() or l[0] in "*!@":
            continue
        resto = l[24:].split()
        linhas.append([l[:6].strip(), l[7:24].strip(), *resto[1:2 + len(COEF_CUL)]])
    df = pd.DataFrame(linhas, columns=["VAR#", "VAR-NAME", "ECO#", *COEF_CUL])
    df[COEF_CUL] = df[COEF_CUL].astype(float)
    df["MG"] = df["VAR-NAME"].map(_mg_do_nome)
    df["GENERICO"] = df["VAR#"].str.startswith("99")
    return df.set_index("VAR#")


def limites_dssat():
    """Linhas MINIMA e MAXIMA do .CUL (faixas de busca usadas pelo GENCALC/GLUE do DSSAT)."""
    cul = ler_cul()
    return cul.loc[["999991", "999992"], COEF_CUL].set_axis(["MINIMA", "MAXIMA"])


def tabela_genotipos():
    """Cultivares do .CUL com os parâmetros do ecótipo de cada uma (NaN se o ecótipo não existe)."""
    cul, eco = ler_cul().drop(["999991", "999992"]), ler_eco()
    df = cul.join(eco, on="ECO#")
    df["ECO_AUSENTE"] = df["ECONAME"].isna()
    df["MG"] = df["MG"].fillna(df["MG_ECO"])          # cultivar sem MG no nome: usa o do ecótipo
    return df
