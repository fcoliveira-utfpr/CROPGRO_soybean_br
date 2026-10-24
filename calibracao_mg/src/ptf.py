"""Pedotransferência e montagem de perfis de solo do DSSAT a partir do MapBiomas Solo + SoilGrids.

Função principal: Medrado & Lima (2014, Geoderma Regional 1:59–66; "Savanna" em Ottoni et al., 2018), que
estima os parâmetros de van Genuchten (θr, θs, α, n) a partir de argila, silte, areia, matéria orgânica e
densidade do solo. Coeficientes do modelo completo da Tabela 2 do artigo (pdf em dados/literatura/pdfs/medrado_lima2014.pdf),
conferidos contra o exemplo numérico da amostra nº 1120 (seção 3.2).

Unidades do artigo: umidades em base de massa (kg/kg), h em kPa, α em 1/kPa. Para o DSSAT as umidades são
convertidas para volume multiplicando pela densidade do solo (ρágua = 1 Mg/m³).
"""
import numpy as np
import pandas as pd

from solo_mapbiomas import CAMADAS_MB, CAMADAS_SG

# Potenciais para os limites do DSSAT (kPa)
H_CC = 33             # capacidade de campo (SUPOSIÇÃO; 10 kPa é alternativa para Latossolos)
H_PMP = 1500          # ponto de murcha permanente
FATOR_MO = 1.724      # carbono orgânico → matéria orgânica
DP = 2.65             # densidade de partículas (Mg/m³), SUPOSIÇÃO: o MapBiomas e o SoilGrids não têm essa variável

# Tabela 2, modelo completo: {parâmetro: ({i: (a_i, b_i)} para X1 … X5, F)}
# X1 = densidade do solo (Mg/m³), X2 = argila (%), X3 = areia total (%), X4 = silte (%), X5 = matéria orgânica (%)
COEF_MEDRADO_LIMA = {
    "K": ({2: (-0.01831805, 0.89935543), 3: (-0.01131157, 1.00134021),
           4: (-0.00684340, 1.11564515), 5: (0.01622120, 0.48555009)}, 2.01973831),
    "theta_r": ({1: (0.13461391, -1.58555722), 2: (0.04883605, 0.19647777), 3: (-0.00949548, 0.56355803),
                 4: (-0.00005212, 1.43345189), 5: (0.01057297, 1.14951659)}, 0.01153532),
    "alpha": ({1: (-0.00351317, 2.42317427), 2: (-0.03460735, 0.58034639), 3: (-0.01845649, 0.85117517),
               4: (-0.04897867, 0.59621217), 5: (0.00001279, 5.38486557)}, 1.58532681),
    "n": ({1: (1.82865464, -0.57582713), 2: (-0.04262312, 0.59924802), 3: (1.51361637, 0.00985830),
           4: (-0.67293635, 0.05966974), 5: (-0.00163216, 2.66119588)}, -0.64209861),
}


def van_genuchten(h, theta_r, theta_s, alpha, n):
    """θ(h) de van Genuchten (1980), com m = 1 − 1/n; h e 1/α na mesma unidade."""
    m = 1 - 1 / n
    return theta_r + (theta_s - theta_r) / (1 + (alpha * np.abs(h)) ** n) ** m


def _potencia(coef, X):
    termos, F = coef
    return sum(a * X[i] ** b for i, (a, b) in termos.items()) + F


def medrado_lima(argila, silte, areia, mo, bd, dp=DP):
    """Parâmetros de van Genuchten em base de massa (kg/kg; α em 1/kPa) pelo modelo completo."""
    X = {1: bd, 2: argila, 3: areia, 4: silte, 5: max(mo, 0.001)}     # nota da Tabela 2: MO > 0,001
    theta_p = 1 / bd - 1 / dp
    K = min(_potencia(COEF_MEDRADO_LIMA["K"], X), 1.0)                # Eq. 3: K ≤ 1
    theta_s = K * theta_p
    theta_r = max(_potencia(COEF_MEDRADO_LIMA["theta_r"], X), 0.0001)
    alpha = max(_potencia(COEF_MEDRADO_LIMA["alpha"], X), 0.0001)
    n = max(_potencia(COEF_MEDRADO_LIMA["n"], X), 1.0001)
    return theta_r, theta_s, alpha, n


def medrado_lima_volume(argila, silte, areia, mo, bd, dp=DP):
    """Igual a medrado_lima, mas θr e θs em volume (cm³/cm³), como pede o DSSAT."""
    tr, ts, a, n = medrado_lima(argila, silte, areia, mo, bd, dp)
    return tr * bd, ts * bd, a, n


def _no_meio(prof_cm):
    """Índice da camada do SoilGrids (0–5, 5–15, 15–30, 30–60, 60–100 cm) que contém a profundidade."""
    limites = [5, 15, 30, 60, 100]
    return next(i for i, b in enumerate(limites) if prof_cm <= b)


def camadas_dssat(perfil, ptf=medrado_lima_volume, h_cc=H_CC, h_pmp=H_PMP):
    """Converte um perfil (linha de perfis_solo_mapbiomas_*.csv) em camadas do DSSAT de 10 em 10 cm até 1 m."""
    linhas = []
    for k, c in enumerate(CAMADAS_MB):
        meio = 10 * k + 5
        j = _no_meio(meio)
        argila, areia, silte = perfil[f"clay_{c}"], perfil[f"sand_{c}"], perfil[f"silt_{c}"]
        bd, oc = perfil[f"bd_{CAMADAS_SG[j]}"], perfil[f"oc_{CAMADAS_SG[j]}"]
        tr, ts, a, n = ptf(argila, silte, areia, oc * FATOR_MO, bd)
        linhas.append(dict(slb=10 * (k + 1), slll=round(van_genuchten(h_pmp, tr, ts, a, n), 3),
                           sdul=round(van_genuchten(h_cc, tr, ts, a, n), 3), ssat=round(ts, 3), sbdm=round(bd, 2),
                           sloc=round(oc, 2), slcl=round(argila, 1), slsi=round(silte, 1)))
    return pd.DataFrame(linhas)


# ── Perfil do DSSAT a partir de um perfil do MapBiomas ───────────────────────────────────────
MO_MAX = 5.75          # máximo de matéria orgânica do conjunto de validação de Medrado & Lima (2014)
# SRGF por camada de 10 cm, tomado da camada de Battisti & Sentelhas (2019, Tabela S2) que contém o meio da camada
SRGF_10CM = [1.00, 1.00, 0.18, 0.09, 0.05, 0.02, 0.01, 0.0, 0.0, 0.0]
SSKS_CLASSE = {"argiloso": 3.4, "argiloarenoso": 6.7, "francoarenoso": 15.0}   # cm/h, camada superficial da Tabela S1
TEXTURA_DSSAT = {"argiloso": "C", "argiloarenoso": "SC", "francoarenoso": "SL"}


def medrado_lima_limitada(argila, silte, areia, mo, bd):
    """Medrado & Lima em volume, com a matéria orgânica limitada à faixa de validação do artigo."""
    return medrado_lima_volume(argila, silte, areia, min(mo, MO_MAX), bd)


def escalar_carbono(perfil):
    """Escala o carbono do SoilGrids (forma do perfil) para que o estoque de 0–30 cm iguale o do MapBiomas.

    O SoilGrids superestima o carbono nas áreas de lavoura (ex.: ~84 contra ~55 t/ha no oeste do PR).
    """
    p = perfil.copy()
    estoque_sg = sum(p[f"oc_{c}"] * p[f"bd_{c}"] * e for c, e in [("0-5cm", 5), ("5-15cm", 10), ("15-30cm", 15)])
    if pd.notna(p.get("cos_mb_0_30")) and p["cos_mb_0_30"] > 0 and estoque_sg > 0:
        fator = p["cos_mb_0_30"] / estoque_sg
        for c in CAMADAS_SG:
            p[f"oc_{c}"] *= fator
    return p


def camadas_perfil_mb(perfil, h_cc=10):
    """Camadas do DSSAT (0–100 cm de 10 em 10 cm + 100–200 cm) com SRGF e limite de pedregosidade."""
    perfil = escalar_carbono(perfil)
    c = camadas_dssat(perfil, ptf=medrado_lima_limitada, h_cc=h_cc)
    c["srgf"] = SRGF_10CM
    c.loc[c.slb > perfil.get("prof_pedreg50", 200), "srgf"] = 0.0          # raiz não passa da pedregosidade
    c = pd.concat([c, c.iloc[[-1]].assign(slb=200, srgf=0.0)], ignore_index=True)
    c["ssks"] = SSKS_CLASSE[perfil.solo]
    c["slhw"] = 5.3
    c["slni"] = (c.sloc / 12).round(3)
    return c


def degenerado(camadas, minimo=0.005):
    return bool((camadas.sdul - camadas.slll < minimo).any())
