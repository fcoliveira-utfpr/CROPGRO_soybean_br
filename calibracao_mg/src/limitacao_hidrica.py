"""Limitação hídrica relativa (decisão de 22/10/2026): Yw/Yp por ponto, região e safra, em vez da lacuna absoluta
Yw − IBGE. Na opção A o Yw fica abaixo do IBGE em 42% dos ponto-safra (crescimento da BRS 284 sem ganho genético;
IBGE +30% em 2001–2024 contra Yw estável), por isso o IBGE entra só para validar a variação entre anos:
anomalias sem tendência (resíduo de uma reta por célula) do Yw simulado × do rendimento do IBGE.

Entrada: resultados/produtividade/lacuna_pontos.csv (de lacuna_produtividade.py).
Saídas em resultados/produtividade/: limitacao_regiao.csv, limitacao_regiao_ano.csv, anomalias_yw_ibge.csv.
Uso: .venv/bin/python src/limitacao_hidrica.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

RAIZ = Path(__file__).resolve().parents[1]
DIR = RAIZ / "resultados" / "produtividade"


def sem_tendencia(g, col):
    """Resíduo de uma reta no ano (mín. 6 safras); relativo à média da célula."""
    if g[col].notna().sum() < 6:
        return pd.Series(np.nan, index=g.index)
    a, b = np.polyfit(g.ano, g[col], 1)
    return (g[col] - (a * g.ano + b)) / g[col].mean()


if __name__ == "__main__":
    d = pd.read_csv(DIR / "lacuna_pontos.csv").dropna(subset=["yp", "yw"])
    d["razao"] = (d.yw / d.yp).clip(upper=1)
    d["w"] = d.area_colhida_ha.fillna(d.area_colhida_ha.median())

    def resumo(g):
        w = g.w / g.w.sum()
        return pd.Series({"n": len(g), "yp": (g.yp * w).sum(), "yw": (g.yw * w).sum(), "razao": (g.razao * w).sum(),
                          "p_razao_lt_0_5": float((g.razao < 0.5).mean()),
                          "cv_yw_entre_anos": float(g.groupby(["lon", "lat"]).yw.agg(lambda x: x.std() / x.mean()).mean()),
                          "tend_razao_dec": float(np.polyfit(g.ano, g.razao, 1)[0] * 10)})

    reg = d.groupby("regiao").apply(resumo)
    reg.loc["Brasil"] = resumo(d)
    reg.round(3).to_csv(DIR / "limitacao_regiao.csv")
    d.groupby(["regiao", "ano"]).apply(lambda g: pd.Series({"razao": np.average(g.razao, weights=g.w),
                                                            "yw": np.average(g.yw, weights=g.w)})) \
        .round(3).to_csv(DIR / "limitacao_regiao_ano.csv")

    # validação da variação entre anos: anomalias sem tendência do Yw e do IBGE por célula
    v = d.dropna(subset=["ibge_kg_ha"]).copy()
    v["anom_yw"] = v.groupby(["lon", "lat"], group_keys=False).apply(lambda g: sem_tendencia(g, "yw"))
    v["anom_ibge"] = v.groupby(["lon", "lat"], group_keys=False).apply(lambda g: sem_tendencia(g, "ibge_kg_ha"))
    v = v.dropna(subset=["anom_yw", "anom_ibge"])
    v[["lon", "lat", "ano", "regiao", "anom_yw", "anom_ibge", "w"]].to_csv(DIR / "anomalias_yw_ibge.csv", index=False)
    # agregado por região × ano (média ponderada pela área): é a escala em que o clima domina a variação
    ra = v.groupby(["regiao", "ano"]).apply(lambda g: pd.Series({"yw": np.average(g.anom_yw, weights=g.w),
                                                                 "ibge": np.average(g.anom_ibge, weights=g.w)}))
    cor = pd.DataFrame({"r_ponto_safra": v.groupby("regiao").apply(lambda g: g.anom_yw.corr(g.anom_ibge)),
                        "r_regiao_ano": ra.groupby("regiao").apply(lambda g: g.yw.corr(g.ibge)),
                        "n_anos": ra.groupby("regiao").size()})
    cor.loc["Brasil"] = [v.anom_yw.corr(v.anom_ibge), ra.yw.corr(ra.ibge), len(ra)]
    cor.round(3).to_csv(DIR / "anomalias_correlacao.csv")
    pd.set_option("display.width", 200)
    print(reg.round(2).to_string())
    print("\nCorrelação das anomalias sem tendência (Yw simulado × IBGE):")
    print(cor.round(2).to_string())
