"""Etapa 4: estatísticas por opção (célula × solo × GMR × decêndio), classes de risco e janelas ótimas.

Perda da safra = β · (1 − ry / ry_normal), com ry = Yw/Yp e ry_normal a mediana de ry do manejo representativo do PR
(semeaduras de 21/set a 1/nov, GMR 6,0 e 6,5 — o mesmo da validação) na célula × solo: o rendimento "normal" do local,
como a tendência do IBGE. β calibra a amplitude: no Sul o CROPGRO simula quebras maiores que as observadas (desvio
das anomalias 0,41 × 0,23 no IBGE), e β é a inclinação da anomalia do IBGE (rendimento / tendência − 1) sobre a
anomalia simulada (ry / mediana − 1) (etapa 5, validacao_resumo.json); sem a validação, β = 1. Assim, no manejo
típico a frequência de perda > 20% fica na escala das quebras do IBGE, e opções piores que o típico acumulam perda.
O modelo dá o padrão (quando, onde e com que manejo a quebra acontece) e o IBGE, a amplitude.
β é calibrado pela frequência: no manejo representativo, P(perda > 20%) = fração de quebras do IBGE (anomalia < −20%,
validacao_resumo.json). A inclinação da regressão (β_reg = r · desvio_IBGE / desvio_ry) encolhe a variância e
subestima as caudas — com ela nenhuma opção passa do nível 20% —, e a razão dos desvios (β_var) iguala a variância mas
inclui o ruído municipal do IBGE; as duas ficam como sensibilidade (BETA_MODO = "freq" | "reg" | "var").

Estatísticas (45 safras): ry e Yw medianos, P20 de ry, CVaR20 de ry e de Yw (média dos 20% piores anos), P(perda > 20/30/40%); probabilidades
de seca forte na floração/enchimento (estresse ≥ SECA_FORTE), geada (Tmin ≤ 2 °C entre emergência e R7 corrigido),
calor na floração (≥ CALOR_DIAS dias com Tmax ≥ 36 °C entre R1 e R5), chuva na colheita (≥ CHUVA_COLH mm em 10 dias após R8
corrigido), colheita a tempo do milho safrinha; data mediana do R8 corrigido.

Nível de risco hídrico (análogo ao ZARC): menor nível em {20, 30, 40%} tal que P(perda > 20%) ≤ nível; acima de 40%:
sem indicação. Restrições: geada > 20%, calor > 20%, chuva na colheita > 40%.
Classe: recomendado (nível 20% e sem restrição), recomendado com restrição (nível ≤ 30%, ou nível 20% com
restrição), não recomendado (o resto). Janela ótima: opções da célula × solo não dominadas (Pareto) em Yw mediano,
CVaR20 de Yw e P(safrinha), entre as recomendadas (com ou sem restrição). Yw absoluto e não ry: Yp cai nas semeaduras
tardias, e o ry subiria mesmo com a produtividade caindo.

Saída: resultados/opcoes.csv (uma linha por opção) e resultados/opcoes_celula.csv (ponderado pelos solos).
Uso: ../calibracao_mg/.venv/bin/python src/etapa4_classificacao.py
"""
import json
import os

import numpy as np
import pandas as pd

import config as cfg
from etapa5_validacao import GMR_REP, SEM_REP

SECA_FORTE = 0.5
CALOR_DIAS = 3
CHUVA_COLH = 80
LIM_GEADA, LIM_CALOR, LIM_CHUVA = 0.2, 0.2, 0.4
NIVEIS = (0.2, 0.3, 0.4)


def cvar(x, q=0.2):
    x = np.sort(x[~np.isnan(x)])
    return x[:max(1, int(np.ceil(q * len(x))))].mean() if len(x) else np.nan


BETA_MODO = os.environ.get("ZN_BETA", "freq")


def beta(anom_rep):
    """Escala da perda; anom_rep = 1 − ry/ry_normal no manejo representativo."""
    arq = cfg.DIR_RES / "validacao_resumo.json"
    if not arq.exists():
        return 1.0
    v = json.loads(arq.read_text())
    if BETA_MODO == "reg":
        return v["beta"]["central"]
    if BETA_MODO == "var":
        return v["beta"]["desvio_ibge"] / v["beta"]["desvio_ry"]
    a, b = 0.01, 2.0                                       # bissecção: P(β·anom > 0,2) cresce com β
    for _ in range(40):
        m = (a + b) / 2
        a, b = (m, b) if (m * anom_rep > 0.2).mean() < v["frac_quebras"] else (a, m)
    return round((a + b) / 2, 3)


BETA = 1.0


def estatisticas(g):
    perda = BETA * (1 - g.ry / g.ry_normal)
    return pd.Series({
        "ry_med": g.ry.median(), "ry_p20": g.ry.quantile(0.2), "ry_cvar20": cvar(g.ry.to_numpy()),
        "yw_med": g.yw.median(), "yw_cvar20": cvar(g.yw.to_numpy().astype(float)),
        "p_perda20": (perda > 0.2).mean(), "p_perda30": (perda > 0.3).mean(), "p_perda40": (perda > 0.4).mean(),
        "p_seca_flor": (g.seca_flor >= SECA_FORTE).mean(), "p_seca_ench": (g.seca_ench >= SECA_FORTE).mean(),
        "p_geada": g.geada.mean(), "p_calor": (g.calor36 >= CALOR_DIAS).mean(),
        "p_chuva_colh": (g.chuva_colh >= CHUVA_COLH).mean(), "p_chuva_colh_sc": (g.chuva_colh_sc >= CHUVA_COLH).mean(),
        "p_safrinha": g.safrinha.mean(), "p_safrinha_sc": g.safrinha_sc.mean(),
        "r8c_mediana_doy": g.r8c.dt.dayofyear.median(), "ciclo_r8c": g.ciclo_r8c.median(), "n": len(g)})


def nivel(p30):
    return next((int(100 * n) for n in NIVEIS if p30 <= n), 0)        # 0 = sem indicação (risco > 40%)


def classificar(o):
    o["nivel_risco"] = o.p_perda20.map(nivel)
    restr = []
    for r in o.itertuples():
        m = [n for n, ok in (("geada", r.p_geada > LIM_GEADA), ("calor", r.p_calor > LIM_CALOR),
                             ("chuva na colheita", r.p_chuva_colh > LIM_CHUVA)) if ok]
        restr.append(", ".join(m))
    o["restricoes"] = restr
    o["classe"] = np.select([(o.nivel_risco == 20) & (o.restricoes == ""),
                             o.nivel_risco.isin([20, 30])], ["recomendado", "recomendado com restrição"],
                            "não recomendado")
    return o


def pareto(g):
    """Marca as opções recomendadas não dominadas em (yw_med, yw_cvar20, p_safrinha)."""
    ok = g.classe != "não recomendado"
    x = g[["yw_med", "yw_cvar20", "p_safrinha"]].to_numpy()
    otima = np.zeros(len(g), bool)
    for i in np.flatnonzero(ok.to_numpy()):
        dom = np.all(x[ok.to_numpy()] >= x[i], axis=1) & np.any(x[ok.to_numpy()] > x[i], axis=1)
        otima[i] = not dom.any()
    return pd.Series(otima, index=g.index)


if __name__ == "__main__":
    ind = pd.read_parquet(cfg.DIR_RES / "indicadores_safra.parquet")
    ind = ind[ind.yp > 0]
    rep = ind[ind.semeadura.isin(SEM_REP) & ind.gmr.isin(GMR_REP)]
    normal = rep.groupby(["id_celula", "solo"]).ry.median().rename("ry_normal")
    ind = ind.merge(normal, on=["id_celula", "solo"])
    r = ind[ind.semeadura.isin(SEM_REP) & ind.gmr.isin(GMR_REP)]
    BETA = beta((1 - r.ry / r.ry_normal).dropna().to_numpy())
    o = ind.groupby(["id_celula", "solo", "gmr", "semeadura"]).apply(estatisticas).reset_index()
    o = classificar(o)
    o["janela_otima"] = o.groupby(["id_celula", "solo"], group_keys=False).apply(pareto)
    solos = pd.read_csv(cfg.DIR_RES / "solos_usados.csv")
    o = o.merge(solos[["id_solo", "solo", "peso", "ad_100cm_mm", "awc_ana"]].rename(
        columns={"id_solo": "solo", "solo": "textura"}), on="solo")
    sufixo = "" if BETA_MODO == "freq" else f"_{BETA_MODO}"
    o.to_csv(cfg.DIR_RES / f"opcoes{sufixo}.csv", index=False)
    num = [c for c in o.columns if c.startswith(("ry_", "yw_", "p_")) or c in ("r8c_mediana_doy", "ciclo_r8c")]
    oc = o.assign(**{c: o[c] * o.peso for c in num}).groupby(["id_celula", "gmr", "semeadura"])[num].sum().reset_index()
    oc = classificar(oc)
    oc.to_csv(cfg.DIR_RES / f"opcoes_celula{sufixo}.csv", index=False)
    (cfg.DIR_RES / f"beta{sufixo}.json").write_text(json.dumps({"modo": BETA_MODO, "beta": BETA}))
    print(f"β ({BETA_MODO}) = {BETA}")
    print(o.nivel_risco.value_counts(normalize=True).round(3).to_string())
    print(o.classe.value_counts(normalize=True).round(3).to_string())
    print(o.groupby("semeadura").nivel_risco.apply(lambda x: (x > 0).mean()).round(2).to_string())
