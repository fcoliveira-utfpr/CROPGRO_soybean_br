"""Lacuna de produtividade (opção A): Yw simulado − produtividade do IBGE, por ponto, região e safra.

IBGE PAM (tabela 1612) agregado à célula de 1° centrada em cada ponto (municípios com centróide a ±0,5°), com peso pela
área colhida. Safra do SBI = ano da semeadura; ano do IBGE = ano da colheita (safra + 1).
Sensibilidade (opção B): Yw comparado ao percentil 90 e 95 das produtividades municipais da região no mesmo ano.

Saídas em resultados/produtividade/: ibge_celulas.csv, lacuna_pontos.csv, lacuna_regiao_ano.csv.
Uso: .venv/bin/python src/lacuna_produtividade.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

RAIZ = Path(__file__).resolve().parents[1]
DIR = RAIZ / "resultados" / "produtividade"
# regiões sojícolas para o resumo (lat/lon do ponto); a primeira condição verdadeira prevalece
REGIOES = [("RS/SC", lambda la, lo: la <= -27), ("PR/SP/MS", lambda la, lo: la <= -20),
           ("MT/RO/Norte", lambda la, lo: (la > -17) & (lo <= -50)), ("MATOPIBA", lambda la, lo: la > -15),
           ("GO/MG/DF/norte MS", lambda la, lo: la > -20)]


def regiao(lat, lon):
    return np.select([f(lat, lon) for _, f in REGIOES], [n for n, _ in REGIOES], default="outra")


def ibge_municipios():
    m = pd.read_csv(RAIZ / "dados" / "brutos" / "ibge" / "sidra_1612_soja_municipios_1980_2024.csv")
    c = pd.read_csv(RAIZ / "dados" / "processados" / "centroides_municipios.csv")
    m = m.merge(c, left_on="cod_municipio", right_on="cod").dropna(subset=["rendimento_kg_ha", "area_colhida_ha"])
    m = m[m.area_colhida_ha > 0]
    m["lon1"], m["lat1"] = np.floor(m.lon) + 0.5, np.floor(m.lat) + 0.5
    return m


def ibge_celulas(m):
    g = m.assign(p=m.rendimento_kg_ha * m.area_colhida_ha).groupby(["lon1", "lat1", "ano"])
    return pd.DataFrame({"ibge_kg_ha": g.p.sum() / g.area_colhida_ha.sum(), "area_colhida_ha": g.area_colhida_ha.sum(),
                         "n_mun": g.size()}).reset_index().rename(columns={"lon1": "lon", "lat1": "lat"})


if __name__ == "__main__":
    m = ibge_municipios()
    cel = ibge_celulas(m)
    cel.to_csv(DIR / "ibge_celulas.csv", index=False)
    y = pd.read_csv(DIR / "yp_yw_pontos.csv")
    y["ano"] = y.safra + 1
    d = y.merge(cel, on=["lon", "lat", "ano"], how="left")
    d["regiao"] = regiao(d.lat.values, d.lon.values)
    d["lacuna_kg_ha"] = d.yw - d.ibge_kg_ha
    d["lacuna_rel"] = d.lacuna_kg_ha / d.yw
    d.to_csv(DIR / "lacuna_pontos.csv", index=False)
    # sensibilidade: P90/P95 municipal por região e ano
    m["regiao"] = regiao(m.lat.values, m.lon.values)
    pct = m.groupby(["regiao", "ano"]).rendimento_kg_ha.quantile([0.9, 0.95]).unstack()
    pct.columns = ["ibge_p90", "ibge_p95"]
    ok = d.dropna(subset=["yw", "ibge_kg_ha"])
    w = ok.area_colhida_ha
    r = ok.assign(w=w, yp_w=ok.yp * w, yw_w=ok.yw * w, ib_w=ok.ibge_kg_ha * w).groupby(["regiao", "ano"])
    res = pd.DataFrame({"n_pontos": r.size(), "yp": r.yp_w.sum() / r.w.sum(), "yw": r.yw_w.sum() / r.w.sum(),
                        "ibge": r.ib_w.sum() / r.w.sum()}).join(pct)
    res["lacuna"] = res.yw - res.ibge
    res["lacuna_rel"] = res.lacuna / res.yw
    res["ibge_p95_sobre_yw"] = res.ibge_p95 / res.yw
    res.round(3).to_csv(DIR / "lacuna_regiao_ano.csv")
    print(f"{len(ok)} ponto-safra com Yw e IBGE; {d.ibge_kg_ha.isna().sum()} sem IBGE na célula")
    print(res.groupby("regiao")[["yp", "yw", "ibge", "ibge_p95", "lacuna", "lacuna_rel"]].mean().round(2).to_string())
    print("\nfração de ponto-safra com Yw < IBGE:", round(float((ok.yw < ok.ibge_kg_ha).mean()), 3))
