"""Etapa 5: validação contra o rendimento municipal do IBGE (PAM, tabela 1612), DSSAT × índice tipo ZARC (ISNA).

1. Município do PR → célula: o centroide do município cai na célula (ou na mais próxima a até 0,25°).
2. IBGE: rendimento por município × ano de colheita (1981–2024); tendência por município (regressão local LOWESS,
   fração 0,5, robusta a anos de quebra); anomalia relativa = rendimento / tendência − 1; quebra = anomalia < −20%.
   Só municípios com ≥ 15 anos de dados e ≥ 1.000 ha plantados em média.
3. Preditores da safra (manejo representativo do PR: semeaduras de 21/set a 31/out, GMR 6,0 e 6,5, solos pesados
   pela área): DSSAT — Yw/Yp (ry); ZARC-proxy — ISNA da fase crítica (zarc_proxy.py, AD do ZARC do solo).
   Ambos como anomalia em relação à mediana da célula.
4. Habilidade: correlação de Pearson e Spearman com a anomalia do IBGE (por município e conjunta); AUC para detectar
   as quebras; comparação pareada DSSAT × ISNA (bootstrap por ano, que preserva a correlação espacial).

Saída: resultados/validacao_municipio_ano.csv, resultados/validacao_resumo.json.
Uso: ../calibracao_mg/.venv/bin/python src/etapa5_validacao.py
"""
import json

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import roc_auc_score

import config as cfg
import zarc_proxy as zp
from etapa2_simular import ARQ_CEL, ARQ_SOLO, arq_clima

SEM_REP = ["09-21", "10-01", "10-11", "10-21", "11-01"]
GMR_REP = [6.0, 6.5]
QUEBRA = -0.20


def lowess(y, x, frac=0.5, it=3):
    """LOWESS (Cleveland, 1979): regressão local linear com pesos tricúbicos e iterações robustas (bisquare)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    r = int(np.ceil(frac * n))
    rob = np.ones(n)
    for _ in range(it + 1):
        ajuste = np.empty(n)
        for i in range(n):
            dist = np.abs(x - x[i])
            h = np.sort(dist)[r - 1]
            w = np.clip(1 - (dist / max(h, 1e-12)) ** 3, 0, None) ** 3 * rob
            A = np.c_[np.ones(n), x - x[i]] * np.sqrt(w)[:, None]
            ajuste[i] = np.linalg.lstsq(A, y * np.sqrt(w), rcond=None)[0][0]
        e = y - ajuste
        s = np.median(np.abs(e))
        rob = np.clip(1 - (e / (6 * max(s, 1e-12))) ** 2, 0, None) ** 2
    return ajuste


def municipios_celulas():
    m = pd.read_csv(cfg.CALIB / "dados" / "processados" / "centroides_municipios.csv")
    m = m[m.cod.astype(str).str.startswith("41")]
    c = pd.read_csv(ARQ_CEL)
    d = np.hypot(m.lon.to_numpy()[:, None] - c.lon.to_numpy()[None], m.lat.to_numpy()[:, None] - c.lat.to_numpy()[None])
    j = d.argmin(1)
    m = m.assign(id_celula=c.id_celula.to_numpy()[j], dist=d.min(1))
    return m[m.dist <= 0.25]


def ibge():
    i = pd.read_csv(cfg.CALIB / "dados" / "brutos" / "ibge" / "sidra_1612_soja_municipios_1980_2024.csv")
    i = i[i.cod_municipio.astype(str).str.startswith("41")].copy()
    i["ano"] = i.ano.astype(int)
    i = i[(i.ano >= 1981) & (i.rendimento_kg_ha > 0)]
    grupos = []
    for cod, g in i.groupby("cod_municipio"):
        if len(g) < 15 or g.area_plantada_ha.mean() < 1000:
            continue
        g = g.sort_values("ano")
        t = lowess(g.rendimento_kg_ha, g.ano)
        grupos.append(g.assign(tendencia=t, anom_ibge=g.rendimento_kg_ha / t - 1))
    return pd.concat(grupos)


def preditores():
    ind = pd.read_parquet(cfg.DIR_RES / "indicadores_safra.parquet",
                          columns=["id_celula", "solo", "gmr", "semeadura", "ano", "ry"])
    solos = pd.read_csv(cfg.DIR_RES / "solos_usados.csv")
    ind = ind.merge(solos[["id_solo", "peso"]], left_on="solo", right_on="id_solo")
    ind = ind[ind.semeadura.isin(SEM_REP) & ind.gmr.isin(GMR_REP)]
    ind["ryp"] = ind.ry * ind.peso
    dssat = ind.groupby(["id_celula", "ano", "gmr", "semeadura"]).ryp.sum().groupby(["id_celula", "ano"]).mean()
    # ISNA: mesmo manejo, solos pesados pela área com a AD do ZARC
    cel = pd.read_csv(ARQ_CEL)
    perf = pd.read_csv(ARQ_SOLO)
    sem = [(int(s[:2]), int(s[3:])) for s in SEM_REP]
    partes = []
    for r in cel.itertuples():
        clima = pd.read_csv(arq_clima(r.lon, r.lat), index_col=0, parse_dates=True).interpolate(limit=5).ffill().bfill()
        p = perf[(perf.id_celula == r.id_celula) & perf.valido]
        p = p[p.area_ha / p.area_ha.sum() >= 0.05]
        acc = None
        for s in p.itertuples():
            x = zp.isna_celula(clima, s.ad_zarc, GMR_REP, sem, range(1980, 2025)).groupby("ano").indice_isna.mean() * s.area_ha
            acc = x if acc is None else acc.add(x, fill_value=0)
        partes.append((acc / p.area_ha.sum()).rename("isna").to_frame().assign(id_celula=r.id_celula))
    isna = pd.concat(partes).reset_index().set_index(["id_celula", "ano"])["isna"]
    d = pd.concat([dssat.rename("ry"), isna], axis=1).reset_index()
    d["ano_colheita"] = d.ano + 1
    for v in ("ry", "isna"):
        d[f"anom_{v}"] = d[v] / d.groupby("id_celula")[v].transform("median") - 1
    return d


def auc(y, x):
    return roc_auc_score(y, -x) if 0 < y.sum() < len(y) else np.nan


if __name__ == "__main__":
    mc = municipios_celulas()
    i = ibge().merge(mc[["cod", "id_celula"]], left_on="cod_municipio", right_on="cod")
    p = preditores()
    p.to_csv(cfg.DIR_RES / "preditores_celula_ano.csv", index=False)
    v = i.merge(p, left_on=["id_celula", "ano"], right_on=["id_celula", "ano_colheita"], suffixes=("", "_sim"))
    v = v.dropna(subset=["anom_ibge", "anom_ry", "anom_isna"])
    v["quebra"] = (v.anom_ibge < QUEBRA).astype(int)
    v.to_csv(cfg.DIR_RES / "validacao_municipio_ano.csv", index=False)
    res = dict(n_municipios=int(v.cod_municipio.nunique()), n_municipio_ano=int(len(v)),
               anos=[int(v.ano.min()), int(v.ano.max())], frac_quebras=round(float(v.quebra.mean()), 3))
    for nome in ("ry", "isna"):
        x = v[f"anom_{nome}"]
        por_mun = v.groupby("cod_municipio").apply(lambda g: pearsonr(g.anom_ibge, g[f"anom_{nome}"])[0] if len(g) > 5 else np.nan)
        res[nome] = dict(pearson=round(float(pearsonr(v.anom_ibge, x)[0]), 3),
                         spearman=round(float(spearmanr(v.anom_ibge, x)[0]), 3),
                         auc_quebra=round(float(auc(v.quebra, x)), 3),
                         pearson_municipio_mediana=round(float(por_mun.median()), 3))
    # amplitude: inclinação (pela origem) da anomalia do IBGE contra a anomalia simulada → β da etapa 4
    def inclinacao(b):
        return float(np.sum(b.anom_ibge * b.anom_ry) / np.sum(b.anom_ry ** 2))
    # bootstrap por ano da diferença DSSAT − ISNA
    rng = np.random.default_rng(1)
    anos = v.ano.unique()
    dif_r, dif_auc, betas = [], [], []
    for _ in range(1000):
        b = pd.concat([v[v.ano == a] for a in rng.choice(anos, len(anos))])
        dif_r.append(pearsonr(b.anom_ibge, b.anom_ry)[0] - pearsonr(b.anom_ibge, b.anom_isna)[0])
        dif_auc.append(auc(b.quebra, b.anom_ry) - auc(b.quebra, b.anom_isna))
        betas.append(inclinacao(b))
    res["dssat_menos_isna"] = dict(pearson_ic90=[round(float(q), 3) for q in np.nanquantile(dif_r, [0.05, 0.95])],
                                   auc_ic90=[round(float(q), 3) for q in np.nanquantile(dif_auc, [0.05, 0.95])])
    res["beta"] = dict(central=round(inclinacao(v), 3), ic90=[round(float(q), 3) for q in np.quantile(betas, [0.05, 0.95])],
                       desvio_ibge=round(float(v.anom_ibge.std()), 3), desvio_ry=round(float(v.anom_ry.std()), 3))
    # por ano (média estadual ponderada pela área)
    w = v.assign(peso=v.area_plantada_ha.fillna(v.area_colhida_ha).fillna(1))
    anual = w.groupby("ano").apply(lambda g: pd.Series({c: np.average(g[c], weights=g.peso)
                                                       for c in ("anom_ibge", "anom_ry", "anom_isna")}))
    anual.to_csv(cfg.DIR_RES / "validacao_estado_ano.csv")
    res["estado_ano"] = {n: round(float(pearsonr(anual.anom_ibge, anual[f"anom_{n}"])[0]), 3) for n in ("ry", "isna")}
    (cfg.DIR_RES / "validacao_resumo.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print(json.dumps(res, indent=1, ensure_ascii=False))
