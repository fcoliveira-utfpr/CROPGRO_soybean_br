"""Produtos de uso público (24/10/2026): coeficientes do CROPGRO-Soybean por ponto e safra, por GMR, e arquivos do DSSAT.

Saídas em produtos/dados/:
  coeficientes_ponto_safra.csv  — por ponto × safra com satélite: GMR efetivo (média, dp, quantis), semeadura, os seis
                                  coeficientes de fenologia na mediana da posterior, Yp, Yw, Yw/Yp e rendimento do IBGE;
  coeficientes_ponto.csv        — por ponto: médias de 2001–2024 e das colheitas 2021–2024, tendências por década;
  theta_rmg.csv                 — coeficientes de 5,0 a 9,5 (passo 0,1): mediana e intervalo de 90% da posterior;
  SBGRO048_BR.CUL / .ECO        — cultivares genéricas brasileiras por GMR (BRG050…BRG095) e uma cultivar por ponto
                                  (BR0001…), com o GMR efetivo das colheitas 2021–2024 e o crescimento da BRS 284.
Uso: DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1 .venv/bin/python src/gerar_produtos.py
"""
import json
import sys
import warnings
from datetime import date, timedelta
from pathlib import Path

import emcee
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))
import bs2019 as bs  # noqa: E402
import calibrar_fenologia as cf  # noqa: E402
import lacuna_produtividade as lp  # noqa: E402
import treino_fenologia as tf  # noqa: E402

RAIZ = tf.RAIZ
OUT = RAIZ / "produtos" / "dados"
GEN = Path(bs.__file__).parent  # só para referência
CUL_ORIG = Path(sys.modules["DSSATTools"].__file__).parent / "dssat-csm-os" / "Data" / "Genotype"
CADEIA = [RAIZ / "resultados" / "emulador_fenologia" / f"cadeia_v13_{c}.h5" for c in "ab"]
FASE = ["csdl", "ppsen", "em-fl", "fl-sd", "sd-pm", "r1ppo"]

# colunas do .CUL: (nome, formatador) — reproduzem exatamente as linhas do SBGRO048.CUL (verificado em verificar_formato)
F_CUL = [("csdl", "{:6.2f}"), ("ppsen", "{:6.3f}"), ("em-fl", "{:6.1f}"), ("fl-sh", "{:6.1f}"), ("fl-sd", "{:6.1f}"),
         ("sd-pm", "{:6.2f}"), ("fl-lf", "{:6.2f}"), ("lfmax", "{:6.3f}"), ("slavr", "{:5.0f}."), ("sizlf", "{:6.1f}"),
         ("xfrt", "{:6.2f}"), ("wtpsd", "{:6.2f}"), ("sfdur", "{:6.1f}"), ("sdpdv", "{:6.2f}"), ("podur", "{:6.1f}"),
         ("thrsh", "{:6.1f}"), ("sdpro", "dot3"), ("sdlip", "dot3")]
ECO_COLS = ["thvar", "pl-em", "em-v1", "v1-ju", "ju-r0", "pm06", "pm09", "lngsh", "r7-r8", "fl-vs", "trifl", "rwdth",
            "rhght", "r1ppo", "optbi", "slobi"]


def fmt_num(f, v):
    if f == "dot3":                       # o DSSAT grava .400 em vez de 0.400 em SDPRO/SDLIP (e R1PPO no .ECO)
        return "{:6.3f}".format(v).replace("0.", " .", 1)
    return f.format(v)


def linha_cul(var, nome, eco, c):
    cab = f"{var:<6} {nome:<16} {'.':>5} {eco:<6}"
    return cab + "".join(fmt_num(f, c[k]) for k, f in F_CUL)


def ler_eco(cod):
    for l in open(CUL_ORIG / "SBGRO048.ECO", encoding="latin-1"):
        if l.startswith(cod):
            return l.rstrip("\n")
    raise KeyError(cod)


def linha_eco(eco, nome, base, r1ppo, fl_vs):
    """Copia a linha do ecótipo genérico e troca FL-VS e R1PPO nas suas posições (larguras de 6 caracteres)."""
    cab, vals = base[:30], base[30:]
    campos = [vals[i:i + 6] for i in range(0, len(vals), 6)]
    campos[ECO_COLS.index("fl-vs")] = f"{fl_vs:6.2f}"
    r = f"{r1ppo:6.3f}".replace("0.", " .", 1)
    campos[ECO_COLS.index("r1ppo")] = r
    return f"{eco:<6} {nome:<17}" + cab[24:] + "".join(campos)


def verificar_formato():
    linhas = {l[:6]: l.rstrip("\n") for l in open(CUL_ORIG / "SBGRO048.CUL", encoding="latin-1") if l[:6].strip()}
    orig = linhas["990005"]
    v = orig[36:]
    vals = [v[i:i + 6] for i in range(0, len(v), 6)]
    c = {k: float(x) for (k, _), x in zip(F_CUL, vals)}
    nova = linha_cul("990005", orig[7:23].strip(), orig[30:36].strip(), c)
    nova = nova[:7] + orig[7:30] + nova[30:]           # nome e EXPNO como no original
    assert nova == orig, f"formato do .CUL não confere:\n{orig}\n{nova}"
    base = ler_eco("SB0501")
    assert linha_eco("SB0501", base[7:24], base, 0.414, float(base[30:][ECO_COLS.index('fl-vs') * 6:][:6])) == base, \
        "formato do .ECO não confere"


def posterior():
    ch = np.concatenate([emcee.backends.HDFBackend(str(f), read_only=True).get_chain()[4500:].reshape(-1, len(cf.PARAMS))
                         for f in CADEIA])
    return ch


def coef(th, g):
    ix = {n: i for i, n in enumerate(cf.NOMES)}
    return {k: np.clip(th[..., ix[f"a_{k}"]] + th[..., ix[f"b_{k}"]] * (np.asarray(g) - 7), *cf.FAIXAS[k]) for k in FASE}


def crescimento(g):
    mg = min(bs.COEF_S2, key=lambda m: abs(m - g))
    c = dict(bs.COEF_S2_COMUNS)
    c["fl-sh"], c["fl-lf"] = bs.COEF_S2[mg][3], bs.COEF_S2[mg][6]
    return c


def cultivar_completa(th_med, g):
    c = {k: float(v) for k, v in coef(th_med, g).items()}
    c.update(crescimento(g))
    return c


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    verificar_formato()
    ch = posterior()
    med = np.median(ch, axis=0)

    # θ(GMR) com intervalo de 90%
    grade = np.round(np.arange(5.0, 9.51, 0.1), 1)
    am = ch[np.random.default_rng(0).integers(0, len(ch), 4000)]
    linhas = []
    for g in grade:
        cm, ca = coef(med, g), coef(am, g)
        linhas.append({"rmg": g, **{f"{k}": round(float(cm[k]), 4) for k in FASE},
                       **{f"{k}_q05": round(float(np.percentile(ca[k], 5)), 4) for k in FASE},
                       **{f"{k}_q95": round(float(np.percentile(ca[k], 95)), 4) for k in FASE}})
    theta = pd.DataFrame(linhas)
    theta.to_csv(OUT / "theta_rmg.csv", index=False)

    # por ponto × safra
    ge = pd.read_csv(RAIZ / "resultados" / "sbi" / "gmr_efetivo_final.csv")
    pts = pd.read_csv(RAIZ / "dados" / "processados" / "pontos_bs2019_digitalizados.csv")[["lon", "lat", "elev"]] \
        .drop_duplicates(["lon", "lat"])
    yl = pd.read_csv(RAIZ / "resultados" / "produtividade" / "lacuna_pontos.csv")[["lon", "lat", "safra", "yp", "yw", "ibge_kg_ha"]]
    d = ge.merge(pts, on=["lon", "lat"], how="left").merge(yl, on=["lon", "lat", "safra"], how="left")
    d["regiao"] = lp.regiao(d.lat.values, d.lon.values)
    d["ano_colheita"] = d.safra + 1
    d["data_semeadura"] = [(date(int(s), 9, 10) + timedelta(days=int(round(x)))).isoformat() for s, x in zip(d.safra, d.s_q50)]
    cc = coef(med, d.g_q50.values)
    for k in FASE:
        d[k] = np.round(cc[k], 4)
    d["ecotipo_mg"] = np.clip(np.round(d.g_q50), 5, 9).astype(int)
    d["yw_yp"] = (d.yw / d.yp).clip(upper=1).round(3)
    cols = ["lon", "lat", "elev", "regiao", "safra", "ano_colheita", "g_media", "g_dp", "g_q05", "g_q50", "g_q95",
            "s_q50", "s_dp", "data_semeadura"] + FASE + ["ecotipo_mg", "yp", "yw", "yw_yp", "ibge_kg_ha"]
    ren = {"g_media": "rmg_mean", "g_dp": "rmg_sd", "g_q05": "rmg_q05", "g_q50": "rmg_q50", "g_q95": "rmg_q95",
           "s_q50": "sowing_days_after_10sep", "s_dp": "sowing_sd_days", "data_semeadura": "sowing_date",
           "regiao": "region", "safra": "season_start_year", "ano_colheita": "harvest_year", "elev": "elevation_m",
           "ecotipo_mg": "ecotype_mg", "yp": "yp_kg_ha", "yw": "yw_kg_ha", "ibge_kg_ha": "ibge_yield_kg_ha"}
    out = d[cols].rename(columns=ren).round(3)
    out.to_csv(OUT / "coefficients_point_season.csv", index=False)

    # por ponto
    def resumo(g):
        rec = g[g.ano_colheita >= 2021]
        tr = (lambda c: 10 * np.polyfit(g.ano_colheita, g[c], 1)[0]) if len(g) >= 8 else (lambda c: np.nan)
        return pd.Series({"elevation_m": g.elev.iloc[0], "region": g.regiao.iloc[0], "n_seasons": len(g),
                          "rmg_mean": g.g_q50.mean(), "rmg_median": g.g_q50.median(),
                          # mediana das colheitas 2021–2024: robusta a safras com curva de NDVI anômala (ex.: seca de 2021/22)
                          "rmg_recent": rec.g_q50.median() if len(rec) else np.nan,
                          "rmg_trend_per_decade": tr("g_q50"), "sowing_days_mean": g.s_q50.mean(),
                          "sowing_trend_days_per_decade": tr("s_q50"), "yw_yp_mean": g.yw_yp.mean(),
                          "yw_cv": g.yw.std() / g.yw.mean()})
    p = d.groupby(["lon", "lat"]).apply(resumo).reset_index()
    fr = pd.read_csv(RAIZ / "dados" / "processados" / "fracoes_solo_mapbiomas_soja2016.csv")
    fr = fr.rename(columns={"f_argiloso": "soil_clay", "f_argiloarenoso": "soil_sandy_clay", "f_francoarenoso": "soil_sandy_loam"})
    p = p.merge(fr[["lon", "lat", "soil_clay", "soil_sandy_clay", "soil_sandy_loam"]], on=["lon", "lat"], how="left")
    p["rmg_for_cul"] = p.rmg_recent.fillna(p.rmg_median)
    p["cultivar_id"] = [f"BR{i + 1:04d}" for i in range(len(p))]
    cp = coef(med, p.rmg_for_cul.values)
    for k in FASE:
        p[k] = np.round(cp[k], 4)
    p.round(3).to_csv(OUT / "coefficients_point.csv", index=False)

    # arquivos do DSSAT
    cab_cul = [l.rstrip("\n") for l in open(CUL_ORIG / "SBGRO048.CUL", encoding="latin-1")]
    i_var = next(i for i, l in enumerate(cab_cul) if l.startswith("@VAR#"))
    cab_eco = [l.rstrip("\n") for l in open(CUL_ORIG / "SBGRO048.ECO", encoding="latin-1")]
    i_eco = next(i for i, l in enumerate(cab_eco) if l.startswith("@ECO#"))
    nota = ["!", "! Brazilian soybean cultivars by relative maturity group (RMG), CROPGRO-Soybean (DSSAT 4.8).",
            "! Phenology: Bayesian calibration as functions of RMG (posterior median, version 13).",
            "! Growth: BRS 284 (Battisti and Sentelhas, 2019); FL-SH and FL-LF from the nearest maturity group of that study.",
            "! BRG050-BRG095: generic cultivars for RMG 5.0-9.5. BR0001-BR0224: one cultivar per grid point, with the",
            "! effective RMG inferred from MODIS phenology (median of the harvests 2021-2024; see coefficients_point.csv).", "!"]
    lc, le = [], []
    for g in np.round(np.arange(5.0, 9.51, 0.5), 1):
        var, eco = f"BRG{int(round(g * 10)):03d}", f"BRE{int(round(g * 10)):03d}"
        c = cultivar_completa(med, g)
        lc.append(linha_cul(var, f"BR RMG {g:.1f}", eco, c))
        mg = int(np.clip(round(g), 5, 9))
        le.append(linha_eco(eco, f"BR RMG {g:.1f}", ler_eco(f"SB0{mg}01"), c["r1ppo"], bs.ECO_S2["fl-vs"]))
    for r in p.itertuples():
        c = cultivar_completa(med, r.rmg_for_cul)
        nome = f"{r.lon:+.1f}{r.lat:+.1f} {r.rmg_for_cul:.1f}"
        eco = "E" + r.cultivar_id[2:] + "B"
        lc.append(linha_cul(r.cultivar_id, nome[:16], eco, c))
        mg = int(np.clip(round(r.rmg_for_cul), 5, 9))
        le.append(linha_eco(eco, nome[:17], ler_eco(f"SB0{mg}01"), c["r1ppo"], bs.ECO_S2["fl-vs"]))
    open(OUT / "SBGRO048_BR.CUL", "w").write("\n".join([cab_cul[0]] + nota + cab_cul[i_var:i_var + 1] + lc) + "\n")
    open(OUT / "SBGRO048_BR.ECO", "w").write("\n".join([cab_eco[0]] + nota + cab_eco[i_eco:i_eco + 1] + le) + "\n")
    meta = {"posterior": "v13, 2 ensembles x 40 walkers x 9000 steps, second half", "n_points": int(len(p)),
            "n_point_seasons": int(len(out)), "growth": "BRS 284 (Battisti & Sentelhas 2019, Table S2)"}
    json.dump(meta, open(OUT / "metadata.json", "w"), indent=1)
    print(meta)
    print(theta.loc[theta.rmg.isin([5.0, 7.0, 9.0]), ["rmg"] + FASE].to_string(index=False))
