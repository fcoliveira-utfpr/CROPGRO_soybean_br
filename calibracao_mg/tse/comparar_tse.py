"""Comparação TSE × MCMC (v13) × GLUE na calibração fenológica do CROPGRO-Soybean.

1. DSSAT real nos 209 tratamentos experimentais (o terreno do TSE: os mesmos FileX/FileA que ele usou), com as
   cultivares por classe de GMR de cada método: MCMC = SBGRO048_BR publicado (θ(GMR) da v13 nas classes);
   GLUE = melhor conjunto de 60.000 sorteios (glue_amostras_60000.npz) nas classes; TSE-grade = ótimo
   multi-tratamento da grade + R1PPO da etapa do ecótipo; TSE-redução = ótimo da redução de faixa (R1PPO da v13).
   Erro simulado − observado por evento (R1, R5, R7; R8 entra como R7 + 11 d), como no FileA.
2. Retas θ(GMR) ajustadas aos ótimos do TSE por classe (mínimos quadrados ponderados pelo nº de tratamentos) e
   avaliadas com o emulador em todos os 746 termos — inclusive médias multilocais e catálogos, que o TSE não usa —
   e na log-posterior da v13 (escalas de erro e τ na mediana da posterior).
3. Identificabilidade: na grade (1.024 combinações por classe), quantas combinações ficam a até 1 dia do critério
   mínimo do TSE e quanto variam os coeficientes entre elas.

Saída: resultados/tse/comparacao_dssat.csv, comparacao_erros.csv, coeficientes_metodos.csv, identificabilidade.csv.
Uso: DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1 .venv/bin/python tse/comparar_tse.py
"""
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import emcee
import numpy as np
import pandas as pd

AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parent
sys.path.insert(0, str(RAIZ / "src"))
import calibrar_fenologia as cf  # noqa: E402
import arquivos_genotipo as ag  # noqa: E402
from DSSATTools import run as dr  # noqa: E402

OUT = RAIZ / "resultados" / "tse"
CADEIAS = [RAIZ / "resultados" / "emulador_fenologia" / f"cadeia_v13_{c}.h5" for c in "ab"]
CUL_BR, ECO_BR = RAIZ / "produtos" / "dados" / "SBGRO048_BR.CUL", RAIZ / "produtos" / "dados" / "SBGRO048_BR.ECO"
COEF_CUL = ["CSDL", "PPSEN", "EM-FL", "FL-SD", "SD-PM"]
CHAVE = {"CSDL": "csdl", "PPSEN": "ppsen", "EM-FL": "em-fl", "FL-SD": "fl-sd", "SD-PM": "sd-pm", "R1PPO": "r1ppo"}
DELTA8 = 11


def posterior():
    ch = np.concatenate([emcee.backends.HDFBackend(str(f), read_only=True).get_chain()[4500:].reshape(-1, len(cf.NOMES))
                         for f in CADEIAS])
    return np.median(ch, axis=0)


def coef_linha(th, g):
    p = dict(zip(cf.NOMES, th))
    return {K: float(np.clip(p[f"a_{k}"] + p[f"b_{k}"] * (g - 7), *cf.FAIXAS[k])) for K, k in CHAVE.items()}


def ler_evaluate(pasta):
    txt = (pasta / "Evaluate.OUT").read_text(encoding="latin-1").splitlines()
    linhas, cab = [], None
    for l in txt:
        if l.startswith("@RUN"):
            cab = l.replace("@", " ").split()
        elif cab and l.strip() and not l.startswith("*"):
            v = l.split()
            if len(v) == len(cab):
                linhas.append(v)
    d = pd.DataFrame(linhas, columns=cab)
    return d[["TN", "ADAPS", "ADAPM", "PDFPS", "PDFPM", "MDAPS", "MDAPM"]].apply(pd.to_numeric, errors="coerce").replace(-99, np.nan)


def simular(coefs, rotulo, t):
    """coefs: {classe: {CSDL..., R1PPO}}; roda cada FileX da classe com a cultivar da classe e devolve sim − obs."""
    res = []
    for classe, g in t.groupby("classe"):
        cod = f"{int(round(classe * 10)):03d}"
        trab = Path(f"/tmp/tsecmp{cod}")
        shutil.rmtree(trab, ignore_errors=True)
        trab.mkdir()
        orig = Path(g.dir.iloc[0])
        for f in orig.glob("B*.SB[XA]"):
            shutil.copy(f, trab)
        shutil.copy(orig.parent / "DSSATPRO.L48", trab)
        v = coefs[classe]
        (trab / "SBGRO048.CUL").write_text(ag.substituir(CUL_BR.read_text(), f"BRG{cod}", {k: v[k] for k in COEF_CUL}))
        (trab / "SBGRO048.ECO").write_text(ag.substituir(ECO_BR.read_text(), f"BRE{cod}", {"R1PPO": v["R1PPO"]}))
        for arq, ga in g.groupby("arquivo"):
            for f in trab.glob("*.OUT"):
                f.unlink()
            subprocess.run([dr.BIN_PATH, "A", arq], cwd=trab, env={"DSSAT_HOME": f"/tmp/tse{cod}/"}, capture_output=True)
            ev = ler_evaluate(trab).rename(columns={"TN": "TRNO"})
            res.append(ga[["classe", "arquivo", "TRNO", "fonte", "local", "gmr"]].merge(ev, on="TRNO", how="left"))
    d = pd.concat(res, ignore_index=True)
    d["metodo"] = rotulo
    for e, s, m in (("r1", "ADAPS", "ADAPM"), ("r5", "PDFPS", "PDFPM"), ("r7", "MDAPS", "MDAPM")):
        d[f"erro_{e}"] = d[s] - d[m]
    return d


def ajustar_retas(tab, n):
    """θ = a + b(GMR − 7) por coeficiente, mínimos quadrados ponderados pelo número de tratamentos da classe."""
    g = tab.index.values - 7
    w = n.reindex(tab.index).values
    out = {}
    for K in tab.columns:
        b, a = np.polyfit(g, tab[K].values, 1, w=np.sqrt(w))
        out[K] = (a, b)
    return out


def theta_de_retas(med, retas):
    th = med.copy()
    ix = {n: i for i, n in enumerate(cf.NOMES)}
    for K, (a, b) in retas.items():
        th[ix[f"a_{CHAVE[K]}"]], th[ix[f"b_{CHAVE[K]}"]] = a, b
    return th


def erros_emulador(c, th, rotulo):
    pr = c.prever(th)
    o = c.obs
    linhas = []
    exp = o.tipo.values == "experimento"
    for e in ("r1", "r5", "r7", "r8"):
        k = exp & (o.estadio.values == e)
        if k.any():
            r = o.dias_ve.values[k] - pr[e][k]
            linhas.append(dict(metodo=rotulo, dado="experimento", estadio=e, n=int(k.sum()), rmse=np.sqrt(np.mean(r ** 2)),
                               vies=-np.mean(r)))
    pred = np.array([pr[e][g].mean() for e, g in zip(c.cat_est, c.cat_grupos)])
    for agr, nome in ((True, "media_multilocal"), (False, "catalogo")):
        for e in np.unique(c.cat_est):
            k = (c.cat_agr == agr) & (c.cat_est == e)
            if k.any():
                r = c.cat_y[k] - pred[k]
                linhas.append(dict(metodo=rotulo, dado=nome, estadio=e, n=int(k.sum()), rmse=np.sqrt(np.mean(r ** 2)),
                                   vies=-np.mean(r)))
    return linhas


def identificabilidade(cod):
    f = OUT / f"grade_{cod}" / "Compute" / "RMSE-Multi-TRT-Avareged.txt"
    if not f.exists():
        return None
    d = pd.read_csv(f, sep=r"\s+", skiprows=1, header=None).iloc[:, [0, 1, 2, 3, 4, -1]]
    d.columns = COEF_CUL + ["criterio"]
    d = d[d.criterio < 9000]
    perto = d[d.criterio <= d.criterio.min() + 1.0]
    linha = dict(classe=int(cod) / 10, n_comb=len(d), crit_min=d.criterio.min(), n_ate_1d=len(perto))
    for K in COEF_CUL:
        lo, hi = cf.FAIXAS[CHAVE[K]]
        linha[f"amplitude_{K}"] = (perto[K].max() - perto[K].min()) / (hi - lo)       # fração da faixa
    return linha


if __name__ == "__main__":
    t = pd.read_csv(AQUI / "tratamentos_tse.csv")
    r = pd.read_csv(OUT / "resumo.csv")
    med = posterior()
    n = t.groupby("classe").size()
    classes = sorted(t.classe.unique())
    cul_br = {cl: coef_linha(med, cl) for cl in classes}     # referência numérica (o CUL publicado é a v13 nas classes)
    tse = {}
    for metodo in ("grade", "reducao"):
        rm = r[r.metodo == metodo].set_index("classe")
        eco = r[r.metodo == "eco"].set_index("classe")
        tse[metodo] = {cl: {**{K: rm.loc[cl, K] for K in COEF_CUL},
                            "R1PPO": eco.loc[cl, "R1PPO"] if metodo == "grade" and cl in eco.index and pd.notna(eco.loc[cl].get("R1PPO"))
                            else cul_br[cl]["R1PPO"]} for cl in classes if cl in rm.index}
    glue = np.load(RAIZ / "resultados" / "glue" / "glue_amostras_60000.npz")
    th_glue = med.copy()
    ix = {nm: i for i, nm in enumerate(cf.NOMES)}
    melhor = glue["theta"][np.nanargmax(glue["logpost"])]
    for nm, v in zip(glue["nomes"], melhor):
        th_glue[ix[str(nm)]] = v
    metodos = {"MCMC v13": cul_br, "GLUE (melhor de 60.000)": {cl: coef_linha(th_glue, cl) for cl in classes},
               "TSE grade": tse["grade"], "TSE redução de faixa": tse["reducao"]}
    # coeficientes por classe
    pd.DataFrame([dict(metodo=m, classe=cl, **v) for m, d in metodos.items() for cl, v in d.items()]).to_csv(
        OUT / "coeficientes_metodos.csv", index=False)
    # 1. DSSAT real nos experimentos
    sim = pd.concat([simular(v, m, t[t.classe.isin(v)]) for m, v in metodos.items()], ignore_index=True)
    sim.to_csv(OUT / "comparacao_dssat.csv", index=False)
    tab = []
    for m, g in sim.groupby("metodo", sort=False):
        for e in ("r1", "r5", "r7"):
            x = g[f"erro_{e}"].dropna()
            tab.append(dict(metodo=m, dado="experimento (DSSAT)", estadio=e, n=len(x), rmse=np.sqrt(np.mean(x ** 2)),
                            vies=x.mean(), mae=x.abs().mean()))
    # 2. retas θ(GMR) e emulador em todos os dados
    c = cf.Calibracao()
    th = {"MCMC v13": med, "GLUE (melhor de 60.000)": th_glue}
    for m in ("grade", "reducao"):
        cc = pd.DataFrame(tse[m]).T
        th[f"TSE {'grade' if m == 'grade' else 'redução de faixa'} (reta)"] = theta_de_retas(med, ajustar_retas(cc, n))
    for m, x in th.items():
        tab += erros_emulador(c, x, m)
    erros = pd.DataFrame(tab)
    erros.to_csv(OUT / "comparacao_erros.csv", index=False)
    lp = {m: c.log_post(x) for m, x in th.items()}
    (OUT / "logpost_metodos.json").write_text(json.dumps(lp, indent=1))
    # 3. identificabilidade na grade
    ident = pd.DataFrame([x for x in (identificabilidade(f"{int(round(cl * 10)):03d}") for cl in classes) if x])
    ident.to_csv(OUT / "identificabilidade.csv", index=False)
    pd.set_option("display.width", 220)
    print(erros.pivot_table(index=["dado", "estadio"], columns="metodo", values="rmse").round(1).to_string())
    print(json.dumps({k: round(v, 1) for k, v in lp.items()}, indent=1))
    print(ident.round(2).to_string(index=False))
