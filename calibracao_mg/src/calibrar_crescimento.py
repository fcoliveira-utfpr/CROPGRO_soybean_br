"""Fase 3: calibração do crescimento com os experimentos de Silva (2021) — protocolo, emenda 2.

1. Amostra: hipercubo latino (N conjuntos dos 7 coeficientes em cs.FAIXAS) + referências (BRS 284, apêndice G por
   cultivar e apêndice G fora da amostra), simulados em todos os tratamentos I-1000N e R-1000N (série diária).
2. Para cada dobra (deixa um experimento I-1000N de fora) e para o ajuste completo: melhor J da amostra nos
   experimentos de treino + Nelder-Mead a partir dos 3 melhores.
3. Critérios G1–G3 fora da amostra.

Saídas em resultados/evandro/crescimento/: amostra.csv (conjunto × tratamento), coeficientes.csv (por dobra),
avaliacao.csv (tratamento × versão), criterios_G.csv.
Uso: DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1 .venv/bin/python src/calibrar_crescimento.py [N=4000] [processos=4]
"""
import os
import sys
from multiprocessing import get_context
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import qmc

sys.path.insert(0, str(Path(__file__).resolve().parent))
import crescimento_silva as cs  # noqa: E402

OUT = cs.OUT
CAL = ["PI-1", "PI-2", "PI-3", "TU-1", "TE-1"]          # I-1000N (calibração)
SIG_Y, SIG_S = 0.10, 0.25
SERIES = ["biomassa_kgha", "iaf"]
T = None
OBS = None


def preparar():
    global T, OBS
    t = cs.tratamentos()
    T = t[t.n == "1000N"].reset_index(drop=True)
    s = cs.series_obs()
    OBS = {k: g for k, g in s[s.n == "1000N"].groupby(["experimento", "tratamento"])}


def metricas(t, cres):
    """Simula um tratamento e resume: produtividade, biomassa final, IAF máximo e nRMSE das séries."""
    r, serie = cs.simular(t, cres, diario=True)
    o = OBS.get((t["experimento"], t["tratamento"]))
    lin = dict(experimento=t["experimento"], tratamento=t["tratamento"], hwam=r.get("hwam", np.nan),
               cwam=r.get("cwam", np.nan), r7=r.get("r7", np.nan))
    lin["laix"] = float(serie.iaf.max()) if serie is not None else np.nan
    for v in SERIES:
        lin[f"nrmse_{v}"] = np.nan
        if serie is None or o is None:
            continue
        ov = o[o.variavel == v]
        if len(ov) < 3:
            continue
        sim = np.interp(ov.das, serie.das, serie[v], right=serie[v].iloc[-1])
        lin[f"nrmse_{v}"] = float(np.sqrt(np.mean((sim - ov.valor) ** 2)) / ov.valor.mean())
    return lin


def _tarefa(args):
    rotulo, cres, linhas = args
    out = []
    for i in linhas:
        lin = metricas(T.iloc[i], cres)
        out.append(dict(conjunto=rotulo, **{f"c_{k}": cres[k] for k in cs.CRES}, **lin))
    return out


def referencias():
    """BRS 284; apêndice G por cultivar; apêndice G fora da amostra (mediana das outras cultivares)."""
    g = pd.read_csv(cs.EV / "coeficientes_tese_apG.csv").set_index("coef")
    apg = {c: {k: float(g.loc[k.upper(), c]) for k in cs.CRES} for c in g.columns}
    nome = {"TMG 7062": "TMG 7062", "TMG 7063": "TMG 7063", "NS 7901": "NS7901", "65i65RSF": "65i65RSF",
            "8579RSF": "8579RSF", "BRS 399": "BRS 399"}
    refs = {"BRS 284": cs.BRS284}
    for cult in set(T.cultivar):
        c = nome[cult]
        refs[f"apG {c}"] = apg[c]
        outras = [x for x in apg if x != c]
        refs[f"apG fora {c}"] = {k: float(np.median([apg[x][k] for x in outras])) for k in cs.CRES}
    return refs


def amostra(n, proc):
    arq = OUT / "amostra.csv"
    if arq.exists():
        return pd.read_csv(arq)
    lo = np.array([cs.FAIXAS[k][0] for k in cs.CRES])
    hi = np.array([cs.FAIXAS[k][1] for k in cs.CRES])
    X = qmc.scale(qmc.LatinHypercube(d=len(cs.CRES), seed=20261108).random(n), lo, hi)
    todos = list(range(len(T)))
    tarefas = [(r, c, todos) for r, c in referencias().items()]
    tarefas += [(f"lhs{j:05d}", dict(zip(cs.CRES, x)), todos) for j, x in enumerate(X)]
    L = []
    with get_context("fork").Pool(proc, initializer=cs._ini) as pool:
        for k, r in enumerate(pool.imap_unordered(_tarefa, tarefas, chunksize=4), 1):
            L += r
            if k % 200 == 0:
                print(f"{k}/{len(tarefas)} conjuntos", flush=True)
    d = pd.DataFrame(L)
    d.to_csv(arq, index=False)
    return d


def objetivo(d, exps):
    """J por conjunto, nos tratamentos I-1000N dos experimentos `exps` (linhas de amostra.csv ou de uma simulação)."""
    x = d[(d.tratamento == "I-1000N") & d.experimento.isin(exps)].merge(
        T[["experimento", "tratamento", "yield_obs_kgha"]], on=["experimento", "tratamento"])
    x = x.assign(ey=((x.hwam - x.yield_obs_kgha) / x.yield_obs_kgha / SIG_Y) ** 2,
                 es=(x[[f"nrmse_{v}" for v in SERIES]] / SIG_S).pow(2).mean(axis=1))
    g = x.groupby("conjunto")
    return (g.ey.mean() + g.es.mean()).fillna(np.inf)


def _norm(c):
    return np.array([(c[k] - cs.FAIXAS[k][0]) / (cs.FAIXAS[k][1] - cs.FAIXAS[k][0]) for k in cs.CRES])


def _desnorm(z):
    z = np.clip(z, 0, 1)
    return {k: cs.FAIXAS[k][0] + z[i] * (cs.FAIXAS[k][1] - cs.FAIXAS[k][0]) for i, k in enumerate(cs.CRES)}


def refinar(args):
    """Nelder-Mead (espaço normalizado, recortado em [0, 1]) a partir de um conjunto, nos experimentos `exps`."""
    c0, exps = args
    idx = [i for i in range(len(T)) if T.experimento[i] in exps and T.tratamento[i] == "I-1000N"]

    def f(z):
        c = _desnorm(z)
        d = pd.DataFrame([dict(conjunto="x", **metricas(T.iloc[i], c)) for i in idx])
        return float(objetivo(d, exps).iloc[0])

    r = minimize(f, _norm(c0), method="Nelder-Mead", options=dict(maxfev=250, xatol=1e-3, fatol=1e-3,
                                                                  initial_simplex=None))
    return r.fun, _desnorm(r.x)


def calibrar(A, proc):
    arq = OUT / "coeficientes.csv"
    if arq.exists():
        return pd.read_csv(arq)
    lhs = A[A.conjunto.str.startswith("lhs")]
    dobras = {"completo": CAL, **{f"sem {e}": [x for x in CAL if x != e] for e in CAL}}
    trabalhos, chaves = [], []
    for nome, exps in dobras.items():
        J = objetivo(lhs, exps).nsmallest(3)
        for cj in J.index:
            c0 = lhs[lhs.conjunto == cj].iloc[0]
            trabalhos.append(({k: c0[f"c_{k}"] for k in cs.CRES}, exps))
            chaves.append((nome, cj, J[cj]))
    with get_context("fork").Pool(proc, initializer=cs._ini) as pool:
        res = pool.map(refinar, trabalhos)
    L = [dict(dobra=n, partida=cj, J_amostra=j0, J=j, **c) for (n, cj, j0), (j, c) in zip(chaves, res)]
    d = pd.DataFrame(L)
    d.to_csv(OUT / "refinamento.csv", index=False)
    best = d.loc[d.groupby("dobra").J.idxmin()].reset_index(drop=True)
    best.to_csv(arq, index=False)
    return best


def avaliar(coefs, proc):
    """Simula as versões (ajuste por dobra, BRS 284, apêndice G) em todos os tratamentos 1000N."""
    refs = referencias()
    cf = {r.dobra: {k: r[k] for k in cs.CRES} for _, r in coefs.iterrows()}
    tarefas = [(n, c, list(range(len(T)))) for n, c in {**{f"silva {k}": v for k, v in cf.items()}, **refs}.items()]
    with get_context("fork").Pool(proc, initializer=cs._ini) as pool:
        L = sum(pool.map(_tarefa, tarefas), [])
    d = pd.DataFrame(L).merge(T[["experimento", "tratamento", "agua", "yield_obs_kgha", "cultivar"]],
                              on=["experimento", "tratamento"])
    s = cs.series_obs()
    s = s[s.n == "1000N"]
    fim = s[s.variavel == "biomassa_kgha"].sort_values("das").groupby(["experimento", "tratamento"]).valor.last()
    lai = s[s.variavel == "iaf"].groupby(["experimento", "tratamento"]).valor.max()
    d = d.merge(fim.rename("cwam_obs"), on=["experimento", "tratamento"], how="left") \
         .merge(lai.rename("laix_obs"), on=["experimento", "tratamento"], how="left")
    d.to_csv(OUT / "avaliacao.csv", index=False)
    return d


def fora_da_amostra(d):
    """Uma linha por tratamento e versão comparada: 'silva' (dobra sem o experimento, ou completo se não calibrado),
    'BRS 284', 'apG fora' (mediana das outras cultivares) e 'apG própria' (referência dentro da amostra)."""
    nome = {"NS 7901": "NS7901"}
    L = []
    for (e, tr), g in d.groupby(["experimento", "tratamento"]):
        cult = nome.get(g.cultivar.iloc[0], g.cultivar.iloc[0])
        escolha = {"silva": f"silva sem {e}" if e in CAL else "silva completo", "BRS 284": "BRS 284",
                   "apG fora": f"apG fora {cult}", "apG própria": f"apG {cult}"}
        for v, cj in escolha.items():
            L.append(g[g.conjunto == cj].assign(versao=v))
    return pd.concat(L)


def criterios(f):
    def est(x, sim, obs):
        x = x.dropna(subset=[sim, obs])
        e = x[sim] - x[obs]
        return len(x), float(np.sqrt((e ** 2).mean())), float(e.mean())
    L = []
    for v, g in f.groupby("versao"):
        for nome, sel, sim, obs in [
                ("G1 produtividade I-1000N", g.tratamento == "I-1000N", "hwam", "yield_obs_kgha"),
                ("G2 produtividade R-1000N", g.tratamento == "R-1000N", "hwam", "yield_obs_kgha"),
                ("G3 biomassa final", g.cwam_obs.notna(), "cwam", "cwam_obs"),
                ("G3 IAF máximo", g.laix_obs.notna(), "laix", "laix_obs")]:
            n, rmse, vies = est(g[sel], sim, obs)
            L.append(dict(criterio=nome, versao=v, n=n, rmse=rmse, vies=vies))
    return pd.DataFrame(L)


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
    proc = int(sys.argv[2]) if len(sys.argv) > 2 else os.cpu_count()
    OUT.mkdir(parents=True, exist_ok=True)
    cs.theta()
    preparar()
    A = amostra(n, proc)
    print("amostra:", A.conjunto.nunique(), "conjuntos", flush=True)
    C = calibrar(A, proc)
    print(C.round(3).to_string(index=False), flush=True)
    D = avaliar(C, proc)
    F = fora_da_amostra(D)
    F.to_csv(OUT / "fora_da_amostra.csv", index=False)
    G = criterios(F)
    G.to_csv(OUT / "criterios_G.csv", index=False)
    pd.set_option("display.width", 200)
    print(G.round(2).to_string(index=False))
    print(F[["experimento", "tratamento", "versao", "yield_obs_kgha", "hwam", "cwam_obs", "cwam", "laix_obs", "laix"]]
          .round(0).to_string(index=False))
