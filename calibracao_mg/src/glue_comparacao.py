"""Comparação com o GLUE (29/10/2026): amostragem independente uniforme × MCMC na mesma posterior.

O GLUE do DSSAT (Beven e Binley, 1992; He et al., 2010) sorteia conjuntos de coeficientes de forma independente e
uniforme nas faixas, roda o modelo para cada um e pondera os conjuntos pela verossimilhança. Aqui o método é reproduzido
com o emulador e a mesma verossimilhança hierárquica da calibração (v13: δ8 fixo, prior da literatura, faixa
fisiológica): os 12 coeficientes (a_k, b_k) variam; as escalas de erro e τ ficam fixas na mediana da posterior, como o
GLUE fixa a variância do erro. Para N = 6.000 (padrão do DSSAT), 60.000 e 240.000 sorteios: tamanho efetivo da amostra
(ESS) dos pesos, melhor log-posterior encontrado e mediana/intervalo de 90% ponderados de θ(GMR), comparados com a MCMC.

Saídas: resultados/glue/glue_resumo.csv, glue_coeficientes.csv e glue_amostras_<N>.npz.
Uso: DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1 .venv/bin/python src/glue_comparacao.py [N_max]
"""
import sys
import time
import warnings
from pathlib import Path

import emcee
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))
import calibrar_fenologia as cf  # noqa: E402

RAIZ = Path(__file__).resolve().parents[1]
OUT = RAIZ / "resultados" / "glue"
CADEIAS = [RAIZ / "resultados" / "emulador_fenologia" / f"cadeia_v13_{c}.h5" for c in "ab"]
COEF = [n for n in cf.NOMES if n[:2] in ("a_", "b_")]
GS = (5, 7, 9)


def wquant(x, w, q):
    i = np.argsort(x)
    cw = np.cumsum(w[i]) / w.sum()
    return np.interp(q, cw, x[i])


def theta(TH, ix, k, g):
    return np.clip(TH[:, ix[f"a_{k}"]] + TH[:, ix[f"b_{k}"]] * (g - 7), *cf.FAIXAS[k])


if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 240_000
    OUT.mkdir(parents=True, exist_ok=True)
    ix = {n: i for i, n in enumerate(cf.NOMES)}
    ch = np.concatenate([emcee.backends.HDFBackend(str(f), read_only=True).get_chain()[4500:].reshape(-1, len(cf.NOMES))
                         for f in CADEIAS])
    lp_mcmc = np.concatenate([emcee.backends.HDFBackend(str(f), read_only=True).get_log_prob()[4500:].ravel() for f in CADEIAS])
    med = np.median(ch, axis=0)
    c = cf.Calibracao()
    lo = np.array([p[1] for p in cf.PARAMS])
    hi = np.array([p[2] for p in cf.PARAMS])
    rng = np.random.default_rng(42)
    TH = np.tile(med, (nmax, 1))
    for n in COEF:
        TH[:, ix[n]] = rng.uniform(lo[ix[n]], hi[ix[n]], nmax)
    lp = np.empty(nmax)
    t0 = time.time()
    lote = 40                                     # como na MCMC (40 caminhantes): lotes grandes estouram a memória
    for i in range(0, nmax, lote):
        lp[i:i + lote] = c.log_post_vet(TH[i:i + lote])
        if (i // lote) % 250 == 0:
            print(f"{i + lote}/{nmax}  {time.time() - t0:.0f} s  melhor {np.nanmax(lp[:i + lote]):.1f}", flush=True)
    dt = time.time() - t0
    np.savez_compressed(OUT / f"glue_amostras_{nmax}.npz", theta=TH[:, [ix[n] for n in COEF]], logpost=lp, nomes=np.array(COEF))
    lp_mcmc_med = c.log_post(med)
    resumo, coefs = [], []
    # referência: MCMC (mesmas escalas de erro fixas não se aplicam; usa-se a posterior completa da v13)
    for k in cf.FAIXAS:
        for g in GS:
            v = theta(ch, ix, k, g)
            coefs.append(dict(metodo="MCMC", N=len(ch), coef=k, rmg=g, q05=np.percentile(v, 5), q50=np.median(v), q95=np.percentile(v, 95)))
    for N in sorted({min(6_000, nmax), min(60_000, nmax), nmax}):
        l = lp[:N]
        ok = np.isfinite(l)
        w = np.zeros(N)
        w[ok] = np.exp(l[ok] - l[ok].max())
        ess = w.sum() ** 2 / (w ** 2).sum()
        resumo.append(dict(N=N, ess=ess, melhor_logpost=l[ok].max(), logpost_mediana_mcmc=lp_mcmc_med,
                           logpost_max_mcmc=lp_mcmc.max(), tempo_s=dt * N / nmax,
                           frac_peso_maior=w.max() / w.sum()))
        for k in cf.FAIXAS:
            for g in GS:
                v = theta(TH[:N], ix, k, g)
                coefs.append(dict(metodo="GLUE", N=N, coef=k, rmg=g, q05=wquant(v, w, 0.05), q50=wquant(v, w, 0.5),
                                  q95=wquant(v, w, 0.95)))
    r = pd.DataFrame(resumo)
    r.to_csv(OUT / "glue_resumo.csv", index=False)
    cc = pd.DataFrame(coefs)
    cc.to_csv(OUT / "glue_coeficientes.csv", index=False)
    pd.set_option("display.width", 200)
    print(r.round(3).to_string(index=False))
    print(cc[cc.rmg == 7].pivot_table(index="coef", columns=["metodo", "N"], values="q50").round(3).to_string())
