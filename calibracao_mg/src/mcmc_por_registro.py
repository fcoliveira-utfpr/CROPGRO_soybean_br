"""MCMC por registro para a inversão do satélite: custo e concordância com o NPE (22/10/2026).

Para cada ponto-safra MODIS sorteado, amostra a posterior de (g, s, a, b) com o mesmo modelo de observação do SBI
(src/sbi_simulador.py): SOS = 71 + s + E + a + ε, EOS50 = 71 + s + R7 + b + ε, ε ~ N(0, 3²), com as mesmas priors.
A incerteza dos coeficientes θ(GMR) é marginalizada como no NPE: a verossimilhança é a média sobre K sorteios fixos
da posterior da calibração (Monte Carlo com números aleatórios comuns). As variáveis de ambiente das 172 datas de
semeadura possíveis são pré-calculadas uma vez por registro (o emulador recebe a semeadura em dias inteiros).

Saída: resultados/sbi/mcmc_por_registro.csv (tempo, ESS, médias e desvios de g e s pela MCMC e pelo NPE).
Uso: DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1 .venv/bin/python src/mcmc_por_registro.py [n_registros] [passos] [K]
"""
import sys
import time
import warnings
from pathlib import Path

import emcee
import joblib
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))
import calibrar_fenologia as cf  # noqa: E402
import emulador_fenologia as ef  # noqa: E402
import sbi_simulador as ss  # noqa: E402

CADEIA = ss.RAIZ / "resultados" / "emulador_fenologia" / "v13_3000passos" / "cadeia_v13_a.h5"
NW = 32
DESCARTE = None      # um terço dos passos (definido no main)


def preparar(rec, loc, phi):
    """Tabela de entradas do emulador para as 172 semeaduras × K sorteios (coeficientes recalculados por g)."""
    amb = ef.variaveis_ambiente(pd.DataFrame({"idx_local": rec.idx_local, "safra": rec.safra, "doy_sem": np.arange(172)}), loc)
    return amb.reset_index(drop=True)


def log_post_fn(rec, amb, phi, mod):
    ix = {n: i for i, n in enumerate(cf.NOMES)}
    K = len(phi)
    x = np.array([rec.sos_p50, rec.eos50_p50])
    lo, hi = np.array([5.0, 0.0]), np.array([9.5, 171.0])

    def lp(TH):                          # TH: (nw, 4) = g, s, a, b
        out = np.full(len(TH), -np.inf)
        ok = np.all((TH[:, :2] >= lo) & (TH[:, :2] <= hi), axis=1)
        if not ok.any():
            return out
        T = TH[ok]
        n = len(T)
        g = np.repeat(T[:, 0], K)
        sd = np.repeat(np.floor(T[:, 1]).astype(int).clip(0, 171), K)
        X = amb.iloc[sd].reset_index(drop=True).copy()
        X["doy_sem"], X["gmr"], X["mg_int"] = sd, g, np.clip(np.round(g), 5, 9)
        P = np.tile(phi, (n, 1))
        for k, (a, b) in cf.FAIXAS.items():
            X[k] = np.clip(P[:, ix[f"a_{k}"]] + P[:, ix[f"b_{k}"]] * (g - 7), a, b)
        X = X[ef.colunas()]
        e, r7 = mod["e"].predict(X).reshape(n, K), mod["r7"].predict(X).reshape(n, K)
        base = 71 + np.floor(T[:, 1])[:, None]
        r1 = x[0] - (base + e + T[:, 2:3])
        r2 = x[1] - (base + r7 + T[:, 3:4])
        ll = -0.5 * (r1 ** 2 + r2 ** 2) / ss.SIGMA_OBS ** 2
        m = ll.max(1, keepdims=True)
        ll = (m + np.log(np.exp(ll - m).mean(1, keepdims=True)))[:, 0]           # média sobre os K sorteios
        lpri = -0.5 * ((T[:, 2] - ss.A_MED) / ss.A_DP) ** 2 - 0.5 * ((T[:, 3] - ss.B_MED) / ss.B_DP) ** 2
        out[ok] = ll + lpri
        return out
    return lp


if __name__ == "__main__":
    n_reg = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    passos = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
    K = int(sys.argv[3]) if len(sys.argv) > 3 else 50
    DESCARTE = passos // 3
    rng = np.random.default_rng(5)
    sat, loc = ss.locais_satelite()
    npe = pd.read_csv(ss.DIR / "gmr_efetivo_final.csv")
    ch = emcee.backends.HDFBackend(str(CADEIA), read_only=True).get_chain()
    ch = ch[len(ch) // 2:].reshape(-1, ch.shape[2])
    phi = ch[rng.integers(0, len(ch), K)]
    mod = {k: joblib.load(ef.DIR / f"modelo_{k}.joblib") for k in ("e", "r7")}
    linhas = []
    for j in rng.choice(len(sat), n_reg, replace=False):
        rec = sat.iloc[j]
        t0 = time.time()
        amb = preparar(rec, loc, phi)
        lp = log_post_fn(rec, amb, phi, mod)
        p0 = np.column_stack([rng.uniform(5, 9.5, NW), rng.uniform(0, 171, NW),
                              rng.normal(ss.A_MED, ss.A_DP, NW), rng.normal(ss.B_MED, ss.B_DP, NW)])
        sam = emcee.EnsembleSampler(NW, 4, lp, vectorize=True,
                                    moves=[(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2)])
        sam.run_mcmc(p0, passos, progress=False)
        dt = time.time() - t0
        c = sam.get_chain(discard=DESCARTE)
        tau = sam.get_autocorr_time(discard=DESCARTE, tol=0)
        ess = c.shape[0] * NW / np.nanmax(tau)
        f = c.reshape(-1, 4)
        ref = npe[(npe.lon == rec.lon) & (npe.lat == rec.lat) & (npe.safra == rec.safra)]
        linhas.append(dict(lon=rec.lon, lat=rec.lat, safra=rec.safra, tempo_s=dt, passos=passos, K=K, tau_max=np.nanmax(tau),
                           ess=ess, s_por_1000_ess=dt / ess * 1000, aceitacao=sam.acceptance_fraction.mean(),
                           g_mcmc=f[:, 0].mean(), g_dp_mcmc=f[:, 0].std(), s_mcmc=f[:, 1].mean(), s_dp_mcmc=f[:, 1].std(),
                           g_npe=ref.g_media.iloc[0] if len(ref) else np.nan, g_dp_npe=ref.g_dp.iloc[0] if len(ref) else np.nan,
                           s_npe=ref.s_media.iloc[0] if len(ref) else np.nan, s_dp_npe=ref.s_dp.iloc[0] if len(ref) else np.nan))
        print({k: round(v, 2) if isinstance(v, float) else v for k, v in linhas[-1].items()}, flush=True)
        pd.DataFrame(linhas).to_csv(ss.DIR / "mcmc_por_registro.csv", index=False)
    d = pd.DataFrame(linhas)
    print(d[["tempo_s", "ess", "s_por_1000_ess", "tau_max", "aceitacao"]].describe().round(2).to_string())
    print("g: MCMC × NPE  r =", round(d.g_mcmc.corr(d.g_npe), 3), " RMSE =", round(np.sqrt(((d.g_mcmc - d.g_npe) ** 2).mean()), 3))
    print("s: MCMC × NPE  r =", round(d.s_mcmc.corr(d.s_npe), 3), " RMSE =", round(np.sqrt(((d.s_mcmc - d.s_npe) ** 2).mean()), 3))
