"""Simulações para a inversão amortizada (SBI) da fenologia por satélite → GMR efetivo.

Para cada ponto × safra com fenologia MODIS, os desconhecidos são:
  g  — GMR efetivo da lavoura (5–9,5);
  s  — data de semeadura (dias após 10/set, 0–171);
  a  — defasagem entre a emergência e o início do verde (SOS, 20% da amplitude do NDVI), dias;
  b  — defasagem entre a maturação fisiológica (R7) e a queda de 50% do NDVI (EOS50), dias.
Os coeficientes θ(GMR) vêm de amostras da posterior da calibração (a incerteza deles é marginalizada).

Observação simulada (dias desde 1/jul da safra, como em fenologia_modis_pontos.csv):
  SOS   = (semeadura − 1/jul) + E(θ, g, s, clima) + a + ε
  EOS50 = (semeadura − 1/jul) + R7(θ, g, s, clima) + b + ε,     ε ~ N(0, σ_obs)
com E e R7 dados pelo emulador. O contexto (sem depender dos desconhecidos) é latitude, altitude e médias mensais de
temperatura média e mínima de setembro a abril da safra.

Saída: resultados/sbi/simulacoes_<nome>.npz (theta, x, contexto) e resultados/sbi/observado.npz (x e contexto reais).
Uso: .venv/bin/python src/sbi_simulador.py N CADEIA.h5 [nome]
"""
import sys
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
import treino_fenologia as tf  # noqa: E402

RAIZ = tf.RAIZ
DIR = RAIZ / "resultados" / "sbi"
PRIOR = {"g": (5.0, 9.5), "s": (0.0, 171.0)}
A_MED, A_DP = 10.0, 4.0          # emergência → SOS (verde visível), dias
B_MED, B_DP = 0.0, 6.0           # R7 → EOS50, dias
SIGMA_OBS = 3.0                  # erro da mediana de ~150 pixels, dias
MESES = [9, 10, 11, 12, 1, 2, 3, 4]


def locais_satelite():
    sat = pd.read_csv(RAIZ / "resultados" / "fenologia_satelite" / "fenologia_modis_pontos.csv").dropna(subset=["sos_p50", "eos50_p50"])
    pts = pd.read_csv(RAIZ / "dados" / "processados" / "pontos_bs2019_digitalizados.csv")
    pts = pts[["lon", "lat", "elev"]].drop_duplicates(["lon", "lat"])
    loc = sat[["lon", "lat"]].drop_duplicates().merge(pts, on=["lon", "lat"], how="left").reset_index(drop=True)
    loc["local"], loc["arq"] = "sat", "no"
    loc.index.name = "idx_local"
    sat = sat.merge(loc.reset_index()[["idx_local", "lon", "lat"]], on=["lon", "lat"])
    sat["safra"] = sat.ano_colheita - 1
    return sat, loc


def contexto(sat, loc):
    """lat, elev e médias mensais (set–abr) de Tmed e Tmin da safra: 2 + 16 colunas, por linha de sat."""
    partes = []
    for i, g in sat.groupby("idx_local"):
        c = pd.read_csv(tf.arq_clima(loc.loc[i]), index_col=0, parse_dates=True)
        c["tmed"] = (c.Tmax + c.Tmin) / 2
        c["safra"] = c.index.year - (c.index.month < 7)
        m = c[c.index.month.isin(MESES)].groupby(["safra", c[c.index.month.isin(MESES)].index.month])[["tmed", "Tmin"]].mean().unstack()
        m.columns = [f"{v}_{mes:02d}" for v, mes in m.columns]
        v = g[["safra"]].join(m, on="safra").drop(columns="safra")
        v.insert(0, "elev", loc.loc[i, "elev"])
        v.insert(0, "lat", loc.loc[i, "lat"])
        partes.append(v)
    return pd.concat(partes).loc[sat.index]


def amostras_theta(cadeia, n, rng, descarte=0.5):
    ch = emcee.backends.HDFBackend(str(cadeia), read_only=True).get_chain()
    ch = ch[int(len(ch) * descarte):].reshape(-1, ch.shape[2])
    # SBI final must use the current calibration layout, including a_r1ppo and b_r1ppo.
    # Silently padding an older chain would replace the calibrated R1PPO with its prior start value.
    if ch.shape[1] != len(cf.PARAMS):
        raise ValueError(
            f"cadeia tem {ch.shape[1]} parâmetros; esperado {len(cf.PARAMS)} "
            f"({', '.join(cf.NOMES)})"
        )
    necessarios = {"a_r1ppo", "b_r1ppo"}
    if not necessarios.issubset(cf.NOMES):
        raise RuntimeError(f"calibração não contém parâmetros R1PPO: {necessarios - set(cf.NOMES)}")
    th = ch[rng.integers(0, len(ch), n)]
    return th


def simular(sat, loc, ctx, theta_post, rng):
    n = len(theta_post)
    linhas = rng.integers(0, len(sat), n)
    s0 = sat.iloc[linhas].reset_index(drop=True)
    g = rng.uniform(*PRIOR["g"], n)
    s = rng.uniform(*PRIOR["s"], n)
    a = rng.normal(A_MED, A_DP, n)
    b = rng.normal(B_MED, B_DP, n)
    amb = ef.variaveis_ambiente(pd.DataFrame({"idx_local": s0.idx_local, "safra": s0.safra, "doy_sem": np.floor(s)}), loc)
    X = amb.copy()
    X["doy_sem"], X["gmr"], X["mg_int"] = np.floor(s), g, np.clip(np.round(g), 5, 9)
    ix = {nome: i for i, nome in enumerate(cf.NOMES)}
    # Includes R1PPO(GMR): the v9 posterior values feed the R7 emulator through ef.colunas().
    for k, (lo, hi) in cf.FAIXAS.items():
        X[k] = np.clip(theta_post[:, ix[f"a_{k}"]] + theta_post[:, ix[f"b_{k}"]] * (g - 7), lo, hi)
    if "r1ppo" not in ef.colunas():
        raise RuntimeError("o emulador de fenologia não aceita R1PPO como entrada")
    X = X[ef.colunas()]
    mod = {k: joblib.load(ef.DIR / f"modelo_{k}.joblib") for k in ("e", "r7")}
    e, r7 = mod["e"].predict(X), mod["r7"].predict(X)
    sem_jul1 = 71 + np.floor(s)                       # 10/set = 71 dias após 1/jul
    sos = sem_jul1 + e + a + rng.normal(0, SIGMA_OBS, n)
    eos = sem_jul1 + r7 + b + rng.normal(0, SIGMA_OBS, n)
    theta = np.column_stack([g, s, a, b])
    x = np.column_stack([sos, eos])
    return theta, x, ctx.iloc[linhas].to_numpy(dtype=float)


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 100_000
    cadeia = Path(sys.argv[2])
    nome = sys.argv[3] if len(sys.argv) > 3 else "base"
    rng = np.random.default_rng(11)
    DIR.mkdir(parents=True, exist_ok=True)
    sat, loc = locais_satelite()
    ctx = contexto(sat, loc)
    ok = ctx.notna().all(axis=1)
    sat, ctx = sat[ok].reset_index(drop=True), ctx[ok].reset_index(drop=True)
    np.savez(DIR / "observado.npz", x=sat[["sos_p50", "eos50_p50"]].to_numpy(float), contexto=ctx.to_numpy(float),
             lon=sat.lon.values, lat=sat.lat.values, safra=sat.safra.values, colunas_contexto=np.array(ctx.columns))
    theta, x, c = simular(sat, loc, ctx, amostras_theta(cadeia, n, rng), rng)
    np.savez(DIR / f"simulacoes_{nome}.npz", theta=theta, x=x, contexto=c, nomes_theta=np.array(["g", "s", "a", "b"]))
    print(f"{n} simulações; {len(sat)} ponto-safra observados; contexto {ctx.shape[1]} colunas", flush=True)
