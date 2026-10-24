"""Sensibilidade da calibração à escolha do emulador (22/10/2026).

Treina emuladores MLP (2 × 128, como em comparar_emuladores.py) para E, R1, R5 e R7 numa subamostra da tabela de
treino e avalia a posterior final (v13) com os dois emuladores: diferença das previsões nas observações, log-posterior
e erro por tipo de dado na mediana, e log-posterior de amostras da posterior.

Saída: resultados/emulador_fenologia/sensibilidade_emulador.txt
Uso: DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1 .venv/bin/python src/sensibilidade_emulador.py [n_linhas] [n_amostras]
"""
import sys
import warnings
from pathlib import Path

import emcee
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))
import calibrar_fenologia as cf  # noqa: E402
import comparar_emuladores as ce  # noqa: E402
import emulador_fenologia as ef  # noqa: E402
import treino_fenologia as tf  # noqa: E402


def erros(c, th):
    pr = c.prever(th)
    o = c.obs.copy()
    o["sim"] = [pr[e][i] for i, e in enumerate(o.estadio)]
    exp = o[o.tipo == "experimento"]
    out = {f"exp_{e}": np.sqrt(((g.dias_ve - g.sim) ** 2).mean()) for e, g in exp.groupby("estadio")}
    reg = o[o.tipo != "experimento"].groupby("id_registro").agg(y=("dias_ve", "first"), s=("sim", "mean"),
                                                                t=("tipo", "first"), e=("estadio", "first"))
    out.update({f"{t}_{e}": np.sqrt(((g.y - g.s) ** 2).mean()) for (t, e), g in reg.groupby(["t", "e"])})
    return pd.Series(out), o.sim.values


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 1_000_000
    n_am = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    d = pd.read_parquet(ef.DIR / "tabela_treino.parquet")
    if d.r1ppo.isna().any():
        d["r1ppo"] = d.r1ppo.fillna(pd.Series(tf.r1ppo_padrao(d.gmr.values), index=d.index))
    d = ef.preparar(d.sample(n, random_state=1).reset_index(drop=True))
    X = d[ef.colunas()]
    mlp = {}
    for s in ef.SAIDAS:
        ok = d[s].notna() & X.notna().all(axis=1)
        mlp[s] = ce.MODELOS["mlp_2x128"]().fit(X[ok], d.loc[ok, s])
        print("MLP treinada:", s, int(ok.sum()), "linhas", flush=True)
    ch = np.concatenate([emcee.backends.HDFBackend(str(ef.DIR / f"v13_3000passos/cadeia_v13_{k}.h5"), read_only=True)
                         .get_chain()[1500:].reshape(-1, len(cf.PARAMS)) for k in "ab"])
    med = np.median(ch, axis=0)
    am = ch[np.random.default_rng(0).integers(0, len(ch), n_am)]
    c = cf.Calibracao()
    linhas = []
    res = {}
    for nome, mods in [("gradient_boosting", dict(c.modelos)), ("mlp", mlp)]:
        c.modelos = mods
        e, sim = erros(c, med)
        lp = c.log_post_vet(am)
        res[nome] = sim
        linhas.append(pd.concat([pd.Series({"emulador": nome, "logpost_mediana": float(c.log_post(med)),
                                            "logpost_amostras_media": float(np.mean(lp))}), e.round(3)]))
    t = pd.DataFrame(linhas)
    dif = res["mlp"] - res["gradient_boosting"]
    o = c.obs
    txt = [t.to_string(index=False), "",
           "Diferença das previsões (MLP − GB) nas observações, na mediana da posterior (dias):",
           pd.DataFrame({"estadio": o.estadio, "tipo": o.tipo, "dif": dif}).groupby(["tipo", "estadio"]).dif
           .agg(["mean", "std", lambda x: np.sqrt((x ** 2).mean())]).rename(columns={"<lambda_0>": "rms"}).round(2).to_string()]
    open(ef.DIR / "sensibilidade_emulador.txt", "w").write("\n".join(txt))
    print("\n".join(txt))
