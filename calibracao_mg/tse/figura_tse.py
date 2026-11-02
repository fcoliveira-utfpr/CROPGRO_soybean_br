"""Figura da comparação TSE × MCMC × GLUE (resultados de comparar_tse.py), no estilo das figuras dos artigos.

(a) RMSE por tipo de dado e estádio: experimentos no DSSAT real (o terreno do TSE) e médias multilocais e catálogos
    pelo emulador (fora da amostra do TSE).
(b) Coeficientes por GMR: retas da MCMC (mediana e IC 90%), ótimos do TSE por classe e reta ajustada a eles.
(c) Identificabilidade na grade do TSE: fração da faixa de cada coeficiente coberta pelas combinações a até 1 dia
    do critério mínimo.

Saída: resultados/figuras_en/10_5_tse_comparacao.png/.pdf
Uso: FIG_LANG=en DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1 .venv/bin/python tse/figura_tse.py
"""
import sys
from pathlib import Path

import emcee
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parent
sys.path.insert(0, str(RAIZ / "src"))
import calibrar_fenologia as cf  # noqa: E402
import figuras as fg  # noqa: E402

OUT = RAIZ / "resultados" / "tse"
CORES = {"MCMC v13": fg.CAT[0], "GLUE (melhor de 60.000)": fg.CAT[1], "TSE grade": fg.CAT[2],
         "TSE redução de faixa": fg.CAT[3]}
ROT = {"MCMC v13": "MCMC (this study)", "GLUE (melhor de 60.000)": "GLUE (best of 60,000)", "TSE grade": "TSE, full grid",
       "TSE redução de faixa": "TSE, range reduction"}
COEFS = [("csdl", "CSDL", "CSDL (h)"), ("ppsen", "PPSEN", "PPSEN (h⁻¹)"), ("em-fl", "EM-FL", "EM-FL (PTD)"),
         ("fl-sd", "FL-SD", "FL-SD (PTD)"), ("sd-pm", "SD-PM", "SD-PM (PTD)"), ("r1ppo", "R1PPO", "R1PPO (h)")]


def metodo_base(m):
    return m.replace(" (reta)", "")


if __name__ == "__main__":
    fg.estilo()
    e = pd.read_csv(OUT / "comparacao_erros.csv")
    coef = pd.read_csv(OUT / "coeficientes_metodos.csv")
    ident = pd.read_csv(OUT / "identificabilidade.csv")
    ch = np.concatenate([emcee.backends.HDFBackend(str(RAIZ / "resultados" / "emulador_fenologia" / f"cadeia_v13_{c}.h5"),
                                                   read_only=True).get_chain()[4500::50].reshape(-1, len(cf.NOMES)) for c in "ab"])
    ix = {n: i for i, n in enumerate(cf.NOMES)}
    fig = plt.figure(figsize=(14, 8.6))
    gs = fig.add_gridspec(2, 6, height_ratios=[1, 1.05], hspace=0.42, wspace=0.8)
    # (a) erros
    ax = fig.add_subplot(gs[0, :4])
    grupos = [("experimento (DSSAT)", "r1", "Exp. R1"), ("experimento (DSSAT)", "r5", "Exp. R5"),
              ("experimento (DSSAT)", "r7", "Exp. R7/R8"), ("media_multilocal", "r1", "Multi-site R1"),
              ("media_multilocal", "r8", "Multi-site R8"), ("catalogo", "r8", "Catalogue cycle")]
    ms = list(CORES)
    w = 0.8 / len(ms)
    for j, m in enumerate(ms):
        v = []
        for dado, est, _ in grupos:
            s = e[(e.dado == dado) & (e.estadio == est) & (e.metodo.map(metodo_base) == m)]
            v.append(s.rmse.iloc[0] if len(s) else np.nan)
        x = np.arange(len(grupos)) + (j - (len(ms) - 1) / 2) * w
        ax.bar(x, v, w * 0.92, color=CORES[m], label=ROT[m])
    ax.set_xticks(range(len(grupos)))
    ax.set_xticklabels([g[2] for g in grupos])
    ax.axvline(2.5, color=fg.GRADE, lw=1)
    ax.text(1, ax.get_ylim()[1] * 0.97, "in-sample for TSE (real DSSAT)", ha="center", va="top", fontsize=8, color=fg.TEXTO_2)
    ax.text(4, ax.get_ylim()[1] * 0.97, "out-of-sample for TSE (emulator)", ha="center", va="top", fontsize=8, color=fg.TEXTO_2)
    ax.set_ylabel("RMSE (days)")
    ax.set_title("(a) Error by data type and stage", loc="left")
    ax.legend(ncol=2, fontsize=8, loc="upper left", bbox_to_anchor=(0, 0.9))
    # (c) identificabilidade
    ax = fig.add_subplot(gs[0, 4:])
    nomes = ["CSDL", "PPSEN", "EM-FL", "FL-SD", "SD-PM"]
    M = ident[[f"amplitude_{n}" for n in nomes]].values
    im = ax.imshow(M, aspect="auto", cmap=fg.CMAP_SEQ, vmin=0, vmax=1)
    ax.set_xticks(range(len(nomes)))
    ax.set_xticklabels(nomes, rotation=30)
    ax.set_yticks(range(len(ident)))
    ax.set_yticklabels([f"{c:.1f} · {'R1' if c >= 8.5 else 'R7'} · {k}" for c, k in zip(ident.classe, ident.n_ate_1d)],
                       fontsize=8)  # evento usado pelo TSE: o observado em todos os tratamentos
    for (i, j), v in np.ndenumerate(M):
        ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7, color="white" if v > 0.55 else fg.TEXTO)
    ax.grid(False)
    ax.set_title("(c) TSE grid: range spanned by combinations\nwithin 1 d of the optimum\n(rows: RMG · event used · n)", loc="left", fontsize=10)
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
    cb.set_label("fraction of the range", color=fg.TEXTO_2)
    cb.outline.set_visible(False)
    # (b) coeficientes por GMR
    g = np.linspace(5, 9.5, 46)
    for k, (c, C, rot) in enumerate(COEFS):
        ax = fig.add_subplot(gs[1, k])
        v = np.clip(ch[:, ix[f"a_{c}"]][:, None] + ch[:, ix[f"b_{c}"]][:, None] * (g - 7), *cf.FAIXAS[c])
        ax.fill_between(g, *np.percentile(v, [5, 95], axis=0), color=fg.SEQ[2], lw=0)
        ax.plot(g, np.median(v, axis=0), color=CORES["MCMC v13"], lw=2)
        for m in ("TSE grade", "TSE redução de faixa", "GLUE (melhor de 60.000)"):
            s = coef[coef.metodo == m]
            if len(s) and C in s:
                ax.plot(s.classe, s[C], "o" if "TSE" in m else "-", color=CORES[m], ms=4.5, lw=1.5,
                        mec="white", mew=0.6)
        lo, hi = cf.FAIXAS[c]
        for y in (lo, hi):
            ax.axhline(y, color=fg.TEXTO_2, lw=0.8, ls="--")
        ax.set_ylim(lo - 0.05 * (hi - lo), hi + 0.05 * (hi - lo))
        ax.set_xlabel("RMG")
        ax.set_title(rot, loc="left", fontsize=9.5)
    fig.text(0.07, 0.47, "(b) Coefficients by RMG: MCMC median and 90% interval (band), GLUE best (line), "
                         "TSE optima per RMG class (points); dashed = calibration range", fontsize=9.5, fontweight="bold")
    fg.salvar(fig, "10_5_tse_comparacao")
