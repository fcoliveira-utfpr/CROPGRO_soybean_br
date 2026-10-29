"""Figuras das curvas observadas × simuladas em Paragominas (cultivar padrão × calibrada) — 29/10/2026.
Usa resultados/curvas/ (de src/curvas_paragominas.py). Respeita FIG_LANG=en (figuras.py traduz e grava em figuras_en/).
Uso: .venv/bin/python src/figuras_curvas.py   (e FIG_LANG=en .venv/bin/python src/figuras_curvas.py)
"""
import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))
import figuras as fg  # noqa: E402

RAIZ = Path(__file__).resolve().parents[1]
D = RAIZ / "resultados" / "curvas"
fg.estilo()
COR = {"padrao": fg.CAT[1], "sampaio": fg.CAT[2], "br": fg.CAT[0]}
NOME = {"padrao": "padrão do DSSAT", "sampaio": "Sampaio et al. (calibrada no experimento)", "br": "nacional por GMR (este trabalho)"}
TRATS = {3: "GMR 7,7 · 30 pl m⁻² · 6 jan", 6: "GMR 8,8 · 20 pl m⁻² · 6 jan"}


def series():
    pg = pd.read_csv(D / "paragominas_series.csv")
    ob = pd.read_csv(D / "paragominas_obs_series.csv")
    ob = ob.replace(-99, np.nan)
    for c in ["CWAD", "LAID", "GWAD", "CDAY"]:
        ob[c] = pd.to_numeric(ob[c], errors="coerce")
    fig, axs = plt.subplots(3, 2, figsize=(11, 9), sharex=True)
    for j, (trno, titulo) in enumerate(TRATS.items()):
        o = ob[ob.TRNO == trno]
        for i, (var, rot, esc) in enumerate([("LAID", "IAF (m² m⁻²)", 1), ("CWAD", "biomassa (t ha⁻¹)", 1000),
                                             ("GWAD", "grãos (t ha⁻¹)", 1000)]):
            ax = axs[i, j]
            for v in ["padrao", "sampaio", "br"]:
                s = pg[(pg.TRNO == trno) & (pg.versao == v)]
                ax.plot(s.DAS - 5, s[var] / esc, color=COR[v], label=NOME[v])          # DAS desde a semeadura
            ax.plot(o.CDAY, o[var] / esc, "o", color=fg.TEXTO, ms=5, mfc="none", label="observado")
            ax.set_ylabel(rot)
            if i == 0:
                ax.set_title(titulo)
        axs[2, j].set_xlabel("dias após a semeadura")
    fig.legend(*axs[0, 0].get_legend_handles_labels(), loc="upper center", ncol=4, bbox_to_anchor=(0.5, 1.04),
               frameon=False, fontsize=9)
    fg.salvar(fig, "13_1_paragominas_series")


def dispersao():
    ev = pd.read_csv(D / "paragominas_avaliacao.csv")
    vars_ = [("ADAP", "florescimento R1 (dias)"), ("MDAP", "maturação R7 (dias)"), ("HWAM", "grãos (kg ha⁻¹)"),
             ("LAIX", "IAF máximo")]
    fig, axs = plt.subplots(1, 4, figsize=(15, 4))
    for ax, (v, rot) in zip(axs, vars_):
        lo, hi = np.inf, -np.inf
        for ver in ["padrao", "sampaio", "br"]:
            e = ev[ev.versao == ver]
            o, s = pd.to_numeric(e[v + "M"], errors="coerce"), pd.to_numeric(e[v + "S"], errors="coerce")
            rmse = np.sqrt(((s - o) ** 2).mean())
            ax.plot(o, s, "o", color=COR[ver], ms=5, label=f"{NOME[ver].split(' (')[0]} · RMSE {rmse:.1f}" if v != "HWAM" else f"{NOME[ver].split(' (')[0]} · RMSE {rmse:.0f}")
            lo, hi = min(lo, o.min(), s.min()), max(hi, o.max(), s.max())
        pad = (hi - lo) * 0.08
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color=fg.TEXTO_2, lw=1)
        ax.set_xlim(lo - pad, hi + pad); ax.set_ylim(lo - pad, hi + pad)
        ax.set_title(rot); ax.set_xlabel("observado"); ax.legend(fontsize=7, loc="upper left")
    axs[0].set_ylabel("simulado")
    fg.salvar(fig, "13_2_paragominas_dispersao")


if __name__ == "__main__":
    series()
    dispersao()
    print("figuras em", fg.DIR_FIG)
