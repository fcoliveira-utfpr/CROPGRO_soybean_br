"""Figuras da validação independente com os experimentos de Silva (2021, tese ESALQ).

- 5_1_mapa_observacoes: observações da calibração v13 + os 5 locais de Silva (validação independente, fora da
  calibração). Lê os arquivos oficiais da v13 (observacoes_fenologia.csv e locais.csv); os da v14 (com Silva) estão em
  emulador_fenologia/emulador_v14/.
- 9_7_validacao_silva: observado × simulado (DSSAT real, mediana da v13) de R1, R3, R5 e R7 nos 9 experimentos de
  avaliação, por local, sobre os 209 tratamentos de calibração (cinza).

Uso: .venv/bin/python src/figuras_validacao_silva.py   (FIG_LANG=en para a versão em inglês, em figuras_en/)
Requer resultados/evandro/v13_fora_amostra.csv (validar_evandro.py) e avaliacoes_experimentos_detalhe.csv
(avaliar_experimentos.py, rótulo "v13").
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parent
sys.path.insert(0, str(AQUI))
import figuras as fg  # noqa: E402

DE = RAIZ / "resultados" / "emulador_fenologia"
EV = RAIZ / "resultados" / "evandro"
LOCAIS = ["Piracicaba-SP", "Sorriso-MT", "Cruz Alta-RS", "Tupanciretã-RS", "Teresina-PI"]   # ordem fixa das cores
FORA = ["TE-100", "TE-50", "SO-2"]


def silva():
    d = pd.read_csv(EV / "v13_fora_amostra.csv")
    return d[(d.cenario_gmr == "gmr_calibracao") & ~d.experimento.isin(FORA)]


def mapa_observacoes():
    CAT = fg.CAT
    L = pd.read_csv(DE / "locais.csv", index_col="idx_local")
    o = pd.read_csv(DE / "observacoes_fenologia.csv")
    sat = pd.read_csv(RAIZ / "resultados/fenologia_satelite/fenologia_modis_pontos.csv").dropna(subset=["sos_p50"])[
        ["lon", "lat"]].drop_duplicates()
    import validar_evandro as ve
    fig, ax = plt.subplots(figsize=(6.6, 6.6))
    fg.mapa_base(ax)
    ax.scatter(sat.lon, sat.lat, s=7, color="#b5b3ad", zorder=2, label="satélite (224 pontos)")
    for k, (tipo, nome, mk) in enumerate([("experimento", "experimentos", "o"), ("agregado", "Alliprandini et al. (2009)", "s"),
                                          ("catalogo", "catálogos", "^")]):
        ii = o[o.tipo == tipo].idx_local.unique()
        ax.scatter(L.loc[ii, "lon"], L.loc[ii, "lat"], marker=mk, s=46, color=CAT[k], edgecolor="white", linewidth=0.8,
                   zorder=4, label=nome)
    xy = np.array([ve.CELULA[l][:2] for l in LOCAIS])
    ax.scatter(xy[:, 0], xy[:, 1], marker="D", s=52, color=CAT[3], edgecolor=fg.TEXTO, linewidth=0.9, zorder=5,
               label="validação independente (Silva 2021)")
    ax.legend(loc="lower left")
    ax.set_title("Observações de fenologia usadas na calibração")
    fg.salvar(fig, "5_1_mapa_observacoes")


def observado_simulado():
    CAT = fg.CAT
    s = silva()
    cal = pd.read_csv(EV / "avaliacoes_experimentos_detalhe.csv")
    cal = cal[cal.metodo == "v13"]
    colcal = {"r1": ("ADAPM", "ADAPS"), "r5": ("PDFPM", "PDFPS"), "r7": ("MDAPM", "MDAPS")}
    fig, axs = plt.subplots(1, 4, figsize=(14, 4.0))
    for ax, e in zip(axs, ["r1", "r3", "r5", "r7"]):
        x, y = s[f"{e}_obs"].astype(float), s[f"{e}_sim"].astype(float)
        todos = [x, y]
        tit = f"{e.upper()}\nRMSE {np.sqrt(((y - x) ** 2).mean()):.1f} d"
        if e in colcal:
            co, cs = colcal[e]
            c = cal[[co, cs]].dropna()
            ax.scatter(c[co], c[cs], s=10, color="#c9c7c1", linewidth=0, zorder=2,
                       label="tratamentos da calibração" if e == "r1" else None)
            tit += f" · calibração {np.sqrt(((c[cs] - c[co]) ** 2).mean()):.1f} d"
            todos += [c[co], c[cs]]
        for k, loc in enumerate(LOCAIS):
            g = s[s.local == loc]
            ax.scatter(g[f"{e}_obs"], g[f"{e}_sim"], s=44, color=CAT[k], edgecolor="white", linewidth=0.9, zorder=4,
                       label=loc if e == "r1" else None)
        v = pd.concat(todos)
        lo, hi = v.min() - 4, v.max() + 4
        ax.plot([lo, hi], [lo, hi], color=fg.TEXTO_2, lw=1, zorder=3)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect("equal")
        ax.set_title(tit, fontsize=10)
        ax.set_xlabel("observado (dias após a semeadura)")
    axs[0].set_ylabel("simulado (dias após a semeadura)")
    fig.legend(*axs[0].get_legend_handles_labels(), loc="lower center", ncol=6, bbox_to_anchor=(0.5, -0.06))
    fig.suptitle("Validação independente — experimentos de Silva (2021), DSSAT com a mediana da posterior",
                 x=0.01, ha="left", weight="bold")
    fg.salvar(fig, "9_7_validacao_silva")


if __name__ == "__main__":
    fg.estilo()
    mapa_observacoes()
    observado_simulado()
    print("figuras em", fg.DIR_FIG)
