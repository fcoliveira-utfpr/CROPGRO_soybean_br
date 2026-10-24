"""Funções de apoio às figuras dos artigos (notebooks/06_figuras_resultados.ipynb).

Paleta: instância de referência validada (categórica em ordem fixa; sequencial azul; divergente azul↔vermelho com cinza
no meio; status bom/atenção/sério/crítico reservado para classificações, sempre com rótulo de texto).
"""
import json
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, ListedColormap
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon

RAIZ = Path(__file__).resolve().parents[1]
# FIG_LANG=en: mesmas figuras com os textos traduzidos (traducao_en.json) em resultados/figuras_en/ (artigos em inglês)
LINGUA = os.environ.get("FIG_LANG", "pt")
DIR_FIG = RAIZ / "resultados" / ("figuras_en" if LINGUA == "en" else "figuras")
ARQ_TRAD = Path(__file__).with_name("traducao_en.json")

SUPERFICIE = "#fcfcfb"
TEXTO = "#0b0b0b"
TEXTO_2 = "#52514e"
GRADE = "#e4e3df"
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
SEQ = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab",
       "#184f95", "#104281", "#0d366b"]
CINZA_MEIO = "#f0efec"
STATUS = {"bom": "#0ca30c", "atenção": "#fab219", "sério": "#ec835a", "crítico": "#d03b3b"}

CMAP_SEQ = LinearSegmentedColormap.from_list("seq_azul", SEQ)
CMAP_DIV = LinearSegmentedColormap.from_list("div_azul_vermelho", ["#0d366b", "#3987e5", "#b7d3f6", CINZA_MEIO,
                                                                   "#f6c0bf", "#e66767", "#a32424"])


def estilo():
    mpl.rcParams.update({
        "figure.facecolor": SUPERFICIE, "axes.facecolor": SUPERFICIE, "savefig.facecolor": SUPERFICIE,
        "axes.edgecolor": GRADE, "axes.labelcolor": TEXTO_2, "xtick.color": TEXTO_2, "ytick.color": TEXTO_2,
        "text.color": TEXTO, "axes.titlecolor": TEXTO, "axes.titleweight": "bold", "axes.titlesize": 11,
        "axes.titlelocation": "left", "axes.labelsize": 9.5, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRADE,
        "grid.linewidth": 0.6, "lines.linewidth": 2.0, "lines.solid_capstyle": "round", "legend.frameon": False,
        "legend.fontsize": 8.5, "font.size": 9.5, "figure.dpi": 110, "savefig.dpi": 300, "savefig.bbox": "tight",
    })


def _tr(texto, trad):
    if not texto or texto in trad["exato"]:
        return trad["exato"].get(texto, texto)
    for pt, en in trad["trecho"]:          # trechos mais longos primeiro (ordenados no JSON)
        texto = texto.replace(pt, en)
    return texto


def traduzir(fig):
    """Traduz todos os textos da figura (títulos, eixos, legendas, anotações, rótulos fixos dos ticks, colorbars)."""
    from matplotlib.text import Text
    from matplotlib.ticker import FixedFormatter
    trad = json.load(open(ARQ_TRAD)) if ARQ_TRAD.exists() else {"exato": {}, "trecho": []}
    trad["trecho"] = sorted(trad["trecho"], key=lambda x: -len(x[0]))
    fig.canvas.draw()
    for ax in fig.axes:
        for eixo in (ax.xaxis, ax.yaxis):
            f = eixo.get_major_formatter()
            if isinstance(f, FixedFormatter):
                f.seq = [_tr(t, trad) for t in f.seq]
    vistos = set()
    for t in fig.findobj(Text):
        txt = t.get_text()
        vistos.add(txt)
        t.set_text(_tr(txt, trad))
    if os.environ.get("FIG_COLETAR"):      # lista os textos para montar o dicionário
        arq = DIR_FIG / "_textos.txt"
        with open(arq, "a") as f:
            for v in sorted(vistos):
                if v.strip() and any(c.isalpha() for c in v):
                    f.write(repr(v) + "\n")


def salvar(fig, nome):
    DIR_FIG.mkdir(parents=True, exist_ok=True)
    if LINGUA == "en":
        traduzir(fig)
    for ext in ("png", "pdf"):
        fig.savefig(DIR_FIG / f"{nome}.{ext}")


# ------------------------------------------------------------------ mapas
def ufs():
    g = json.load(open(RAIZ / "dados" / "processados" / "brasil_ufs_gaul.geojson"))
    poligonos = []
    def coletar(geo):
        if geo["type"] == "GeometryCollection":
            for g2 in geo["geometries"]:
                coletar(g2)
        elif geo["type"] == "MultiPolygon":
            for p in geo["coordinates"]:
                poligonos.append(np.array(p[0]))
        elif geo["type"] == "Polygon":
            poligonos.append(np.array(geo["coordinates"][0]))

    for f in g["features"]:
        coletar(f["geometry"])
    return poligonos


def mapa_base(ax, fundo=True, contorno=TEXTO_2):
    for p in ufs():
        ax.add_patch(Polygon(p, closed=True, facecolor="#f4f3f0" if fundo else "none", edgecolor=contorno,
                             linewidth=0.45, zorder=1))
    ax.set_xlim(-74.5, -34)
    ax.set_ylim(-34.5, 6)
    ax.set_aspect("equal")
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])


# ------------------------------------------------------------------ fluxogramas
def caixa(ax, x, y, texto, cor=SEQ[1], borda=SEQ[8], w=0.2, h=0.09, fs=8.5, peso="normal"):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.008,rounding_size=0.015",
                                facecolor=cor, edgecolor=borda, linewidth=1.2, zorder=2))
    ax.text(x, y, texto, ha="center", va="center", fontsize=fs, color=TEXTO, weight=peso, zorder=3, wrap=True)


def seta(ax, x0, y0, x1, y1, texto=None, cor=TEXTO_2, curva=0.0):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=11, linewidth=1.2, color=cor,
                                 connectionstyle=f"arc3,rad={curva}", zorder=1))
    if texto:
        ax.text((x0 + x1) / 2, (y0 + y1) / 2, texto, fontsize=7.5, color=TEXTO_2, ha="center", va="center",
                backgroundcolor=SUPERFICIE, zorder=3)


def tela(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")


# ------------------------------------------------------------------ heatmap classificado
def heatmap_classes(ax, M, linhas, colunas, limites, rotulos, cores, fmt="{:.2f}", texto_classe=False, textos=None):
    """Heatmap em que a cor é a CLASSE (faixas `limites`) e o número fica escrito na célula; legenda com os rótulos.
    limites: bordas crescentes (len = n_classes + 1); cores: uma por classe; rotulos: nome de cada classe."""
    cmap = ListedColormap(cores)
    norm = BoundaryNorm(limites, cmap.N)
    M = np.asarray(M, dtype=float)
    ax.imshow(np.ma.masked_invalid(M), cmap=cmap, norm=norm, aspect="auto")
    ax.set_xticks(range(len(colunas)), colunas)
    ax.set_yticks(range(len(linhas)), linhas)
    ax.tick_params(length=0)
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks(np.arange(-0.5, len(colunas)), minor=True)
    ax.set_yticks(np.arange(-0.5, len(linhas)), minor=True)
    ax.grid(which="minor", color=SUPERFICIE, linewidth=2)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            if np.isnan(v):
                ax.text(j, i, "—", ha="center", va="center", fontsize=8, color=TEXTO_2)
                continue
            k = int(np.clip(np.searchsorted(limites, v, side="right") - 1, 0, len(cores) - 1))
            claro = _luminancia(cores[k]) < 0.45
            txt = (textos[i][j] if textos is not None else fmt.format(v)) + (f"\n{rotulos[k]}" if texto_classe else "")
            ax.text(j, i, txt, ha="center", va="center", fontsize=7.8, color="white" if claro else TEXTO)
    alcas = [mpl.patches.Patch(facecolor=c, label=r) for c, r in zip(cores, rotulos)]
    ax.legend(handles=alcas, loc="upper left", bbox_to_anchor=(1.01, 1), title="Classe", title_fontsize=8.5)


def _luminancia(hexcor):
    r, g, b = mpl.colors.to_rgb(hexcor)
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def rhat(amostras):
    n = min(len(a) for a in amostras)
    x = np.stack([a[-n:] for a in amostras])
    W = x.var(axis=1, ddof=1).mean()
    B = n * x.mean(axis=1).var(ddof=1)
    return float(np.sqrt(((n - 1) / n * W + B / n) / W))
