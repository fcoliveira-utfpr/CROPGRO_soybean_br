"""Inventário dos dados do projeto: quantidade de registros por local, em três mapas (30/10/2026).

(a) Observações de fenologia da calibração (746 termos de verossimilhança): experimentos de campo, médias multilocais
    de Alliprandini et al. (2009) e catálogos de cultivares; tamanho do símbolo = registros no local.
(b) Satélite e emulador: pontos com fenologia MODIS (cor = safras invertidas pelo NPE) e os 73 locais do treino do
    emulador.
(c) Produtividade e zoneamento: municípios com soja no IBGE (PAM, 1980–2024), as 270 células de 0,25° do zoneamento
    do PR e os ensaios de época de semeadura usados no viés fenológico do zoneamento (subconjunto dos experimentos).

Saída: resultados/figuras_en/0_1_inventario_dados.png/.pdf (textos já em inglês) e resultados/inventario_dados.csv.
Uso: FIG_LANG=en .venv/bin/python src/figura_inventario.py
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

import figuras as fg

RAIZ = fg.RAIZ
DE = RAIZ / "resultados" / "emulador_fenologia"
ZONE = RAIZ.parent / "dssat_zone"
ESCALA = 2.2                     # área do símbolo (pt²) por registro


def tamanho(n):
    return 12 + ESCALA * np.asarray(n, float)


if __name__ == "__main__":
    fg.estilo()
    o = pd.read_csv(DE / "observacoes_fenologia.csv")
    loc = pd.read_csv(DE / "locais.csv")
    exp = o[o.tipo == "experimento"].groupby(["local", "lon", "lat"]).size().rename("n").reset_index()
    agr = o[o.tipo == "agregado"].groupby(["local", "lon", "lat"]).id_registro.nunique().rename("n").reset_index()
    c = o[o.tipo == "catalogo"]              # cada registro de catálogo se repete nos pontos representativos da região
    cat = c.groupby("local").agg(lon=("lon", "mean"), lat=("lat", "mean"), n=("id_registro", "nunique")).reset_index()
    n_cat = c.id_registro.nunique()
    n_agr = o[o.tipo == "agregado"].id_registro.nunique()
    sat = pd.read_csv(RAIZ / "produtos" / "dados" / "coefficients_point_season.csv").groupby(["lon", "lat"]).size() \
        .rename("n").reset_index()
    ib = pd.read_csv(RAIZ / "dados" / "brutos" / "ibge" / "sidra_1612_soja_municipios_1980_2024.csv")
    ib = ib[ib.rendimento_kg_ha > 0]
    cen = pd.read_csv(RAIZ / "dados" / "processados" / "centroides_municipios.csv")
    mun = ib.groupby("cod_municipio").size().rename("n").reset_index().merge(cen, left_on="cod_municipio", right_on="cod")
    cel = pd.read_csv(ZONE / "dados" / "processados" / "celulas.csv")
    epoca = pd.read_csv(ZONE / "resultados" / "etapa0_fenologia.csv").groupby("local").size().rename("n").reset_index() \
        .merge(exp[["local", "lon", "lat"]], on="local")

    fig, ax = plt.subplots(1, 3, figsize=(16, 7.4))
    for a in ax:
        fg.mapa_base(a)
        a.set_xlim(-66, -34)
        a.set_ylim(-34.5, 1)

    # (a) calibração
    a = ax[0]
    a.scatter(agr.lon, agr.lat, s=tamanho(agr.n), marker="s", color=fg.CAT[1], edgecolor="white", lw=0.8, zorder=3,
              alpha=0.9)
    a.scatter(cat.lon, cat.lat, s=tamanho(cat.n), marker="^", color=fg.CAT[2], edgecolor="white", lw=0.8, zorder=4,
              alpha=0.9)
    a.scatter(exp.lon, exp.lat, s=tamanho(exp.n), marker="o", color=fg.CAT[0], edgecolor="white", lw=0.8, zorder=5,
              alpha=0.9)
    a.set_title("(a) Phenology observations in the calibration", fontsize=10.5)
    leg = [Line2D([], [], ls="", marker="o", color=fg.CAT[0], ms=8,
                  label=f"field experiments: {exp.n.sum()} records, {len(exp)} sites"),
           Line2D([], [], ls="", marker="s", color=fg.CAT[1], ms=8,
                  label=f"multi-site means (Alliprandini 2009): {n_agr} records, {len(agr)} sites"),
           Line2D([], [], ls="", marker="^", color=fg.CAT[2], ms=8,
                  label=f"cultivar catalogues: {n_cat} records, {len(cat)} regions")]
    l1 = a.legend(handles=leg, loc="upper left", bbox_to_anchor=(0, 0), fontsize=7.8, title="746 likelihood terms",
                  title_fontsize=8.5)
    l1._legend_box.align = "left"
    a.add_artist(l1)
    ref = [10, 50, 100]
    a.legend(handles=[Line2D([], [], ls="", marker="o", color=fg.TEXTO_2, mfc="none",
                             ms=np.sqrt(tamanho(r)), label=f"{r}") for r in ref],
             loc="upper right", fontsize=7.8, title="records per site", title_fontsize=8, labelspacing=1.1)
    desloc = {"Passo Fundo-RS": (10, -12), "Ijuí-RS": (-12, 8), "Jaboticabal-SP": (9, -3)}
    for r in exp.itertuples():                                   # rótulos só dos locais com mais registros
        if r.local in desloc:
            dx, dy = desloc[r.local]
            a.annotate(f"{r.local.split('-')[0]} ({r.n})", (r.lon, r.lat), xytext=(dx, dy), textcoords="offset points",
                       fontsize=7.5, color=fg.TEXTO_2, ha="left" if dx > 0 else "right")

    # (b) satélite e emulador
    a = ax[1]
    sc = a.scatter(sat.lon, sat.lat, c=sat.n, cmap=fg.CMAP_SEQ, s=26, vmin=0, vmax=24, zorder=3, edgecolor="none")
    a.scatter(loc.lon, loc.lat, s=34, facecolor="none", edgecolor=fg.TEXTO, lw=0.9, zorder=4)
    a.set_title("(b) Satellite phenology and emulator training", fontsize=10.5)
    a.legend(handles=[Line2D([], [], ls="", marker="o", color=fg.SEQ[7], ms=6,
                             label=f"MODIS points: {len(sat)} points, {sat.n.sum():,} point-seasons (2000–2024)"),
                      Line2D([], [], ls="", marker="o", mfc="none", mec=fg.TEXTO, ms=7,
                             label=f"emulator training: {len(loc)} sites, 2.81 M simulated seasons")],
             loc="upper left", bbox_to_anchor=(0, 0), fontsize=7.8)
    cb = fig.colorbar(sc, cax=a.inset_axes([0.6, 0.2, 0.36, 0.025]), orientation="horizontal")
    cb.set_label("seasons inverted by the NPE per point", color=fg.TEXTO_2, fontsize=8)
    cb.outline.set_visible(False)

    # (c) produtividade e zoneamento
    a = ax[2]
    a.scatter(mun.lon, mun.lat, s=2.5, color="#b5b3ad", zorder=2, edgecolor="none")
    a.scatter(cel.lon, cel.lat, s=7, marker="s", color=fg.CAT[0], zorder=3, edgecolor="none")
    a.scatter(epoca.lon, epoca.lat, s=tamanho(epoca.n), marker="D", facecolor="none", edgecolor=fg.CAT[1], lw=1.4,
              zorder=4)
    a.set_title("(c) Yield data and Paraná zoning", fontsize=10.5)
    a.legend(handles=[Line2D([], [], ls="", marker="o", color="#b5b3ad", ms=5,
                             label=f"IBGE municipal yield: {len(mun):,} municipalities, {len(ib):,} municipality-years"),
                      Line2D([], [], ls="", marker="s", color=fg.CAT[0], ms=5,
                             label=f"zoning: {len(cel)} cells of 0.25°, 1.69 M simulated seasons"),
                      Line2D([], [], ls="", marker="D", mfc="none", mec=fg.CAT[1], ms=7,
                             label=f"sowing-date trials (zoning bias): {epoca.n.sum()} records, {len(epoca)} sites")],
             loc="upper left", bbox_to_anchor=(0, 0), fontsize=7.8)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.2)
    fg.salvar(fig, "0_1_inventario_dados")      # bbox "tight" corta legendas fora dos eixos: salva também com folga
    for ext in ("png", "pdf"):
        fig.savefig(fg.DIR_FIG / f"0_1_inventario_dados.{ext}", bbox_inches="tight", pad_inches=0.3,
                    bbox_extra_artists=[l1] + [a.get_legend() for a in ax[1:]])

    inv = pd.DataFrame([
        ("calibração", "experimentos de campo (R1, R5, R7, R8)", exp.n.sum(), f"{len(exp)} locais, 8 estudos"),
        ("calibração", "médias multilocais (Alliprandini et al., 2009; R1, R8)", n_agr, f"{len(agr)} locais, 88 cultivares"),
        ("calibração", "catálogos de cultivares (ciclo)", n_cat, f"{len(cat)} regiões, 7 fontes"),
        ("satélite", "pontos-safra MODIS invertidos pelo NPE", sat.n.sum(), f"{len(sat)} pontos"),
        ("emulador", "locais de treino", len(loc), "2,81 milhões de safras simuladas"),
        ("IBGE", "município-anos com produtividade de soja (PAM)", len(ib), f"{len(mun)} municípios"),
        ("zoneamento", "células de 0,25° no PR", len(cel), "1,69 milhão de safras simuladas"),
        ("zoneamento", "registros dos ensaios de época (viés fenológico)", epoca.n.sum(), f"{len(epoca)} locais"),
    ], columns=["bloco", "dado", "n", "detalhe"])
    inv.to_csv(RAIZ / "resultados" / "inventario_dados.csv", index=False)
    print(inv.to_string(index=False))
