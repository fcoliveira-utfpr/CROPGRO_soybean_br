"""Etapa 6: figuras do relatório, tabela municipal no formato do ZARC e números-chave.

Figuras (resultados/figuras/*.png): contexto (soja e classe AD), janela recomendada por célula, matrizes GMR ×
decêndio em quatro células representativas, riscos por fase, milho safrinha (limite e efeito da correção do R8) e
validação (séries estaduais e habilidade DSSAT × ISNA).
Tabela: resultados/tabela_municipal.csv — município × classe AD × grupo de cultivar × decêndio (numeração do ZARC,
1/set = 25) com nível de risco, classe, restrições, Yw mediano e CVaR20, Yw/Yp mediano, P(safrinha) e R8 corrigido mediano.
Números: resultados/numeros_relatorio.json.
Uso: ../calibracao_mg/.venv/bin/python src/etapa6_figuras.py
"""
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, ListedColormap  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

import config as cfg  # noqa: E402
import etapa5_validacao as e5  # noqa: E402
from etapa2_simular import ARQ_CEL, ARQ_SOLO  # noqa: E402

FIG = cfg.DIR_RES / "figuras"
LADO = 0.25
# paleta de referência (dataviz): superfície, textos, sequencial azul, divergente azul↔vermelho, séries
SUP, TXT, TXT2, GRADE = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
AZUL = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
SEQ = LinearSegmentedColormap.from_list("seq", AZUL)
DIV = LinearSegmentedColormap.from_list("div", ["#184f95", "#6da7ec", "#f0efec", "#ec8c8b", "#b8302f"])
SERIES = {"IBGE": "#2a78d6", "DSSAT (Yw/Yp)": "#eb6834", "Índice tipo ZARC (ISNA)": "#1baf7a"}
CLASSES = ["não recomendado", "recomendado com restrição", "recomendado"]
COR_CLASSE = ["#e4e3df", "#86b6ef", "#1c5cab"]
LETRA = {"não recomendado": "–", "recomendado com restrição": "r", "recomendado": "R"}
REPRESENTATIVAS = {"Oeste (Cascavel)": (-53.46, -24.96), "Norte (Londrina)": (-51.16, -23.31),
                   "Campos Gerais (Ponta Grossa)": (-50.16, -25.09), "Sudoeste (Pato Branco)": (-52.67, -26.23)}
DEC = [f"{m:02d}-{d:02d}" for m, d in cfg.DECENDIOS]
ROT_DEC = {s: f"{int(s[3:])}/{['set', 'out', 'nov', 'dez'][int(s[:2]) - 9]}" for s in DEC}

plt.rcParams.update({"font.size": 9, "axes.edgecolor": GRADE, "axes.labelcolor": TXT2, "xtick.color": TXT2,
                     "ytick.color": TXT2, "text.color": TXT, "figure.facecolor": SUP, "axes.facecolor": SUP,
                     "savefig.facecolor": SUP, "axes.spines.top": False, "axes.spines.right": False})


def contorno_pr(ax):
    g = json.load(open(cfg.CALIB / "dados" / "processados" / "brasil_ufs_gaul.geojson"))
    f = next(f for f in g["features"] if f["properties"]["ADM1_NAME"] in ("Paraná", "Parana"))
    geo = f["geometry"]
    partes = geo["geometries"] if geo["type"] == "GeometryCollection" else [geo]
    for g in partes:                         # GAUL: coleção de polígonos, linhas e pontos
        if g["type"] == "Polygon":
            linhas = [g["coordinates"][0]]
        elif g["type"] == "MultiPolygon":
            linhas = [p[0] for p in g["coordinates"]]
        elif g["type"] == "LineString":
            linhas = [g["coordinates"]]
        else:
            continue
        for l in linhas:
            x, y = np.array(l).T
            ax.plot(x, y, color=TXT2, lw=0.6)


def mapa(ax, df, col, cmap, norm=None, titulo="", rotulo="", vmin=None, vmax=None, ticks=None, ticklabels=None):
    c = pd.read_csv(ARQ_CEL).merge(df, on="id_celula", how="left")
    norm = norm or plt.Normalize(vmin if vmin is not None else np.nanmin(c[col]),
                                 vmax if vmax is not None else np.nanmax(c[col]))
    for r in c.itertuples():
        v = getattr(r, col)
        cor = GRADE if pd.isna(v) else cmap(norm(v))
        ax.add_patch(Rectangle((r.lon - LADO / 2, r.lat - LADO / 2), LADO, LADO, facecolor=cor, edgecolor=SUP, lw=0.4))
    contorno_pr(ax)
    ax.set_xlim(-54.8, -47.9)
    ax.set_ylim(-26.8, -22.4)
    ax.set_aspect(1 / np.cos(np.radians(24.5)))
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(titulo, loc="left", fontsize=9.5, color=TXT)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = plt.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.05, pad=0.02, aspect=30)
    cb.outline.set_visible(False)
    cb.set_label(rotulo, color=TXT2)
    if ticks is not None:
        cb.set_ticks(ticks)
        if ticklabels is not None:
            cb.set_ticklabels(ticklabels)
    return c


def celula_mais_proxima(lon, lat, ids=None):
    c = pd.read_csv(ARQ_CEL)
    if ids is not None:
        c = c[c.id_celula.isin(ids)]
    return c.id_celula.iloc[int(np.argmin(np.hypot(c.lon - lon, c.lat - lat)))]


def fig_contexto(num):
    p = pd.read_csv(ARQ_SOLO)
    p = p[p.valido]
    dom = p.sort_values("area_ha").groupby("id_celula").tail(1)[["id_celula", "classe_ad"]]
    dom["ad_num"] = dom.classe_ad.str[2:].astype(int)
    c = pd.read_csv(ARQ_CEL)[["id_celula", "soja_ha"]].assign(soja_kha=lambda d: d.soja_ha / 1e3)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.9))
    mapa(ax[0], c, "soja_kha", SEQ, titulo="Área de soja por célula de 0,25° (MapBiomas 2023)", rotulo="mil ha")
    cm = ListedColormap(["#9ec5f4", "#3987e5", "#184f95"])
    mapa(ax[1], dom, "ad_num", cm, norm=BoundaryNorm([2.5, 3.5, 4.5, 5.5], 3),
         titulo="Classe de água disponível do ZARC do solo dominante", rotulo="", ticks=[3, 4, 5],
         ticklabels=["AD3", "AD4", "AD5"])
    fig.tight_layout()
    fig.savefig(FIG / "01_contexto.png", dpi=160)
    num["classe_ad_area"] = p.groupby("classe_ad").area_ha.sum().div(p.area_ha.sum()).round(3).to_dict()


def fig_janela(oc, num):
    rec = oc[oc.classe != "não recomendado"]
    jan = rec.groupby(["id_celula", "gmr"]).semeadura.agg(["count", "min", "max"]).reset_index()
    g6 = jan[jan.gmr == 6.0].copy()
    d6 = oc[oc.gmr == 6.0]
    sem = d6[d6.classe == "recomendado"].groupby("id_celula").size().rename("n").reset_index()
    pico = d6.loc[d6.groupby("id_celula").yw_med.idxmax(), ["id_celula", "semeadura"]]
    pico["dec"] = pico.semeadura.map(DEC.index)
    melhor = oc[oc.classe == "recomendado"].groupby("id_celula").size().rename("n_rec").reset_index()
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    mapa(ax[0], sem, "n", SEQ, titulo="GMR 6,0: decêndios recomendados sem restrição", rotulo="decêndios (de 12)",
         vmin=0, vmax=12)
    mapa(ax[1], pico, "dec", SEQ, titulo="GMR 6,0: decêndio de maior Yw mediano", rotulo="", vmin=0, vmax=11,
         ticks=[0, 3, 6, 9, 11], ticklabels=[ROT_DEC[DEC[i]] for i in (0, 3, 6, 9, 11)])
    mapa(ax[2], melhor, "n_rec", SEQ, titulo="Opções GMR × decêndio sem restrição (de 72)", rotulo="opções",
         vmin=0, vmax=72)
    fig.tight_layout()
    fig.subplots_adjust(top=0.9)
    num["pico_yw_gmr6"] = pico.semeadura.value_counts(normalize=True).sort_index().round(3).to_dict()
    fig.savefig(FIG / "02_janela.png", dpi=160)
    num["celulas_sem_janela_gmr6"] = int(270 - len(g6))
    num["frac_opcoes_por_classe"] = oc.classe.value_counts(normalize=True).round(3).to_dict()


def fig_matrizes(oc, num, o_otima):
    fig, ax = plt.subplots(2, 4, figsize=(14, 6.4), sharey=True)
    cm = ListedColormap(COR_CLASSE)
    num["representativas"] = {}
    for j, (nome, (lon, lat)) in enumerate(REPRESENTATIVAS.items()):
        cid = celula_mais_proxima(lon, lat, oc.id_celula.unique())
        d = oc[oc.id_celula == cid]
        cls = d.pivot(index="semeadura", columns="gmr", values="classe").reindex(DEC)
        ry = d.pivot(index="semeadura", columns="gmr", values="yw_med").reindex(DEC)
        ry = ry / np.nanmax(ry.to_numpy(float))               # relativo ao melhor manejo da célula
        saf = d.pivot(index="semeadura", columns="gmr", values="p_safrinha").reindex(DEC)
        oti = o_otima[o_otima.id_celula == cid].pivot(index="semeadura", columns="gmr", values="janela_otima") \
            .reindex(DEC).fillna(False)
        im0 = ax[0, j].imshow(ry.to_numpy(float), cmap=SEQ, vmin=0.5, vmax=1.0, aspect="auto")
        for (i, k), v in np.ndenumerate(cls.to_numpy()):
            escuro = ry.iloc[i, k] > 0.78
            ax[0, j].text(k, i, LETRA[v], ha="center", va="center", fontsize=7, color=SUP if escuro else TXT)
            if oti.iloc[i, k]:
                ax[0, j].add_patch(Rectangle((k - 0.5, i - 0.5), 1, 1, fill=False, edgecolor=TXT, lw=1.6))
        im = ax[1, j].imshow(saf.to_numpy(float), cmap=SEQ, vmin=0, vmax=1, aspect="auto")
        ax[0, j].set_title(f"{nome} ({cid})", loc="left", fontsize=9)
        for a in ax[:, j]:
            a.set_xticks(range(len(cfg.GMR)))
            a.set_xticklabels([f"{g:.1f}" for g in cfg.GMR])
            a.set_yticks(range(len(DEC)))
            a.set_yticklabels([ROT_DEC[s] for s in DEC])
            for s in a.spines.values():
                s.set_visible(False)
        ax[1, j].set_xlabel("GMR")
        num["representativas"][nome] = dict(
            id=cid, recomendadas=int((d.classe == "recomendado").sum()),
            com_restricao=int((d.classe == "recomendado com restrição").sum()),
            melhor=d.sort_values("yw_med", ascending=False).iloc[0][["gmr", "semeadura", "yw_med", "ry_med",
                                                                     "p_safrinha", "classe"]].to_dict())
    ax[0, 0].set_ylabel("Semeadura (Yw / Yw máx. da célula)")
    ax[1, 0].set_ylabel("Semeadura (P colher a tempo da safrinha)")
    fig.text(0.5, 0.975, "Letras: R recomendado · r recomendado com restrição · – não recomendado.   "
             "Borda preta: janela ótima (Pareto em Yw mediano, CVaR20 de Yw e safrinha) em algum solo da célula.",
             ha="center", va="top", color=TXT2, fontsize=8.5)
    fig.subplots_adjust(top=0.9, hspace=0.25, bottom=0.21, right=0.9)
    cb0 = fig.colorbar(im0, cax=fig.add_axes([0.915, 0.57, 0.008, 0.31]))
    cb0.set_label("Yw mediano / melhor opção", color=TXT2)
    cb0.outline.set_visible(False)
    cb = fig.colorbar(im, cax=fig.add_axes([0.3, 0.085, 0.4, 0.015]), orientation="horizontal")
    cb.set_label("probabilidade de colher a tempo do milho safrinha", color=TXT2)
    cb.outline.set_visible(False)
    fig.savefig(FIG / "03_matrizes.png", dpi=160)


def fig_riscos(o_cel, num):
    d = o_cel[(o_cel.gmr == 6.0) & (o_cel.semeadura == "10-01")]
    fig, ax = plt.subplots(1, 4, figsize=(16, 3.9))
    for a, (col, tit) in zip(ax, [("p_seca_flor", "Seca forte na floração"), ("p_geada", "Geada (Tmin ≤ 2 °C) até R7"),
                                  ("p_calor", "≥ 3 dias com Tmax ≥ 36 °C em R1–R5"),
                                  ("p_chuva_colh", "≥ 80 mm nos 10 dias após R8")]):
        mapa(a, d, col, SEQ, titulo=tit, rotulo="probabilidade (GMR 6,0, semeadura 1/out)", vmin=0,
             vmax=max(0.05, float(np.nanmax(d[col]))))
        num[f"{col}_gmr6_1out"] = dict(mediana=round(float(d[col].median()), 3), max=round(float(d[col].max()), 3))
    fig.tight_layout()
    fig.savefig(FIG / "04_riscos.png", dpi=160)


def fig_safrinha(o_cel, num):
    lim = pd.read_csv(cfg.DIR_RES / "safrinha_limite.csv")
    d = o_cel[(o_cel.gmr == 6.0) & (o_cel.semeadura == "09-21")].copy()
    d["dif"] = d.p_safrinha - d.p_safrinha_sc
    fig, ax = plt.subplots(1, 3, figsize=(14, 3.9))
    sem_data = lim.safrinha_limite_doy <= 1                 # nenhuma data entre 1/jan e 15/mar com geada ≤ 20%: em branco
    num["safrinha_sem_data"] = int(sem_data.sum())
    mapa(ax[0], lim[~sem_data], "safrinha_limite_doy", SEQ, titulo="Última semeadura do milho safrinha (geada ≤ 20%)",
         rotulo="", vmin=32, vmax=74, ticks=[32, 46, 60, 74], ticklabels=["1/fev", "15/fev", "1/mar", "15/mar"])
    mapa(ax[1], d, "p_safrinha", SEQ, titulo="P(colher a tempo) — GMR 6,0, semeadura 21/set, R8 corrigido",
         rotulo="probabilidade", vmin=0, vmax=1)
    mapa(ax[2], d, "dif", DIV, norm=plt.Normalize(-0.6, 0.6),
         titulo="Efeito da correção do R8 (corrigido − simulado)", rotulo="diferença de probabilidade")
    fig.tight_layout()
    fig.savefig(FIG / "05_safrinha.png", dpi=160)
    lim = lim[~sem_data]
    num["safrinha_limite"] = dict(min=lim.safrinha_limite.iloc[lim.safrinha_limite_doy.argmin()],
                                  max=lim.safrinha_limite.iloc[lim.safrinha_limite_doy.argmax()],
                                  mediana_doy=float(lim.safrinha_limite_doy.median()))
    num["efeito_correcao_safrinha_21set"] = dict(mediana=round(float(d.dif.median()), 3), min=round(float(d.dif.min()), 3))


def fig_validacao(num):
    res = json.loads((cfg.DIR_RES / "validacao_resumo.json").read_text())
    a = pd.read_csv(cfg.DIR_RES / "validacao_estado_ano.csv")
    fig, ax = plt.subplots(1, 2, figsize=(13, 4), gridspec_kw={"width_ratios": [2.4, 1]})
    for (nome, cor), col in zip(SERIES.items(), ["anom_ibge", "anom_ry", "anom_isna"]):
        z = (a[col] - a[col].mean()) / a[col].std()
        ax[0].plot(a.ano, z, color=cor, lw=2, label=nome)
    ax[0].axhline(0, color=GRADE, lw=1)
    ax[0].set_ylabel("anomalia padronizada (desvios-padrão)")
    ax[0].set_title("Paraná, média ponderada pela área: anomalias da safra (ano de colheita)", loc="left", fontsize=9.5)
    ax[0].legend(frameon=False, loc="lower left", ncol=3)
    ax[0].grid(axis="y", color=GRADE, lw=0.5)
    met = ["pearson", "spearman", "auc_quebra"]
    x = np.arange(len(met))
    for k, (nome, chave) in enumerate([("DSSAT (Yw/Yp)", "ry"), ("Índice tipo ZARC (ISNA)", "isna")]):
        v = [res[chave][m] for m in met]
        ax[1].bar(x + (k - 0.5) * 0.38, v, 0.36, color=SERIES[nome], label=nome)
        for xi, vi in zip(x + (k - 0.5) * 0.38, v):
            ax[1].text(xi, vi + 0.01, f"{vi:.2f}", ha="center", fontsize=7.5, color=TXT2)
    ax[1].set_xticks(x)
    ax[1].set_xticklabels(["r de Pearson", "r de Spearman", "AUC das quebras"])
    ax[1].set_ylim(0, 1)
    ax[1].legend(frameon=False, loc="upper left", fontsize=7.5)
    ax[1].set_title("Habilidade: município × ano", loc="left", fontsize=9.5)
    fig.tight_layout()
    fig.savefig(FIG / "06_validacao.png", dpi=160)
    num["validacao"] = res


def tabela_municipal(o):
    mc = e5.municipios_celulas()
    nomes = pd.read_csv(cfg.CALIB / "dados" / "brutos" / "ibge" / "sidra_1612_soja_municipios_1980_2024.csv",
                        usecols=["cod_municipio", "municipio"]).drop_duplicates("cod_municipio")
    mc = mc.merge(nomes, left_on="cod", right_on="cod_municipio", how="left")
    solos = pd.read_csv(cfg.DIR_RES / "solos_usados.csv")[["id_solo", "classe_ad"]]
    t = o.merge(solos, left_on="solo", right_on="id_solo").merge(mc[["cod", "municipio", "id_celula"]], on="id_celula")
    t["grupo"] = t.gmr.map(lambda g: "I" if g <= 5.5 else "II" if g <= 6.5 else "III")
    t["decendio_zarc"] = t.semeadura.map(lambda s: 24 + DEC.index(s) + 1)
    cols = ["cod", "municipio", "id_celula", "textura", "classe_ad", "grupo", "gmr", "decendio_zarc", "semeadura",
            "nivel_risco", "classe", "restricoes", "yw_med", "yw_cvar20", "ry_med", "ry_cvar20", "p_perda30", "p_safrinha", "p_geada",
            "p_calor", "p_chuva_colh", "r8c_mediana_doy", "janela_otima"]
    t[cols].sort_values(["municipio", "classe_ad", "gmr", "decendio_zarc"]).round(3).to_csv(
        cfg.DIR_RES / "tabela_municipal.csv", index=False)
    return t


if __name__ == "__main__":
    FIG.mkdir(parents=True, exist_ok=True)
    o = pd.read_csv(cfg.DIR_RES / "opcoes.csv")
    oc = pd.read_csv(cfg.DIR_RES / "opcoes_celula.csv")
    num = {}
    fig_contexto(num)
    fig_janela(oc, num)
    otima = o.groupby(["id_celula", "gmr", "semeadura"]).janela_otima.any().reset_index()
    fig_matrizes(oc, num, otima)
    fig_riscos(oc, num)
    fig_safrinha(oc, num)
    if (cfg.DIR_RES / "validacao_resumo.json").exists():
        fig_validacao(num)
    t = tabela_municipal(o)
    num["tabela_municipal"] = dict(municipios=int(t.cod.nunique()), linhas=int(len(t)))
    (cfg.DIR_RES / "numeros_relatorio.json").write_text(json.dumps(num, indent=1, ensure_ascii=False, default=str))
    print(json.dumps(num, indent=1, ensure_ascii=False, default=str)[:3000])
