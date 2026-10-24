"""Análise da replicação de Battisti & Sentelhas (2019): estatísticas por ponto, zonas e comparação com o artigo."""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score

RAIZ = Path(__file__).resolve().parents[1]
DIR_SIM = RAIZ / "resultados" / "bs2019" / "sim"
PONTOS = RAIZ / "dados" / "processados" / "pontos_bs2019_digitalizados.csv"
ALVOS = RAIZ / "dados" / "literatura" / "bs2019_alvos.csv"
ORDEM_SEMEADURA = ["15-set", "15-out", "15-nov", "15-dez", "15-jan", "15-fev"]


def _data(col):
    return pd.to_datetime(col.astype("Int64").astype(str), format="%Y%j", errors="coerce")


def carregar():
    """Todas as simulações, com dias após a semeadura até E, R1 e R7 e a zona digitalizada de cada ponto."""
    sim = pd.concat([pd.read_csv(f) for f in sorted(DIR_SIM.glob("*.csv"))], ignore_index=True)
    sem = _data(sim.PDAT)
    for nome, col in [("E_das", "EDAT"), ("R1_das", "ADAT"), ("R7_das", "MDAT")]:
        sim[nome] = (_data(sim[col]) - sem).dt.days
    sim["semeadura"] = pd.Categorical(sim.semeadura, ORDEM_SEMEADURA, ordered=True)
    pts = pd.read_csv(PONTOS)
    return sim.merge(pts[["lon", "lat", "elev", "zona_bs"]], on=["lon", "lat"], how="left")


def tabela_largura(sim, cenario="atingivel"):
    """Matriz pontos × (semeadura, solo, safra) da produtividade — a mesma usada no agrupamento do artigo."""
    m = sim[sim.cenario == cenario].pivot_table(index=["lon", "lat"], columns=["semeadura", "solo", "safra"],
                                                 values="HWAM", observed=True)
    return m.dropna(axis=1, how="any")


def estatisticas_ponto(sim):
    """Média, CV, P25 e P75 da produtividade atingível por ponto (Fig. 2) e as médias potencial/atingível."""
    a = sim[sim.cenario == "atingivel"].groupby(["lon", "lat"]).HWAM
    p = sim[sim.cenario == "potencial"].groupby(["lon", "lat"]).HWAM
    est = pd.DataFrame({"media": a.mean(), "cv": 100 * a.std() / a.mean(), "p25": a.quantile(.25),
                        "p75": a.quantile(.75), "yp": p.mean()})
    est["perda_pct"] = 100 * (1 - est.media / est.yp)
    return est.reset_index()


def agrupar(sim, k=10):
    """Agrupamento hierárquico de Ward (distância euclidiana) dos pontos pela produtividade atingível."""
    m = tabela_largura(sim)
    grupos = fcluster(linkage(m.values, "ward"), k, "maxclust")
    return pd.Series(grupos, index=m.index, name="zona_rep").reset_index()


def casar_zonas(z, col_ref="zona_bs", col_rep="zona_rep"):
    """Renumera as zonas replicadas pela correspondência de maior sobreposição com as do artigo (Húngaro)."""
    ok = z.dropna(subset=[col_ref])
    tab = pd.crosstab(ok[col_rep], ok[col_ref])
    lin, col = linear_sum_assignment(-tab.values)
    mapa = {tab.index[i]: int(tab.columns[j]) for i, j in zip(lin, col)}
    z = z.assign(zona_rep_casada=z[col_rep].map(mapa))
    ari = adjusted_rand_score(ok[col_ref], ok[col_rep])
    acerto = (z.dropna(subset=[col_ref]).eval(f"{col_ref} == zona_rep_casada")).mean()
    return z, ari, acerto


def resumo_por_zona(sim, col="zona_bs"):
    """Médias por zona de Yp e Yw (todas as semeaduras, solos e safras), como no texto e na Fig. 5."""
    g = sim.dropna(subset=[col]).groupby([col, "cenario"]).HWAM
    r = g.agg(["mean", "min", "max"]).unstack("cenario")
    r.columns = [f"{c}_{e}" for c, e in r.columns]
    r["perda_pct"] = 100 * (1 - r.mean_atingivel / r.mean_potencial)
    return r.round(0)


def comparar_alvos(resumo):
    """Coloca lado a lado os números do texto do artigo e os replicados."""
    alvos = pd.read_csv(ALVOS)
    mapa = {"produtividade_atingivel_media": "mean_atingivel", "produtividade_potencial_media": "mean_potencial",
            "produtividade_potencial_min": "min_potencial", "produtividade_potencial_max": "max_potencial",
            "perda_por_deficit_hidrico": "perda_pct", "perda_por_deficit_hidrico_max": "perda_pct",
            "produtividade_atingivel_min": "min_atingivel", "produtividade_atingivel_max": "max_atingivel"}
    linhas = []
    for a in alvos.itertuples():
        col = mapa.get(a.variavel)
        if col is None:
            continue
        if a.zona == "todas":
            rep = resumo[col].min() if col.startswith("min") else resumo[col].max()
        elif "-" in str(a.zona):
            zs = [float(z) for z in range(int(a.zona.split("-")[0]), int(a.zona.split("-")[1]) + 1)]
            sub = resumo.loc[resumo.index.isin(zs), col]
            rep = sub.min() if col.startswith("min") else sub.max()
        else:
            rep = resumo[col].get(float(a.zona), np.nan)
        linhas.append(dict(variavel=a.variavel, zona=a.zona, artigo=a.valor, replicado=rep, unidade=a.unidade))
    df = pd.DataFrame(linhas)
    df["dif_pct"] = (100 * (df.replicado - df.artigo) / df.artigo).round(1)
    return df


def comparar_ibge_microrregiao(sim, est=None):
    """Fig. 6: média simulada (atingível, 1980–2013) × média do IBGE (2007–2016) por microrregião.

    Cada município com soja em 2016 é ligado ao ponto de clima mais próximo; a média da microrregião
    pondera os pontos pela área de soja dos municípios. A zona da microrregião é a de maior área.
    """
    proc = RAIZ / "dados" / "processados"
    est = est if est is not None else estatisticas_ponto(sim)
    pts = est.merge(pd.read_csv(PONTOS)[["lon", "lat", "zona_bs"]], on=["lon", "lat"])
    mun = pd.read_csv(proc / "sidra2016_soja_municipios.csv")[["codigo_ibge", "area_ha", "lon", "lat"]]
    mic = pd.read_csv(proc / "municipio_microrregiao.csv").rename(columns={"cod": "codigo_ibge"})
    mun = mun.merge(mic, on="codigo_ibge", how="left")
    mun["nome_micro"] = mun.microrregiao.str.strip() + " - " + mun.uf_sigla.str.strip()
    d = np.hypot(mun.lon.values[:, None] - pts.lon.values[None], mun.lat.values[:, None] - pts.lat.values[None])
    j = d.argmin(1)
    mun = mun.assign(sim=pts.media.values[j], zona_bs=pts.zona_bs.values[j], dist=d.min(1))
    mun = mun[mun.dist <= 0.75]                       # município a até ~0,75° do ponto
    g = mun.groupby("nome_micro")
    micro = pd.DataFrame({"sim": g.apply(lambda x: np.average(x.sim, weights=x.area_ha)),
                          "zona_bs": g.apply(lambda x: x.groupby("zona_bs").area_ha.sum().idxmax() if x.zona_bs.notna().any() else np.nan),
                          "area_ha": g.area_ha.sum()})
    ibge = pd.read_csv(proc / "sidra_soja_microrregiao_2007_2016.csv")
    obs = ibge.groupby("microrregiao").rend_kg_ha.mean().rename("ibge")
    micro = micro.join(obs, how="inner").reset_index().rename(columns={"index": "microrregiao"})
    r = micro.dropna(subset=["zona_bs"]).groupby("zona_bs").apply(
        lambda x: pd.Series({"n": len(x), "r": x.sim.corr(x.ibge) if len(x) > 2 else np.nan,
                             "sim_media": x.sim.mean(), "ibge_media": x.ibge.mean()}))
    return micro, r.round(2)


def estatisticas_ponto_ponderadas(sim, fracoes):
    """Como estatisticas_ponto, mas a média atingível de cada ponto pondera os 3 solos pela área de soja
    em cada classe de textura do MapBiomas (em vez de peso igual). Pontos sem soja ficam de fora."""
    f = fracoes.melt(id_vars=["lon", "lat"], value_vars=[f"f_{s}" for s in ("francoarenoso", "argiloarenoso", "argiloso")],
                     var_name="solo", value_name="peso")
    f["solo"] = f.solo.str[2:]
    a = sim[sim.cenario == "atingivel"].groupby(["lon", "lat", "solo"]).HWAM.mean().reset_index().merge(f, on=["lon", "lat", "solo"])
    media = a.groupby(["lon", "lat"]).apply(lambda g: np.average(g.HWAM, weights=g.peso) if g.peso.sum() > 0 else np.nan)
    est = estatisticas_ponto(sim).set_index(["lon", "lat"])
    est["media_igual"] = est.media
    est["media"] = media
    return est.dropna(subset=["media"]).reset_index()


def estatisticas_solo_mb(pasta):
    """Média atingível por célula com os perfis do MapBiomas, ponderada pela área de soja de cada perfil."""
    s = pd.concat([pd.read_csv(f) for f in sorted(Path(pasta).glob("*.csv"))], ignore_index=True)
    m = s.groupby(["lon", "lat", "solo"]).agg(HWAM=("HWAM", "mean"), area_ha=("area_ha", "first"), ad80=("ad80", "first")).reset_index()
    est = m.groupby(["lon", "lat"]).apply(lambda g: pd.Series({"media": np.average(g.HWAM, weights=g.area_ha),
                                                               "ad80": np.average(g.ad80, weights=g.area_ha)})).reset_index()
    return est, s
