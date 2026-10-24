"""Monta a tabela de observações de fenologia (GMR conhecido) para a calibração.

Cada linha: local de simulação (um dos locais do treino de fenologia), safra(s), data de semeadura, GMR e a
observação (dias até R1, R7 ou R8, contados da EMERGÊNCIA ou, se o trabalho só informa isso, da SEMEADURA: coluna
`base`; a coluna `dias_ve` guarda o valor nos dois casos). Tipos:
- "experimento": registros de dados/literatura/fenologia_observada.csv (safra e data conhecidas);
- "catalogo": ciclo (emergência → maturação de colheita, R8) dos catálogos da Embrapa, atribuído aos pontos
  representativos da região e da faixa de altitude, semeadura na época recomendada e safras 2005/06–2017/18;
- "agregado": médias por cultivar de Alliprandini et al. (2009) em vários locais e duas safras (R1 e R8 desde a semeadura).
"""
from pathlib import Path

import numpy as np
import pandas as pd

RAIZ = Path(__file__).resolve().parents[1]
DIR_EMU = RAIZ / "resultados" / "emulador_fenologia"
SAFRAS_CATALOGO = range(2005, 2018, 3)      # 5 safras: o catálogo é uma média, mais safras não acrescentam informação
SEMEADURA_CATALOGO = {"Centro-Sul (MR 1–3)": "11-01", "GO": "10-25", "MT": "10-20", "MT/RO": "10-20"}
FAIXA = {"até 500 m": "CS<500", "500 a 700 m": "CS500-700", "acima de 700 m": "CS>700", "até 700 m": "CS<500"}


def experimentos(locais):
    f = pd.read_csv(RAIZ / "dados" / "literatura" / "fenologia_observada.csv")
    f = f[f.data_semeadura.notna() & (f.data_semeadura != "")]
    if "usar" in f:
        f = f[f.usar.fillna(1).astype(int) == 1]          # registros marcados como inconsistentes ficam de fora
    idx = {nome: i for i, nome in zip(locais.index, locais.local)}
    linhas = []
    for r in f.itertuples():
        ds = pd.Timestamp(r.data_semeadura)
        safra = ds.year - (ds.month < 7)
        for est in ("r1", "r5", "r7", "r8"):
            # contagem desde a emergência quando existe; senão desde a semeadura (base = "semeadura")
            v, base = getattr(r, f"ve_{est}"), "emergencia"
            if pd.isna(v):
                v, base = getattr(r, f"das_{est}"), "semeadura"
            if pd.notna(v):
                linhas.append(dict(tipo="experimento", fonte=r.fonte, idx_local=idx[r.local], safra=safra,
                                   data_semeadura=ds.strftime("%m-%d"), gmr=r.gmr, estadio=est, dias_ve=float(v), peso=1.0,
                                   base=base))
    return pd.DataFrame(linhas)


def catalogos(locais):
    c = pd.read_csv(RAIZ / "dados" / "literatura" / "catalogos_embrapa_ciclo.csv").dropna(subset=["ciclo_medio"])
    linhas = []
    for r in c.itertuples():
        if r.regiao.startswith("Centro-Sul"):
            grupo = FAIXA.get(r.altitude) if isinstance(r.altitude, str) else None
            sem = SEMEADURA_CATALOGO["Centro-Sul (MR 1–3)"]
        elif r.regiao == "GO":
            grupo, sem = "GO", SEMEADURA_CATALOGO["GO"]
        elif r.regiao in ("MT", "MT/RO"):
            grupo, sem = "MT/RO", SEMEADURA_CATALOGO[r.regiao]
        else:
            continue
        if grupo is None:
            continue
        idxs = locais.index[locais.local == grupo]
        # um registro de catálogo vale como uma observação: divide o peso entre locais e safras
        w = 1.0 / (len(idxs) * len(SAFRAS_CATALOGO))
        for i in idxs:
            for s in SAFRAS_CATALOGO:
                linhas.append(dict(tipo="catalogo", fonte=r.arquivo, idx_local=i, safra=s, data_semeadura=sem, gmr=r.gmr,
                                   estadio="r8", dias_ve=float(r.ciclo_medio), peso=w, id_registro=f"{r.arquivo}:{r.linha}:{r.altitude}",
                                   base="emergencia"))
    return pd.DataFrame(linhas)


LOCAIS_ALLIPRANDINI = {
    "centro_oeste": ["Sorriso-MT", "Primavera do Leste-MT", "Campo Verde-MT", "Goiatuba-GO", "Morrinhos-GO", "Uberlândia-MG",
                     "Brasília-DF", "São Gotardo-MG"],
    "sul": ["Cruz Alta-RS", "Santa Tereza do Oeste-PR", "Cascavel-PR", "Palotina-PR", "Cruzeiro do Sul-PR", "Ponta Grossa-PR",
            "Morro Agudo-SP"]}


def agregados(locais):
    """Alliprandini et al. (2009, Crop Sci. 49:801, Tabs. 2 e 3): dias da semeadura até R1 (NDF) e R8 (NDM), média de
    7 (Sul) ou 8 (Centro-Oeste) locais × 2 safras (2002/03 e 2003/04), semeadura nas duas primeiras semanas de novembro
    (aqui 07/11). Cada valor é uma observação; a previsão é a média sobre os locais e safras do grupo."""
    a = pd.read_csv(RAIZ / "dados" / "literatura" / "alliprandini2009_ciclo_regioes.csv")
    idx = {nome: i for i, nome in zip(locais.index, locais.local)}
    linhas = []
    for r in a.itertuples():
        ii = [idx[n] for n in LOCAIS_ALLIPRANDINI[r.regiao]]
        for est, v in (("r1", r.ndf), ("r8", r.ndm)):
            for i in ii:
                for s in (2002, 2003):
                    linhas.append(dict(tipo="agregado", fonte="Alliprandini et al. (2009)", idx_local=i, safra=s,
                                       data_semeadura="11-07", gmr=r.mg, estadio=est, dias_ve=float(v), peso=1 / (2 * len(ii)),
                                       id_registro=f"alli:{r.regiao}:{r.cultivar}:{est}", base="semeadura"))
    return pd.DataFrame(linhas)


def catalogos_extra(locais):
    """Registros tipo catálogo com locais, safras e semeadura explícitos (dados/literatura/catalogos_extra_ciclo.csv):
    catálogo da Embrapa para MS 2010/11 e Santana et al. (2022). Locais por nome ou por índice ("#35")."""
    c = pd.read_csv(RAIZ / "dados" / "literatura" / "catalogos_extra_ciclo.csv")
    idx = {nome: i for i, nome in zip(locais.index, locais.local)}
    linhas = []
    for r in c.itertuples():
        ii = [int(n[1:]) if n.startswith("#") else idx[n] for n in r.locais.split(";")]
        ss = [int(x) for x in str(r.safras).split(";")]
        for i in ii:
            for s in ss:
                linhas.append(dict(tipo="catalogo", fonte=r.fonte, idx_local=i, safra=s, data_semeadura=r.data_semeadura, gmr=r.gmr,
                                   estadio="r8", dias_ve=float(r.ciclo_ve_r8), peso=1 / (len(ii) * len(ss)),
                                   id_registro=f"extra:{r.fonte[:20]}:{r.cultivar}:{r.regiao}", base="emergencia"))
    return pd.DataFrame(linhas)


def montar():
    locais = pd.read_csv(DIR_EMU / "locais.csv", index_col="idx_local")
    obs = pd.concat([experimentos(locais), catalogos(locais), catalogos_extra(locais), agregados(locais)], ignore_index=True)
    obs["doy_sem"] = obs.data_semeadura.map(lambda md: (pd.Timestamp(f"2001-{md}") - pd.Timestamp("2001-09-10")).days % 365)
    obs = obs.join(locais[["local", "lon", "lat", "elev"]], on="idx_local")
    obs.to_csv(DIR_EMU / "observacoes_fenologia.csv", index=False)
    return obs


if __name__ == "__main__":
    o = montar()
    print(o.groupby(["tipo", "estadio"]).agg(n=("dias_ve", "size"), peso=("peso", "sum"), media=("dias_ve", "mean")).round(1))
