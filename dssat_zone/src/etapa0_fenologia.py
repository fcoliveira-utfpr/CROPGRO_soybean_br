"""Etapa 0: viés fenológico por época de semeadura nos experimentos do Sul e do PR, com as cultivares nacionais.

Simula, sem estresse hídrico, os experimentos de Umuarama-PR (27/out e 10/nov, 2005–2007), Passo Fundo-RS (24/set,
14/nov e 1/fev, 2012/13), Ijuí-RS (dez–jan) e Jaboticabal-SP (out e nov, controle tropical) com a cultivar do GMR exato
de cada registro (θ(GMR) na mediana da posterior v13, gerada com gerar_produtos) e compara R1 e R8 com o observado.
R8 simulado = R7 (MDAT) + δ8 (11 d, como na calibração). Dias contados da emergência ou da semeadura, como no registro.

Saída: resultados/etapa0_fenologia.csv e resumo por local × época na tela.
Uso: DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1 ../calibracao_mg/.venv/bin/python src/etapa0_fenologia.py
"""
import numpy as np
import pandas as pd

import config as cfg
import bs2019 as bs
import dssat_direto as dd
import gerar_produtos as gp
import ptf
import treino_fenologia as tf

LOCAIS = ["Umuarama-PR", "Passo Fundo-RS", "Ijuí-RS", "Jaboticabal-SP"]
DELTA8 = 11


def data(yyddd):
    return pd.to_datetime(yyddd.astype(int).astype(str), format="%Y%j")


def cultivares_exatas(pasta, gmrs):
    """Acrescenta ao .CUL/.ECO da pasta as cultivares BRGxxx/BRExxx dos GMR que não estão nos arquivos publicados."""
    th = np.median(gp.posterior(), axis=0)
    cul, eco = (pasta / "SBGRO048.CUL").read_text(encoding="latin-1"), (pasta / "SBGRO048.ECO").read_text(encoding="latin-1")
    for g in gmrs:
        var = cfg.cultivar_id(g)
        if f"\n{var} " in cul:
            continue
        c = gp.cultivar_completa(th, g)
        ecoid = "BRE" + var[3:]
        cul += gp.linha_cul(var, f"BR RMG {g:.1f}", ecoid, c) + "\n"
        eco += gp.linha_eco(ecoid, f"BR RMG {g:.1f}", gp.ler_eco(f"SB0{int(np.clip(round(g), 5, 9))}01"), c["r1ppo"],
                            bs.ECO_S2["fl-vs"]) + "\n"
    (pasta / "SBGRO048.CUL").write_text(cul, encoding="latin-1")
    (pasta / "SBGRO048.ECO").write_text(eco, encoding="latin-1")


if __name__ == "__main__":
    obs = pd.read_csv(cfg.CALIB / "resultados" / "emulador_fenologia" / "observacoes_fenologia.csv")
    obs = obs[(obs.tipo == "experimento") & obs.local.isin(LOCAIS) & obs.estadio.isin(["r1", "r8"])].copy()
    obs["mes"] = obs.data_semeadura.str[:2].astype(int)
    obs["dia"] = obs.data_semeadura.str[3:].astype(int)
    locais = pd.read_csv(cfg.CALIB / "resultados" / "emulador_fenologia" / "locais.csv", index_col="idx_local")
    perfis = pd.read_csv(cfg.CALIB / "dados" / "processados" / "perfis_solo_mapbiomas_soja2016.csv")
    p = perfis[(perfis.lon == -53.5) & (perfis.lat == -24.5) & (perfis.solo == "argiloso")].iloc[0]
    solos = {"ZN00000000": (ptf.camadas_perfil_mb(p), "C")}      # sem estresse hídrico: o solo não altera a fenologia
    pasta = cfg.DIR_RUN / "fenologia"
    res = []
    for (il, safra), o in obs.groupby(["idx_local", "safra"]):
        loc = locais.loc[il]
        dd.preparar_pasta(pasta)
        cultivares_exatas(pasta, sorted(o.gmr.unique()))
        clima = pd.read_csv(tf.arq_clima(dict(loc, idx=il)), index_col=0, parse_dates=True).interpolate(limit=5).ffill().bfill()
        wsta = dd.escrever_clima(pasta, clima, loc.lat, loc.lon, loc.elev)
        dd.escrever_solos(pasta, solos, loc.lat, loc.lon)
        cfg.ANO_INI, cfg.NYERS = int(safra), 1
        trat = [dict(solo="ZN00000000", gmr=g, semeadura=(m, d), agua=False)
                for g, m, d in o[["gmr", "mes", "dia"]].drop_duplicates().itertuples(index=False)]
        s = dd.rodar_tratamentos(pasta, wsta, solos, trat)
        for r in o.itertuples():
            x = s[(s.gmr == r.gmr) & (s.semeadura == (r.mes, r.dia))].iloc[0]
            ref = x.EDAT if r.base == "emergencia" else x.PDAT
            est = x.ADAT if r.estadio == "r1" else x.MDAT
            sim = (data(pd.Series([est]))[0] - data(pd.Series([ref]))[0]).days + (DELTA8 if r.estadio == "r8" else 0)
            res.append(dict(local=r.local, safra=r.safra, semeadura=r.data_semeadura, mes=r.mes, gmr=r.gmr,
                            estadio=r.estadio, obs=r.dias_ve, sim=sim, erro=sim - r.dias_ve))
    d = pd.DataFrame(res)
    d.to_csv(cfg.DIR_RES / "etapa0_fenologia.csv", index=False)
    print(d.groupby(["local", "estadio", "semeadura"]).agg(n=("erro", "size"), obs=("obs", "mean"), sim=("sim", "mean"),
                                                         vies=("erro", "mean"), rmse=("erro", lambda e: np.sqrt(np.mean(e ** 2))))
          .round(1).to_string())
