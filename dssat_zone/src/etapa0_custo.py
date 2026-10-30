"""Etapa 0: custo da execução direta do DSSAT num ponto do oeste do PR (−53,5; −24,5).

Um FileX com os solos MapBiomas da célula × 6 GMR × 12 decêndios (com água) + 6 GMR × 12 decêndios sem estresse
hídrico (Yp), 45 safras cada. Mede o tempo e confere as saídas (número de linhas, datas, Yw ≤ Yp).
Uso: ../calibracao_mg/.venv/bin/python src/etapa0_custo.py
"""
import time

import pandas as pd

import config as cfg
import dssat_direto as dd
import ptf

LON, LAT, ELEV = -53.5, -24.5, 434

if __name__ == "__main__":
    pasta = cfg.DIR_RUN / "custo"
    dd.preparar_pasta(pasta)
    clima = pd.read_csv(cfg.CALIB / "dados" / "brutos" / "xavier" / f"xavier_no_{LON:+.2f}_{LAT:+.2f}_1980_2025.csv.gz",
                        index_col=0, parse_dates=True).interpolate(limit=5).ffill().bfill()
    wsta = dd.escrever_clima(pasta, clima, LAT, LON, ELEV)
    perfis = pd.read_csv(cfg.CALIB / "dados" / "processados" / "perfis_solo_mapbiomas_soja2016.csv")
    perfis = perfis[(perfis.lon == LON) & (perfis.lat == LAT)]
    solos = {f"ZN{k:08d}": (ptf.camadas_perfil_mb(p), ptf.TEXTURA_DSSAT[p.solo])
             for k, (_, p) in enumerate(perfis.iterrows())}
    dd.escrever_solos(pasta, solos, LAT, LON)
    s0 = next(iter(solos))
    trat = [dict(solo=s, gmr=g, semeadura=sem, agua=True) for s in solos for g in cfg.GMR for sem in cfg.DECENDIOS]
    trat += [dict(solo=s0, gmr=g, semeadura=sem, agua=False) for g in cfg.GMR for sem in cfg.DECENDIOS]
    t0 = time.time()
    s = dd.rodar_tratamentos(pasta, wsta, solos, trat)
    dt = time.time() - t0
    print(f"{len(trat)} tratamentos × {cfg.NYERS} safras = {len(trat) * cfg.NYERS} rodadas; {len(s)} no Summary")
    print(f"tempo {dt:.1f} s → {1000 * dt / len(s):.1f} ms por safra")
    s.to_csv(cfg.DIR_RES / "etapa0_custo_summary.csv", index=False)
    print(s.groupby(["agua", "gmr"]).HWAM.median().unstack(0).round(0))
