"""Etapa 2: simulações fatoriais por célula de 0,25° do PR.

Por célula: solos MapBiomas com fração ≥ F_MIN da soja da célula × 6 GMR × 12 decêndios, com água (Yw), mais o
solo de maior área × 6 GMR × 12 decêndios sem estresse hídrico (Yp, que não depende do solo); 45 safras
(1980/81–2024/25) cada. Saídas do Summary (datas, HWAM, chuva, ET) e, por fase, do OVERVIEW (estresse hídrico,
temperaturas, dias de Tmin < 0/2 °C e Tmax > 32/34 °C).

Solos (chave SOLO): "bs2019" (padrão) = os 3 solos de Battisti & Sentelhas (2019, Tab. S1) das classes texturais
presentes na soja da célula (MapBiomas Solo), pesados pela área — a configuração dos artigos, que na calibração
acompanhou melhor o IBGE (r 0,44) que os perfis por PTF (r 0,27; calibracao_mg/resultados/solo_mapbiomas/
comparacao_variantes_solo.csv); "ptf" = perfis MapBiomas + SoilGrids pela PTF de Medrado & Lima (sensibilidade).
Chave REESCALAR_AWC (só com "ptf"; padrão desligada): reescala DUL − LL para que a água disponível de 0–100 cm iguale a AWC da
ANA × 1000 mm (mesmo fator abaixo de 100 cm; DUL limitado a SAT − 0,01).

Saída: resultados/sim/<id_celula>.parquet (ou sim_ptf/ com ZN_SOLO=ptf; uma linha por safra × tratamento) e
resultados/solos_usados.csv.
Uso: ../calibracao_mg/.venv/bin/python src/etapa2_simular.py [n_processos] [n_celulas]   (retoma: pula as feitas)
"""
import os
import sys
import time
from multiprocessing import current_process, get_context

import numpy as np
import pandas as pd

import config as cfg
import bs2019 as bs
import dssat_direto as dd
import gee
import ptf

F_MIN = 0.05
SOLO = os.environ.get("ZN_SOLO", "bs2019")
REESCALAR_AWC = False
DIR_SIM = cfg.DIR_RES / ("sim" if SOLO == "bs2019" else f"sim_{SOLO}")
ARQ_CEL = cfg.DIR_DADOS / "processados" / "celulas.csv"
ARQ_SOLO = cfg.DIR_DADOS / "processados" / "perfis_solo.csv"


def arq_clima(lon, lat):
    x, y = gee.centro_celula_xavier(lon, lat)
    return gee.DIR_XAVIER / f"xavier_px_{float(x):+.2f}_{float(y):+.2f}_1980_2025.csv.gz"


def agua_disponivel(cam, ate=100):
    """Água disponível (mm) de 0 a `ate` cm."""
    esp = np.diff(np.r_[0, cam.slb.to_numpy()])
    m = cam.slb.to_numpy() <= ate
    return float(np.sum((cam.sdul - cam.slll).to_numpy()[m] * esp[m] * 10))


def reescalar(cam, awc):
    f = awc * 1000 / agua_disponivel(cam)
    c = cam.copy()
    c["sdul"] = np.minimum(c.slll + f * (c.sdul - c.slll), c.ssat - 0.01).round(3)
    return c


def solos_celula(perfis):
    perfis = perfis[perfis.valido]
    perfis = perfis[perfis.area_ha / perfis.area_ha.sum() >= F_MIN]
    solos, meta = {}, []
    for k, (_, p) in enumerate(perfis.sort_values("area_ha", ascending=False).iterrows()):
        cam = bs.camadas_solo(p.solo) if SOLO == "bs2019" else ptf.camadas_perfil_mb(p)
        if SOLO == "ptf" and REESCALAR_AWC and pd.notna(p.get("awc_ana")):
            cam = reescalar(cam, p.awc_ana)
        if ptf.degenerado(cam):
            continue
        sid = f"Z{p.id_celula[1:]}{int(p.classe)}{k:05d}"[:10]
        solos[sid] = (cam, ptf.TEXTURA_DSSAT[p.solo])
        meta.append(dict(id_celula=p.id_celula, id_solo=sid, solo=p.solo, area_ha=p.area_ha,
                         ad_100cm_mm=agua_disponivel(cam), awc_ana=p.get("awc_ana"),
                         ad_zarc=p.ad_zarc, classe_ad=p.classe_ad))
    m = pd.DataFrame(meta)
    m["peso"] = m.area_ha / m.area_ha.sum()
    return solos, m


def rodar_celula(cel):
    try:
        return _rodar_celula(pd.Series(cel))
    except Exception as erro:              # uma célula com problema não derruba a rodada
        return cel["id_celula"], None, f"ERRO: {str(erro)[:300]}"


def _rodar_celula(cel):
    saida = DIR_SIM / f"{cel.id_celula}.parquet"
    if saida.exists():
        return cel.id_celula, None, 0.0
    t0 = time.time()
    pasta = cfg.DIR_RUN / f"p{current_process()._identity[0]:02d}"
    dd.preparar_pasta(pasta)
    clima = pd.read_csv(arq_clima(cel.lon, cel.lat), index_col=0, parse_dates=True).interpolate(limit=5).ffill().bfill()
    wsta = dd.escrever_clima(pasta, clima, cel.lat, cel.lon, cel.elev)
    perfis = pd.read_csv(ARQ_SOLO)
    solos, meta = solos_celula(perfis[perfis.id_celula == cel.id_celula])
    dd.escrever_solos(pasta, solos, cel.lat, cel.lon)
    s0 = meta.id_solo.iloc[0]
    trat = [dict(solo=s, gmr=g, semeadura=sem, agua=True) for s in solos for g in cfg.GMR for sem in cfg.DECENDIOS]
    trat += [dict(solo=s0, gmr=g, semeadura=sem, agua=False) for g in cfg.GMR for sem in cfg.DECENDIOS]
    s = dd.rodar_tratamentos(pasta, wsta, solos, trat)
    s["semeadura"] = s.semeadura.map(lambda t: f"{t[0]:02d}-{t[1]:02d}")
    s.insert(0, "id_celula", cel.id_celula)
    s.to_parquet(saida, index=False)
    return cel.id_celula, meta, time.time() - t0


if __name__ == "__main__":
    n_proc = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    DIR_SIM.mkdir(parents=True, exist_ok=True)
    c = pd.read_csv(ARQ_CEL)
    if len(sys.argv) > 2:
        c = c.iloc[:int(sys.argv[2])]
    c = c[[arq_clima(r.lon, r.lat).exists() for r in c.itertuples()]]
    print(f"{len(c)} células com clima; {sum((DIR_SIM / f'{i}.parquet').exists() for i in c.id_celula)} já feitas", flush=True)
    arq = cfg.DIR_RES / ("solos_usados.csv" if SOLO == "bs2019" else f"solos_usados_{SOLO}.csv")
    t0 = time.time()
    with get_context("fork").Pool(n_proc) as pool:
        for n, (cid, meta, dt) in enumerate(pool.imap_unordered(rodar_celula, c.to_dict("records")), 1):
            if isinstance(dt, str):
                print(f"{n}/{len(c)} {cid} {dt}", flush=True)
            elif meta is not None:          # grava a cada célula: os pesos dos solos não se perdem se a rodada cair
                antigo = pd.read_csv(arq) if arq.exists() else pd.DataFrame()
                pd.concat([antigo, meta]).drop_duplicates("id_solo", keep="last").to_csv(arq, index=False)
                print(f"{n}/{len(c)} {cid} {dt / 60:.1f} min ({len(meta)} solos); decorrido {(time.time() - t0) / 3600:.2f} h",
                      flush=True)
