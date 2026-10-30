"""Etapa 0: confere a execução direta contra o DSSATTools (caminho dos artigos) num tratamento da etapa0_custo.

Mesmo clima, solo (perfil MapBiomas de maior área), cultivar (rodar_produtividade.cultivar ≡ BRG do .CUL), semeadura,
início e condição inicial; 45 safras. As datas e a produtividade devem coincidir.
Uso: DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1 ../calibracao_mg/.venv/bin/python src/etapa0_validar.py
"""
from datetime import date

import numpy as np
import pandas as pd

import config as cfg
import bs2019 as bs
import dssat_direto as dd
import gerar_produtos as gp
import ptf
import rodar_produtividade as rp
from DSSATTools.filex import Field, Planting, SCGeneral, SCManagement, SCMethods, SCOptions, SCOutputs, SimulationControls
from DSSATTools.run import DSSAT
from DSSATTools.soil import SoilLayer, SoilProfile
from DSSATTools.filex import InitialConditions, InitialConditionsLayer
from etapa0_custo import LAT, LON, ELEV

GMR, SEM, SOLO = 6.5, (10, 11), "argiloso"

if __name__ == "__main__":
    rp.TH = dict(zip(gp.cf.NOMES, np.median(gp.posterior(), axis=0)))
    d = pd.read_csv(cfg.DIR_RES / "etapa0_custo_summary.csv")
    perfis = pd.read_csv(cfg.CALIB / "dados" / "processados" / "perfis_solo_mapbiomas_soja2016.csv")
    perfis = perfis[(perfis.lon == LON) & (perfis.lat == LAT)].reset_index(drop=True)
    k = int(perfis.index[perfis.solo == SOLO][0])
    cam = ptf.camadas_perfil_mb(perfis.loc[k])
    a = d[(d.solo == f"ZN{k:08d}") & (d.gmr == GMR) & (d.semeadura == str(SEM)) & d.agua].sort_values("PDAT")

    clima = pd.read_csv(cfg.CALIB / "dados" / "brutos" / "xavier" / f"xavier_no_{LON:+.2f}_{LAT:+.2f}_1980_2025.csv.gz",
                        index_col=0, parse_dates=True).interpolate(limit=5).ffill().bfill()
    ws = bs.estacao(clima, LAT, LON, elev=ELEV)
    sp = SoilProfile(table=[SoilLayer(**r) for r in cam[dd.COLS_SOLO].to_dict("records")], name="ZNVAL00001",
                     soil_data_source="MapBiomas", soil_clasification=ptf.TEXTURA_DSSAT[SOLO], soil_series_name="v",
                     site="Parana", country="Brazil", lat=LAT, long=LON, scom="R", salb=bs.ALBEDO, slu1=bs.SLU1,
                     sldr=bs.SLDR, slro=bs.CN, slnf=1.0, slpf=1.0)
    sdate = date(cfg.ANO_INI, *cfg.INICIO)
    ic = InitialConditions(pcr="SB", icdat=sdate, icres=0, icren=0, icrip=0, icrid=0, icrn=1, icre=1,
                           table=[InitialConditionsLayer(icbl=b, sh2o=round(l + bs.FRACAO_AGUA_INICIAL * (u - l), 3),
                                                         snh4=0.5, sno3=2.0) for b, l, u in zip(cam.slb, cam.slll, cam.sdul)])
    ctrl = SimulationControls(general=SCGeneral(sdate=sdate, nyers=cfg.NYERS),
                              options=SCOptions(water="Y", nitro="N", symbi="N", co2="M"),
                              methods=SCMethods(evapo=bs.EVAPO), management=SCManagement(irrig="N", ferti="N", harvs="M"),
                              outputs=SCOutputs(grout="N", waout="N", caout="N"))
    campo = Field(id_field="BSXV0001", wsta=ws, id_soil=sp, flob=0, fldt="DR000", fldd=0, flds=0,
                  sltx=ptf.TEXTURA_DSSAT[SOLO], sldp=200)
    plant = Planting(pdate=date(cfg.ANO_INI, *SEM), ppop=bs.POPULACAO, ppoe=bs.POPULACAO, plrs=bs.ESPACAMENTO,
                     pldp=bs.PROF_SEMEADURA)
    D = DSSAT(str(cfg.DIR_RUN / "validar"))
    D.run_treatment(field=campo, cultivar=rp.cultivar(GMR), planting=plant, initial_conditions=ic,
                    simulation_controls=ctrl, verbose=False)
    b = bs._ler_summary(D.output_files["Summary"]).sort_values("PDAT")
    for c in ["PDAT", "ADAT", "MDAT", "HWAM"]:
        x, y = a[c].to_numpy(), b[c].to_numpy()
        print(f"{c}: {len(x)} × {len(y)} safras; iguais {np.sum(x == y)}; máx. dif. {np.nanmax(np.abs(x - y)):.0f}; "
              f"máx. dif. relativa {np.nanmax(np.abs(x - y) / np.maximum(np.abs(y), 1)):.4f}")
