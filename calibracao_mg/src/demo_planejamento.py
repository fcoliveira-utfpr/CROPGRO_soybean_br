"""Demonstração de uso das cultivares por ponto no planejamento (24/10/2026).

Em 5 pontos representativos, simula 7 datas de semeadura (15/set a 15/dez) × 3 GMR (o efetivo das colheitas
2021–2024 arredondado a 0,5, e ±1) × 24 safras (2000/01–2023/24), com a cultivar brasileira do GMR (fenologia
calibrada + crescimento da BRS 284) e os solos do ponto ponderados pela área (≥ 5%). Saídas por safra: datas de R1,
R7 e colheita (R8), Yp e Yw.

Saída: resultados/planejamento/demo_planejamento.csv
Uso: DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1 .venv/bin/python src/demo_planejamento.py [n_processos]
"""
import os
import shutil
import sys
import warnings
from datetime import date, timedelta
from multiprocessing import current_process, get_context
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))
import bs2019 as bs  # noqa: E402
import gerar_produtos as gp  # noqa: E402
import rodar_produtividade as rp  # noqa: E402
import treino_fenologia as tf  # noqa: E402
from DSSATTools.filex import (Field, Planting, SCGeneral, SCManagement, SCMethods, SCOptions,  # noqa: E402
                              SCOutputs, SimulationControls)
from DSSATTools.run import DSSAT  # noqa: E402

RAIZ = tf.RAIZ
OUT = RAIZ / "resultados" / "planejamento"
PONTOS = {"South (RS)": (-53.5, -28.5), "West Paraná": (-53.5, -24.5), "Goiás": (-51.5, -17.5),
          "Mato Grosso": (-55.5, -12.5), "MATOPIBA (BA)": (-46.5, -11.5)}
SEMEADURAS = [(9, 15), (10, 1), (10, 15), (11, 1), (11, 15), (12, 1), (12, 15)]
SAFRAS = range(2000, 2024)
_D = None


def _ini():
    global _D
    pasta = rp.rb.DIR_RUN / f"pl{current_process()._identity[0]:02d}"
    shutil.rmtree(pasta, ignore_errors=True)
    _D = DSSAT(str(pasta))


def resumo(txt):
    linhas = txt.splitlines()
    i = next(k for k, l in enumerate(linhas) if l.startswith("@"))
    cab = linhas[i].replace("@", " ").split()
    dados = [l.split() for l in linhas[i + 1:] if l.strip()]
    df = pd.DataFrame([x for x in dados if len(x) == len(cab)], columns=cab)
    return df[["PDAT", "ADAT", "MDAT", "HDAT", "HWAM"]].apply(pd.to_numeric, errors="coerce").replace(-99, np.nan).iloc[0]


def dias(yyddd, p):
    if pd.isna(yyddd):
        return np.nan
    return (pd.to_datetime(str(int(yyddd)), format="%Y%j") - pd.Timestamp(p)).days


def rodar(args):
    nome, lon, lat, elev, solos, g, mes, dia = args
    loc = dict(lon=lon, lat=lat, elev=elev, arq="no")
    clima = pd.read_csv(tf.arq_clima(loc), index_col=0, parse_dates=True).interpolate(limit=5).ffill().bfill()
    ws = bs.estacao(clima, lat, lon, elev=elev)
    cul = rp.cultivar(g)
    out = []
    for safra in SAFRAS:
        pdate = date(safra + (mes < 7), mes, dia)
        sdate = pdate - timedelta(days=30)
        dom = max(solos, key=solos.get)
        for solo, agua in [(dom, False)] + [(k, True) for k in solos]:
            ctrl = SimulationControls(general=SCGeneral(sdate=sdate, nyers=1),
                                      options=SCOptions(water="Y" if agua else "N", nitro="N", symbi="N", co2="M"),
                                      methods=SCMethods(evapo=bs.EVAPO),
                                      management=SCManagement(irrig="N", ferti="N", harvs="M"),
                                      outputs=SCOutputs(grout="N", waout="N", caout="N"))
            campo = Field(id_field="BSXV0001", wsta=ws, id_soil=bs.perfil_solo(solo, lat, lon), flob=0, fldt="DR000",
                          fldd=0, flds=0, sltx=bs.TEXTURA_DSSAT[solo], sldp=200)
            plant = Planting(pdate=pdate, ppop=bs.POPULACAO, ppoe=bs.POPULACAO, plrs=bs.ESPACAMENTO, pldp=bs.PROF_SEMEADURA)
            try:
                _D.run_treatment(field=campo, cultivar=cul, planting=plant,
                                 initial_conditions=bs.condicoes_iniciais(solo, sdate), simulation_controls=ctrl, verbose=False)
                s = resumo(_D.output_files["Summary"])
            except Exception:  # noqa: BLE001
                continue
            out.append(dict(ponto=nome, lon=lon, lat=lat, rmg=g, sowing=f"{dia:02d}/{mes:02d}", sowing_doy=pdate.timetuple().tm_yday,
                            season=safra, sowing_date=pdate.isoformat(), soil=solo if agua else "all",
                            weight=solos[solo] if agua else 1.0, scenario="yw" if agua else "yp",
                            r1_dap=dias(s.ADAT, pdate), r7_dap=dias(s.MDAT, pdate), harvest_dap=dias(s.HDAT, pdate),
                            harvest_date=None if pd.isna(s.HDAT) else pd.to_datetime(str(int(s.HDAT)), format="%Y%j").date().isoformat(),
                            yield_kg_ha=s.HWAM))
    return out


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    OUT.mkdir(parents=True, exist_ok=True)
    rp.TH = dict(zip(gp.cf.NOMES, np.median(gp.posterior(), axis=0)))
    pts = pd.read_csv(RAIZ / "produtos" / "dados" / "coefficients_point.csv")
    fr = pd.read_csv(RAIZ / "dados" / "processados" / "fracoes_solo_mapbiomas_soja2016.csv")
    tarefas = []
    for nome, (lon, lat) in PONTOS.items():
        r = pts[(pts.lon == lon) & (pts.lat == lat)].iloc[0]
        f = fr[(fr.lon == lon) & (fr.lat == lat)].iloc[0]
        solos = {k: float(f[f"f_{k}"]) for k in rp.SOLOS if f[f"f_{k}"] >= rp.F_MIN}
        tot = sum(solos.values())
        solos = {k: v / tot for k, v in solos.items()}
        g0 = round(r.rmg_recent * 2) / 2
        for g in (g0 - 1, g0, g0 + 1):
            for mes, dia in SEMEADURAS:
                tarefas.append((nome, lon, lat, float(r.elevation_m), solos, float(np.clip(g, 5, 9.5)), mes, dia))
    print(len(tarefas), "tarefas", flush=True)
    linhas = []
    with get_context("fork").Pool(n, initializer=_ini) as pool:
        for k, res in enumerate(pool.imap_unordered(rodar, tarefas), 1):
            linhas += res
            if k % 10 == 0:
                print(f"{k}/{len(tarefas)}", flush=True)
                pd.DataFrame(linhas).to_csv(OUT / "demo_planejamento.csv", index=False)
    pd.DataFrame(linhas).to_csv(OUT / "demo_planejamento.csv", index=False)
    print("fim", flush=True)
