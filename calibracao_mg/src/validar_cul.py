"""Valida os arquivos SBGRO048_BR.CUL/.ECO: a mesma safra simulada (i) com a cultivar montada pelo pipeline
(rodar_produtividade.cultivar, via DSSATTools) e (ii) com o código BR do ponto lido dos arquivos gerados, chamando o
executável do DSSAT diretamente. As datas e a produtividade devem coincidir.
Uso: DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1 .venv/bin/python src/validar_cul.py [n_pontos]
"""
import re
import shutil
import subprocess
import sys
import warnings
from datetime import date, timedelta
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))
import bs2019 as bs  # noqa: E402
import gerar_produtos as gp  # noqa: E402
import rodar_produtividade as rp  # noqa: E402
import treino_fenologia as tf  # noqa: E402
from DSSATTools import run as dr  # noqa: E402
from DSSATTools.filex import (Field, Planting, SCGeneral, SCManagement, SCMethods, SCOptions,  # noqa: E402
                              SCOutputs, SimulationControls)
from DSSATTools.run import DSSAT  # noqa: E402

PROD = tf.RAIZ / "produtos" / "dados"

if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    pts = pd.read_csv(PROD / "coefficients_point.csv")
    rp.TH = dict(zip(gp.cf.NOMES, np.median(gp.posterior(), axis=0)))
    pasta = rp.rb.DIR_RUN / "valcul"
    shutil.rmtree(pasta, ignore_errors=True)
    D = DSSAT(str(pasta))
    res = []
    for r in pts.sample(n, random_state=3).itertuples():
        loc = dict(lon=r.lon, lat=r.lat, elev=r.elevation_m, arq="no")
        clima = pd.read_csv(tf.arq_clima(loc), index_col=0, parse_dates=True).interpolate(limit=5).ffill().bfill()
        ws = bs.estacao(clima, r.lat, r.lon, elev=r.elevation_m)
        pdate = date(2019, 9, 10) + timedelta(days=int(round(r.sowing_days_mean)))
        sdate = pdate - timedelta(days=30)
        ctrl = SimulationControls(general=SCGeneral(sdate=sdate, nyers=1),
                                  options=SCOptions(water="Y", nitro="N", symbi="N", co2="M"),
                                  methods=SCMethods(evapo=bs.EVAPO), management=SCManagement(irrig="N", ferti="N", harvs="M"),
                                  outputs=SCOutputs(grout="N", waout="N", caout="N"))
        campo = Field(id_field="BSXV0001", wsta=ws, id_soil=bs.perfil_solo("argiloarenoso", r.lat, r.lon), flob=0,
                      fldt="DR000", fldd=0, flds=0, sltx="SC", sldp=200)
        plant = Planting(pdate=pdate, ppop=bs.POPULACAO, ppoe=bs.POPULACAO, plrs=bs.ESPACAMENTO, pldp=bs.PROF_SEMEADURA)
        D.run_treatment(field=campo, cultivar=rp.cultivar(r.rmg_for_cul), planting=plant,
                        initial_conditions=bs.condicoes_iniciais("argiloarenoso", sdate), simulation_controls=ctrl,
                        verbose=False)
        a = bs._ler_summary(D.output_files["Summary"]).iloc[0]
        # mesma rodada com os arquivos gerados: troca CUL/ECO e o código da cultivar no arquivo X
        shutil.copy(PROD / "SBGRO048_BR.CUL", pasta / "SBGRO048.CUL")
        shutil.copy(PROD / "SBGRO048_BR.ECO", pasta / "SBGRO048.ECO")
        fx = next(pasta.glob("*.SBX"))
        txt = fx.read_text()
        txt = re.sub(r"(\*CULTIVARS\n@C CR INGENO CNAME\n 1 SB )\S+", lambda m: m.group(1) + r.cultivar_id, txt)
        fx.write_text(txt)
        subprocess.run([dr.BIN_PATH, "C", fx.name, "1"], cwd=pasta, capture_output=True)
        b = bs._ler_summary((pasta / "Summary.OUT").read_text()).iloc[0]
        res.append(dict(id=r.cultivar_id, rmg=r.rmg_for_cul, adat_pipe=a.ADAT, adat_cul=b.ADAT, mdat_pipe=a.MDAT,
                        mdat_cul=b.MDAT, hwam_pipe=a.HWAM, hwam_cul=b.HWAM))
    d = pd.DataFrame(res)
    print(d.to_string(index=False))
    ok = (d.adat_pipe == d.adat_cul).all() and (d.mdat_pipe == d.mdat_cul).all() and \
        (np.abs(d.hwam_pipe - d.hwam_cul) <= 0.01 * d.hwam_pipe + 1).all()
    print("VALIDADO" if ok else "DIFERENÇAS ENCONTRADAS")
