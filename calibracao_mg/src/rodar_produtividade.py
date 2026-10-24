"""Produtividade potencial (Yp) e limitada por água (Yw) por ponto × safra 2000/01–2023/24 (opção A, 18/10/2026).

Por ponto-safra com fenologia MODIS invertida pelo SBI:
  - GMR efetivo e semeadura: mediana da posterior do SBI (resultados/sbi/gmr_efetivo_<sbi>.csv; s = dias após 10/set);
  - fenologia: θ(GMR) na mediana da posterior da calibração (cadeia; δ8 não entra no DSSAT);
  - crescimento: BRS 284 de Battisti & Sentelhas (2019, Tab. S2; bs2019.COEF_S2_COMUNS e ECO_S2);
  - solo: os 3 solos da Tab. S1 com peso pela fração de textura do MapBiomas na célula de 1° (base v1);
  - manejo, população, ET (Priestley-Taylor) e condição inicial: configuração validada da replicação (bs2019).
Uma safra por rodada do DSSAT (início 30 d antes da semeadura).

Saída: resultados/produtividade/sim/<lon>_<lat>.csv (uma linha por safra × solo × cenário) e, ao final,
resultados/produtividade/yp_yw_pontos.csv (média ponderada pelos solos).

Uso: DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1 .venv/bin/python src/rodar_produtividade.py CADEIA.h5 [sbi=final] [n_processos]\n(FAIXA_FISIO precisa estar definido para a v12: os coeficientes são recortados em cf.FAIXAS)
"""
import os
import shutil
import sys
import warnings
from datetime import date, timedelta
from multiprocessing import current_process, get_context
from pathlib import Path

import emcee
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))
import bs2019 as bs  # noqa: E402
import calibrar_fenologia as cf  # noqa: E402
import rodar_bs2019 as rb  # noqa: E402
import treino_fenologia as tf  # noqa: E402
from DSSATTools.filex import (Field, Planting, SCGeneral, SCManagement, SCMethods, SCOptions,  # noqa: E402
                              SCOutputs, SimulationControls)
from DSSATTools.run import DSSAT  # noqa: E402

RAIZ = tf.RAIZ
DIR = RAIZ / "resultados" / "produtividade"
SOLOS = ["argiloso", "argiloarenoso", "francoarenoso"]
F_MIN = 0.05        # fração mínima de um solo na célula para entrar no Yw (21/10: corta ~45% das rodadas)
_D = None
TH = None        # θ(GMR) na mediana da posterior (definido no processo principal e herdado pelo fork)


def theta_mediana(cadeia, descarte=0.5):
    ch = emcee.backends.HDFBackend(str(cadeia), read_only=True).get_chain()
    ch = ch[int(ch.shape[0] * descarte):].reshape(-1, ch.shape[2])
    return dict(zip(cf.NOMES, np.median(ch, axis=0)))


def coeficientes(gmr):
    return {k: float(np.clip(TH[f"a_{k}"] + TH[f"b_{k}"] * (gmr - 7), *cf.FAIXAS[k])) for k in cf.FAIXAS}


def cultivar(gmr):
    """Fenologia calibrada + crescimento da BRS 284 (B&S 2019). FL-SH e FL-LF da Tab. S2 (MG mais próximo): não alteram
    R1/R5/R7, mas o FL-SH genérico (10 d em MG ≥ 7) maior que FL-SD impede o enchimento de grão (Yp = 0)."""
    c = tf.cultivar(dict(gmr=gmr, **coeficientes(gmr)))
    mg = min(bs.COEF_S2, key=lambda m: abs(m - gmr))
    c["fl-sh"], c["fl-lf"] = bs.COEF_S2[mg][3], bs.COEF_S2[mg][6]
    return c


def _ini():
    global _D
    pasta = rb.DIR_RUN / f"y{current_process()._identity[0]:02d}"
    shutil.rmtree(pasta, ignore_errors=True)
    _D = DSSAT(str(pasta))


def rodar_ponto(args):
    lon, lat, elev, linhas = args
    saida = DIR / "sim" / f"{lon:+.1f}_{lat:+.1f}.csv"
    if saida.exists():
        return "pulado"
    loc = dict(lon=lon, lat=lat, elev=elev, arq="no")
    clima = pd.read_csv(tf.arq_clima(loc), index_col=0, parse_dates=True).interpolate(limit=5).ffill().bfill()
    ws = bs.estacao(clima, lat, lon, elev=elev)
    out = []
    for r in linhas:
        pdate = date(int(r["safra"]), 9, 10) + timedelta(days=int(round(r["s_q50"])))
        sdate = pdate - timedelta(days=30)
        cul = cultivar(r["g_q50"])
        # Yp (sem estresse hídrico) não depende do solo: uma rodada por safra, no solo de maior fração (peso 1).
        # Yw: só os solos com fração ≥ F_MIN; os pesos são renormalizados em consolidar().
        dom = max(SOLOS, key=lambda k: r[f"f_{k}"])
        casos = [(dom, False)] + [(k, True) for k in SOLOS if r[f"f_{k}"] >= F_MIN or k == dom]
        for solo, agua in casos:
            ctrl = SimulationControls(
                general=SCGeneral(sdate=sdate, nyers=1),
                options=SCOptions(water="Y" if agua else "N", nitro="N", symbi="N", co2="M"),
                methods=SCMethods(evapo=bs.EVAPO), management=SCManagement(irrig="N", ferti="N", harvs="M"),
                outputs=SCOutputs(grout="N", waout="N", caout="N"))
            campo = Field(id_field="BSXV0001", wsta=ws, id_soil=bs.perfil_solo(solo, lat, lon), flob=0,
                          fldt="DR000", fldd=0, flds=0, sltx=bs.TEXTURA_DSSAT[solo], sldp=200)
            plant = Planting(pdate=pdate, ppop=bs.POPULACAO, ppoe=bs.POPULACAO, plrs=bs.ESPACAMENTO,
                             pldp=bs.PROF_SEMEADURA)
            try:
                _D.run_treatment(field=campo, cultivar=cul, planting=plant,
                                 initial_conditions=bs.condicoes_iniciais(solo, sdate),
                                 simulation_controls=ctrl, verbose=False)
                s = bs._ler_summary(_D.output_files["Summary"]).iloc[0]
                hw, erro = s.get("HWAM", np.nan), ""
            except Exception as e:      # noqa: BLE001
                hw, erro, s = np.nan, str(e)[:60], pd.Series(dtype=float)
            out.append(dict(lon=lon, lat=lat, safra=int(r["safra"]), gmr=r["g_q50"], sem_dias=r["s_q50"],
                            data_semeadura=pdate.isoformat(), solo=solo if agua else "todos",
                            peso=r[f"f_{solo}"] if agua else 1.0,
                            cenario="yw" if agua else "yp", hwam=hw, mdat=s.get("MDAT", np.nan),
                            prcp=s.get("PRCP", np.nan), etcp=s.get("ETCP", np.nan), erro=erro))
    pd.DataFrame(out).to_csv(saida, index=False)
    return "ok"


def tarefas(sbi):
    ge = pd.read_csv(RAIZ / "resultados" / "sbi" / f"gmr_efetivo_{sbi}.csv")
    pts = pd.read_csv(RAIZ / "dados" / "processados" / "pontos_bs2019_digitalizados.csv")[["lon", "lat", "elev"]]
    fr = pd.read_csv(RAIZ / "dados" / "processados" / "fracoes_solo_mapbiomas_soja2016.csv")
    ge = ge.merge(pts.drop_duplicates(["lon", "lat"]), on=["lon", "lat"], how="left")
    ge = ge.merge(fr[["lon", "lat"] + [f"f_{s}" for s in SOLOS]], on=["lon", "lat"], how="left")
    sem_solo = ge[f"f_{SOLOS[0]}"].isna()
    ge.loc[sem_solo, [f"f_{s}" for s in SOLOS]] = [0.0, 1.0, 0.0]      # sem fração do MapBiomas: argiloarenoso
    return [(lon, lat, float(g.elev.iloc[0]), g.to_dict("records")) for (lon, lat), g in ge.groupby(["lon", "lat"])], int(sem_solo.sum())


def consolidar():
    d = pd.concat([pd.read_csv(f) for f in sorted((DIR / "sim").glob("*.csv"))])
    d.to_csv(DIR / "yp_yw_simulacoes.csv", index=False)
    d = d.dropna(subset=["hwam"])
    d["peso_n"] = d.peso / d.groupby(["lon", "lat", "safra", "cenario"]).peso.transform("sum")
    m = (d.assign(v=d.hwam * d.peso_n).groupby(["lon", "lat", "safra", "gmr", "sem_dias", "cenario"]).v.sum()
         .unstack("cenario").reset_index())
    m.to_csv(DIR / "yp_yw_pontos.csv", index=False)
    return m


if __name__ == "__main__":
    cadeia = Path(sys.argv[1])
    sbi = sys.argv[2] if len(sys.argv) > 2 else "final"
    n = int(sys.argv[3]) if len(sys.argv) > 3 else os.cpu_count()
    (DIR / "sim").mkdir(parents=True, exist_ok=True)
    TH = theta_mediana(cadeia)
    print({g: {k: round(v, 3) for k, v in coeficientes(g).items()} for g in (5, 7, 9)}, flush=True)
    lista, sem_solo = tarefas(sbi)
    print(f"{len(lista)} pontos; {sum(len(x[3]) for x in lista)} ponto-safra; {sem_solo} sem fração de solo", flush=True)
    with get_context("fork").Pool(n, initializer=_ini) as pool:
        for k, r in enumerate(pool.imap_unordered(rodar_ponto, lista), 1):
            if k % 10 == 0:
                print(f"{k}/{len(lista)} pontos", flush=True)
    m = consolidar()
    print(m.describe().round(0).to_string(), flush=True)
