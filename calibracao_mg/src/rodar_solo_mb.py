"""Roda a configuração validada de Battisti & Sentelhas (2019, v2) trocando os 3 solos fixos pelos perfis do
MapBiomas Solo (+ SoilGrids) com pedotransferência de Medrado & Lima (2014), só para a produtividade atingível.

Cada célula tem até 3 perfis (classes de argila), com a área de soja de cada um. Saída por célula em
resultados/solo_mapbiomas/sim_cc<h>kPa/<lon>_<lat>.csv.

Uso: .venv/bin/python src/rodar_solo_mb.py <h_cc em kPa> [n_processos]
"""
import os
import shutil
import sys
import warnings
from multiprocessing import current_process, get_context
from pathlib import Path

import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))
import bs2019 as bs  # noqa: E402
import ptf  # noqa: E402
import rodar_bs2019 as rb  # noqa: E402
from DSSATTools.filex import InitialConditions, InitialConditionsLayer  # noqa: E402
from DSSATTools.run import DSSAT  # noqa: E402
from DSSATTools.soil import SoilLayer, SoilProfile  # noqa: E402

RAIZ = rb.RAIZ
PERFIS = RAIZ / "dados" / "processados" / "perfis_solo_mapbiomas_soja2016.csv"
H_CC = float(sys.argv[1]) if len(sys.argv) > 1 else 10.0
DIR_SAIDA = RAIZ / "resultados" / "solo_mapbiomas" / f"sim_cc{H_CC:g}kPa"
_D = None


def _ini():
    global _D
    pasta = rb.DIR_RUN / f"m{current_process()._identity[0]:02d}"
    shutil.rmtree(pasta, ignore_errors=True)
    _D = DSSAT(str(pasta))


def _perfil_dssat(camadas, nome, lat, lon, classe):
    cols = ["slb", "slll", "sdul", "ssat", "sbdm", "sloc", "slcl", "slsi", "slhw", "srgf", "ssks", "slni"]
    return SoilProfile(table=[SoilLayer(**r) for r in camadas[cols].to_dict("records")], name=nome,
                       soil_data_source="MapBiomas+SG", soil_clasification=ptf.TEXTURA_DSSAT[classe],
                       soil_series_name=classe, site="Brasil", country="Brazil", lat=lat, long=lon, scom="R",
                       salb=bs.ALBEDO, slu1=bs.SLU1, sldr=bs.SLDR, slro=bs.CN, slnf=1.0, slpf=1.0)


def rodar_celula(args):
    lon, lat, elev, perfis = args
    saida = DIR_SAIDA / f"{lon:+.1f}_{lat:+.1f}.csv"
    if saida.exists():
        return "pulado"
    clima = pd.read_csv(rb.arq_clima(lon, lat), index_col=0, parse_dates=True).interpolate(limit=5).ffill().bfill()
    ws = bs.estacao(clima, lat, lon, elev=elev)
    partes = []
    for k, p in enumerate(perfis):
        p = pd.Series(p)
        cam = ptf.camadas_perfil_mb(p, h_cc=H_CC)
        if ptf.degenerado(cam):
            continue
        nome = f"MB{k}{abs(int(lon * 10)):04d}{abs(int(lat * 10)):03d}"[:10]
        # substitui temporariamente o perfil e as condições iniciais do bs2019 pelos do MapBiomas
        bs.perfil_solo = lambda _n, la, lo, cam=cam, nome=nome, cl=p.solo: _perfil_dssat(cam, nome, la, lo, cl)
        bs.condicoes_iniciais = lambda _n, dia, cam=cam: InitialConditions(
            pcr="SB", icdat=dia, icres=0, icren=0, icrip=0, icrid=0, icrn=1, icre=1,
            table=[InitialConditionsLayer(icbl=b, sh2o=round(l + bs.FRACAO_AGUA_INICIAL * (u - l), 3), snh4=0.5, sno3=2.0)
                   for b, l, u in zip(cam.slb, cam.slll, cam.sdul)])
        bs.TEXTURA_DSSAT[p.solo] = ptf.TEXTURA_DSSAT[p.solo]
        for sem in bs.SEMEADURAS:
            s = bs.simular(_D, ws, p.solo, lat, lon, sem, agua=True)
            partes.append(s.assign(area_ha=p.area_ha, ad80=float(((cam.sdul - cam.slll) * 100)[:8].sum())))
    if partes:
        pd.concat(partes).to_csv(saida, index=False)
    return "ok"


if __name__ == "__main__":
    n = int(sys.argv[2]) if len(sys.argv) > 2 else os.cpu_count()
    DIR_SAIDA.mkdir(parents=True, exist_ok=True)
    pf = pd.read_csv(PERFIS)
    pts = pd.read_csv(RAIZ / "dados" / "processados" / "pontos_bs2019_digitalizados.csv")
    pf = pf.merge(pts[["lon", "lat", "elev"]], on=["lon", "lat"])
    tarefas = [(lon, lat, g.elev.iloc[0], g.to_dict("records")) for (lon, lat), g in pf.groupby(["lon", "lat"])]
    print(f"CC = {H_CC:g} kPa | {len(tarefas)} células | {len(pf)} perfis | {n} processos", flush=True)
    with get_context("fork").Pool(n, initializer=_ini) as pool:
        for k, _ in enumerate(pool.imap_unordered(rodar_celula, tarefas), 1):
            if k % 25 == 0 or k == len(tarefas):
                print(f"{k}/{len(tarefas)}", flush=True)
