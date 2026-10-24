"""Conjunto de treino do emulador de FENOLOGIA do CROPGRO-Soybean.

- Coeficientes: CSDL, PPSEN, EM-FL, FL-SD, SD-PM (faixas do Morris) e GMR contínuo (4,8–9,3), com o ecótipo do
  M GROUP inteiro mais próximo (5 a 9). Crescimento fixo nos valores de Battisti & Sentelhas (2019, Tabela S2).
- Semeadura sorteada entre 10/set e 28/fev; safras 1980/81 a 2024/25; sem balanço hídrico (fenologia potencial).
- Locais: 3 experimentos (Jaboticabal, Ijuí, Paragominas), 18 pontos representativos dos catálogos da Embrapa
  (Centro-Sul por faixa de altitude, GO, MT/RO) e 30 pontos de validação do satélite (estratificados por zona).

Saída por safra: dias após a semeadura até emergência, R1 e R7.
Uso: .venv/bin/python src/treino_fenologia.py [n_amostras_por_local] [n_processos] [lote]
Lotes > 1 acrescentam amostras (ids a partir de lote × 100000, outra semente) sem apagar as dos lotes anteriores.
"""
import os
import shutil
import sys
import time
import warnings
from datetime import date, timedelta
from multiprocessing import current_process, get_context
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import qmc

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))
import bs2019 as bs  # noqa: E402
import gee  # noqa: E402
import morris as mo  # noqa: E402
import rodar_bs2019 as rb  # noqa: E402
from DSSATTools.crop import Soybean  # noqa: E402
from DSSATTools.filex import (Field, Planting, SCGeneral, SCManagement, SCMethods, SCOptions,  # noqa: E402
                              SCOutputs, SimulationControls)
from DSSATTools.run import DSSAT  # noqa: E402

RAIZ = rb.RAIZ
DIR = RAIZ / "resultados" / "emulador_fenologia"
DIR_TREINO = DIR / "treino"
# Lote 1: PPSEN 0,20–0,385 e EM-FL 16,5–30. Lote 2 (16/10/2026): faixas ampliadas, porque a calibração encostou
# em PPSEN mínimo e EM-FL máximo; PPSEN até o MINIMA do DSSAT (0,129) e EM-FL além do MAXIMA (30) até 36.
# Lote 6 (17/10/2026, v6): faixas ampliadas outra vez (a v5 encostou em CSDL 13, FL-SD 18,5, SD-PM 25 e EM-FL 36 nos
# extremos de GMR) e R1PPO do ecótipo calibrado (a resposta à época de semeadura ficou fraca: semeaduras antecipadas
# 10–34 d mais longas que o simulado). Nos lotes 1–5 o R1PPO é o do ecótipo do MG inteiro (r1ppo_padrao).
PARAM = [("csdl", 11.5, 14.0), ("ppsen", 0.129, 0.385), ("em-fl", 16.5, 40.0), ("fl-sd", 9.0, 22.0),
         ("sd-pm", 20.0, 38.0), ("gmr", 4.8, 9.3), ("doy_sem", 0, 171), ("r1ppo", 0.2, 1.2)]   # doy_sem: dias após 10/set
ANO_INI, NYERS = 1980, 45
EXPERIMENTOS = {"Jaboticabal-SP": (-48.25, -21.25, 600), "Ijuí-RS": (-53.85, -28.45, 370), "Paragominas-PA": (-46.95, -3.05, 122)}
# Lote 3: ensaios do MATOPIBA de Silva et al. (2022, PLOS ONE; dados em dados/brutos/dataverse_C7KUE1) e outros, centro da célula
# do Xavier que contém o ensaio. Ficam no fim da lista para não mudar os índices dos 47 locais dos lotes 1 e 2.
EXPERIMENTOS_LOTE3 = {"Tasso Fragoso-MA": (-45.75, -8.45, 242), "São Raimundo das Mangabeiras-MA": (-45.45, -7.05, 234),
                      "Uruçuí-PI": (-44.55, -7.25, 167), "Baixa Grande do Ribeiro-PI": (-45.25, -7.85, 325),
                      "Balsas-MA": (-46.05, -7.55, 247), "Chapadinha-MA": (-43.35, -3.75, 105),
                      # Barbosa et al. (2013, Semina) e Cruz et al. (2010, Sci. Agrar.)
                      "Umuarama-PR": (-53.35, -23.75, 480), "São Desidério-BA": (-45.95, -12.75, 728),
                      # Alliprandini et al. (2009, Crop Sci.): 8 locais do Centro-Oeste e 7 do Sul (altitude da Fig. 1)
                      "Sorriso-MT": (-55.75, -12.55, 400), "Primavera do Leste-MT": (-54.25, -15.55, 500),
                      "Campo Verde-MT": (-55.15, -15.55, 800), "Goiatuba-GO": (-49.35, -18.05, 800),
                      "Morrinhos-GO": (-49.15, -17.75, 900), "Uberlândia-MG": (-48.25, -18.95, 900),
                      "Brasília-DF": (-47.85, -15.75, 1000), "São Gotardo-MG": (-46.05, -19.35, 1000),
                      "Cruz Alta-RS": (-53.65, -28.65, 400), "Santa Tereza do Oeste-PR": (-53.65, -25.05, 800),
                      "Cascavel-PR": (-53.45, -24.95, 800), "Palotina-PR": (-53.85, -24.25, 400),
                      "Cruzeiro do Sul-PR": (-52.15, -22.95, 500), "Ponta Grossa-PR": (-50.15, -25.05, 1000),
                      "Morro Agudo-SP": (-48.05, -20.75, 600)}
# Lote 4 (17/10/2026): Strieder et al. (2013, Passo Fundo) e Centro-Norte de MS (catálogo MS 2010/11; Santana et al. 2022)
EXPERIMENTOS_LOTE4 = {"Passo Fundo-RS": (-52.45, -28.25, 687), "São Gabriel do Oeste-MS": (-54.55, -19.35, 658),
                      "Chapadão do Sul-MS": (-52.65, -18.75, 810)}
SEMENTE = 7
_D, _WS = None, {}


def locais():
    pts = pd.read_csv(RAIZ / "dados" / "processados" / "pontos_bs2019_digitalizados.csv")
    fr = pd.read_csv(RAIZ / "dados" / "processados" / "fracoes_solo_mapbiomas_soja2016.csv")
    sat = pd.read_csv(RAIZ / "resultados" / "fenologia_satelite" / "fenologia_modis_pontos.csv").dropna(subset=["sos_p50"])
    pts = pts.merge(fr[["lon", "lat", "soja_ha"]], on=["lon", "lat"])
    rng = np.random.default_rng(SEMENTE)
    L = [dict(local=k, lon=v[0], lat=v[1], elev=v[2], tipo="experimento", arq="px") for k, v in EXPERIMENTOS.items()]
    # Catálogos: Centro-Sul (MR 1–3) por faixa de altitude; GO; MT/RO
    cs = pts[pts.lat.between(-31, -22) & pts.lon.between(-57, -47)]
    for nome, (a, b) in {"CS<500": (0, 500), "CS500-700": (500, 700), "CS>700": (700, 2000)}.items():
        s = cs[cs.elev.between(a, b)].nlargest(4, "soja_ha")
        L += [dict(local=nome, lon=r.lon, lat=r.lat, elev=r.elev, tipo="catalogo", arq="no") for r in s.itertuples()]
    for nome, (la, lo) in {"GO": ((-19.5, -13), (-53, -46)), "MT/RO": ((-17, -8), (-62, -52))}.items():
        s = pts[pts.lat.between(*la) & pts.lon.between(*lo)].nlargest(3, "soja_ha")
        L += [dict(local=nome, lon=r.lon, lat=r.lat, elev=r.elev, tipo="catalogo", arq="no") for r in s.itertuples()]
    # Validação com satélite: 3 pontos por zona, com mais safras observadas
    n_obs = sat.groupby(["lon", "lat"]).size().rename("n_obs").reset_index()
    v = pts.merge(n_obs, on=["lon", "lat"]).dropna(subset=["zona_bs"])
    ja = {(x["lon"], x["lat"]) for x in L}
    v = v[[(a, b) not in ja for a, b in zip(v.lon, v.lat)]]
    v = v.sort_values("n_obs", ascending=False).groupby("zona_bs").head(3)
    L += [dict(local=f"sat_z{int(r.zona_bs)}", lon=r.lon, lat=r.lat, elev=r.elev, tipo="satelite", arq="no") for r in v.itertuples()]
    for d in L:
        d["lote_local"] = 1
    L += [dict(local=k, lon=v[0], lat=v[1], elev=v[2], tipo="experimento", arq="px", lote_local=3) for k, v in EXPERIMENTOS_LOTE3.items()]
    L += [dict(local=k, lon=v[0], lat=v[1], elev=v[2], tipo="experimento", arq="px", lote_local=4) for k, v in EXPERIMENTOS_LOTE4.items()]
    return pd.DataFrame(L).reset_index(drop=True)


def arq_clima(loc):
    if loc["arq"] == "px":
        return gee.DIR_XAVIER / f"xavier_px_{loc['lon']:+.2f}_{loc['lat']:+.2f}_1980_2025.csv.gz"
    return gee.DIR_XAVIER / f"xavier_no_{loc['lon']:+.2f}_{loc['lat']:+.2f}_1980_2025.csv.gz"


def tarefas(n, lote=1):
    L = locais()
    lo, hi = np.array([p[1] for p in PARAM]), np.array([p[2] for p in PARAM])
    lista = []
    for j, loc in L.iterrows():
        # lotes 1 e 2: 47 locais originais; lotes 3 e 4: só os seus locais; lote ≥ 5: todos (ressimulação com R3/R5/R8)
        if lote < 5 and loc["lote_local"] != (1 if lote <= 2 else lote):
            continue
        sem = SEMENTE + j if lote == 1 else SEMENTE + 1000 * lote + j
        X = qmc.scale(qmc.LatinHypercube(d=len(PARAM), seed=sem).random(n), lo, hi)
        base = 0 if lote == 1 else lote * 100000
        for i, x in enumerate(X):
            lista.append(dict(id=base + j * n + i, amostra=i, lote=lote, idx_local=j, **loc.to_dict(), **dict(zip([p[0] for p in PARAM], x))))
    return L, lista


def _ini():
    global _D
    pasta = rb.DIR_RUN / f"fe{current_process()._identity[0]:02d}"
    shutil.rmtree(pasta, ignore_errors=True)
    _D = DSSAT(str(pasta))


def _estacao(t):
    k = (t["lon"], t["lat"])
    if k not in _WS:
        c = pd.read_csv(arq_clima(t), index_col=0, parse_dates=True).interpolate(limit=5).ffill().bfill()
        ws = bs.estacao(c, t["lat"], t["lon"], elev=t["elev"])
        texto = ws._write_wth()
        ws._write_wth = lambda texto=texto: texto
        _WS[k] = ws
    return _WS[k]


def r1ppo_padrao(gmr):
    """R1PPO dos ecótipos genéricos SB0501–SB0801 (0,414 + 0,045 por MG, até 0,549), usado nos lotes 1–5."""
    return 0.414 + 0.045 * (np.clip(np.round(gmr), 5, 8) - 5)


def cultivar(t):
    mg_int = int(np.clip(round(t["gmr"]), 5, 9))
    c = Soybean(f"99000{mg_int}")
    for k in ["csdl", "ppsen", "em-fl", "fl-sd", "sd-pm"]:
        c[k] = float(t[k])
    for k, v in bs.COEF_S2_COMUNS.items():
        c[k] = v
    for k, v in bs.ECO_S2.items():
        c["eco#"][k] = v
    if "r1ppo" in t and pd.notna(t["r1ppo"]):
        c["eco#"]["r1ppo"] = float(t["r1ppo"])
    return c


def rodar(t):
    pdate = date(2001, 9, 10) + timedelta(days=int(t["doy_sem"]))
    pdate = date(ANO_INI + (pdate.month < 7), pdate.month, pdate.day)
    sdate = pdate - timedelta(days=30)
    ctrl = SimulationControls(general=SCGeneral(sdate=sdate, nyers=NYERS),
                              options=SCOptions(water="N", nitro="N", symbi="N", co2="M"),
                              methods=SCMethods(evapo=bs.EVAPO), management=SCManagement(irrig="N", ferti="N", harvs="M"),
                              outputs=SCOutputs(grout="N", waout="N", caout="N"))
    campo = Field(id_field="BSXV0001", wsta=_estacao(t), id_soil=bs.perfil_solo("argiloarenoso", t["lat"], t["lon"]),
                  flob=0, fldt="DR000", fldd=0, flds=0, sltx="SC", sldp=200)
    plant = Planting(pdate=pdate, ppop=bs.POPULACAO, ppoe=bs.POPULACAO, plrs=bs.ESPACAMENTO, pldp=bs.PROF_SEMEADURA)
    try:
        _D.run_treatment(field=campo, cultivar=cultivar(t), planting=plant,
                         initial_conditions=bs.condicoes_iniciais("argiloarenoso", sdate), simulation_controls=ctrl, verbose=False)
        s = bs._ler_summary(_D.output_files["Summary"])
        ev = ler_evaluate(_D.output_files.get("Evaluate", ""))
    except Exception as erro:
        return [dict(id=t["id"], erro=str(erro)[:80])]
    pd_ = mo.bs_data(s.PDAT)
    out = pd.DataFrame({"id": t["id"], "safra": ANO_INI + np.arange(len(s)), "data_semeadura": pd_.dt.strftime("%Y-%m-%d"),
                        "e": (mo.bs_data(s.EDAT) - pd_).dt.days, "r1": (mo.bs_data(s.ADAT) - pd_).dt.days,
                        "r7": (mo.bs_data(s.MDAT) - pd_).dt.days})
    if len(ev) == len(out):          # a partir do lote 4: R3, R5 e R8 do próprio DSSAT (dias após a semeadura)
        out["r3"], out["r5"], out["r8_dssat"] = ev.PD1PS.values, ev.PDFPS.values, ev.R8APS.values
    return out.to_dict("records")


def ler_evaluate(txt):
    """Evaluate.OUT: uma linha por safra; PD1PS (R3), PDFPS (R5) e R8APS (R8) em dias após a semeadura (-99 → NaN)."""
    if not txt:
        return pd.DataFrame()
    linhas = txt.splitlines()
    i = next(k for k, l in enumerate(linhas) if l.startswith("@"))
    cab = linhas[i].replace("@", " ").split()
    dados = [l.split() for l in linhas[i + 1:] if l.strip()]
    df = pd.DataFrame([d for d in dados if len(d) == len(cab)], columns=cab)
    return df[["PD1PS", "PDFPS", "R8APS"]].apply(pd.to_numeric, errors="coerce").replace(-99, np.nan)


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else os.cpu_count()
    lote = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    DIR_TREINO.mkdir(parents=True, exist_ok=True)
    L, lista = tarefas(n, lote)
    tar = pd.DataFrame(lista)
    if lote == 1:
        L.to_csv(DIR / "locais.csv", index_label="idx_local")
    else:
        L0 = pd.read_csv(DIR / "locais.csv", index_col="idx_local")
        n0 = (L.lote_local == 1).sum()
        assert np.allclose(L0[["lon", "lat"]].values[:n0], L[["lon", "lat"]].values[:n0]), "locais diferem do lote 1"
        L.to_csv(DIR / "locais.csv", index_label="idx_local")
        antigas = pd.read_csv(DIR / "tarefas.csv")
        antigas = antigas[antigas.get("lote", pd.Series(1, index=antigas.index)).fillna(1) != lote]
        tar = pd.concat([antigas, tar], ignore_index=True)
    tar.to_csv(DIR / "tarefas.csv", index=False)
    feitos = set()
    for f in DIR_TREINO.glob("parte_*.csv"):
        feitos |= set(pd.read_csv(f, usecols=["id"]).id)
    pend = sorted([t for t in lista if t["id"] not in feitos], key=lambda t: (t["amostra"], t["idx_local"]))
    parte = len(list(DIR_TREINO.glob("parte_*.csv")))
    print(f"{len(L)} locais | {len(lista)} tarefas | {len(pend)} pendentes | {nproc} processos", flush=True)
    t0, buf = time.time(), []
    with get_context("fork").Pool(nproc, initializer=_ini) as pool:
        for k, out in enumerate(pool.imap_unordered(rodar, pend, chunksize=8), 1):
            buf += out
            if k % 1000 == 0 or k == len(pend):
                pd.DataFrame(buf).to_csv(DIR_TREINO / f"parte_{parte:04d}.csv", index=False)
                parte, buf = parte + 1, []
                taxa = k / (time.time() - t0)
                print(f"{k}/{len(pend)}  {time.time() - t0:.0f} s  faltam ~{(len(pend) - k) / taxa / 3600:.1f} h", flush=True)
