"""Produtividade potencial (Yp) com a fenologia de cada método de calibração: o teste a jusante do MLCSE aplicado ao TSE.

Mesma configuração de src/rodar_produtividade.py (crescimento da BRS 284 de Battisti & Sentelhas 2019, GMR efetivo e
semeadura do SBI por ponto × safra, solo de maior fração, sem estresse hídrico e de N), trocando só a fenologia:
  MCMC v13            θ(GMR) na mediana da posterior (retas contínuas);
  GLUE                melhor conjunto de 60.000 (retas contínuas);
  TSE grade/redução   cultivar da classe de GMR mais próxima (ótimo multi-tratamento + R1PPO do ecótipo na grade;
                      R1PPO da v13 na redução), como o usuário do TSE faria com uma cultivar por classe.
Saída: resultados/tse/produtividade/<metodo>.csv e resumo_produtividade.csv (Yp, índice de colheita, ciclo).
Uso: DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1 .venv/bin/python tse/produtividade_tse.py [n_processos] [n_safras_por_ponto]
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
AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parent
sys.path.insert(0, str(RAIZ / "src"))
sys.path.insert(0, str(AQUI))
import bs2019 as bs  # noqa: E402
import calibrar_fenologia as cf  # noqa: E402
import comparar_tse as ct  # noqa: E402
import rodar_bs2019 as rb  # noqa: E402
import rodar_produtividade as rp  # noqa: E402
import treino_fenologia as tf  # noqa: E402
from DSSATTools.filex import (Field, Planting, SCGeneral, SCManagement, SCMethods, SCOptions,  # noqa: E402
                              SCOutputs, SimulationControls)
from DSSATTools.run import DSSAT  # noqa: E402

OUT = RAIZ / "resultados" / "tse" / "produtividade"
METODOS = ["MCMC v13", "GLUE (melhor de 60.000)", "TSE grade", "TSE redução de faixa"]
_D = None
COEF = None      # {metodo: função gmr -> {csdl, ppsen, em-fl, fl-sd, sd-pm, r1ppo}} (herdado pelo fork)


def funcoes_coeficientes():
    med = ct.posterior()
    glue = np.load(RAIZ / "resultados" / "glue" / "glue_amostras_60000.npz")
    th_glue = med.copy()
    ix = {n: i for i, n in enumerate(cf.NOMES)}
    for n, v in zip(glue["nomes"], glue["theta"][np.nanargmax(glue["logpost"])]):
        th_glue[ix[str(n)]] = v
    tab = pd.read_csv(RAIZ / "resultados" / "tse" / "coeficientes_metodos.csv")

    def linha(th):
        return lambda g: {ct.CHAVE[K]: v for K, v in ct.coef_linha(th, g).items()}

    def classe(m):
        t = tab[tab.metodo == m].set_index("classe")
        return lambda g: {ct.CHAVE[K]: float(t.loc[min(t.index, key=lambda c: abs(c - g)), K]) for K in ct.CHAVE}

    return {"MCMC v13": linha(med), "GLUE (melhor de 60.000)": linha(th_glue),
            "TSE grade": classe("TSE grade"), "TSE redução de faixa": classe("TSE redução de faixa")}


def cultivar(metodo, gmr):
    c = tf.cultivar(dict(gmr=gmr, **COEF[metodo](gmr)))
    mg = min(bs.COEF_S2, key=lambda m: abs(m - gmr))
    c["fl-sh"], c["fl-lf"] = bs.COEF_S2[mg][3], bs.COEF_S2[mg][6]
    return c


def _ini():
    global _D
    pasta = rb.DIR_RUN / f"t{current_process()._identity[0]:02d}"
    shutil.rmtree(pasta, ignore_errors=True)
    _D = DSSAT(str(pasta))


def rodar_ponto(args):
    lon, lat, elev, linhas = args
    saida = OUT / "sim" / f"{lon:+.1f}_{lat:+.1f}.csv"
    if saida.exists():
        return "pulado"
    loc = dict(lon=lon, lat=lat, elev=elev, arq="no")
    clima = pd.read_csv(tf.arq_clima(loc), index_col=0, parse_dates=True).interpolate(limit=5).ffill().bfill()
    ws = bs.estacao(clima, lat, lon, elev=elev)
    out = []
    for r in linhas:
        pdate = date(int(r["safra"]), 9, 10) + timedelta(days=int(round(r["s_q50"])))
        sdate = pdate - timedelta(days=30)
        solo = max(rp.SOLOS, key=lambda k: r[f"f_{k}"])
        for m in METODOS:
            ctrl = SimulationControls(
                general=SCGeneral(sdate=sdate, nyers=1),
                options=SCOptions(water="N", nitro="N", symbi="N", co2="M"),
                methods=SCMethods(evapo=bs.EVAPO), management=SCManagement(irrig="N", ferti="N", harvs="M"),
                outputs=SCOutputs(grout="N", waout="N", caout="N"))
            campo = Field(id_field="BSXV0001", wsta=ws, id_soil=bs.perfil_solo(solo, lat, lon), flob=0,
                          fldt="DR000", fldd=0, flds=0, sltx=bs.TEXTURA_DSSAT[solo], sldp=200)
            plant = Planting(pdate=pdate, ppop=bs.POPULACAO, ppoe=bs.POPULACAO, plrs=bs.ESPACAMENTO,
                             pldp=bs.PROF_SEMEADURA)
            try:
                _D.run_treatment(field=campo, cultivar=cultivar(m, r["g_q50"]), planting=plant,
                                 initial_conditions=bs.condicoes_iniciais(solo, sdate),
                                 simulation_controls=ctrl, verbose=False)
                s = bs._ler_summary(_D.output_files["Summary"]).iloc[0]
                erro = ""
            except Exception as e:      # noqa: BLE001
                s, erro = pd.Series(dtype=float), str(e)[:60]
            out.append(dict(lon=lon, lat=lat, safra=int(r["safra"]), gmr=r["g_q50"], metodo=m,
                            data_semeadura=pdate.isoformat(), hwam=s.get("HWAM", np.nan), cwam=s.get("CWAM", np.nan),
                            adat=s.get("ADAT", np.nan), mdat=s.get("MDAT", np.nan), erro=erro))
    pd.DataFrame(out).to_csv(saida, index=False)
    return "ok"


def ciclo(d, col):
    """Dias da semeadura até a data YYYYDDD do DSSAT."""
    v = pd.to_numeric(d[col], errors="coerce")
    ok = v > 0
    dt = pd.Series(pd.NaT, index=d.index)
    dt[ok] = pd.to_datetime((v[ok] // 1000).astype(int).astype(str), format="%Y") + pd.to_timedelta(v[ok] % 1000 - 1, "D")
    return (dt - pd.to_datetime(d.data_semeadura)).dt.days


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    k_safras = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    (OUT / "sim").mkdir(parents=True, exist_ok=True)
    COEF = funcoes_coeficientes()
    for g in (5, 7, 9.5):
        print(g, {m: {k: round(v, 2) for k, v in f(g).items()} for m, f in COEF.items()}, flush=True)
    lista, _ = rp.tarefas("final")
    if k_safras:      # subamostra: k safras por ponto, espaçadas
        lista = [(lo, la, el, ls[:: max(1, len(ls) // k_safras)][:k_safras]) for lo, la, el, ls in lista]
    print(f"{len(lista)} pontos; {sum(len(x[3]) for x in lista)} ponto-safra × {len(METODOS)} métodos", flush=True)
    with get_context("fork").Pool(n, initializer=_ini) as pool:
        for k, r in enumerate(pool.imap_unordered(rodar_ponto, lista), 1):
            if k % 10 == 0:
                print(f"{k}/{len(lista)} pontos", flush=True)
    d = pd.concat([pd.read_csv(f) for f in sorted((OUT / "sim").glob("*.csv"))], ignore_index=True)
    d["ic"] = d.hwam / d.cwam
    d["dias_r1"], d["dias_mat"] = ciclo(d, "adat"), ciclo(d, "mdat")
    d.to_csv(OUT / "yp_metodos.csv", index=False)
    reg = pd.read_csv(RAIZ / "resultados" / "produtividade" / "lacuna_pontos.csv")[["lon", "lat", "regiao"]].drop_duplicates()
    d = d.merge(reg, on=["lon", "lat"], how="left")
    linhas = []
    for (m, rg), g in pd.concat([d, d.assign(regiao="Brasil")]).groupby(["metodo", "regiao"]):
        linhas.append(dict(metodo=m, regiao=rg, n=len(g), yp=g.hwam.mean(), yp_p10=g.hwam.quantile(0.1),
                           frac_yp_zero=(g.hwam.fillna(0) < 100).mean(), ic=g.ic.median(), dias_r1=g.dias_r1.mean(),
                           dias_mat=g.dias_mat.mean(), erros=(g.erro.fillna("") != "").sum()))
    r = pd.DataFrame(linhas)
    r.to_csv(OUT / "resumo_produtividade.csv", index=False)
    pd.set_option("display.width", 200)
    print(r.round(2).to_string(index=False))
