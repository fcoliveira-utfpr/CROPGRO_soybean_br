"""Correção empírica do R7/R8 pela data de semeadura (decisão de 29/10/2026).

Nas semeaduras de set–nov o CROPGRO calibrado capta só ~50–60% do encurtamento observado do ciclo (etapa0_fenologia).
Correção: R7/R8 corrigido = simulado + k · max(0, d_ref − d), com d o dia da semeadura e d_ref = 10/nov (as médias
multilocais de Alliprandini, semeadas em 7/nov, são reproduzidas sem viés). O intercepto de cada local não entra:
só a parte do erro que depende da época. k é a inclinação dentro de cada local (efeito fixo de local × safra), ajustada
nos registros de R8 de set–nov; o intervalo vem de bootstrap por local.

Só os indicadores de calendário usam a correção (janela da safrinha, chuva na colheita, geada tardia); a
produtividade não é corrigida. Uso: ../calibracao_mg/.venv/bin/python src/fenologia_correcao.py (ajusta e grava).
"""
import json

import numpy as np
import pandas as pd

import config as cfg

ARQ = cfg.DIR_RES / "fenologia_correcao.json"
D_REF = pd.Timestamp("2001-11-10").dayofyear


def doy(mes_dia):
    return pd.to_datetime("2001-" + pd.Series(mes_dia), format="%Y-%m-%d").dt.dayofyear.to_numpy()


def ajustar(d):
    """Inclinação comum do erro do R8 contra a antecipação da semeadura, com efeito fixo de local × safra."""
    x = np.maximum(0, D_REF - d.doy.to_numpy()).astype(float)
    y = d.erro.to_numpy().astype(float)
    g = d.grupo.to_numpy()
    xc = x - pd.Series(x).groupby(g).transform("mean").to_numpy()
    yc = y - pd.Series(y).groupby(g).transform("mean").to_numpy()
    return -np.sum(xc * yc) / np.sum(xc ** 2)          # erro = sim − obs fica mais negativo quanto mais cedo


def carregar():
    return json.loads(ARQ.read_text())


def corrigir(data_estadio, data_semeadura, k=None):
    """Datas (Timestamp ou Series) de R7/R8 corrigidas; k padrão = estimativa central."""
    k = carregar()["k"] if k is None else k
    d = pd.to_datetime(data_semeadura)
    antecip = np.maximum(0, D_REF - getattr(d, "dt", d).dayofyear)
    return pd.to_datetime(data_estadio) + pd.to_timedelta(np.round(k * antecip), unit="D")


if __name__ == "__main__":
    d = pd.read_csv(cfg.DIR_RES / "etapa0_fenologia.csv")
    d = d[(d.estadio == "r8") & d.mes.isin([9, 10, 11])].copy()
    d["doy"] = doy(d.semeadura)
    d["grupo"] = d.local + d.safra.astype(str) + d.gmr.astype(str)
    d = d[d.groupby("grupo").doy.transform("nunique") > 1]        # só quem tem mais de uma época
    k = ajustar(d)
    rng = np.random.default_rng(1)
    locais = d.local.unique()
    boot = []
    for _ in range(2000):
        amostra = rng.choice(locais, len(locais))
        boot.append(ajustar(pd.concat([d[d.local == l].assign(grupo=d.grupo + str(i)) for i, l in enumerate(amostra)])))
    por_local = {l: round(float(ajustar(d[d.local == l])), 3) for l in locais}
    out = dict(k=round(float(k), 3), ic90=[round(float(v), 3) for v in np.quantile(boot, [0.05, 0.95])],
               d_ref="10/nov", n=int(len(d)), por_local=por_local)
    ARQ.write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print(json.dumps(out, indent=1, ensure_ascii=False))
    for nome, dia in [("1/set", "09-01"), ("21/set", "09-21"), ("11/out", "10-11"), ("1/nov", "11-01")]:
        a = max(0, D_REF - doy([dia])[0])
        print(f"{nome}: +{k * a:.1f} d (IC90 {out['ic90'][0] * a:.1f}–{out['ic90'][1] * a:.1f})")
