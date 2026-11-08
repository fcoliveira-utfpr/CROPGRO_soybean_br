"""Critérios E2 e E3 da emenda 1 do protocolo (resultados/linha_base_v13/protocolo.md): RMSE por local (DSSAT real),
média entre locais com peso igual, para R1, R5, R7 (e R3 só em Silva 2021).

Locais: os dos 209 tratamentos (resultados/evandro/avaliacoes_experimentos_detalhe.csv, coluna metodo) + os 5 de
Silva (9 experimentos de avaliação; sem TE-100, TE-50, SO-2). Silva: v13 = fora da amostra (v13_fora_amostra.csv,
GMR de calibração); demais = "<rótulo> completa" em avaliacoes_cadeias.csv.

Saída: resultados/evandro/criterios_por_local.csv (por local) e criterios_E2_E3.csv (médias).
Uso: .venv/bin/python src/avaliar_por_local.py [rótulos...]   (padrão: todos os disponíveis)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

RAIZ = Path(__file__).resolve().parent.parent
EV = RAIZ / "resultados" / "evandro"
FORA = ["TE-100", "TE-50", "SO-2"]


def silva(rotulo):
    if rotulo == "v13":
        d = pd.read_csv(EV / "v13_fora_amostra.csv")
        d = d[d.cenario_gmr == "gmr_calibracao"]
    else:
        d = pd.read_csv(EV / "avaliacoes_cadeias.csv")
        d = d[d.conjunto == f"{rotulo} completa"]
    d = d[~d.experimento.isin(FORA)]
    if d.empty:
        return None
    d = d.rename(columns={f"{e}_erro": f"erro_{e}" for e in ("r1", "r3", "r5", "r7")})
    return d.assign(local=d.local + " (Silva)", grupo="Silva 2021")


def por_local(rotulo):
    a = pd.read_csv(EV / "avaliacoes_experimentos_detalhe.csv")
    a = a[a.metodo == rotulo].assign(grupo="experimentos atuais")
    s = silva(rotulo)
    if a.empty or s is None:
        return None
    d = pd.concat([a, s], ignore_index=True)
    L = []
    for (g, loc), x in d.groupby(["grupo", "local"]):
        lin = dict(versao=rotulo, grupo=g, local=loc)
        for e in ("r1", "r3", "r5", "r7"):
            v = x[f"erro_{e}"].dropna() if f"erro_{e}" in x else pd.Series(dtype=float)
            lin[f"n_{e}"], lin[f"rmse_{e}"] = len(v), (float(np.sqrt((v ** 2).mean())) if len(v) else np.nan)
        L.append(lin)
    return pd.DataFrame(L)


def medias(t):
    L = []
    for (v, nome), x in [((v, "todos os locais (E2)"), x) for v, x in t.groupby("versao")] + \
                       [((v, f"só {g} (E3)" if g == "experimentos atuais" else f"só {g}"), x)
                        for (v, g), x in t.groupby(["versao", "grupo"])]:
        lin = dict(versao=v, conjunto=nome)
        for e in ("r1", "r3", "r5", "r7"):
            r = x[f"rmse_{e}"].dropna()
            lin[f"locais_{e}"], lin[f"media_rmse_{e}"] = len(r), r.mean()
        L.append(lin)
    return pd.DataFrame(L)


if __name__ == "__main__":
    rot = sys.argv[1:] or list(pd.read_csv(EV / "avaliacoes_experimentos_detalhe.csv").metodo.unique())
    t = pd.concat([x for x in (por_local(r) for r in rot) if x is not None], ignore_index=True)
    t.to_csv(EV / "criterios_por_local.csv", index=False)
    m = medias(t)
    m.to_csv(EV / "criterios_E2_E3.csv", index=False)
    pd.set_option("display.width", 200)
    print(m.round(2).to_string(index=False))
