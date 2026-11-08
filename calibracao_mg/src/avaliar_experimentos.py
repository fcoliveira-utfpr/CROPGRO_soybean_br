"""Critérios C1 e C3 do protocolo (resultados/linha_base_v13/protocolo.md): DSSAT real nos 209 tratamentos
experimentais (os FileX/FileA do TSE), com as cultivares por classe de GMR da mediana de uma posterior.

Mesma simulação de tse/comparar_tse.py (linha "MCMC v13", experimento (DSSAT)), mas para qualquer cadeia
(19 parâmetros da v13 ou 20 da v14). Acrescenta ao fim de resultados/evandro/avaliacoes_experimentos.csv.

Uso: DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1 .venv/bin/python src/avaliar_experimentos.py <rótulo> <burn> <cadeias...>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parent
sys.path.insert(0, str(AQUI))
sys.path.insert(0, str(RAIZ / "tse"))
import comparar_tse as ct  # noqa: E402
import validar_evandro as ve  # noqa: E402

OUT = RAIZ / "resultados" / "evandro"


def avaliar(arquivos, burn, rotulo):
    t = pd.read_csv(RAIZ / "tse" / "tratamentos_tse.csv")
    th = ve.posterior(arquivos, burn)
    coefs = {cl: ct.coef_linha(th, cl) for cl in sorted(t.classe.unique())}
    for cl in coefs:                       # DSSAT_HOME curto do TSE (/tmp/tse<cod> → tse/raiz_<cod>/DSSAT48)
        cod = f"{int(round(cl * 10)):03d}"
        curto = Path(f"/tmp/tse{cod}")
        if not curto.exists():
            curto.unlink(missing_ok=True)
            curto.symlink_to(RAIZ / "tse" / f"raiz_{cod}" / "DSSAT48")
    sim = ct.simular(coefs, rotulo, t)
    L = []
    for e in ("r1", "r5", "r7"):
        x = sim[f"erro_{e}"].dropna()
        L.append(dict(conjunto=rotulo, estadio=e, n=len(x), rmse=float(np.sqrt(np.mean(x ** 2))), vies=float(x.mean()),
                      mae=float(x.abs().mean())))
    E = pd.DataFrame(L)
    OUT.mkdir(parents=True, exist_ok=True)
    for arq, d in ((OUT / "avaliacoes_experimentos.csv", E), (OUT / "avaliacoes_experimentos_detalhe.csv", sim)):
        if arq.exists():
            ant = pd.read_csv(arq)
            col = "conjunto" if "conjunto" in d else "metodo"
            d = pd.concat([ant[ant[col] != rotulo], d])
        d.to_csv(arq, index=False)
    return E


if __name__ == "__main__":
    a = sys.argv[1:]
    print(avaliar(a[2:], int(a[1]), a[0]).round(2).to_string(index=False))
