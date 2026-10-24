"""Valida o emulador de produtividade em duas etapas: fenologia prevista (fora do treino) → clima e balanço
hídrico por fase (S–R1, R1–R5, R5–R7) → produtividade. Mesma validação por pontos (GroupKFold, 5 grupos).

Uso: .venv/bin/python src/validar_emulador_fen.py
"""
import sys
import time
import warnings
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import GroupKFold

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))
import emulador as em  # noqa: E402

DIR = em.te.DIR


def modelo():
    return HistGradientBoostingRegressor(max_iter=800, learning_rate=0.06, max_leaf_nodes=63,
                                         l2_regularization=1.0, random_state=0)


if __name__ == "__main__":
    t = time.time()
    d = pd.read_parquet(DIR / "tabela_treino.parquet")
    d = d.join(pd.read_parquet(DIR / "validacao_predicoes.parquet")[["r1_pred", "r7_pred"]])
    d = d.dropna(subset=["r1_pred", "r7_pred"])
    d = em.variaveis_fenologicas(d)
    d.to_parquet(DIR / "tabela_treino_fen.parquet")
    print(f"variáveis por fase: {d.shape} em {time.time() - t:.0f} s", flush=True)
    X, y, g = d[em.colunas_entrada_fen(d)], d.hwam, d.ponto
    pred = pd.Series(np.nan, index=y.index)
    for k, (tr, tt) in enumerate(GroupKFold(5).split(X, y, g)):
        m = modelo().fit(X.iloc[tr], y.iloc[tr])
        pred.iloc[tt] = m.predict(X.iloc[tt])
        print(f"fold {k}: {time.time() - t:.0f} s", flush=True)
    e = pred - y
    met = dict(saida="hwam_fen", n=len(y), rmse=float(np.sqrt((e ** 2).mean())), vies=float(e.mean()),
               r2=float(1 - (e ** 2).sum() / ((y - y.mean()) ** 2).sum()))
    print(met, flush=True)
    pd.DataFrame([met]).to_csv(DIR / "validacao_metricas_fen.csv", index=False)
    d[["id", "safra", "ponto", "zona_bs", "amostra", "hwam"]].assign(hwam_pred=pred).to_parquet(DIR / "validacao_predicoes_fen.parquet")
    joblib.dump(modelo().fit(X, y), DIR / "modelo_hwam_fen.joblib")
    print("modelo final salvo", flush=True)
