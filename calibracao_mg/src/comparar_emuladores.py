"""Compara o emulador de fenologia (gradient boosting) com alternativas, na mesma validação por local (22/10/2026).

Modelos: HistGradientBoosting (configuração do emulador), rede neural densa (MLP, 2 × 128, entradas padronizadas,
Adam, parada antecipada) e regressão ridge (referência linear). Mesma subamostra de safras e as mesmas 5 dobras com
locais inteiros fora do treino (GroupKFold), para R1 e R7. Mede também o tempo de treino e de predição.

Saída: resultados/emulador_fenologia/comparacao_emuladores.csv
Uso: .venv/bin/python src/comparar_emuladores.py [n_linhas]
"""
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold
from sklearn.neural_network import MLPRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))
import emulador_fenologia as ef  # noqa: E402
import treino_fenologia as tf  # noqa: E402

MODELOS = {
    "gradient_boosting": ef.modelo,
    "mlp_2x128": lambda: make_pipeline(StandardScaler(), MLPRegressor(hidden_layer_sizes=(128, 128), learning_rate_init=1e-3, batch_size=512,
                                                                      max_iter=60, early_stopping=True, n_iter_no_change=5,
                                                                      random_state=0)),
    "ridge": lambda: make_pipeline(StandardScaler(), Ridge(alpha=1.0)),
}

if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 400_000
    d = pd.read_parquet(ef.DIR / "tabela_treino.parquet")
    if d.r1ppo.isna().any():
        d["r1ppo"] = d.r1ppo.fillna(pd.Series(tf.r1ppo_padrao(d.gmr.values), index=d.index))
    d = ef.preparar(d.sample(n, random_state=0).reset_index(drop=True))
    X, g = d[ef.colunas()], d.idx_local
    linhas = []
    for saida in ("r1", "r7"):
        ok = d[saida].notna() & X.notna().all(axis=1)          # MLP e ridge não aceitam NaN: mesmas linhas para todos
        Xs, y, gs = X[ok].reset_index(drop=True), d.loc[ok, saida].reset_index(drop=True), g[ok].reset_index(drop=True)
        for nome, f in MODELOS.items():
            pred = pd.Series(np.nan, index=y.index)
            t_tr = t_pr = 0.0
            for tr, tt in GroupKFold(5).split(Xs, y, gs):
                t0 = time.time()
                m = f().fit(Xs.iloc[tr], y.iloc[tr])
                t_tr += time.time() - t0
                t0 = time.time()
                pred.iloc[tt] = m.predict(Xs.iloc[tt])
                t_pr += time.time() - t0
            e = pred - y
            linhas.append(dict(saida=saida, modelo=nome, n=len(y), rmse=float(np.sqrt((e ** 2).mean())), mae=float(e.abs().mean()),
                               vies=float(e.mean()), r2=float(1 - (e ** 2).sum() / ((y - y.mean()) ** 2).sum()),
                               treino_s_por_dobra=t_tr / 5, predicao_us_por_linha=t_pr / len(y) * 1e6))
            print(linhas[-1], flush=True)
            pd.DataFrame(linhas).to_csv(ef.DIR / "comparacao_emuladores.csv", index=False)
    print(pd.DataFrame(linhas).round(3).to_string(index=False))
