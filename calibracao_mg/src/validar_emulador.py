"""Valida o emulador deixando pontos inteiros de fora (GroupKFold, 5 grupos de 8 pontos) e treina os modelos finais.

Saídas: resultados/emulador/validacao_metricas.csv, validacao_predicoes.parquet, modelo_<saida>.joblib
Uso: .venv/bin/python src/validar_emulador.py
"""
import sys
import time
import warnings
from pathlib import Path

import joblib
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))
import emulador as em  # noqa: E402

DIR = em.te.DIR

if __name__ == "__main__":
    d = pd.read_parquet(DIR / "tabela_treino.parquet")
    metricas, preds = [], d[["id", "safra", "ponto", "zona_bs", "amostra"] + em.SAIDAS].copy()
    for s in em.SAIDAS:
        t = time.time()
        pred, m = em.validar(d, s)
        preds[f"{s}_pred"] = pred
        metricas.append(m)
        print(f"{s}: RMSE {m['rmse']:.1f} | viés {m['vies']:.1f} | R² {m['r2']:.3f} | {time.time() - t:.0f} s", flush=True)
    pd.DataFrame(metricas).to_csv(DIR / "validacao_metricas.csv", index=False)
    preds.to_parquet(DIR / "validacao_predicoes.parquet")
    for s in em.SAIDAS:
        joblib.dump(em.treinar(d, s), DIR / f"modelo_{s}.joblib")
        print(f"modelo final {s} salvo", flush=True)
