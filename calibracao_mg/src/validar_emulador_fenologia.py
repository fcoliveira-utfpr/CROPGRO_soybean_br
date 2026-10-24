"""Monta a tabela, valida (locais fora do treino) e treina o emulador de fenologia; salva os modelos.

Uso: .venv/bin/python src/validar_emulador_fenologia.py
"""
import sys
import time
import warnings
from pathlib import Path

import joblib

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))
import emulador_fenologia as ef  # noqa: E402

if __name__ == "__main__":
    t = time.time()
    d = ef.montar_tabela()
    d.to_parquet(ef.DIR / "tabela_treino.parquet")
    print(f"tabela {d.shape} em {time.time() - t:.0f} s | locais {d.idx_local.nunique()}", flush=True)
    met, pred = ef.validar(d)
    met.to_csv(ef.DIR / "validacao_metricas.csv", index=False)
    d[["id", "idx_local", "safra"] + ef.SAIDAS].join(pred.add_suffix("_pred")).to_parquet(ef.DIR / "validacao_predicoes.parquet")
    print(met.round(2).to_string(index=False), f"\n{time.time() - t:.0f} s", flush=True)
    for s, m in ef.treinar(d).items():
        joblib.dump(m, ef.DIR / f"modelo_{s}.joblib")
    print("modelos salvos", flush=True)
