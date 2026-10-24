"""Estende o clima do Xavier (BR-DWGD) até 2025 nos 462 pontos de Battisti & Sentelhas (2019) e junta com 1980–2013.

Gera dados/brutos/xavier/xavier_no_<lon>_<lat>_1980_2025.csv.gz (safras 1980/81 a 2024/25).
Uso: .venv/bin/python src/baixar_xavier_extensao.py
"""
import pandas as pd

import gee

if __name__ == "__main__":
    gee.iniciar()
    pontos = pd.read_csv(gee.RAIZ / "dados" / "processados" / "pontos_bs2019_digitalizados.csv")
    falhas = gee.baixar_series(pontos, 2014, 2025, no=True, n_threads=24)
    print("falhas:", falhas, flush=True)
    for r in pontos.itertuples():
        base = gee.DIR_XAVIER / f"xavier_no_{r.lon:+.2f}_{r.lat:+.2f}"
        a, b = f"{base}_1980_2013.csv.gz", f"{base}_2014_2025.csv.gz"
        try:
            d = pd.concat([pd.read_csv(a, index_col=0, parse_dates=True), pd.read_csv(b, index_col=0, parse_dates=True)])
        except FileNotFoundError:
            continue
        d[~d.index.duplicated()].sort_index().to_csv(f"{base}_1980_2025.csv.gz")
    print("junção concluída", flush=True)
