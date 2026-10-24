"""Baixa o clima do Xavier (BR-DWGD, 0,1°) nos 462 pontos digitalizados de Battisti & Sentelhas (2019).

Cada ponto é a média dos 4 pixels em volta do nó x,5°. Período 1980–2013 (33 safras).
Uso: .venv/bin/python src/baixar_xavier_bs2019.py
"""
import pandas as pd

import gee

if __name__ == "__main__":
    gee.iniciar()
    pontos = pd.read_csv(gee.RAIZ / "dados" / "processados" / "pontos_bs2019_digitalizados.csv")
    falhas = gee.baixar_series(pontos, 1980, 2013, no=True, n_threads=24)
    print("falhas:", falhas)
