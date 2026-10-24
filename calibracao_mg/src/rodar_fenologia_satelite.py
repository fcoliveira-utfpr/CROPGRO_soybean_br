"""Roda a fenologia MODIS nos 462 pontos de Battisti & Sentelhas (2019), safras 2000/01 a 2023/24.

Uso: .venv/bin/python src/rodar_fenologia_satelite.py
"""
import pandas as pd

import fenologia_satelite as fs
import gee

if __name__ == "__main__":
    gee.iniciar()
    pts = pd.read_csv(gee.RAIZ / "dados" / "processados" / "pontos_bs2019_digitalizados.csv")
    fr = pd.read_csv(gee.RAIZ / "dados" / "processados" / "fracoes_solo_mapbiomas_soja2016.csv")
    pts = pts.merge(fr[["lon", "lat"]], on=["lon", "lat"])        # só células com soja
    fs.rodar(pts, range(2001, 2025), n=150, threads=16)
