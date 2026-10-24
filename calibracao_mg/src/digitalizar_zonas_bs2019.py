"""Atribui a cada ponto digitalizado a zona ambiental da Fig. 3 de Battisti & Sentelhas (2019).

A cor do mapa em volta de cada nó (janela de 19 × 19 px, só pixels coloridos) vota na cor mais
próxima da legenda. Pontos fora da área colorida ficam sem zona (NaN).

Uso: .venv/bin/python src/digitalizar_zonas_bs2019.py
"""
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

RAIZ = Path(__file__).resolve().parents[1]
PDF = RAIZ / "dados" / "literatura" / "pdfs" / "1-s2.0-S0378429018315715-main.pdf"
PONTOS = RAIZ / "dados" / "processados" / "pontos_bs2019_digitalizados.csv"

# Imagem da pág. 5 (1798 × 1743 px): ticks da moldura e cores das caixas da legenda (zonas 1–10)
TICK_X_PX = [120, 281, 442, 604, 764, 926, 1087, 1248, 1409, 1570]   # 72°W … 36°W
TICK_Y_PX = [156, 317, 479, 640, 802, 964, 1126, 1287, 1449, 1611]   # 4°N … 32°S
CORES = np.array([[53, 100, 22], [69, 129, 29], [89, 172, 32], [127, 246, 40], [167, 246, 43],
                  [214, 253, 146], [230, 247, 47], [227, 205, 42], [214, 155, 39], [204, 101, 34]], float)


def extrair_imagem():
    tmp = Path(tempfile.mkdtemp())
    subprocess.run(["pdfimages", "-f", "5", "-l", "5", "-j", str(PDF), str(tmp / "p")], check=True)
    return np.asarray(Image.open(next(tmp.glob("p-000.*"))).convert("RGB")).astype(float)


def zona(im, x, y):
    w = im[y - 9:y + 10, x - 9:x + 10].reshape(-1, 3)
    w = w[w.std(1) > 20]                                    # só pixels coloridos
    if len(w) < 20:
        return np.nan
    d = np.linalg.norm(w[:, None, :] - CORES[None], axis=2)
    votos = np.bincount(d.argmin(1)[d.min(1) < 40], minlength=10)
    return int(votos.argmax()) + 1 if votos.max() >= 15 else np.nan


if __name__ == "__main__":
    im = extrair_imagem()
    px = np.polyfit(np.arange(-72, -35, 4), TICK_X_PX, 1)
    py = np.polyfit(np.arange(4, -33, -4), TICK_Y_PX, 1)
    p = pd.read_csv(PONTOS)
    p["zona_bs"] = [zona(im, int(round(np.polyval(px, r.lon))), int(round(np.polyval(py, r.lat)))) for r in p.itertuples()]
    p.to_csv(PONTOS, index=False)
    print(p.zona_bs.value_counts(dropna=False).sort_index().to_dict())
