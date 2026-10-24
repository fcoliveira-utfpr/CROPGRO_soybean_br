"""Digitaliza os pontos de clima da Fig. 1 de Battisti & Sentelhas (2019, Field Crops Res. 240:95–105).

Os pontos formam uma grade regular de 1° com centros em x,5°. A imagem é extraída do PDF, a escala
pixel → grau é calibrada pelos ticks da moldura, e cada nó da grade é testado quanto à presença de
um ponto (núcleo escuro e "grosso"; linhas de divisa têm 1–2 px). Nós isolados (< 2 vizinhos) são
descartados como falsos positivos. Resultado: 462 pontos (o artigo informa 453).

Uso: .venv/bin/python src/digitalizar_pontos_bs2019.py
"""
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from scipy import ndimage as nd

RAIZ = Path(__file__).resolve().parents[1]
PDF = RAIZ / "dados" / "literatura" / "pdfs" / "1-s2.0-S0378429018315715-main.pdf"
SAIDA = RAIZ / "dados" / "processados" / "pontos_bs2019_digitalizados.csv"

# Ticks da moldura detectados na imagem da pág. 3 (1498 × 1451 px, 300 ppi)
TICK_X_PX = [110, 246, 382, 518, 656, 792, 928, 1064, 1202, 1338]   # 72°W … 36°W
TICK_Y_PX = [110, 246, 383, 520, 656, 792, 928, 1065, 1202, 1338]   # 4°N … 32°S
LIMIAR_ESCURO = 110


def extrair_imagem():
    tmp = Path(tempfile.mkdtemp())
    subprocess.run(["pdfimages", "-f", "3", "-l", "3", "-j", str(PDF), str(tmp / "p")], check=True)
    return np.asarray(Image.open(next(tmp.glob("p-000.*"))).convert("L")).astype(float)


def eh_ponto(g, x, y):
    xi, yi = int(round(x)), int(round(y))
    lab, _ = nd.label(g[yi - 6:yi + 7, xi - 6:xi + 7] < LIMIAR_ESCURO)
    k = lab[5:8, 5:8]
    k = k[k > 0]
    if not len(k):
        return False
    comp = lab == np.bincount(k).argmax()
    nucleo = (g[yi - 1:yi + 2, xi - 1:xi + 2] < 90).sum()
    return nucleo >= 4 and nd.binary_erosion(comp, np.ones((2, 2))).sum() >= 2


def digitalizar():
    g = extrair_imagem()
    px = np.polyfit(np.arange(-72, -35, 4), TICK_X_PX, 1)
    py = np.polyfit(np.arange(4, -33, -4), TICK_Y_PX, 1)
    achados = set()
    for lon in np.arange(-74, -34) + 0.5:
        for lat in np.arange(-34, 6) + 0.5:
            x, y = np.polyval(px, lon), np.polyval(py, lat)
            legenda = x < 520 and y > 1000
            rosa = x > 1200 and y < 280
            if 8 < x < 1405 and 40 < y < 1442 and not (legenda or rosa) and eh_ponto(g, x, y):
                achados.add((lon, lat))
    viz = lambda k: sum((k[0] + dx, k[1] + dy) in achados for dx in (-1, 0, 1) for dy in (-1, 0, 1)) - 1
    pontos = sorted(k for k in achados if viz(k) >= 2)
    return pd.DataFrame(pontos, columns=["lon", "lat"])


if __name__ == "__main__":
    df = digitalizar()
    df.to_csv(SAIDA, index=False)
    print(f"{len(df)} pontos → {SAIDA.relative_to(RAIZ)}")
