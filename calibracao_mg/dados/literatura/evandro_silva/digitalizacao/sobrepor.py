"""Desenha os pontos digitalizados sobre a figura, para conferência visual."""
import sys

import pandas as pd
from PIL import Image, ImageDraw

sys.path.insert(0, sys.path[0])
import numpy as np
from grade import calibrar
from series import EXP, MANUAL, PADRAO

COR = {"R-0N": (0, 0, 255), "R-1000N": (255, 0, 255), "I-0N": (0, 200, 255), "I-1000N": (0, 0, 0)}
k = int(sys.argv[2])
d = pd.read_csv(sys.argv[1])
d = d[d.experimento == EXP[k]]
im = Image.open(f"{sys.argv[3]}/t-{k:03d}.png").convert("RGB")
a = np.asarray(im).astype(int)
E = calibrar(a)
dr = ImageDraw.Draw(im)
inv = {v[0]: p for p, v in PADRAO.items()}
for _, r in d.iterrows():
    p = inv[r.variavel]
    e = E[p]
    m = MANUAL.get((k, p), {})
    yb = m.get("yb", e["yb"])
    x = e["xe"] + r.das / 20 * r.px_20d
    y = yb - r.valor / r.passo * r.px_passo
    dr.ellipse([x - 9, y - 9, x + 9, y + 9], outline=COR[r.tratamento], width=3)
im.resize((im.width * 6 // 10, im.height * 6 // 10)).save(sys.argv[4])
