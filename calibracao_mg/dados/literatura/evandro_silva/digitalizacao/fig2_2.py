"""Digitaliza a Fig. 2.2 da tese de E. H. F. M. Silva: grão (círculos) e biomassa total (triângulos), pretos, 4 painéis."""
import sys

import numpy as np
import pandas as pd
from PIL import Image
from scipy import ndimage


a = np.asarray(Image.open(sys.argv[1]).convert("RGB")).astype(int)
H, W, _ = a.shape
CAIXAS = {"PI1-100": (0, 0, W // 2 + 5, H // 2 + 10), "PI2-100": (W // 2 - 25, 0, W, H // 2 + 10),
          "TE-100": (0, H // 2 - 10, W // 2 + 5, H), "TE-50": (W // 2 - 25, H // 2 - 10, W, H)}
preto = a.sum(2) < 200
# calibração medida nos ticks (px): zero do eixo x, zero do eixo y, px por 20 d, px por 2000 kg/ha
CAL = {"PI1-100": (104, 432, 66.9, 80.0), "PI2-100": (528, 432, 67.5, 80.0),
       "TE-100": (102, 859, 66.9, 80.7), "TE-50": (530, 859, 67.0, 80.7)}
out = []
for exp, cx in CAIXAS.items():
    x0, y0, x1, y1 = cx
    xe, yb, dx, py = CAL[exp]
    e = dict(xe=xe, yb=yb)
    if exp == "TE-50":                            # legenda no canto superior direito
        x1 = xe + int(4.6 * dx)
    sub = np.zeros_like(preto)
    sub[y0:int(e["yb"]) - 2, e["xe"] + 3:x1] = preto[y0:int(e["yb"]) - 2, e["xe"] + 3:x1]
    miolo = ndimage.distance_transform_edt(sub) >= 2.5        # barras de erro e texto não têm miolo
    lab, n = ndimage.label(ndimage.binary_dilation(miolo, iterations=2))
    for i in range(1, n + 1):
        yy, xx = np.where((lab == i) & sub)
        if len(yy) < 25:
            continue
        h, w = np.ptp(yy) + 1, np.ptp(xx) + 1
        cheio = len(yy) / (h * w)
        forma = "graos" if len(yy) >= 60 else "biomassa"     # círculo (~70–90 px) maior que o triângulo (~30–45 px)
        cy = (yy.min() + yy.max()) / 2 if forma == "graos" else yy.min() + 2 * h / 3   # centroide do triângulo
        out.append(dict(experimento=exp, variavel=forma + "_kgha", das=round((xx.mean() - e["xe"]) / dx * 20, 1),
                        valor=round((e["yb"] - cy) / py * 2000), preenchimento=round(cheio, 2), area_px=len(yy)))
d = pd.DataFrame(out)
d = d[~((d.experimento == "TE-50") & (d.valor > 6000))]               # marcadores da legenda
# círculos de grão parcialmente cobertos (barra de erro/curva) que a área classificou como triângulo (conferido na figura)
for exp, das in (("PI1-100", 114.1), ("TE-100", 56.2)):
    k = (d.experimento == exp) & (d.das.sub(das).abs() < 0.2) & (d.variavel == "biomassa_kgha")
    k &= d.valor == d[k].valor.min()
    d.loc[k, "variavel"] = "graos_kgha"
d = d.sort_values(["experimento", "variavel", "das"])
d["fonte"] = "Silva (2021) tese ESALQ, Fig. 2.2 (digitalizada; média ± dp no original)"
d.to_csv(sys.argv[2], index=False)
print(d.groupby(["experimento", "variavel"]).size().unstack().to_string())
print(d.to_string(index=False))
