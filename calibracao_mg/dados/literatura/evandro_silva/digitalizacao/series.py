"""Digitaliza os pontos observados (IAF, folha, haste, grãos, biomassa) das Figs. 3.6–3.13 da tese de E. H. F. M. Silva.

Eixos e ticks pela grade fixa (grade.py); valor de y = passo × nº do tick a partir de 0 (passos lidos nas figuras);
marcadores por cor (R-0N vermelho, R-1000N laranja, I-0N oliva, I-1000N verde), linhas removidas por abertura 5×5.
"""
import sys

import numpy as np
from PIL import Image
from scipy import ndimage

from grade import caixas, calibrar

EXP = {19: "PI-1", 20: "PI-2", 21: "PI-3", 22: "SO-1", 23: "SO-2", 24: "CA-1", 25: "TU-1", 26: "TE-1"}
PADRAO = {(0, 0): ("iaf", 1), (0, 1): ("graos_kgha", 1000), (1, 0): ("folha_kgha", 500),
          (1, 1): ("biomassa_kgha", 2000), (2, 0): ("haste_kgha", 1000)}
PASSOS = {22: {(0, 1): 2000, (2, 0): 500}, 23: {(0, 1): 2000}, 24: {(0, 0): 2, (0, 1): 2000},
          26: {(2, 0): 500}}
SEM = {25: [(2, 0)]}      # TU-1: sem painel de haste (fenologia no lugar)
# correções manuais (pixels lidos nas figuras): yb = zero do eixo y; py = px por passo; dx = px por 20 d
MANUAL = {(24, (0, 1)): dict(yb=470, py=111.3), (26, (0, 1)): dict(dx=89.9), (26, (1, 1)): dict(dx=89.9)}


def classes(a):
    r, g, b = a[..., 0], a[..., 1], a[..., 2]
    return {"R-0N": (r > 90) & (g < 70) & (b < 70),
            "R-1000N": (r > 200) & (g > 85) & (g < 175) & (b < 110),
            "I-0N": (r > 140) & (g > 140) & (b < 110) & (np.abs(r - g) < 45),
            "I-1000N": (g > 80) & (r < 90) & (b < 90) & (g - r > 40)}


def ticks_regulares(ty, yb):
    ty = sorted([t for t in ty if t <= yb + 4], reverse=True)
    if not ty or abs(ty[0] - yb) > 6:
        ty = [yb] + ty
    d0 = ty[0] - ty[1]
    reg = ty[:2]
    for t in ty[2:]:
        if abs((reg[-1] - t) - d0) <= 3:
            reg.append(t)
        else:
            break
    return reg


def digitalizar(a, k):
    out = []
    mascaras = classes(a)
    E = calibrar(a)
    for pos, cx in caixas(a).items():
        if pos not in PADRAO or pos in SEM.get(k, []):
            continue
        var, passo = PADRAO[pos]
        passo = PASSOS.get(k, {}).get(pos, passo)
        e = E[pos]
        m_ = MANUAL.get((k, pos), {})
        yb = m_.get("yb", e["yb"])
        py = m_.get("py", float(np.mean(-np.diff(e["ty"]))) if len(e["ty"]) > 1 else np.nan)   # px por passo
        dx = m_.get("dx", float(np.mean(np.diff(e["tx"]))) if len(e["tx"]) > 1 else np.nan)    # px por 20 d
        n_reg = len(e["ty"]) if "py" not in m_ else 4
        y_top = max(0, yb - (n_reg - 0.1) * py)              # acima disso há quebra de eixo ou tabela
        e = dict(e, yb=int(yb))
        x0, y0, x1, y1 = cx
        for trat, m in mascaras.items():
            sub = np.zeros_like(m)
            sub[int(y_top):e["yb"] + 8, e["xe"] + 2:x1] = m[int(y_top):e["yb"] + 8, e["xe"] + 2:x1]
            sub = ndimage.binary_closing(sub, structure=np.ones((3, 3)))      # tapa o contorno preto/JPEG
            miolo = ndimage.distance_transform_edt(sub) >= 2.0                  # linhas (2–3 px) não têm miolo
            lab, n = ndimage.label(miolo)
            for i in range(1, n + 1):
                yy, xx = np.where(lab == i)
                if len(yy) < 6:
                    continue
                cy, cxp = yy.mean(), xx.mean()
                out.append(dict(experimento=EXP[k], variavel=var, tratamento=trat,
                                passo=passo, px_passo=round(py, 1), px_20d=round(dx, 1),
                                das=round((cxp - e["xe"]) / dx * 20, 1),
                                valor=round((e["yb"] - cy) / py * passo, 3 if var == "iaf" else 0),
                                area_px=len(yy)))
    return out


if __name__ == "__main__":
    import pandas as pd
    linhas = []
    for k in EXP:
        a = np.asarray(Image.open(f"{sys.argv[1]}/t-{k:03d}.png").convert("RGB")).astype(int)
        linhas += digitalizar(a, k)
    d = pd.DataFrame(linhas).sort_values(["experimento", "variavel", "tratamento", "das"])
    # o mesmo marcador partido em dois (contorno/oclusão): mesmo tratamento a < 2,5 d -> fica o de maior área
    grupo = (d.groupby(["experimento", "variavel", "tratamento"]).das.diff().fillna(99) >= 2.5).cumsum()
    d = d.loc[d.groupby(grupo).area_px.idxmax()]
    d.to_csv(sys.argv[2], index=False)
    print(d.groupby(["experimento", "variavel", "tratamento"]).size().unstack().to_string())
    print("área dos marcadores:", d.area_px.describe().round(0).to_dict())
