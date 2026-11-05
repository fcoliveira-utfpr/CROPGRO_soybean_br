"""Digitaliza as barras de fenologia (painel inferior direito) das Figs. 3.6–3.13 da tese de E. H. F. M. Silva."""
import sys
import numpy as np
from PIL import Image

EXP = {19: "PI-1", 20: "PI-2", 21: "PI-3", 22: "SO-1", 23: "SO-2", 24: "CA-1", 25: "TU-1", 26: "TE-1"}
LAV = np.array([214, 217, 234])


def grupos(xs, gap=2):
    g = []
    for x in xs:
        if g and x - g[-1][-1] <= gap:
            g[-1].append(x)
        else:
            g.append([x])
    return [float(np.mean(k)) for k in g]


def analisar(arq):
    a = np.asarray(Image.open(arq).convert("RGB")).astype(int)
    H, W, _ = a.shape
    dark = a.sum(2) < 300
    # painel de fenologia: onde estão as barras lavanda (planting → R1); eixo = 1ª linha escura longa abaixo delas
    lav_mask = (np.abs(a - LAV).sum(2) < 30) & (a[..., 2] - a[..., 0] > 12)
    lav_mask[lav_mask.sum(1) < 70] = False        # tira o quadradinho da legenda
    lav_mask[:int(0.4 * H)] = False               # painel fica na metade de baixo
    yy, xx = np.where(lav_mask)
    xa, xb = np.percentile(xx, 1), np.percentile(xx, 99)
    y0 = int(yy.max()) + 8                                    # pula a borda da última barra
    faixa = dark[y0:y0 + 150, int(xa):int(xb)].sum(1)
    y_eixo = y0 + int(np.argmax(faixa))                       # linha do eixo x
    marca = (a[y_eixo + 2:y_eixo + 9, :].sum(2) < 380).sum(0) >= 5          # ticks abaixo do eixo
    ticks = [t for t in grupos(np.where(marca)[0]) if xa - 30 < t]
    ticks = [t for i, t in enumerate(ticks) if i == 0 or t - ticks[i - 1] > 40]
    ticks = ticks[:next((i for i in range(1, len(ticks)) if abs((ticks[i] - ticks[i - 1]) - (ticks[1] - ticks[0])) > 5), len(ticks))]
    assert 6 <= len(ticks) <= 8, ticks     # passo de 20 d
    x0, x1 = ticks[0], ticks[-1]
    xmax = 20 * (len(ticks) - 1)
    dia = lambda x: (x - x0) / (x1 - x0) * xmax
    # barras: faixas lavanda na coluna logo após x=0
    col = a[:y_eixo, int(x0) + 12]
    lav = np.where(np.abs(col - LAV).sum(1) < 40)[0]
    linhas = []
    for y in lav:
        if linhas and y - linhas[-1][-1] <= 2:
            linhas[-1].append(y)
        else:
            linhas.append([y])
    linhas = [int(np.mean(l)) for l in linhas if len(l) > 15]
    out = []
    for y in linhas:
        s = a[y - 8:y + 9, int(x0) - 5:int(x1) + 5].mean(0).sum(1)
        bordas = [dia(x + int(x0) - 5) for x in grupos(np.where(s < 450)[0] + 0)]
        bordas = [b for b in bordas if 1 < b < xmax - 0.5]   # tira o quadro
        out.append(bordas)
    return out


if __name__ == "__main__":
    print("experimento,linha,r1,r3,r5,r7")
    for k, e in EXP.items():
        for i, b in enumerate(analisar(f"{sys.argv[1]}/t-{k:03d}.png")):
            print(e, "observado" if i == 0 else f"sim{i}", *[round(v, 1) for v in b], sep=",")
