"""Detecta os painéis (eixos x e y) e os ticks nas Figs. 3.6–3.13 da tese de E. H. F. M. Silva."""
import numpy as np
from PIL import Image


def grupos(xs, gap=2):
    g = []
    for x in xs:
        if g and x - g[-1][-1] <= gap:
            g[-1].append(x)
        else:
            g.append([x])
    return [float(np.mean(k)) for k in g]


def carregar(arq):
    return np.asarray(Image.open(arq).convert("RGB")).astype(int)


def corridas(lin, gap=4, minimo=300):
    xs = np.where(lin)[0]
    if not len(xs):
        return []
    c = np.where(np.diff(xs) > gap)[0]
    ini, fim = np.r_[xs[0], xs[c + 1]], np.r_[xs[c], xs[-1]]
    return [(i0, i1) for i0, i1 in zip(ini, fim) if i1 - i0 > minimo]


def paineis(a, altura=470):
    """Painéis a partir das linhas de eixo x (contínuas); o eixo y sobe a partir da ponta esquerda."""
    H, W, _ = a.shape
    dark = a.sum(2) < 300
    cand = []
    for y in range(H):
        for x0, x1 in corridas(dark[y]):
            # é eixo x se tem ticks logo abaixo
            tx = grupos(np.where((a[y + 2:y + 9, x0 - 2:x1 + 3].sum(2) < 380).sum(0) >= 5)[0] + x0 - 2)
            if 5 <= len(tx) <= 9:
                cand.append((y, x0, x1, tx))
    res = []
    for y, x0, x1, tx in cand:
        yt = max(0, y - altura)
        cols = dark[yt:y, max(0, x0 - 40):x0 + 10].sum(0)
        xe = max(0, x0 - 40) + int(np.argmax(cols))          # eixo y = coluna mais escura perto da ponta
        if cols.max() < 300:
            continue                                          # borda de barra, não painel
        if any(abs(y - r["yb"]) < 6 and abs(xe - r["xe"]) < 6 for r in res):
            continue
        seg = (a[yt:y + 3, xe - 9:xe - 1].sum(2) < 380).sum(1) >= 5
        ty = grupos(np.where(seg)[0] + yt)
        res.append(dict(xe=xe, xr=x1, yb=y, ticks_x=tx, ticks_y=ty))
    return res


if __name__ == "__main__":
    import sys
    for k in range(19, 27):
        for p in paineis(carregar(f"{sys.argv[1]}/t-{k:03d}.png")):
            print(k, p["xe"], p["yb"], p["xr"], "tx", len(p["ticks_x"]), [round(t) for t in p["ticks_x"]][:2],
                  "ty", len(p["ticks_y"]), [round(t) for t in p["ticks_y"]])
