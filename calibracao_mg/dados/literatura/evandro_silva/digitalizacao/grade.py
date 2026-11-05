"""Eixos e ticks de cada painel numa grade fixa 3 × 2 (Figs. 3.6–3.13 da tese de E. H. F. M. Silva)."""
import numpy as np

from paineis import grupos


def caixas(a):
    H, W, _ = a.shape
    linhas = [(0.0, 0.31), (0.31, 0.605), (0.605, 0.905)]
    cols = [(0.05, 0.50), (0.53, 0.99)]
    return {(i, j): (int(c0 * W), int(l0 * H), int(c1 * W), int(l1 * H))
            for i, (l0, l1) in enumerate(linhas) for j, (c0, c1) in enumerate(cols)}


def maior_corrida(v):
    best = cur = 0
    for x in v:
        cur = cur + 1 if x else 0
        best = max(best, cur)
    return best


def regulares(ts, tol=3):
    """Ticks igualmente espaçados a partir do primeiro da lista."""
    if len(ts) < 3:
        return ts
    d0 = ts[1] - ts[0]
    reg = ts[:2]
    for t in ts[2:]:
        if abs((t - reg[-1]) - d0) <= tol:
            reg.append(t)
        else:
            break
    return reg


def eixos(a, cx, xe=None):
    x0, y0, x1, y1 = cx
    dark = a[y0:y1, x0:x1].sum(2) < 300
    claro = a[y0:y1, x0:x1].sum(2) < 450
    corr = [maior_corrida(claro[:, j]) for j in range(35, min(160, dark.shape[1]))]
    if xe is None:
        xe = x0 + 35 + int(np.argmax(corr))                          # eixo y: maior linha vertical contínua
        if max(corr) < 250:
            return dict(xe=None, yb=None, ty=[], tx=[])
    sub = a[y0:y1, xe - 9:xe - 1].sum(2) < 380
    ty = [t + y0 for t in grupos(np.where(sub.sum(1) >= 5)[0])]
    ts = sorted(ty, reverse=True)                                    # de baixo para cima
    reg = max((regulares(ts[i:]) for i in range(min(3, len(ts)))), key=len, default=[])   # ignora tick espúrio embaixo
    yb = reg[0] if reg else None
    tx = []
    if yb is not None:
        sub = a[int(yb) + 2:int(yb) + 9, xe - 3:x1].sum(2) < 380
        tx = [t + xe - 3 for t in grupos(np.where(sub.sum(0) >= 5)[0])]
        tx = [t for t in tx if t > xe + 20]
        tx = regulares([float(xe)] + tx, tol=4) if tx else []
    return dict(xe=xe, yb=yb, ty=reg, tx=tx)


def calibrar(a):
    """Eixos de todos os painéis; quem não tem ticks x regulares herda do painel da mesma coluna com o mesmo eixo y."""
    C = caixas(a)
    E = {p: eixos(a, cx) for p, cx in C.items()}
    for p, e in E.items():          # eixo y apagado demais: herda a coluna do vizinho
        if e["xe"] is None:
            viz = [f["xe"] for q, f in E.items() if q[1] == p[1] and f["xe"] is not None]
            if viz:
                E[p] = eixos(a, C[p], xe=viz[0])
    for p, e in E.items():
        if len(e["tx"]) < 5:
            viz = [q for q, f in E.items() if q[1] == p[1] and q != p and len(f["tx"]) >= 5 and f["xe"] and abs(f["xe"] - e["xe"]) <= 4]
            if viz:
                e["tx"] = E[viz[0]]["tx"]
                e["tx_herdado"] = viz[0]
    return E


if __name__ == "__main__":
    import sys
    from PIL import Image
    for k in range(19, 27):
        a = np.asarray(Image.open(f"{sys.argv[1]}/t-{k:03d}.png").convert("RGB")).astype(int)
        for p, e in calibrar(a).items():
            if p == (2, 1):
                continue
            print(k, p, e["xe"], e["yb"], "tx", len(e["tx"]), round(float(np.mean(np.diff(e["tx"]))), 1) if len(e["tx"]) > 1 else None,
                  e.get("tx_herdado", ""), "ty", len(e["ty"]), round(float(-np.mean(np.diff(e["ty"]))), 1) if len(e["ty"]) > 1 else None)
