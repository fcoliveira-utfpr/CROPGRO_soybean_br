"""Compara cadeias independentes da calibração de fenologia: R-hat (Gelman-Rubin, cada conjunto de caminhantes
tratado como uma cadeia), medianas por cadeia, coeficientes θ(GMR) e resíduos por fonte e estádio na mediana conjunta.

Uso: .venv/bin/python src/comparar_cadeias.py cadeia1.h5 cadeia2.h5 [...] [--descarte 0.5]
"""
import argparse
import sys
import warnings
from pathlib import Path

import emcee
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))
import calibrar_fenologia as cf  # noqa: E402


def rhat(amostras):
    """amostras: lista de vetores 1-D (uma cadeia cada), mesmo tamanho."""
    n = min(len(a) for a in amostras)
    x = np.stack([a[-n:] for a in amostras])
    W = x.var(axis=1, ddof=1).mean()
    B = n * x.mean(axis=1).var(ddof=1)
    return float(np.sqrt(((n - 1) / n * W + B / n) / W))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cadeias", nargs="+", type=Path)
    ap.add_argument("--descarte", type=float, default=0.5)
    args = ap.parse_args()
    ch = []
    for f in args.cadeias:
        c = emcee.backends.HDFBackend(str(f), read_only=True).get_chain()
        ch.append(c[int(c.shape[0] * args.descarte):])
        print(f"{f.name}: {c.shape[0]} passos")
    linhas = []
    npar = min(c.shape[2] for c in ch)          # cadeias antigas podem ter menos parâmetros (os novos ficam no fim)
    for i, (nome, lo, hi, _) in enumerate(cf.PARAMS[:npar]):
        meds = [np.median(c[:, :, i]) for c in ch]
        perto = np.mean([np.mean((c[:, :, i] - lo < 0.02 * (hi - lo)) | (hi - c[:, :, i] < 0.02 * (hi - lo))) for c in ch])
        linhas.append(dict(parametro=nome, rhat=rhat([c[:, :, i].ravel() for c in ch]), **{f"med_{k}": m for k, m in enumerate(meds)},
                           faixa=f"[{lo:.3g}, {hi:.3g}]", frac_limite=perto))
    tab = pd.DataFrame(linhas)
    pd.set_option("display.width", 200)
    print(tab.round(3).to_string(index=False))
    th = np.median(np.concatenate([c[:, :, :npar].reshape(-1, npar) for c in ch]), axis=0)
    th = np.concatenate([th, [p[3] for p in cf.PARAMS[npar:]]])
    p = dict(zip(cf.NOMES, th))
    print("\nθ(GMR) na mediana conjunta:")
    print(pd.DataFrame({g: {k: float(np.clip(p[f"a_{k}"] + p[f"b_{k}"] * (g - 7), *cf.FAIXAS[k])) for k in cf.FAIXAS}
                        for g in (5, 6, 7, 8, 9)}).round(3).to_string())
    c = cf.Calibracao()
    pr = c.prever(th)
    o = c.obs.copy()
    o["res"] = o.dias_ve - np.array([pr[e][i] for i, e in enumerate(o.estadio)])
    ex = o[o.tipo == "experimento"].assign(fonte=lambda x: x.fonte.str.slice(0, 22) + " " + x.local.str.slice(0, 14) + " " + x.data_semeadura)
    print("\nResíduos dos experimentos (obs − sim, dias):")
    print(ex.groupby(["fonte", "estadio"]).res.agg(["count", "mean", "std"]).round(1).to_string())
    for tipo in ("catalogo", "agregado"):
        cat = o[o.tipo == tipo].assign(pred=lambda x: x.dias_ve - x.res)
        g = cat.groupby("id_registro").agg(y=("dias_ve", "first"), p=("pred", "mean"), gmr=("gmr", "first"), est=("estadio", "first"))
        g["res"] = g.y - g.p
        for est, gg in g.groupby("est"):
            print(f"\n{tipo} {est}: n={len(gg)}  viés={gg.res.mean():.1f}  RMSE={np.sqrt((gg.res ** 2).mean()):.1f}  "
                  f"r={np.corrcoef(gg.y, gg.p)[0, 1]:.2f}")
            print(gg.groupby(pd.cut(gg.gmr, [0, 6, 7, 8, 9, 12])).res.agg(["count", "mean", "std"]).round(1).to_string())
    print(f"\nescalas de erro (dias): σ_R1={np.exp(p['ln_s_r1']):.1f}  σ_R7/R8={np.exp(p['ln_s_r7']):.1f}  "
          f"σ_R5={np.exp(p.get('ln_s_r5', np.nan)):.1f}  σ_cat={np.exp(p['ln_s_cat']):.1f}  σ_agr={np.exp(p['ln_s_agr']):.1f}  δ8={p['delta8']:.1f}")
