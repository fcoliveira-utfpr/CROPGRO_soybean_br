"""Treina a estimação neural da posterior (NPE) para a inversão do satélite e aplica aos ponto-safra observados.

Roda no ambiente .venv_sbi (torch 2.2.2 + sbi, numpy < 2); as simulações vêm de src/sbi_simulador.py.
Parâmetros θ = (g, s, a, b): GMR efetivo, semeadura (dias após 10/set) e as defasagens satélite↔estádio.
Condicionamento: x = (SOS, EOS50) + contexto (latitude, altitude, Tmed e Tmin mensais de set a abr).

Saídas em resultados/sbi/:
  npe_<nome>.pt                — posterior treinada;
  validacao_<nome>.csv         — por simulação de teste: verdadeiro, média e quantis da posterior, rank (SBC);
  gmr_efetivo_<nome>.csv       — por ponto-safra observado: média, desvio e quantis 5/50/95% de g e de s.

Uso: .venv_sbi/bin/python src/sbi_treinar.py [nome] [n_teste] [n_amostras]
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch.distributions import Normal, Uniform

from sbi.utils import MultipleIndependent

try:                                  # sbi ≥ 0.23
    from sbi.inference import NPE as Inferencia
except ImportError:                   # versões antigas
    from sbi.inference import SNPE as Inferencia

RAIZ = Path(__file__).resolve().parents[1]
DIR = RAIZ / "resultados" / "sbi"
torch.manual_seed(0)


def prior():
    return MultipleIndependent([
        Uniform(torch.tensor([5.0]), torch.tensor([9.5])),         # g
        Uniform(torch.tensor([0.0]), torch.tensor([171.0])),       # s
        Normal(torch.tensor([10.0]), torch.tensor([4.0])),         # a
        Normal(torch.tensor([0.0]), torch.tensor([6.0])),          # b
    ], validate_args=False)


def amostrar(post, xs, n):
    """Amostras da posterior para várias observações; usa sample_batched quando existe."""
    xs = torch.as_tensor(xs, dtype=torch.float32)
    if hasattr(post, "sample_batched"):
        return post.sample_batched((n,), x=xs, show_progress_bars=False).numpy()      # (n, lote, 4)
    return np.stack([post.sample((n,), x=x, show_progress_bars=False).numpy() for x in xs], axis=1)


if __name__ == "__main__":
    nome = sys.argv[1] if len(sys.argv) > 1 else "base"
    n_teste = int(sys.argv[2]) if len(sys.argv) > 2 else 500
    n_am = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
    sim = np.load(DIR / f"simulacoes_{nome}.npz", allow_pickle=True)
    theta = torch.as_tensor(sim["theta"], dtype=torch.float32)
    xc = torch.as_tensor(np.hstack([sim["x"], sim["contexto"]]), dtype=torch.float32)
    th_tr, x_tr, th_te, x_te = theta[n_teste:], xc[n_teste:], theta[:n_teste], xc[:n_teste]

    t0 = time.time()
    inf = Inferencia(prior=prior(), density_estimator="nsf")
    est = inf.append_simulations(th_tr, x_tr).train(show_train_summary=True)
    post = inf.build_posterior(est)
    t_treino = time.time() - t0
    torch.save(est.state_dict(), DIR / f"npe_{nome}.pt")
    print(f"treino: {t_treino:.0f} s com {len(th_tr)} simulações", flush=True)

    # validação em simulações não usadas no treino: erro, cobertura de 90% e rank (SBC)
    t0 = time.time()
    am = amostrar(post, x_te, n_am)
    t_inf = (time.time() - t0) / n_teste
    nomes = ["g", "s", "a", "b"]
    linhas = {}
    for j, p in enumerate(nomes):
        v, a = th_te[:, j].numpy(), am[:, :, j]
        linhas |= {f"{p}_true": v, f"{p}_media": a.mean(0), f"{p}_q05": np.quantile(a, 0.05, 0), f"{p}_q95": np.quantile(a, 0.95, 0),
                   f"{p}_rank": (a < v).sum(0)}
    val = pd.DataFrame(linhas)
    val.to_csv(DIR / f"validacao_{nome}.csv", index=False)
    for p in nomes:
        e = val[f"{p}_media"] - val[f"{p}_true"]
        cob = ((val[f"{p}_q05"] <= val[f"{p}_true"]) & (val[f"{p}_true"] <= val[f"{p}_q95"])).mean()
        print(f"{p}: RMSE {np.sqrt((e ** 2).mean()):.3f}  viés {e.mean():.3f}  cobertura 90% {cob:.3f}", flush=True)
    print(f"tempo por inferência: {1000 * t_inf:.1f} ms ({n_am} amostras)", flush=True)

    # aplicação aos ponto-safra observados
    obs = np.load(DIR / "observado.npz", allow_pickle=True)
    xo = np.hstack([obs["x"], obs["contexto"]])
    t0 = time.time()
    am = np.concatenate([amostrar(post, xo[i:i + 500], n_am) for i in range(0, len(xo), 500)], axis=1)
    print(f"{len(xo)} ponto-safra invertidos em {time.time() - t0:.0f} s", flush=True)
    res = pd.DataFrame({"lon": obs["lon"], "lat": obs["lat"], "safra": obs["safra"]})
    for j, p in enumerate(["g", "s"]):
        a = am[:, :, j]
        res[f"{p}_media"], res[f"{p}_dp"] = a.mean(0), a.std(0)
        for q in (5, 50, 95):
            res[f"{p}_q{q:02d}"] = np.percentile(a, q, axis=0)
    res.to_csv(DIR / f"gmr_efetivo_{nome}.csv", index=False)
    print(res[["g_media", "g_dp", "s_media"]].describe().round(2).to_string(), flush=True)
