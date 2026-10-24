"""Resume a cadeia MCMC da calibração de fenologia.

Uso: .venv/bin/python src/analisar_calibracao_fenologia.py
Gera `resumo_posterior.csv` e `diagnostico_cadeia.txt` em resultados/emulador_fenologia.
"""
import argparse
import sys
from pathlib import Path

import emcee
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import calibrar_fenologia as cf  # noqa: E402


def analisar():
    parser = argparse.ArgumentParser(description="Resume uma cadeia MCMC da calibração de fenologia.")
    parser.add_argument("--chain", type=Path, default=cf.DIR / "cadeia_fenologia.h5", help="arquivo HDF5 da cadeia")
    args = parser.parse_args()
    arq = args.chain.expanduser().resolve()
    backend = emcee.backends.HDFBackend(str(arq), read_only=True)
    cadeia = backend.get_chain()
    if cadeia.ndim != 3 or cadeia.shape[0] < 2:
        raise ValueError(f"Cadeia vazia ou inválida: {cadeia.shape}")

    # Descarta metade inicial para o resumo exploratório. A convergência ainda
    # precisa ser julgada pela autocorrelação e por cadeias independentes.
    pos_burn = cadeia[cadeia.shape[0] // 2:]
    amostras = pos_burn.reshape(-1, pos_burn.shape[-1])
    q = np.quantile(amostras, [0.05, 0.5, 0.95], axis=0)
    resumo = pd.DataFrame({
        "parametro": cf.NOMES,
        "media": amostras.mean(axis=0),
        "desvio_padrao": amostras.std(axis=0, ddof=1),
        "q05": q[0], "mediana": q[1], "q95": q[2],
    })
    limites = {p[0]: (p[1], p[2]) for p in cf.PARAMS}
    resumo["fracao_perto_limite"] = [
        np.mean((amostras[:, i] - limites[n][0] < 0.01 * (limites[n][1] - limites[n][0])) |
                (limites[n][1] - amostras[:, i] < 0.01 * (limites[n][1] - limites[n][0])))
        for i, n in enumerate(cf.NOMES)
    ]
    sufixo = arq.stem.removeprefix("cadeia_fenologia")
    if not arq.stem.startswith("cadeia_fenologia"):
        sufixo = "_" + arq.stem

    try:
        tau = backend.get_autocorr_time(tol=0)
        tau_str = ", ".join(f"{n}={v:.1f}" for n, v in zip(cf.NOMES, tau))
        tau_min, tau_max = float(tau.min()), float(tau.max())
    except Exception as exc:  # cadeia curta para estimar a autocorrelação
        tau = np.full(len(cf.NOMES), np.nan)
        tau_str = f"indisponível ({type(exc).__name__}: {exc})"
        tau_min = tau_max = float("nan")
    ace = np.asarray(backend.accepted)
    ace_media = float(ace.mean() / backend.iteration)
    resumo["tau_integrado"] = tau if np.isfinite(tau).all() else np.nan
    resumo["ess_estimado"] = pos_burn.shape[1] * pos_burn.shape[0] / resumo.tau_integrado
    resumo.to_csv(cf.DIR / f"resumo_posterior{sufixo}.csv", index=False)
    texto = (
        f"Arquivo: {arq.name}\nPassos: {backend.iteration}\nCaminhantes: {cadeia.shape[1]}\n"
        f"Parâmetros: {cadeia.shape[2]}\nAceitação média: {ace_media:.4f}\n"
        f"Autocorrelação por parâmetro: {tau_str}\n"
        f"Razão passos/tau: {backend.iteration / tau_max:.1f} (pior parâmetro)\n"
        f"ESS estimado mínimo entre parâmetros: {resumo.ess_estimado.min():.0f}\n"
        f"Amostras após descarte exploratório de 50%: {len(amostras)}\n"
        "Status: diagnóstico exploratório; não comprova convergência. Rode cadeias independentes "
        "e avalie R-hat/ESS antes de usar as estimativas como resultado final.\n"
    )
    (cf.DIR / f"diagnostico_cadeia{sufixo}.txt").write_text(texto)
    print(texto)
    print(resumo.round(3).to_string(index=False))
    return resumo, texto


if __name__ == "__main__":
    analisar()
