"""Validação fora da amostra: fenologia da v13 (DSSAT real) nos experimentos de Evandro H. F. M. da Silva.

Mesma configuração do treino do emulador e da comparação com o TSE (água e N desligados, solo argiloarenoso de B&S,
clima do Xavier no pixel do experimento, cultivar 99000<MG> + θ(GMR) da mediana da posterior da v13 + R1PPO(GMR)).
A fenologia observada é comum aos tratamentos: uma linha por experimento (dias após a semeadura até R1, R3, R5, R7).

Cenários de GMR: 'calibracao' (comercial conferido; TMG 7063 em Sorriso = 7,0), 'comercial' (TMG 7063 = 6,3) e
'tese' (grupo do DSSAT usado na tese).

Saída: resultados/evandro/v13_fora_amostra.csv (por experimento × cenário) e v13_fora_amostra_erros.csv.
Uso: DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1 .venv/bin/python src/validar_evandro.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parent
sys.path.insert(0, str(AQUI))
sys.path.insert(0, str(RAIZ / "tse"))
import comparar_tse as ct  # noqa: E402
import treino_fenologia as tf  # noqa: E402
from DSSATTools.run import DSSAT  # noqa: E402

EV = RAIZ / "dados" / "literatura" / "evandro_silva"
OUT = RAIZ / "resultados" / "evandro"
# Células do Xavier (resultados/evandro/celulas_xavier.csv). Piracicaba: as duas coordenadas da tese são a área
# experimental da ESALQ (22°42'S 47°38'W); a do cap. 3 (47°30') está arredondada → mesma célula.
CELULA = {"Piracicaba-SP": (-47.65, -22.75, 546), "Teresina-PI": (-42.75, -5.05, 78), "Sorriso-MT": (-55.75, -12.75, 375),
          "Cruz Alta-RS": (-53.55, -28.65, 476), "Tupanciretã-RS": (-53.75, -28.85, 466)}
ESTADIOS = ["r1", "r3", "r5", "r7"]


def observacoes():
    exp = pd.read_csv(EV / "experimentos_tese.csv")
    f2 = pd.read_csv(EV / "fenologia_tese_cap2.csv")
    f3 = pd.read_csv(EV / "fenologia_tese_cap3.csv")
    return exp, f2, f3


def tabela_obs():
    """Uma linha por experimento com R1/R3/R5/R7 observados (dias após a semeadura)."""
    exp, f2, f3 = observacoes()
    f3 = f3[f3.linha == "observado"]
    obs = pd.concat([f2, f3])[["experimento"] + [f"das_{s}" for s in ESTADIOS]]
    obs = obs.rename(columns={f"das_{s}": s for s in ESTADIOS}).drop_duplicates("experimento")
    return exp.merge(obs, on="experimento", how="inner")


def rodar(obs, coefs_por_gmr, coluna_gmr):
    res = []
    for r in obs.itertuples():
        lon, lat, elev = CELULA[r.local]
        g = float(getattr(r, coluna_gmr))
        sem = pd.Timestamp(r.data_semeadura)
        # O DSSAT repete o dia juliano da 1ª safra (1980, bissexto): escolher o dia/mês de 1980 com o mesmo dia
        # juliano da semeadura real, e convertê-lo no doy_sem de tf.rodar (dias após 10/set, base 2001/02).
        ano80 = 1980 if sem.month >= 7 else 1981
        d80 = pd.Timestamp(ano80, 1, 1) + pd.Timedelta(days=sem.dayofyear - 1)
        doy = (pd.Timestamp(2001 if d80.month >= 7 else 2002, d80.month, d80.day) - pd.Timestamp(2001, 9, 10)).days
        t = dict(id=0, lon=lon, lat=lat, elev=elev, arq="px", gmr=g, doy_sem=doy, **coefs_por_gmr(g))
        out = pd.DataFrame(tf.rodar(t))
        if "erro" in out:
            raise RuntimeError(f"{r.experimento}: {out.erro.iloc[0]}")
        s = out[out.data_semeadura == sem.strftime("%Y-%m-%d")]
        assert len(s) == 1, (r.experimento, sem)
        s = s.iloc[0]
        lin = dict(experimento=r.experimento, local=r.local, cultivar=r.cultivar, gmr=g, cenario_gmr=coluna_gmr)
        for e in ESTADIOS:
            lin[f"{e}_obs"] = getattr(r, e)
            lin[f"{e}_sim"] = float(s[e])
            lin[f"{e}_erro"] = float(s[e]) - getattr(r, e)
        res.append(lin)
    return pd.DataFrame(res)


def erros(d, rotulo):
    L = []
    for e in ESTADIOS:
        x = d[f"{e}_erro"].dropna()
        L.append(dict(conjunto=rotulo, estadio=e, n=len(x), rmse=float(np.sqrt((x ** 2).mean())), vies=float(x.mean()),
                      mae=float(x.abs().mean())))
    return L


def posterior(arquivos, burn):
    """Mediana da posterior; aceita cadeias com 19 (v13) ou 20 parâmetros (v14, σ_R3 no fim, ignorado aqui)."""
    import emcee
    ch = np.concatenate([emcee.backends.HDFBackend(str(f), read_only=True).get_chain()[burn:] for f in arquivos], axis=1)
    return np.median(ch.reshape(-1, ch.shape[-1]), axis=0)


def avaliar_cadeia(arquivos, burn, rotulo, experimentos=None):
    """DSSAT real com a mediana de outra posterior (v14, cadeias da validação deixando um local de fora)."""
    OUT.mkdir(parents=True, exist_ok=True)
    tf._D = DSSAT(str(tf.rb.DIR_RUN / "evandro"))
    th = posterior(arquivos, burn)
    coefs = lambda g: {ct.CHAVE[K]: v for K, v in ct.coef_linha(th, g).items()}  # noqa: E731
    obs = tabela_obs()
    if experimentos:
        obs = obs[obs.experimento.isin(experimentos)]
    d = rodar(obs, coefs, "gmr_calibracao").assign(conjunto=rotulo)
    arq = OUT / "avaliacoes_cadeias.csv"
    if arq.exists():
        ant = pd.read_csv(arq)
        d = pd.concat([ant[ant.conjunto != rotulo], d])
    d.to_csv(arq, index=False)
    return d


if __name__ == "__main__" and len(sys.argv) > 1:
    # uso: validar_evandro.py <rótulo> <burn> <cadeia1.h5> [cadeia2.h5 ...] [--exp PI-1,PI-2]
    args = sys.argv[1:]
    exps = args[args.index("--exp") + 1].split(",") if "--exp" in args else None
    if "--exp" in args:
        i = args.index("--exp"); args = args[:i] + args[i + 2:]
    d = avaliar_cadeia(args[2:], int(args[1]), args[0], exps)
    d = d[d.conjunto == args[0]]
    print(d.round(1).to_string(index=False))
    print(pd.DataFrame(erros(d, args[0])).round(2).to_string(index=False))
elif __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    tf._D = DSSAT(str(tf.rb.DIR_RUN / "evandro"))
    th = ct.posterior()
    coefs = lambda g: {ct.CHAVE[K]: v for K, v in ct.coef_linha(th, g).items()}  # noqa: E731
    obs = tabela_obs()
    obs["gmr_comercial_7063"] = obs.gmr_comercial
    print(obs[["experimento", "cultivar", "gmr_tese", "gmr_calibracao"] + ESTADIOS].to_string(index=False))
    tudo = pd.concat([rodar(obs, coefs, c) for c in ["gmr_calibracao", "gmr_comercial", "gmr_tese"]])
    tudo.to_csv(OUT / "v13_fora_amostra.csv", index=False)
    # Conjunto de avaliação do protocolo: um registro por experimento; Teresina (TE-100, TE-50, TE-1) = TE-1; sem SO-2.
    aval = ~tudo.experimento.isin(["TE-100", "TE-50", "SO-2"])
    L = []
    for c, g in tudo[aval].groupby("cenario_gmr"):
        L += erros(g, f"v13 | GMR {c} | avaliação (9 exp.)")
    for c, g in tudo.groupby("cenario_gmr"):
        L += erros(g, f"v13 | GMR {c} | todos (12 exp.)")
    E = pd.DataFrame(L)
    E.to_csv(OUT / "v13_fora_amostra_erros.csv", index=False)
    print(tudo[tudo.cenario_gmr == "gmr_calibracao"].round(1).to_string(index=False))
    print(E.round(2).to_string(index=False))
