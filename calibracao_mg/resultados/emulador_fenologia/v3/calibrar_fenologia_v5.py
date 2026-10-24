"""Calibração bayesiana dos coeficientes de fenologia do CROPGRO-Soybean como funções do GMR.

θ_k(GMR) = a_k + b_k (GMR − 7), para CSDL, PPSEN, EM-FL, FL-SD e SD-PM (limitados às faixas do treino).
Observações com GMR conhecido: experimentos (dias após a emergência — ou após a semeadura, coluna `base` — até R1, R7
e R8) e catálogos da Embrapa
(ciclo até R8, média de locais × safras representativos). R8 = R7 + δ8.
Erros: t de Student (ν = 4) nos experimentos, com escalas σ_R1 e σ_R7 (também para R8), robusta a um
experimento atípico; gaussiano nos catálogos (σ_cat) e nas médias multilocais de Alliprandini et al. 2009 (σ_agr).
O modelo é o emulador de fenologia (resultados/emulador_fenologia/modelo_*.joblib). Amostrador: emcee.

Uso: .venv/bin/python src/calibrar_fenologia.py [passos] [caminhantes]
"""
import argparse
import sys
import time
import warnings
from pathlib import Path

import emcee
import joblib
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))
import emulador_fenologia as ef  # noqa: E402

DIR = ef.DIR
# Faixas do emulador (lote 2 do treino): PPSEN até o MINIMA do DSSAT e EM-FL além do MAXIMA (30), porque a v2 encostou
# em PPSEN = 0,20 e EM-FL = 30
FAIXAS = {"csdl": (11.5, 13.0), "ppsen": (0.129, 0.385), "em-fl": (16.5, 36.0), "fl-sd": (9.0, 18.5), "sd-pm": (25.0, 38.0)}
# (nome, mínimo, máximo, valor inicial): intercepto em GMR 7 e inclinação por unidade de GMR
PARAMS = ([(f"a_{k}", lo, hi, ini) for (k, (lo, hi)), ini in zip(FAIXAS.items(), [12.3, 0.32, 21.0, 14.0, 30.0])] +
          [("b_csdl", -0.6, 0.2, -0.25), ("b_ppsen", -0.03, 0.05, 0.01), ("b_em-fl", -2.0, 4.0, 0.7),
           ("b_fl-sd", -3.0, 3.0, 0.0), ("b_sd-pm", -3.0, 4.0, 0.5),
           ("delta8", 3.0, 20.0, 10.0), ("ln_s_r1", np.log(0.5), np.log(20), np.log(4)),
           ("ln_s_r7", np.log(0.5), np.log(25), np.log(6)), ("ln_s_cat", np.log(0.5), np.log(25), np.log(6)),
           ("ln_s_agr", np.log(0.5), np.log(25), np.log(4)),          # médias multilocais (Alliprandini et al., 2009)
           ("ln_s_r5", np.log(0.5), np.log(20), np.log(5))])          # R5 dos experimentos (v5); sempre por último
NOMES = [p[0] for p in PARAMS]
NU = 4.0              # graus de liberdade da t de Student nos experimentos


class Calibracao:
    def __init__(self):
        locais = pd.read_csv(DIR / "locais.csv", index_col="idx_local")
        self.obs = pd.read_csv(DIR / "observacoes_fenologia.csv")
        amb = ef.variaveis_ambiente(self.obs[["idx_local", "safra", "doy_sem"]], locais)
        self.X = pd.concat([self.obs[["gmr", "doy_sem"]], amb.drop(columns=[c for c in ["doy_sem"] if c in amb])], axis=1)
        for k in FAIXAS:
            self.X[k] = 0.0
        self.X["mg_int"] = np.clip(np.round(self.X.gmr), 5, 9)
        self.X = self.X[ef.colunas()]
        # linhas com as mesmas entradas (R1/R7/R8 do mesmo registro) são previstas uma vez só e depois expandidas
        chave = self.X.round(6).astype(str).agg("|".join, axis=1)
        codigos, unicos = pd.factorize(chave)
        self.inv = codigos
        self.X = self.X.iloc[pd.Series(range(len(chave))).groupby(codigos).first().values].reset_index(drop=True)
        self.modelos = {s: joblib.load(DIR / f"modelo_{s}.joblib") for s in ef.SAIDAS if (DIR / f"modelo_{s}.joblib").exists()}
        o = self.obs
        self.exp = (o.tipo == "experimento").values
        self.cat = ~self.exp
        if "r5" not in self.modelos:        # emulador anterior à v5: sem R5, as observações de R5 ficam de fora
            o.loc[o.estadio == "r5", "tipo"] = "sem_r5"
            self.exp = (o.tipo == "experimento").values
        self.est = o.estadio.values
        self.y = o.dias_ve.values
        # observações contadas da semeadura: soma de volta os dias até a emergência simulados
        self.da_semeadura = (o.get("base", pd.Series("emergencia", index=o.index)) == "semeadura").values
        cat = o[self.cat].reset_index()
        self.cat_grupos = cat.groupby("id_registro")["index"].apply(np.array).tolist()
        self.cat_y = np.array([o.dias_ve.values[g[0]] for g in self.cat_grupos])
        self.cat_est = np.array([o.estadio.values[g[0]] for g in self.cat_grupos])
        self.cat_agr = np.array([o.tipo.values[g[0]] == "agregado" for g in self.cat_grupos])

    def coeficientes(self, th):
        p = dict(zip(NOMES, th))
        dg = self.X.gmr.values - 7
        return {k: np.clip(p[f"a_{k}"] + p[f"b_{k}"] * dg, lo, hi) for k, (lo, hi) in FAIXAS.items()}

    def prever(self, th):
        X = self.X.copy()
        for k, v in self.coeficientes(th).items():
            X[k] = v
        e, r1, r5, r7 = (self.modelos[s].predict(X)[self.inv] if s in self.modelos else np.zeros(len(self.inv))
                         for s in ("e", "r1", "r5", "r7"))
        d8 = dict(zip(NOMES, th))["delta8"]
        e = np.where(self.da_semeadura, 0.0, e)
        return {"r1": r1 - e, "r5": r5 - e, "r7": r7 - e, "r8": r7 - e + d8}

    def log_post(self, th):
        for v, (_, lo, hi, _) in zip(th, PARAMS):
            if not lo <= v <= hi:
                return -np.inf
        p = dict(zip(NOMES, th))
        pr = self.prever(th)
        ll = 0.0
        for est, s in [("r1", np.exp(p["ln_s_r1"])), ("r5", np.exp(p["ln_s_r5"])), ("r7", np.exp(p["ln_s_r7"])),
                       ("r8", np.exp(p["ln_s_r7"]))]:
            m = self.exp & (self.est == est)
            r = self.y[m] - pr[est][m]
            ll += np.sum(-(NU + 1) / 2 * np.log1p((r / s) ** 2 / NU)) - m.sum() * np.log(s)
        s = np.where(self.cat_agr, np.exp(p["ln_s_agr"]), np.exp(p["ln_s_cat"]))
        cat_pred = np.array([pr[e][g].mean() for e, g in zip(self.cat_est, self.cat_grupos)])
        r = self.cat_y - cat_pred
        ll += np.sum(-0.5 * (r / s) ** 2 - np.log(s))
        return ll + (-0.5 * ((p["delta8"] - 11) / 2) ** 2)          # δ8 ~ N(11, 2): R8 − R7 observado em Paragominas (11 ± 3,5 d); com N(10, 4) δ8 e FL-SD não eram identificáveis

    def log_post_vet(self, TH):
        """log_post para várias cadeias de uma vez: uma única predição do emulador com as linhas empilhadas."""
        TH = np.atleast_2d(TH)
        lo, hi = np.array([p[1] for p in PARAMS]), np.array([p[2] for p in PARAMS])
        ok = np.all((TH >= lo) & (TH <= hi), axis=1)
        out = np.full(len(TH), -np.inf)
        if not ok.any():
            return out
        TH_ok = TH[ok]
        n, m = len(TH_ok), len(self.X)
        X = pd.concat([self.X] * n, ignore_index=True)
        dg = np.tile(self.X.gmr.values - 7, n)
        ix = {nome: i for i, nome in enumerate(NOMES)}
        for k, (a, b) in FAIXAS.items():
            X[k] = np.clip(np.repeat(TH_ok[:, ix[f"a_{k}"]], m) + np.repeat(TH_ok[:, ix[f"b_{k}"]], m) * dg, a, b)
        e, r1, r5, r7 = (self.modelos[s].predict(X).reshape(n, m)[:, self.inv] if s in self.modelos
                         else np.zeros((n, len(self.inv))) for s in ("e", "r1", "r5", "r7"))
        d8 = TH_ok[:, ix["delta8"]][:, None]
        e = np.where(self.da_semeadura[None, :], 0.0, e)
        pr = {"r1": r1 - e, "r5": r5 - e, "r7": r7 - e, "r8": r7 - e + d8}
        s_r1, s_r5, s_r7, s_cat, s_agr = (np.exp(TH_ok[:, ix[c]])[:, None]
                                          for c in ("ln_s_r1", "ln_s_r5", "ln_s_r7", "ln_s_cat", "ln_s_agr"))
        ll = np.zeros(n)
        for est, sg in [("r1", s_r1), ("r5", s_r5), ("r7", s_r7), ("r8", s_r7)]:
            msk = self.exp & (self.est == est)
            r = self.y[msk][None, :] - pr[est][:, msk]
            ll += (-(NU + 1) / 2 * np.log1p((r / sg) ** 2 / NU) - np.log(sg)).sum(1)
        cat_pred = np.stack([pr[e][:, g].mean(1) for e, g in zip(self.cat_est, self.cat_grupos)], axis=1)
        r = self.cat_y[None, :] - cat_pred
        sg = np.where(self.cat_agr[None, :], s_agr, s_cat)
        ll += (-0.5 * (r / sg) ** 2 - np.log(sg)).sum(1)
        ll += -0.5 * ((TH_ok[:, ix["delta8"]] - 11) / 2) ** 2
        out[ok] = ll
        return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Amostra a posterior da calibração de fenologia.")
    parser.add_argument("passos", nargs="?", type=int, default=3000, help="iterações a executar")
    parser.add_argument("caminhantes", nargs="?", type=int, default=40)
    parser.add_argument("--output", type=Path, help="arquivo HDF5 da cadeia (padrão: cadeia_fenologia.h5)")
    parser.add_argument("--resume", action="store_true", help="continua uma cadeia já existente no arquivo escolhido")
    parser.add_argument("--seed", type=int, default=0, help="semente aleatória para inicialização (padrão: 0)")
    parser.add_argument("--dispersao", type=float, default=1e-3,
                        help="desvio inicial dos caminhantes como fração da faixa de cada parâmetro (padrão 0,001)")
    parser.add_argument("--initial-chain", type=Path, help="inicializa os caminhantes pelo último estado de outra cadeia")
    parser.add_argument("--move", choices=["mix", "stretch", "de", "de_snooker"], default="mix",
                        help="movimento: mix (padrão atual), stretch, de ou de_snooker")
    args = parser.parse_args()
    passos, nw = args.passos, args.caminhantes
    c = Calibracao()
    t = time.time(); c.log_post(np.array([p[3] for p in PARAMS])); print(f"1 avaliação: {1000 * (time.time() - t):.0f} ms", flush=True)
    arq = (args.output or (DIR / "cadeia_fenologia.h5")).expanduser().resolve()
    arq.parent.mkdir(parents=True, exist_ok=True)
    existe = arq.exists()
    if existe and not args.resume:
        parser.error(f"{arq} já existe; use --resume para continuar ou --output para iniciar outra cadeia")
    if args.resume and not existe:
        parser.error(f"não existe cadeia para retomar: {arq}")
    if args.resume and args.initial_chain:
        parser.error("--resume e --initial-chain não podem ser usados juntos")

    backend = emcee.backends.HDFBackend(str(arq))
    if args.resume:
        cadeia = backend.get_chain()
        if cadeia.shape[1:] != (nw, len(PARAMS)):
            parser.error(f"dimensão da cadeia existente {cadeia.shape[1:]} difere de {(nw, len(PARAMS))}")
        estado_inicial = backend.get_last_sample()
        print(f"Retomando {arq} a partir do passo {backend.iteration}", flush=True)
    else:
        if args.initial_chain:
            origem = args.initial_chain.expanduser().resolve()
            if origem == arq:
                parser.error("--initial-chain e --output devem apontar para arquivos diferentes")
            if not origem.exists():
                parser.error(f"cadeia inicial não encontrada: {origem}")
            inicial = emcee.backends.HDFBackend(str(origem), read_only=True)
            estado = inicial.get_last_sample()
            c0 = estado.coords
            if c0.shape[0] != nw or c0.shape[1] > len(PARAMS):
                parser.error(f"dimensão da cadeia inicial {c0.shape} incompatível com {(nw, len(PARAMS))}")
            if c0.shape[1] < len(PARAMS):      # parâmetros novos (sempre no fim de PARAMS): valor inicial + ruído
                rng = np.random.default_rng(args.seed)
                novos = PARAMS[c0.shape[1]:]
                extra = np.array([p[3] for p in novos]) + args.dispersao * rng.standard_normal((nw, len(novos))) * \
                    np.array([p[2] - p[1] for p in novos])
                c0 = np.hstack([c0, extra])
                print(f"cadeia inicial com {estado.coords.shape[1]} parâmetros; acrescentados {[p[0] for p in novos]}", flush=True)
            estado_inicial = c0
        else:
            rng = np.random.default_rng(args.seed)
            lo, hi = np.array([p[1] for p in PARAMS]), np.array([p[2] for p in PARAMS])
            estado_inicial = np.array([p[3] for p in PARAMS]) + args.dispersao * rng.standard_normal((nw, len(PARAMS))) * (hi - lo)
            estado_inicial = np.clip(estado_inicial, lo + 1e-6 * (hi - lo), hi - 1e-6 * (hi - lo))
        backend.reset(nw, len(PARAMS))
        origem_msg = f"; inicialização de {args.initial_chain}" if args.initial_chain else ""
        print(f"Nova cadeia ({args.move}): {arq}{origem_msg}", flush=True)

    movimentos = {
        "mix": [(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2)],
        "stretch": emcee.moves.StretchMove(),
        "de": emcee.moves.DEMove(),
        "de_snooker": emcee.moves.DESnookerMove(),
    }
    amostrador = emcee.EnsembleSampler(nw, len(PARAMS), c.log_post_vet, vectorize=True, backend=backend,
                                       moves=movimentos[args.move])
    t0 = time.time()
    for k, _ in enumerate(amostrador.sample(estado_inicial, iterations=passos), 1):
        if k % 250 == 0:
            lp = amostrador.get_log_prob()[-1]
            print(f"{backend.iteration} passos totais (+{k}/{passos})  {time.time() - t0:.0f} s  log-post máx {lp.max():.1f}  aceitação acumulada {amostrador.acceptance_fraction.mean():.3f}", flush=True)
