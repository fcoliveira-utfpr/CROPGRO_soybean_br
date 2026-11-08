"""Por que a BRS 284 com a fenologia v13 erra a produtividade de Silva (2021) e a tese não: decomposição trocando um
fator por vez (tratamentos 1000N; clima Xavier, solo do apêndice F, densidade/semeadura reais em todos os nossos casos).

  A. tese        — simulado por Silva (Tab. 3.4): estação no local, cultivar própria calibrada, N e FAO-56 ligados.
  B. apG         — cultivar completa do apêndice G (fenologia + crescimento), na nossa configuração → A−B: clima e
                   configuração (Xavier × estação, Priestley-Taylor, N desligado).
  C. v13 + apG   — fenologia v13 θ(GMR) + crescimento do apêndice G → B−C: fenologia genérica por GMR.
  D. v13 + BRS 284 (cadeia) → C−D: cultivar (crescimento da BRS 284 × cultivares modernas).

Saída: resultados/evandro/crescimento/diagnostico.csv (tratamento × caso) e diagnostico_resumo.csv.
Uso: DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1 .venv/bin/python src/diagnostico_crescimento.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import crescimento_silva as cs  # noqa: E402

NOME = {"NS 7901": "NS7901"}


def main():
    g = pd.read_csv(cs.EV / "coeficientes_tese_apG.csv").set_index("coef")
    t = cs.tratamentos()
    t = t[t.n == "1000N"]
    L = []
    for _, r in t.iterrows():
        c = NOME.get(r.cultivar, r.cultivar)
        apg = {k: float(g.loc[k.upper(), c]) for k in cs.CRES}
        casos = {"B apG": dict(cres=None, cul=cs.cultivar_apg(c)),
                 "C v13 + apG": dict(cres=apg, cul=None),
                 "D v13 + BRS 284": dict(cres=cs.BRS284, cul=None)}
        base = dict(experimento=r.experimento, tratamento=r.tratamento, cultivar=c, obs=r.yield_obs_kgha)
        L.append(dict(base, caso="A tese", hwam=r.yield_sim_kgha, r7=np.nan))
        for nome, kw in casos.items():
            s, _ = cs.simular(r, kw["cres"], cul=kw["cul"])
            L.append(dict(base, caso=nome, hwam=s.get("hwam"), r7=s.get("r7")))
    d = pd.DataFrame(L)
    obs = cs.ve.tabela_obs()[["experimento", "r7"]].rename(columns={"r7": "r7_obs"})
    d = d.merge(obs, on="experimento", how="left")
    d["erro"] = d.hwam - d.obs
    cs.OUT.mkdir(parents=True, exist_ok=True)
    d.to_csv(cs.OUT / "diagnostico.csv", index=False)
    R = []
    for (agua, caso), x in d.assign(agua=d.tratamento.str[0]).groupby(["agua", "caso"]):
        R.append(dict(agua=agua, caso=caso, n=len(x), rmse=np.sqrt((x.erro ** 2).mean()), vies=x.erro.mean(),
                      erro_r7=(x.r7 - x.r7_obs).mean()))
    R = pd.DataFrame(R)
    R.to_csv(cs.OUT / "diagnostico_resumo.csv", index=False)
    pd.set_option("display.width", 200)
    print(R.round(0).to_string(index=False))
    print(d.pivot_table(index=["experimento", "tratamento", "obs", "r7_obs"], columns="caso", values=["hwam", "r7"])
          .round(0).to_string())


if __name__ == "__main__":
    main()
