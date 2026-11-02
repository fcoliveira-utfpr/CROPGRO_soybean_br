"""Roda o TSE original (tse_mac.py) em todas as classes de GMR, em paralelo, e guarda os resultados.

Por classe (instalação tse/raiz_<classe>, preparada por preparar_tse.py), em sequência:
  grade    Run Range setup com N_GRADE valores por coeficiente + Multi TRT comb (caminho documentado para vários
           experimentos; manual do TSE, apêndice A.2);
  reducao  Run Range Reduction (3 valores, o padrão da interface; 4 fases) + Multi TRT comb;
  eco      depois da grade: cultivar fixada no ótimo da grade e R1PPO (arquivo .ECO, caixa "Ecotype") em grade de
           N_ECO valores entre 0,2 e 1,2 + Multi TRT comb — a sequência cultivar → ecótipo do próprio TSE.
Coeficientes da cultivar: CSDL, PPSEN, EM-FL, FL-SD e SD-PM, entre os limites fisiológicos da calibração bayesiana
(cf.FAIXAS), para os dois métodos buscarem na mesma caixa. Eventos: os 4 do PhenologicalEventsDAP.txt padrão.

Saída: resultados/tse/<metodo>_<classe>/ (log do TSE, pasta Compute, CUL do workspace) e resultados/tse/resumo.csv.
Uso: .venv/bin/python tse/rodar_tse.py [n_processos] [classes...]
"""
import json
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

import arquivos_genotipo as ag

AQUI = Path(__file__).resolve().parent
OUT = AQUI.parent / "resultados" / "tse"
PY37 = os.environ.get("PY37", "/private/tmp/claude-501/-Users-coted-sh-Python-DSSAT/6dd8fd9f-ecb9-480e-bb79-14ef9f28138c/"
                              "scratchpad/tools/py37/bin/python3.7")
FAIXAS = {"CSDL": (11.5, 14.0), "PPSEN": (0.129, 0.385), "EM-FL": (16.5, 40.0), "FL-SD": (10.5, 22.0), "SD-PM": (25.0, 38.0)}
N_GRADE, N_ECO = 4, 9
METODOS = {"grade": ("faixa", N_GRADE), "reducao": ("reducao", 3), "eco": ("faixa", N_ECO)}


def otimo_multi(log):
    """Coeficientes de 'Selected optimals' do bloco Multi TRT do log do TSE."""
    bloco = log.split("Multi TRT optimum processing:")[-1].split("Multi TRT optimum finished")[0]
    crit = dict(re.findall(r"\n([A-Za-z_()]+):\t([-\d.]+)", bloco))
    crit.update({"todos": v for v in re.findall(r"All selected: \|Avg\(Sim-obs\)\|:\t([-\d.]+)", bloco)[-1:]})
    coef = {k: float(v) for k, v in re.findall(r"\n(CSDL|PPSEN|EM-FL|FL-SD|SD-PM|R1PPO)\s+\t?([-\d.]+)", bloco)}
    return coef, crit


def rodar_classe(classe, t):
    cod = f"{int(round(classe * 10)):03d}"
    raiz, curto = AQUI / f"raiz_{cod}", f"/tmp/tse{cod}"
    g = t[t.classe == classe]
    linhas = []
    for nome, (metodo, n) in METODOS.items():
        destino = OUT / f"{nome}_{cod}"
        if (destino / "log_tse.txt").exists():
            log = (destino / "log_tse.txt").read_text()
        else:
            ws = raiz / "DSSAT48" / "TSE_workspace"
            shutil.rmtree(ws, ignore_errors=True)
            cul = raiz / "DSSAT48" / "GENOTYPE" / "SBGRO048.CUL"
            original = AQUI.parent / "produtos" / "dados" / "SBGRO048_BR.CUL"
            coefs = {k: [a, b, n] for k, (a, b) in FAIXAS.items()}
            if nome == "eco":                      # cultivar fixada no ótimo multi-tratamento da grade
                otimo, _ = otimo_multi((OUT / f"grade_{cod}" / "log_tse.txt").read_text())
                cul.write_text(ag.substituir(original.read_text(), f"BRG{cod}", otimo))
                coefs = {"R1PPO": [0.2, 1.2, n]}
            else:
                shutil.copy(original, cul)
            cfg = {"cultura": "Soybean", "modelo": "CRGRO", "cultivar": f"BRG{cod}", "arquivo": "ECO" if nome == "eco" else "CUL",
                   "tratamentos": [[a, int(k)] for a, k in zip(g.arquivo, g.TRNO)], "modo": "fenologia",
                   "variaveis": "todos", "coeficientes": coefs, "metodo": metodo}
            arq = raiz / f"cfg_{nome}.json"
            arq.write_text(json.dumps(cfg))
            t0 = time.time()
            with open(raiz / f"saida_{nome}.txt", "w") as f:
                subprocess.run([PY37, str(AQUI / "tse_mac.py"), str(arq)], cwd=AQUI, stdout=f, stderr=subprocess.STDOUT,
                               env={**os.environ, "TSE_RAIZ": str(raiz), "TSE_CURTO": curto, "QT_QPA_PLATFORM": "offscreen"})
            destino.mkdir(parents=True, exist_ok=True)
            shutil.copy(raiz / "log_tse.txt", destino)
            shutil.copy(raiz / f"saida_{nome}.txt", destino)
            if (ws / "Compute").exists():
                shutil.copytree(ws / "Compute", destino / "Compute", dirs_exist_ok=True)
            for f in ("SBGRO048.CUL", "TemporalFiles/combination_check.txt"):
                if (ws / f).exists():
                    shutil.copy(ws / f, destino)
            log = (destino / "log_tse.txt").read_text()
            (destino / "tempo_s.txt").write_text(f"{time.time() - t0:.0f}\n")
        coef, crit = otimo_multi(log)
        n_dssat = log.count("\nDSSAT ")
        linhas.append(dict(metodo=nome, classe=classe, n_trat=len(g), n_simulacoes=n_dssat, **coef,
                           criterio=json.dumps(crit), erros=log.count("TSE encounterd very important error")))
        print(f"{nome} {classe}: {coef} {crit} ({n_dssat} simulações)", flush=True)
    return linhas


if __name__ == "__main__":
    n_proc = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    t = pd.read_csv(AQUI / "tratamentos_tse.csv")
    classes = [float(c) for c in sys.argv[2:]] or sorted(t.classe.unique(), key=lambda c: -(t.classe == c).sum())
    OUT.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(n_proc) as ex:
        res = [l for ls in ex.map(lambda c: rodar_classe(c, t), classes) for l in ls]
    r = pd.DataFrame(res).sort_values(["metodo", "classe"])
    antigo = OUT / "resumo.csv"
    if antigo.exists():
        r = pd.concat([pd.read_csv(antigo), r]).drop_duplicates(["metodo", "classe"], keep="last")
    r.to_csv(antigo, index=False)
    print(r.to_string(index=False))
