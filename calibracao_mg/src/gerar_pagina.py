"""Gera a página estática de exploração dos produtos (produtos/pagina/): mapa do GMR efetivo, séries por ponto,
coeficientes da cultivar e curvas θ(GMR). Todos os dados vão embutidos no HTML (sem servidor).

Saídas: produtos/pagina/index.html (documento completo, para o GitHub Pages) e produtos/pagina/artifact.html
(mesmo conteúdo sem <html>/<head>/<body>, para publicar como artifact).
Uso: .venv/bin/python src/gerar_pagina.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import figuras as fg  # noqa: E402

RAIZ = Path(__file__).resolve().parents[1]
DADOS = RAIZ / "produtos" / "dados"
OUT = RAIZ / "produtos" / "pagina"
LON0, LAT0, K = -74.5, 6.0, 10.0          # projeção equirretangular simples: 10 unidades por grau


def xy(lon, lat):
    return (lon - LON0) * K * np.cos(np.radians(15)), (LAT0 - lat) * K


def caminhos_estados():
    paths = []
    for poly in fg.ufs():
        pts, last = [], None
        for lon, lat in poly:
            p = tuple(round(v, 1) for v in xy(lon, lat))
            if p != last:
                pts.append(p)
                last = p
        if len(pts) > 3:
            paths.append("M" + "L".join(f"{x:g},{y:g}" for x, y in pts) + "Z")
    return paths


def main():
    p = pd.read_csv(DADOS / "coefficients_point.csv")
    s = pd.read_csv(DADOS / "coefficients_point_season.csv")
    th = pd.read_csv(DADOS / "theta_rmg.csv")
    cul = {l[:6]: l for l in open(DADOS / "SBGRO048_BR.CUL") if l[:2] == "BR"}
    eco = {l[:6]: l for l in open(DADOS / "SBGRO048_BR.ECO") if l[:1] in "EB"}
    pts = []
    for r in p.itertuples():
        x, y = xy(r.lon, r.lat)
        ss = s[(s.lon == r.lon) & (s.lat == r.lat)].sort_values("harvest_year")
        e_id = "E" + r.cultivar_id[2:] + "B"
        pts.append({
            "id": r.cultivar_id, "lon": r.lon, "lat": r.lat, "x": round(x, 2), "y": round(y, 2), "elev": round(r.elevation_m),
            "region": r.region.replace("Norte", "North").replace("norte MS", "northern MS"), "n": int(r.n_seasons),
            "rmg": round(r.rmg_mean, 2), "rmgRecent": round(r.rmg_recent, 2) if pd.notna(r.rmg_recent) else None,
            "rmgCul": round(r.rmg_for_cul, 2),
            "trend": round(r.rmg_trend_per_decade, 2) if pd.notna(r.rmg_trend_per_decade) else None,
            "sow": round(r.sowing_days_mean, 1),
            "sowTrend": round(r.sowing_trend_days_per_decade, 1) if pd.notna(r.sowing_trend_days_per_decade) else None,
            "ywyp": round(r.yw_yp_mean, 3) if pd.notna(r.yw_yp_mean) else None,
            "cv": round(r.yw_cv, 3) if pd.notna(r.yw_cv) else None,
            "soil": [round(v, 2) if pd.notna(v) else 0 for v in (r.soil_clay, r.soil_sandy_clay, r.soil_sandy_loam)],
            "coef": [round(float(p.loc[r.Index, k]), 4) for k in ["csdl", "ppsen", "em-fl", "fl-sd", "sd-pm", "r1ppo"]],
            "cul": cul.get(r.cultivar_id, "").rstrip("\n"), "eco": eco.get(e_id, "").rstrip("\n"),
            "series": [[int(a.harvest_year), round(a.rmg_q05, 2), round(a.rmg_q50, 2), round(a.rmg_q95, 2),
                        round(a.sowing_days_after_10sep, 1)] for a in ss.itertuples()],
        })
    theta = {k: th[k].round(4).tolist() for k in th.columns}
    dados = {"points": pts, "theta": theta, "states": caminhos_estados()}
    corpo = (Path(__file__).parent / "pagina_modelo.html").read_text()
    corpo = corpo.replace("/*__DATA__*/null", json.dumps(dados, separators=(",", ":")))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "artifact.html").write_text(corpo.replace("<!--__DOWNLOADS__-->", ""))
    # versão do GitHub Pages: pasta autossuficiente, com os arquivos para baixar (links funcionam fora do visualizador)
    import shutil
    (OUT / "dados").mkdir(exist_ok=True)
    arqs = ["coefficients_point.csv", "coefficients_point_season.csv", "theta_rmg.csv", "SBGRO048_BR.CUL",
            "SBGRO048_BR.ECO", "metadata.json"]
    for a in arqs:
        shutil.copy(DADOS / a, OUT / "dados" / a)
    shutil.copy(RAIZ / "produtos" / "notebooks" / "usage_example.ipynb", OUT / "dados" / "usage_example.ipynb")
    desc = {"coefficients_point.csv": "one row per grid point", "coefficients_point_season.csv": "one row per point and season",
            "theta_rmg.csv": "coefficients by RMG (5.0–9.5)", "SBGRO048_BR.CUL": "DSSAT cultivar file",
            "SBGRO048_BR.ECO": "DSSAT ecotype file", "metadata.json": "provenance",
            "usage_example.ipynb": "notebook: run DSSAT with these cultivars (Jupyter or Google Colab)"}
    lista = "".join(f'<li><a href="dados/{a}" download>{a}</a> <span class="lead">{d}</span></li>' for a, d in desc.items())
    dl = ('<section class="panel note"><h2>Download</h2><ul class="dl">' + lista + "</ul></section>")
    corpo_gh = corpo.replace("<!--__DOWNLOADS__-->", dl)
    doc = ("<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
           "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1, viewport-fit=cover\">\n</head>\n<body>\n"
           + corpo_gh + "\n</body>\n</html>\n")
    (OUT / "index.html").write_text(doc)
    # cópia para o GitHub Pages: Settings → Pages → branch main, pasta /docs (o Pages só serve / ou /docs)
    docs = RAIZ.parent / "docs"
    shutil.rmtree(docs, ignore_errors=True)
    shutil.copytree(OUT, docs, ignore=shutil.ignore_patterns("artifact.html"))
    (docs / ".nojekyll").write_text("")          # serve os arquivos como estão, sem processamento do Jekyll
    print(OUT / "index.html", f"{len(doc) / 1024:.0f} KB", len(pts), "pontos")


if __name__ == "__main__":
    main()
