#!/bin/bash
# Refaz os perfis (com o carbono do MapBiomas) e roda as duas variantes de capacidade de campo.
set -e
cd "$(dirname "$0")/.."
.venv/bin/python -W ignore - <<'PY' 2>&1 | grep -v Survey
import sys; sys.path.insert(0, "src"); import gee, solo_mapbiomas as sm, pandas as pd
gee.iniciar()
p = pd.read_csv("dados/processados/pontos_bs2019_digitalizados.csv")
fr = pd.read_csv("dados/processados/fracoes_solo_mapbiomas_soja2016.csv")
f = sm.perfis_celulas(p.merge(fr[["lon", "lat"]], on=["lon", "lat"]))
f.to_csv("dados/processados/perfis_solo_mapbiomas_soja2016.csv", index=False)
print("perfis:", f.shape, "| carbono MapBiomas 0–30 (t/ha) p10/p50/p90:", f.cos_mb_0_30.quantile([.1, .5, .9]).round(0).tolist())
PY
rm -rf resultados/solo_mapbiomas/sim_cc10kPa resultados/solo_mapbiomas/sim_cc33kPa
for h in 10 33; do .venv/bin/python -W ignore -u src/rodar_solo_mb.py $h 4 2>&1 | grep -v -i -E "warn|created"; done
echo FIM
