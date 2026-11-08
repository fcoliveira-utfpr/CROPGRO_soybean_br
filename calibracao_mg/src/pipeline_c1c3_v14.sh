#!/bin/bash
# C1/C3 do protocolo: DSSAT real nos 209 tratamentos com a mediana da v14 (logo que a v14 termina; ~30 s).
cd "$(dirname "$0")/.."
PY=.venv/bin/python; D=resultados/emulador_fenologia
until grep -q "fim" resultados/logs/log_pipeline_fenologia_v14.txt 2>/dev/null; do sleep 300; done
DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1 $PY src/avaliar_experimentos.py "v14" 4500 $D/cadeia_v14_a.h5 $D/cadeia_v14_b.h5 \
  > resultados/logs/log_c1c3_v14.txt 2>&1
echo "$(date) fim" >> resultados/logs/log_c1c3_v14.txt
