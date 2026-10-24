#!/bin/bash
# Extensão da v13 (22/10/2026): +6.000 passos por cadeia (total 9.000), para N/τ ≥ 25 na metade final (τ ≈ 170).
# Backup das cadeias de 3.000 passos em resultados/emulador_fenologia/v13_3000passos/.
set -e
export DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1
PASSOS=${1:-6000}
PY=.venv/bin/python
D=resultados/emulador_fenologia
L=resultados/log_pipeline_fenologia_v13_ext.txt
echo "$(date) início extensão (+$PASSOS)" >> $L
for s in a b; do
  OMP_NUM_THREADS=2 $PY src/calibrar_fenologia.py $PASSOS 40 --move mix --output $D/cadeia_v13_$s.h5 --resume >> resultados/log_cadeia_v13_$s.txt 2>&1 &
  sleep 30
done
wait
$PY src/comparar_cadeias.py $D/cadeia_v13_a.h5 $D/cadeia_v13_b.h5 > $D/comparacao_cadeias_v13.txt 2>&1
echo "$(date) fim" >> $L
