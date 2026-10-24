#!/bin/bash
# Fenologia v9 (19/10/2026): igual à v8 (δ8 fixo em 11 d), com movimento DE (mix: DEMove 0,8 + DESnooker 0,2),
# porque as cadeias v8 (stretch) não ficaram estacionárias (log-post ainda subindo; aceitação 5–7%; a em −1633, b em −1613).
set -e
export DELTA8_FIXO=11
PASSOS=${1:-4000}
PY=.venv/bin/python
D=resultados/emulador_fenologia
L=resultados/log_pipeline_fenologia_v9.txt
echo "$(date) início (δ8 fixo = $DELTA8_FIXO)" >> $L
for s in a b; do
  OMP_NUM_THREADS=2 $PY src/calibrar_fenologia.py $PASSOS 40 --move mix --initial-chain $D/cadeia_v8_$s.h5 --output $D/cadeia_v9_$s.h5 > resultados/log_cadeia_v9_$s.txt 2>&1 &
  sleep 30
done
wait
$PY src/comparar_cadeias.py $D/cadeia_v9_a.h5 $D/cadeia_v9_b.h5 > $D/comparacao_cadeias_v9.txt 2>&1
echo "$(date) fim" >> $L
