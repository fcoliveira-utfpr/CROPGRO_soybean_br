#!/bin/bash
# Fenologia v8 (18/10/2026): igual à v7, com δ8 = R8 − R7 fixo em 11 d (Paragominas; DELTA8_FIXO), porque δ8 e SD-PM
# não são identificáveis só com R8. Cadeias iniciadas nos estados finais da v7; sem retreino do emulador.
set -e
export DELTA8_FIXO=11
PASSOS=${1:-3500}
PY=.venv/bin/python
D=resultados/emulador_fenologia
L=resultados/log_pipeline_fenologia_v8.txt
echo "$(date) início (δ8 fixo = $DELTA8_FIXO)" >> $L
for s in a b; do
  OMP_NUM_THREADS=2 $PY src/calibrar_fenologia.py $PASSOS 40 --move stretch --initial-chain $D/cadeia_v7_$s.h5 --output $D/cadeia_v8_$s.h5 > resultados/log_cadeia_v8_$s.txt 2>&1 &
  sleep 30
done
wait
$PY src/comparar_cadeias.py $D/cadeia_v8_a.h5 $D/cadeia_v8_b.h5 > $D/comparacao_cadeias_v8.txt 2>&1
echo "$(date) fim" >> $L
