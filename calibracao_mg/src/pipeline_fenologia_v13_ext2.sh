#!/bin/bash
# Retomada da extensão da v13 após a parada de 22/10 (a: 4.343, b: 4.285 passos) até 9.000 passos cada.
set -e
export DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1
PY=.venv/bin/python
D=resultados/emulador_fenologia
L=resultados/log_pipeline_fenologia_v13_ext.txt
echo "$(date) retomada (a +4657, b +4715)" >> $L
OMP_NUM_THREADS=2 $PY src/calibrar_fenologia.py 4657 40 --move mix --output $D/cadeia_v13_a.h5 --resume >> resultados/log_cadeia_v13_a.txt 2>&1 &
sleep 30
OMP_NUM_THREADS=2 $PY src/calibrar_fenologia.py 4715 40 --move mix --output $D/cadeia_v13_b.h5 --resume >> resultados/log_cadeia_v13_b.txt 2>&1 &
wait
$PY src/comparar_cadeias.py $D/cadeia_v13_a.h5 $D/cadeia_v13_b.h5 > $D/comparacao_cadeias_v13.txt 2>&1
echo "$(date) fim" >> $L
