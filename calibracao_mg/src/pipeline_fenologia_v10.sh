#!/bin/bash
# Fenologia v10 (21/10/2026): igual à v9 (δ8 fixo em 11 d, movimento DE mix) + prior da literatura (PRIOR_LIT=1):
# FL-SD ~ N(13,9; 2,75) e SD-PM ~ N(29,95; 3,77) em GMR 5 e 9, R1PPO ~ N(0,5; 0,15). Cadeias iniciadas nos estados finais da v9 a/b (dois modos).
set -e
export DELTA8_FIXO=11
export PRIOR_LIT=1
PASSOS=${1:-4000}
PY=.venv/bin/python
D=resultados/emulador_fenologia
L=resultados/log_pipeline_fenologia_v10.txt
echo "$(date) início (δ8 fixo = $DELTA8_FIXO; prior da literatura)" >> $L
for s in a b; do
  OMP_NUM_THREADS=2 $PY src/calibrar_fenologia.py $PASSOS 40 --move mix --initial-chain $D/cadeia_v9_$s.h5 --output $D/cadeia_v10_$s.h5 > resultados/log_cadeia_v10_$s.txt 2>&1 &
  sleep 30
done
wait
$PY src/comparar_cadeias.py $D/cadeia_v10_a.h5 $D/cadeia_v10_b.h5 > $D/comparacao_cadeias_v10.txt 2>&1
echo "$(date) fim" >> $L
