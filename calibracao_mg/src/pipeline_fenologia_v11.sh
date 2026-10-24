#!/bin/bash
# Fenologia v11 (21/10/2026): mesma posterior da v10 (δ8 fixo em 11 d, PRIOR_LIT=1, movimento DE mix). Na v10 as
# cadeias ficaram em dois modos; o modo b (R1PPO ≈ 1, FL-SD/SD-PM curtos) é melhor em quase todas as fontes
# (decomposição em v10_interrompida/) e a cadeia a era um modo local. As duas cadeias partem do estado final da
# v10 b, cada uma com ruído independente (2% da faixa), para testar a convergência dentro do modo dominante.
set -e
export DELTA8_FIXO=11
export PRIOR_LIT=1
PASSOS=${1:-3000}
PY=.venv/bin/python
D=resultados/emulador_fenologia
L=resultados/log_pipeline_fenologia_v11.txt
echo "$(date) início (δ8 fixo = $DELTA8_FIXO; prior da literatura; início no modo b da v10)" >> $L
i=1
for s in a b; do
  OMP_NUM_THREADS=2 $PY src/calibrar_fenologia.py $PASSOS 40 --move mix --initial-chain $D/v10_interrompida/cadeia_v10_b.h5 \
    --perturbar --dispersao 0.02 --seed $i --output $D/cadeia_v11_$s.h5 > resultados/log_cadeia_v11_$s.txt 2>&1 &
  i=$((i+1))
  sleep 30
done
wait
$PY src/comparar_cadeias.py $D/cadeia_v11_a.h5 $D/cadeia_v11_b.h5 > $D/comparacao_cadeias_v11.txt 2>&1
echo "$(date) fim" >> $L
