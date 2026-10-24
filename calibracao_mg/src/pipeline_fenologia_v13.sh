#!/bin/bash
# Fenologia v13 (21/10/2026): igual à v12 = v11 (δ8 fixo em 11 d, PRIOR_LIT=1, DE mix) + faixa fisiológica (FAIXA_FISIO=1:
# FL-SD ≥ 10,5; SD-PM ≥ 25). A v11 convergiu com FL-SD/SD-PM nos limites do emulador e derrubou o Yp da BRS 284.
# v12 interrompida (~1.500 passos) em dois modos (a −1596, R1PPO ≈ 0,86; b −1610, R1PPO ≈ 0,54); as duas cadeias partem do modo a.
set -e
export DELTA8_FIXO=11
export PRIOR_LIT=1
export FAIXA_FISIO=1
PASSOS=${1:-3000}
PY=.venv/bin/python
D=resultados/emulador_fenologia
L=resultados/log_pipeline_fenologia_v13.txt
echo "$(date) início (δ8 fixo; prior da literatura; faixa fisiológica)" >> $L
i=1
for s in a b; do
  OMP_NUM_THREADS=2 $PY src/calibrar_fenologia.py $PASSOS 40 --move mix --initial-chain $D/v12_interrompida/cadeia_v12_a.h5 \
    --perturbar --dispersao 0.02 --seed $i --output $D/cadeia_v13_$s.h5 > resultados/log_cadeia_v13_$s.txt 2>&1 &
  i=$((i+1))
  sleep 30
done
wait
$PY src/comparar_cadeias.py $D/cadeia_v13_a.h5 $D/cadeia_v13_b.h5 > $D/comparacao_cadeias_v13.txt 2>&1
echo "$(date) fim" >> $L
