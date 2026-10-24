#!/bin/bash
# Fenologia v12 (21/10/2026): v11 (δ8 fixo em 11 d, PRIOR_LIT=1, DE mix) + faixa fisiológica (FAIXA_FISIO=1:
# FL-SD ≥ 10,5; SD-PM ≥ 25). A v11 convergiu com FL-SD/SD-PM nos limites do emulador e derrubou o Yp da BRS 284.
# Inícios distintos: a ← v11 b (R1PPO alto); b ← v10 a (modo de R1PPO ≈ 0,55, FL-SD ≈ 15).
set -e
export DELTA8_FIXO=11
export PRIOR_LIT=1
export FAIXA_FISIO=1
PASSOS=${1:-4000}
PY=.venv/bin/python
D=resultados/emulador_fenologia
L=resultados/log_pipeline_fenologia_v12.txt
echo "$(date) início (δ8 fixo; prior da literatura; faixa fisiológica)" >> $L
OMP_NUM_THREADS=2 $PY src/calibrar_fenologia.py $PASSOS 40 --move mix --initial-chain $D/cadeia_v11_b.h5 --perturbar \
  --dispersao 0.01 --seed 1 --output $D/cadeia_v12_a.h5 > resultados/log_cadeia_v12_a.txt 2>&1 &
sleep 30
OMP_NUM_THREADS=2 $PY src/calibrar_fenologia.py $PASSOS 40 --move mix --initial-chain $D/v10_interrompida/cadeia_v10_a.h5 --perturbar \
  --dispersao 0.01 --seed 2 --output $D/cadeia_v12_b.h5 > resultados/log_cadeia_v12_b.txt 2>&1 &
wait
$PY src/comparar_cadeias.py $D/cadeia_v12_a.h5 $D/cadeia_v12_b.h5 > $D/comparacao_cadeias_v12.txt 2>&1
echo "$(date) fim" >> $L
