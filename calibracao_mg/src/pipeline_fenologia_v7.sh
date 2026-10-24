#!/bin/bash
# Fenologia v7 (18/10/2026): igual à v6, mas o efeito aleatório agrupa por estudo × local × safra × estádio (sem a data de
# semeadura), para que a resposta à época informe o R1PPO. Sem retreino do emulador; cadeias iniciadas nas v6.
set -e
PASSOS=${1:-4000}
PY=.venv/bin/python
D=resultados/emulador_fenologia
L=resultados/log_pipeline_fenologia_v7.txt
echo "$(date) início" >> $L
for s in a b; do
  OMP_NUM_THREADS=2 $PY src/calibrar_fenologia.py $PASSOS 40 --move stretch --initial-chain $D/cadeia_v6_$s.h5 --output $D/cadeia_v7_$s.h5 > resultados/log_cadeia_v7_$s.txt 2>&1 &
  sleep 30
done
wait
$PY src/comparar_cadeias.py $D/cadeia_v7_a.h5 $D/cadeia_v7_b.h5 > $D/comparacao_cadeias_v7.txt 2>&1
echo "$(date) fim" >> $L
