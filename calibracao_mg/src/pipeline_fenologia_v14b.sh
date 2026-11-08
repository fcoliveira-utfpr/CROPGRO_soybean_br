#!/bin/bash
# Fenologia v14b (07/11/2026; emenda 1 do protocolo, resultados/linha_base_v13/protocolo.md): v14 com a escala do R3
# igual à do R5 (R3=compartilhado; 19 parâmetros) e cadeias mais longas, a partir do último estado da cadeia v14 b.
# Ao fim: convergência, DSSAT real nos 209 tratamentos e nos 12 experimentos de Silva, critérios E2/E3 por local.
set -e
cd "$(dirname "$0")/.."
export DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1 R3=compartilhado
PASSOS=${1:-12000}; BURN=$((PASSOS / 2))
PY=.venv/bin/python
D=resultados/emulador_fenologia
L=resultados/logs/log_pipeline_fenologia_v14b.txt
echo "$(date) início v14b ($PASSOS passos)" >> $L
i=1
for s in a b; do
  OMP_NUM_THREADS=2 $PY src/calibrar_fenologia.py $PASSOS 40 --move mix --initial-chain $D/cadeia_v14_b.h5 \
    --perturbar --dispersao 0.02 --seed $((50 + i)) --output $D/cadeia_v14b_$s.h5 > resultados/logs/log_cadeia_v14b_$s.txt 2>&1 &
  i=$((i+1))
  sleep 30
done
wait
$PY src/comparar_cadeias.py $D/cadeia_v14b_a.h5 $D/cadeia_v14b_b.h5 > $D/comparacao_cadeias_v14b.txt 2>&1
$PY src/avaliar_experimentos.py "v14b" $BURN $D/cadeia_v14b_a.h5 $D/cadeia_v14b_b.h5 >> $L 2>&1
$PY src/validar_evandro.py "v14b completa" $BURN $D/cadeia_v14b_a.h5 $D/cadeia_v14b_b.h5 >> $L 2>&1
$PY src/avaliar_por_local.py >> $L 2>&1 || true
echo "$(date) fim" >> $L
