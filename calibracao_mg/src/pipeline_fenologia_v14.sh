#!/bin/bash
# Fenologia v14 (06/11/2026): v13 (δ8 = 11 d, prior da literatura, faixa fisiológica) + dados de Silva (2021, tese ESALQ:
# 9 experimentos com R1/R3/R5/R7 em Piracicaba, Teresina, Sorriso, Cruz Alta e Tupanciretã) + R3 na verossimilhança
# (σ_R3) + emulador retreinado com o lote 7 (locais do Evandro e semeadura de inverno nos trópicos) e a saída R3.
# Protocolo de comparação: resultados/linha_base_v13/protocolo.md. Partida: último estado da v13 (cadeia a), perturbado.
set -e
export DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1 R3=1
PASSOS=${1:-9000}
PY=.venv/bin/python
D=resultados/emulador_fenologia
L=resultados/logs/log_pipeline_fenologia_v14.txt
echo "$(date) início v14 ($PASSOS passos)" >> $L
i=1
for s in a b; do
  OMP_NUM_THREADS=2 $PY src/calibrar_fenologia.py $PASSOS 40 --move mix --initial-chain $D/cadeia_v13_a.h5 \
    --perturbar --dispersao 0.02 --seed $((40 + i)) --output $D/cadeia_v14_$s.h5 > resultados/logs/log_cadeia_v14_$s.txt 2>&1 &
  i=$((i+1))
  sleep 30
done
wait
$PY src/comparar_cadeias.py $D/cadeia_v14_a.h5 $D/cadeia_v14_b.h5 > $D/comparacao_cadeias_v14.txt 2>&1
echo "$(date) fim" >> $L
