#!/bin/bash
# SBI final: simulações com a posterior da fenologia final → NPE → inversão dos ponto-safra MODIS.
# 21/10/2026: fenologia final = v12 (faixa fisiológica); a tentativa com a v11 foi interrompida.
# Mesma configuração do protótipo (100 mil simulações) para comparar com a posterior exploratória v3b.
set -e
export DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1
N=${1:-100000}
CADEIA=${2:-resultados/emulador_fenologia/cadeia_v12_a.h5}
L=resultados/sbi/log_final.txt
echo "$(date) início simulações ($N, $CADEIA)" > $L
.venv/bin/python src/sbi_simulador.py $N $CADEIA final >> $L 2>&1
echo "$(date) início treino" >> $L
.venv_sbi/bin/python src/sbi_treinar.py final >> $L 2>&1
echo "$(date) fim" >> $L
