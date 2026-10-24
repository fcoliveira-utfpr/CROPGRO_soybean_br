#!/bin/bash
# Depois da v12: espera o fim da calibração → SBI final → Yp/Yw por ponto × safra (opção A).
set -e
export DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1
until grep -q "fim" resultados/log_pipeline_fenologia_v12.txt 2>/dev/null; do sleep 120; done
src/pipeline_sbi_final.sh 100000 resultados/emulador_fenologia/cadeia_v12_a.h5
.venv/bin/python src/rodar_produtividade.py resultados/emulador_fenologia/cadeia_v12_a.h5 final 4 > resultados/log_produtividade.txt 2>&1
echo "$(date) fim" >> resultados/log_produtividade.txt
