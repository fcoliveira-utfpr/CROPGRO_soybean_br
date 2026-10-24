#!/bin/bash
# Fenologia v3 (16/10/2026): lote 3 do treino (23 locais novos: MATOPIBA, PR, BA e Alliprandini 2009) → observações → emulador → 2 cadeias independentes → comparação.
# Uso (na pasta calibracao_mg): bash src/pipeline_fenologia_v3.sh [passos]
set -e
PASSOS=${1:-4000}
PY=.venv/bin/python
D=resultados/emulador_fenologia
L=resultados/log_pipeline_fenologia_v3.txt
echo "$(date) início" >> $L
$PY src/treino_fenologia.py 300 4 3 >> $L 2>&1
$PY src/observacoes_fenologia.py >> $L 2>&1
$PY src/validar_emulador_fenologia.py >> $L 2>&1
for s in 1 2; do
  OMP_NUM_THREADS=2 $PY src/calibrar_fenologia.py $PASSOS 40 --move stretch --seed $s --dispersao 0.05 --output $D/cadeia_v3_s$s.h5 > resultados/log_cadeia_v3_s$s.txt 2>&1 &
  sleep 30   # evita corrida na criação dos links do DSSATTools ao importar
done
wait
$PY src/comparar_cadeias.py $D/cadeia_v3_s1.h5 $D/cadeia_v3_s2.h5 > $D/comparacao_cadeias_v3.txt 2>&1
echo "$(date) fim" >> $L
