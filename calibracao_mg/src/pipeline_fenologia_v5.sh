#!/bin/bash
# Fenologia v5 (17/10/2026): junta a v4 (Passo Fundo, catálogo MS, Santana 2022) com o R5.
# Lote 4 (3 locais novos) + lote 5 (ressimulação de 150 amostras × 73 locais guardando R3/R5/R8 do Evaluate.OUT)
# → observações (com R5) → emulador (e, r1, r5, r7) → 2 cadeias iniciadas nas v3c (ln_s_r5 acrescentado) → comparação.
# Uso (na pasta calibracao_mg): bash src/pipeline_fenologia_v5.sh [passos]
set -e
PASSOS=${1:-3000}
PY=.venv/bin/python
D=resultados/emulador_fenologia
L=resultados/log_pipeline_fenologia_v5.txt
echo "$(date) início" >> $L
mkdir -p $D/v3; cp $D/modelo_*.joblib $D/validacao_metricas.csv $D/observacoes_fenologia.csv $D/v3/
$PY src/treino_fenologia.py 300 4 4 >> $L 2>&1
$PY src/treino_fenologia.py 150 4 5 >> $L 2>&1
$PY src/observacoes_fenologia.py >> $L 2>&1
$PY src/validar_emulador_fenologia.py >> $L 2>&1
for s in a b; do
  OMP_NUM_THREADS=2 $PY src/calibrar_fenologia.py $PASSOS 40 --move stretch --seed 7 --initial-chain $D/cadeia_v3c_$s.h5 --output $D/cadeia_v5_$s.h5 > resultados/log_cadeia_v5_$s.txt 2>&1 &
  sleep 30   # evita corrida na criação dos links do DSSATTools ao importar
done
wait
$PY src/comparar_cadeias.py $D/cadeia_v5_a.h5 $D/cadeia_v5_b.h5 > $D/comparacao_cadeias_v5.txt 2>&1
echo "$(date) fim" >> $L
