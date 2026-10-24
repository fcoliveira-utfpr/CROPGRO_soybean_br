#!/bin/bash
# Fenologia v6 (17/10/2026, noite): R1PPO calibrado, faixas ampliadas e efeito aleatório por ensaio × estádio.
# Lote 6 (180 amostras × 73 locais, novas faixas + R1PPO) → observações → emulador (R1PPO como entrada; lotes 1–5
# com o R1PPO do ecótipo) → 2 cadeias iniciadas nas v5 (a_r1ppo, b_r1ppo, ln_tau acrescentados) → comparação.
# Uso (na pasta calibracao_mg): bash src/pipeline_fenologia_v6.sh [passos]
set -e
PASSOS=${1:-4000}
PY=.venv/bin/python
D=resultados/emulador_fenologia
L=resultados/log_pipeline_fenologia_v6.txt
echo "$(date) início" >> $L
mkdir -p $D/v5; cp $D/modelo_*.joblib $D/validacao_metricas.csv $D/observacoes_fenologia.csv $D/comparacao_cadeias_v5.txt $D/v5/
$PY src/treino_fenologia.py 180 4 6 >> $L 2>&1
$PY src/observacoes_fenologia.py >> $L 2>&1
$PY src/validar_emulador_fenologia.py >> $L 2>&1
for s in a b; do
  OMP_NUM_THREADS=2 $PY src/calibrar_fenologia.py $PASSOS 40 --move stretch --seed 8 --dispersao 0.05 --initial-chain $D/cadeia_v5_$s.h5 --output $D/cadeia_v6_$s.h5 > resultados/log_cadeia_v6_$s.txt 2>&1 &
  sleep 30   # evita corrida na criação dos links do DSSATTools ao importar
done
wait
$PY src/comparar_cadeias.py $D/cadeia_v6_a.h5 $D/cadeia_v6_b.h5 > $D/comparacao_cadeias_v6.txt 2>&1
echo "$(date) fim" >> $L
