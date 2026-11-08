#!/bin/bash
# Validação da v14 deixando um local do Evandro de fora (critério C2, protocolo.md): espera a v14 terminar; para cada
# local, cadeia curta (2.000 passos, 40 caminhantes) partindo do último estado da v14 (cadeia a), sem as observações
# do local; descarte de 1.000 passos; DSSAT real nos experimentos do local com a mediana. 2 cadeias por vez.
# Obs.: Cruz Alta-RS de Alliprandini (2009, médias multilocais) é outro estudo e fica na calibração.
set -e
export DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1 R3=1
PY=.venv/bin/python
D=resultados/emulador_fenologia
L=resultados/logs/log_loso_v14.txt
until grep -q "fim" resultados/logs/log_pipeline_fenologia_v14.txt 2>/dev/null; do sleep 300; done
echo "$(date) início" >> $L
declare -a NOMES=("Piracicaba-SP" "Teresina-PI" "Sorriso-MT (Silva 2021)" "Cruz Alta-RS (Silva 2021)" "Tupanciretã-RS")
declare -a COD=(pira tere sorr cruz tupa)
declare -a EXPS=("PI1-100,PI2-100,PI-1,PI-2,PI-3" "TE-1,TE-100,TE-50" "SO-1,SO-2" "CA-1" "TU-1")
rodar() {
  k=$1
  EXCLUIR_LOCAL="${NOMES[$k]}" OMP_NUM_THREADS=2 $PY src/calibrar_fenologia.py 2000 40 --move mix \
    --initial-chain $D/cadeia_v14_a.h5 --perturbar --dispersao 0.01 --seed $((60 + k)) \
    --output $D/loso/cadeia_v14_sem_${COD[$k]}.h5 > resultados/logs/log_loso_v14_${COD[$k]}.txt 2>&1
  $PY src/validar_evandro.py "v14 sem ${COD[$k]}" 1000 $D/loso/cadeia_v14_sem_${COD[$k]}.h5 --exp "${EXPS[$k]}" \
    >> resultados/logs/log_loso_v14_${COD[$k]}.txt 2>&1
  echo "$(date) ${COD[$k]} ok" >> $L
}
mkdir -p $D/loso
rodar 0 & sleep 20; rodar 1 & wait
rodar 2 & sleep 20; rodar 3 & wait
rodar 4
$PY src/validar_evandro.py "v14 completa" 4500 $D/cadeia_v14_a.h5 $D/cadeia_v14_b.h5 >> $L 2>&1
echo "$(date) fim" >> $L
