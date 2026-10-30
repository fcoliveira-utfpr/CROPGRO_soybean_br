#!/bin/bash
# Etapas 2 a 5 em sequência. A etapa 2 é repetida até as 270 células terem simulação (o clima chega aos poucos
# durante a etapa 1); cada rodada pula as células feitas. Células que falham são registradas no log e puladas;
# uma rodada sem nenhuma célula nova encerra a etapa 2.
# Uso: nohup caffeinate -i src/pipeline.sh > resultados/logs/pipeline.txt 2>&1 &
cd "$(dirname "$0")/.."
PY=../calibracao_mg/.venv/bin/python
export DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1
total=$(($(wc -l < dados/processados/celulas.csv) - 1))
while true; do
    feitas=$(ls resultados/sim/*.parquet 2>/dev/null | wc -l | tr -d ' ')
    [ "$feitas" -ge "$total" ] && break
    echo "$(date +%H:%M) etapa 2: $feitas/$total feitas"
    $PY src/etapa2_simular.py 4 2>&1 | grep --line-buffered -v -i warn
    novas=$(ls resultados/sim/*.parquet 2>/dev/null | wc -l | tr -d ' ')
    if [ "$novas" -eq "$feitas" ]; then             # rodada sem célula nova: as que faltam falharam (ver log)
        echo "$(date +%H:%M) etapa 2: nenhuma célula nova; seguindo com $novas/$total"; break
    fi
done
echo "$(date +%H:%M) etapa 3"; $PY src/etapa3_indicadores.py 2>&1 | grep --line-buffered -v -i warn
echo "$(date +%H:%M) etapa 4"; $PY src/etapa4_classificacao.py 2>&1 | grep --line-buffered -v -i warn
echo "$(date +%H:%M) etapa 5"; $PY src/etapa5_validacao.py 2>&1 | grep --line-buffered -v -i warn
echo "$(date +%H:%M) fim"
