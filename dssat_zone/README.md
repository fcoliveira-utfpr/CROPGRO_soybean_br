# dssat_zone — zoneamento da soja no Paraná com o CROPGRO-Soybean calibrado

Zoneamento de risco e de janelas de semeadura da soja no Paraná com o CSM-CROPGRO-Soybean (DSSAT 4.8) e as cultivares
brasileiras por grupo de maturidade relativa (GMR) calibradas em [`../calibracao_mg`](../calibracao_mg) (fenologia
θ(GMR) v13 + crescimento da BRS 284). O objetivo é ir além das classes de 20, 30 e 40% de risco do ZARC: distribuição
completa da produtividade relativa, riscos por fase fenológica, janela do milho safrinha e incerteza.

Resultado principal: ver [`relatorio_tecnico.md`](relatorio_tecnico.md). Em resumo:
- 270 células de 0,25° (5,8 Mha de soja), 1,69 milhão de safras simuladas;
- na detecção das quebras do IBGE, a produtividade relativa simulada supera um índice do tipo ZARC: AUC 0,83 contra 0,79 (IC90 da diferença 0,013–0,067); série estadual r 0,83 contra 0,77;
- 96% das opções no nível de 20%, como no ZARC. O ganho está nas métricas contínuas, nos riscos por fase, na safrinha e na janela ótima multicritério.

## Decisões (29/10/2026; revistas em 30/10)

| Item | Decisão |
|---|---|
| Região | Paraná: 270 células de 0,25° com mais de 1.000 ha de soja (MapBiomas col. 10, 2023) |
| Grade | 0,25°; clima do Xavier (BR-DWGD, 0,1°) no pixel do centro da célula |
| Safras | 1980/81 a 2024/25 (45) |
| Cultivares | `BRG050`–`BRG075` (GMR 5,0 a 7,5, de 0,5 em 0,5) de `SBGRO048_BR.CUL/.ECO`, como publicadas |
| Semeadura | 12 decêndios de 1/set a 21/dez |
| Solos | 3 solos de Battisti & Sentelhas (2019) pela classe textural do MapBiomas Solo sob a soja da célula, pesados pela área (≥ 5%). Perfis MapBiomas + SoilGrids pela PTF de Medrado & Lima ficam como sensibilidade (`ZN_SOLO=ptf`), porque deram Yw/Yp ≈ 0,3 e acompanharam pior o IBGE na calibração. A AWC da ANA não entra (o DSSAT calcula a água ao alcance das raízes); fica como chave opcional de reescalonamento do DUL, desligada |
| Simulação | Executável do DSSAT chamado direto (FileX com até 99 tratamentos, 45 safras por chamada) |
| Produtividade | Só valores relativos: Yw/Yp entre locais, Yw relativo ao melhor manejo dentro da célula |
| Fenologia | R7/R8 corrigidos pela data de semeadura (k = 0,287 d/d antes de 10/nov) só nos indicadores de calendário |
| Perda | β · (1 − ry/ry_normal); β = 0,395 pela frequência de quebras do IBGE (sensibilidade: regressão 0,259, variância 0,482) |
| Milho safrinha | Probabilidade de colher a soja a tempo da data-limite do milho (geada ≤ 20%) na célula |
| Camada econômica | Fica para depois |

## Indicadores por célula × solo × GMR × decêndio

1. Yw/Yp e Yw: mediana, P20, CVaR20 e P(perda > 20/30/40%).
2. Seca forte (estresse ≥ 0,5) na floração e no enchimento (OVERVIEW do DSSAT).
3. Geada entre a emergência e o R7 corrigido (Tmin ≤ 2 °C).
4. Calor na floração (≥ 3 dias com Tmax ≥ 36 °C em R1–R5).
5. Chuva ≥ 80 mm nos 10 dias após o R8 corrigido.
6. Data do R8, duração do ciclo e probabilidade de colher a tempo do milho safrinha.

## Classificação

Nível de risco hídrico como no ZARC (menor de 20/30/40% com P(perda > 20%) ≤ nível). Restrições: geada > 20%, calor > 20%, chuva na colheita > 40%. Classes: recomendado / recomendado com restrição (motivo indicado) / não recomendado. Janela ótima por célula × solo: fronteira de Pareto entre Yw mediano, CVaR20 de Yw e P(safrinha). Tabela municipal no formato do ZARC (município × textura/AD × grupo × decêndio).

## Validação

Anomalias da produtividade municipal do IBGE (PAM; tendência LOWESS), 336 municípios, 1981–2024; comparação com um índice ISNA do tipo ZARC (IN SPA/MAPA nº 1/2022) nos mesmos solos e clima. Pendente: portaria do ZARC e sinistros do seguro (download manual).

## Etapas

| Etapa | Conteúdo | Situação |
|---|---|---|
| 0 | Custo (35 ms por safra; idêntico ao DSSATTools); viés fenológico por época e correção do R7/R8 | concluída |
| 1 | Grade, solos (com limpeza da textura e AD do ZARC), AWC da ANA, clima | concluída |
| 2 | Simulações fatoriais (523 célula × solo × 6 GMR × 12 decêndios × 45 safras; 9,1 h) | concluída |
| 3 | Indicadores por safra e data-limite da safrinha | concluída |
| 4 | Estatísticas, nível de risco, classes e janela ótima | concluída |
| 5 | Validação com o IBGE e comparação com o ISNA | concluída |
| 6 | Figuras, tabela municipal, relatório técnico | concluída |

## Limitações conhecidas

- O modelo capta só 50–60% do alongamento do ciclo nas semeaduras de set–out (etapa 0): R7/R8 corrigidos pela data de semeadura nos indicadores de calendário; a produtividade das semeaduras precoces fica subestimada e os decêndios de setembro têm menor confiança.
- O modelo exagera a amplitude das quebras (desvio 2× o do IBGE): a perda é recalibrada por β, e o nível de risco depende dessa escolha.
- Crescimento da BRS 284 sem ganho genético: produtividades absolutas não são interpretadas.
- Doenças, pragas e excesso de chuva não são simulados.

## Estrutura

- `src/config.py` — caminhos e constantes; importa os módulos de `../calibracao_mg/src`.
- `src/dssat_direto.py` — escrita do FileX com vários tratamentos, chamada do executável, leitura do Summary e do OVERVIEW.
- `src/etapa0_custo.py`, `src/etapa0_validar.py`, `src/etapa0_fenologia.py`, `src/fenologia_correcao.py` — etapa 0.
- `src/etapa1_grade.py` — células, solos, AWC, limpeza, clima.
- `src/etapa2_simular.py` — simulações (`resultados/sim/`, fora do git).
- `src/etapa3_indicadores.py` — indicadores por safra.
- `src/etapa4_classificacao.py` — estatísticas, β, classes, janela ótima.
- `src/etapa5_validacao.py` — IBGE, ISNA, métricas.
- `src/etapa6_figuras.py` — figuras, tabela municipal, números do relatório.
- `src/zarc_proxy.py` — AD e ISNA do tipo ZARC.
- `src/pipeline.sh` — etapas 2 a 5 em sequência.
- `resultados/` — tabelas (algumas compactadas em `.csv.gz`), figuras e relatório.

Ambiente: `../calibracao_mg/.venv` (mesmos pacotes).
