# Dados de Evandro H. F. M. da Silva (ESALQ/USP, orient. F. R. Marin)

Guardados para **validação independente** (não entram na calibração v13). Fonte principal: a tese de doutorado
(Silva, 2021, *Crop modeling for understanding yield-gap causes and the potential for sustainable intensification
of soybean in Brazil*, ESALQ/USP, doi 10.11606/T.11.2021.tde-12042022-143612). PDFs em `docs/artigos/` (fora do git).

## Experimentos (12, em 5 locais)
| Experimento | Local | Safra | Cultivar (GMR) | Semeadura | Tratamentos |
|---|---|---|---|---|---|
| PI1-100, PI2-100 | Piracicaba-SP | 2016/17, 2017/18 | BRS 399 (6.0) | 14/11/2016, 16/12/2017 | irrigado |
| TE-100, TE-50 | Teresina-PI | 2019 (inverno) | 8579RSF (8.0) | 26/07/2019 | irrigação 100% e 50% |
| PI-1, PI-2, PI-3 | Piracicaba-SP | 2017/18 a 2019/20 | TMG 7062 (6.0) | 26/12, 11/11, 27/11 | sequeiro/irrigado × 0/1000 kg N |
| SO-1, SO-2 | Sorriso-MT | 2018/19, 2019/20 | TMG 7063 (7.0), NS 7901 (8.0) | 20/11/2018, 01/11/2019 | sequeiro × 0/1000 N |
| CA-1 | Cruz Alta-RS | 2018/19 | 65i65RSF (6.0) | 28/11/2018 | sequeiro × 0/1000 N |
| TU-1 | Tupanciretã-RS | 2019/20 | 65i65RSF (6.0) | 03/12/2019 | sequeiro/irrigado × 0/1000 N |
| TE-1 | Teresina-PI | 2019 (inverno) | 8579RSF (8.0) | 23/07/2019 | sequeiro/irrigado × 0/1000 N |

## Arquivos
| Arquivo | Conteúdo | Origem |
|---|---|---|
| `experimentos_tese.csv` | local, coordenadas, safra, cultivar, semeadura, colheita, tratamentos | Tabs. 2.1 e 3.1 e texto |
| `fenologia_tese_cap2.csv` | R1/R3/R5/R7 (dias após semeadura), 4 experimentos | Tab. 2.2 (texto) |
| `fenologia_tese_cap3.csv` | R1/R3/R5/R7 observados + simulados por tratamento, 8 experimentos | Figs. 3.6–3.13 (digitalizadas) |
| `crescimento_tese_cap3.csv` | IAF, folha, haste, grãos e biomassa por tratamento ao longo do ciclo (884 pontos) | Figs. 3.6–3.13 (digitalizadas) |
| `crescimento_tese_cap2.csv` | grãos e biomassa total, 4 experimentos (55 pontos) | Fig. 2.2 (digitalizada) |
| `produtividade_tese_tab3_4.csv` | produtividade, óleo e proteína observados e simulados, 26 tratamentos | Tab. 3.4 (texto) |
| `coeficientes_tese_apG.csv` | coeficientes CROPGRO calibrados das 6 cultivares | Apêndice G (texto) |
| `digitalizacao/` | código da digitalização (cores dos marcadores, eixos pelos ticks) | — |

## Qualidade da digitalização
- **Fenologia (barras):** precisão ±1 d. Desvio fixo de 0,45 d descontado, medido em TE-1 contra a Tab. 2.2. As
  diferenças simulado − observado por fase reproduzem os rótulos (B) das figuras.
- **Séries:** grão máximo digitalizado = produtividade da Tab. 3.4 a ±3% em 18 de 26 tratamentos. Nos demais, há
  marcadores escondidos atrás de outros (PI-3 R) ou diferença real entre a última coleta e a colheita.
  Pontos coincidentes de tratamentos diferentes podem ter sido perdidos (o de baixo fica coberto).
- **Cuidados:**
  - SO-1 e SO-2 têm barras de fenologia quase idênticas na tese (possível painel repetido).
  - As tabelas de RMSE dentro das figuras se repetem entre experimentos (TE-1 = PI-1; TU-1 = PI-3); os pontos não.
  - O GMR das cultivares é o informado (arredondado).

## Outras publicações
- **Dissertação de mestrado (2018):** E4 e E5 em Piracicaba (BRS 399, 2015/16 e 2016/17) só em figuras, sem
  fenologia tabelada; E1–E3 são dados de outros autores (Viçosa, Sto. Antônio de Goiás, Paragominas). Não digitalizada.
- **Agric. Water Manag. 252:106929 (2021):** mesmos experimentos de Teresina do cap. 2.
- **Eur. J. Agron. 129:126329 (2021):** fechado. A transcrição anterior da Tab. 5 tinha cultivares e safras que não
  batem com a tese e foi descartada.
- **J. Agric. Sci. 161:279 (2023)** cita 13 experimentos com o mesmo protocolo: além dos 12 acima, há pelo menos mais um.
- **Fonte dos dados brutos:** Evandro (ehfmsilva@usp.br).
