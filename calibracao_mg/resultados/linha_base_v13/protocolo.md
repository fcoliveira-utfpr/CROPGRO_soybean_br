# Linha de base v13 e protocolo de comparação com a v14

Fixado em 06/11/2026, **antes** de qualquer rodada com os dados de Evandro H. F. M. da Silva.
Código da linha de base: commit `f32efc7`, tag local `v13-linha-base`. Métricas: `metricas.csv`.

## Decisões tomadas (06/11/2026)
1. Teresina (semeadura em julho) entra na calibração: o emulador será estendido (lote 7).
2. SO-2 fica de fora da calibração por ter fenologia idêntica a SO-1 na tese (possível painel repetido); entra só
   como verificação.
3. O crescimento calibrado com os dados do Evandro substitui o da BRS 284 **se** vencer na validação abaixo.
4. TSE: só se recalculam as métricas contra a v14; não se relança o lote (46 h).

## Observações
- Uma linha por experimento (a fenologia observada é comum aos tratamentos). TE-100, TE-50 e TE-1 = um registro
  (mesmo local, semeaduras a 3 d, fenologia idêntica).
- GMR: `gmr_calibracao` de `dados/literatura/evandro_silva/experimentos_tese.csv` (GMR comercial conferido; a tese
  usa o grupo do DSSAT arredondado). TMG 7063 em Sorriso: 7,0 na calibração, 6,3 no teste de sensibilidade.

## Critérios (definidos antes de rodar)
**Fenologia** — no DSSAT real, não no emulador:
- C1. RMSE de R1, R5 e R7 nos experimentos atuais (n = 100, 60, 193) não piora mais de 0,5 d em relação à v13.
- C2. Validação deixando um local de fora (5 locais do Evandro: Piracicaba, Teresina, Sorriso, Cruz Alta,
  Tupanciretã): RMSE de R1/R3/R5/R7 menor que o da v13 aplicada fora da amostra aos mesmos experimentos.
- C3. Viés de R7 nos experimentos atuais (v13: −8,7 d) reduzido em módulo.
- C4. Convergência: R-hat ≤ 1,1 e ESS ≥ 500.
- C5. Restrições: a v14 sem faixa fisiológica (e com δ8 livre, se possível) só é adotada se cumprir C1–C4 e não
  derrubar o Yp (FL-SD e SD-PM fora dos limites do emulador = reprovada).

**Crescimento** — deixando um experimento de fora, no DSSAT real:
- G1. Produtividade dos tratamentos irrigados + 1000 N (Tab. 3.4): RMSE menor que o da BRS 284 e que o dos
  coeficientes do apêndice G aplicados fora da amostra.
- G2. Tratamentos de sequeiro (Yw): RMSE e viés menores que os da BRS 284.
- G3. IAF máximo e biomassa final das séries digitalizadas: sem piora em relação à BRS 284 (dado secundário, σ maior).

**Cadeia (depois da adoção)** — mesmas métricas de `metricas.csv`:
- r das anomalias Yw × IBGE por região-ano não cai mais que 0,03 (v13: 0,757 no Brasil).
- A fração de ponto-safras com Yw < IBGE (v13: 42%) é reportada; espera-se queda se o novo crescimento for adotado.

Resultados que não cumprirem os critérios são reportados como tais, sem novo ajuste até cumprir.

## Emenda 1 (07/11/2026, 16h) — antes dos resultados finais da v14 (C2 completo e "v14 completa")
**O que já se sabia ao emendar:** v14 reprovada em C4 (R-hat de ln_s_r3 1,44; ln_s_cat 1,12), C1 (R7 +1,43 d) e
C3 (viés R7 −9,35 d); C2 desfavorável em Piracicaba e Teresina. Sorriso, Cruz Alta e Tupanciretã ainda não avaliados.

**Motivo:** os critérios C1 e C3 medem só os experimentos atuais, que concentram tratamentos em poucas regiões, e por
construção favorecem o modelo ajustado só a eles. Os dados de Silva (2021) entram pela representatividade: Teresina
(baixa latitude, semeadura de inverno), Sorriso (Cerrado de MT), Cruz Alta e Tupanciretã (RS), Piracicaba, além do R3.

**Mudança no modelo (v14b):** só o que causou a não convergência. A escala do R3 é a do R5 (`R3=compartilhado`; sai
ln_s_r3, 19 parâmetros como na v13). 2 cadeias × 40 caminhantes × 12.000 passos a partir do último estado da cadeia
v14 b (melhor log-post), perturbadas (dispersão 0,02). Mesmas DELTA8_FIXO=11, PRIOR_LIT=1, FAIXA_FISIO=1.

**Critérios para adotar a v14b no lugar da v13** (DSSAT real; descarte = metade da cadeia):
- E1. Convergência obrigatória: R-hat ≤ 1,1 em todos os parâmetros (antigo C4).
- E2. Critério principal: média, entre locais (peso igual por local), do RMSE por local de R1, R5 e R7; locais =
  os dos 209 tratamentos + os 5 de Silva (9 experimentos de avaliação). A v14b não pode ser pior que a v13 em nenhum
  estádio. Para a v13 os experimentos de Silva estão fora da amostra; para a v14b, dentro — é a comparação entre
  "não incluir" e "incluir" os dados.
- E3. Salvaguarda nos locais antigos: a mesma média por local, só nos locais dos 209 tratamentos, não piora mais de
  1,0 d em nenhum estádio. (Afrouxa o C1 — 0,5 d no RMSE agregado por tratamento — e passa a pesar locais, não
  tratamentos. Afrouxamento declarado.)
- Diagnósticos, sempre reportados e sem veto: C1, C3 (viés R7), R3 e C2. Se a v14b for adotada, roda-se a validação
  deixando um local de Silva de fora com ela (como na v14) e o resultado é reportado mesmo que desfavorável.
- Se a v14b não cumprir E1–E3, fica a v13, e os dados de Silva entram no artigo só como validação fora da amostra
  (resultados/evandro/v13_fora_amostra*.csv). Não há v14c para fazer cumprir.

## Emenda 2 (08/11/2026) — detalhamento da Fase 3 (crescimento), antes de qualquer rodada de calibração
Fenologia: v13 (decisão de 07/11). Só o teste da BRS 284 em PI-2 e TE-1 (configuração) foi rodado antes desta emenda.

**Simulação** (`src/crescimento_silva.py`): cultivar da cadeia (θ(GMR) da v13 no GMR de calibração + FL-SH/FL-LF da
Tab. S2); solo do apêndice F da tese; clima Xavier na célula; semeadura e densidade reais (Tab. 3.1; espaçamento 0,5 m);
início 30 d antes com 50% da água disponível. Irrigado = água e N desligados (≈ Yp); sequeiro = água ligada, N desligado.

**Calibração** (`src/calibrar_crescimento.py`): 7 coeficientes (LFMAX, SLAVR, SIZLF, XFRT, WTPSD, SFDUR, SDPDV), um
conjunto único (sem dependência do GMR), demais da BRS 284. Dados: tratamentos I-1000N (PI-1, PI-2, PI-3, TU-1, TE-1).
Função objetivo (média sobre os experimentos de treino):
J = média[((Ysim − Yobs)/Yobs / 0,10)²] + média[(nRMSE da série / 0,25)²], séries de biomassa e IAF do mesmo tratamento
(nRMSE = RMSE / média observada). Busca: hipercubo latino (4.000 conjuntos nas faixas de `FAIXAS`) + Nelder-Mead a
partir dos 3 melhores. Deixando um experimento de fora: 5 dobras.

**Comparações fora da amostra:**
- BRS 284 (B&S 2019; não ajustada a Silva).
- Apêndice G fora da amostra: mediana dos coeficientes de crescimento das cultivares do apêndice G diferentes da
  cultivar do experimento deixado de fora (com a fenologia da v13, como as demais). Apêndice G na própria cultivar:
  só referência (dentro da amostra da tese).

**Critérios** (como G1–G3, detalhados):
- G1. RMSE da produtividade dos 5 I-1000N fora da amostra < BRS 284 e < apêndice G fora da amostra.
- G2. Sequeiro R-1000N (8 experimentos, SO-2 incluído): RMSE e |viés| < BRS 284. Coeficientes: os da dobra que excluiu o
  experimento (PI-1, PI-2, PI-3, TU-1, TE-1) ou os do ajuste completo (SO-1, SO-2, CA-1, fora da calibração).
- G3. IAF máximo e biomassa final das séries (tratamentos com série, mesma regra de coeficientes): RMSE não maior que
  1,1 × o da BRS 284.
- Adotado: o ajuste completo (5 experimentos), se G1–G3 forem cumpridos; depois, critérios da cadeia (Fase 4).

## Resultado da Fase 3 (08/11/2026) — crescimento de Silva reprovado; fica a BRS 284
`resultados/evandro/crescimento/criterios_G.csv` (RMSE fora da amostra; kg/ha, IAF em m²/m²):
| Critério | Silva (dobras) | BRS 284 | Ap. G fora | Ap. G própria (ref.) | Resultado |
|---|---|---|---|---|---|
| G1 produtividade I-1000N (n = 5) | 1.326 (viés −497) | 1.132 (−611) | 1.231 | 1.059 | reprova |
| G2 produtividade R-1000N (n = 8) | 1.263 (−950) | 1.356 (−1.003) | 1.447 | 1.323 | cumpre |
| G3 biomassa final (n = 13) | 2.384 | 2.282 (limite 2.511) | 2.459 | 2.213 | cumpre |
| G3 IAF máximo (n = 13) | 1,55 | 1,31 (limite 1,44) | 1,55 | 1,78 | reprova |
Coeficientes instáveis entre dobras (LFMAX 1,16–1,45; WTPSD 0,155–0,216; SDPDV 1,81–2,46) e J ≈ 5–7 mesmo dentro
da amostra (erro relativo de produtividade ~20%): Piracicaba subestimada com qualquer coeficiente (PI-2 4,2–4,6 ×
6,4 t/ha), TU-1 e TE-1 irrigados superestimados. O erro não é de crescimento: vem do ciclo (R7 da v13 cedo nesses
locais) e do clima em ponto de grade. Fica a BRS 284; como fenologia e crescimento não mudam, a cadeia (Fase 4) não
precisa ser refeita.

## Arquivos (08/11/2026)
O estado oficial do repositório é a v13: `resultados/emulador_fenologia/` (modelos, locais, tarefas, observações) e
`dados/literatura/fenologia_observada.csv` (284 linhas, sem Silva). A v14 e a v14b (reprovadas) ficam em
`resultados/emulador_fenologia/emulador_v14/` (observações, locais e métricas no git; modelos e tarefas do lote 7 fora,
pelo tamanho) — `fenologia_observada_v14.csv` é a tabela com as 12 linhas de Silva (usar = 1 em 9). Os dados de Silva
usados na validação estão em `dados/literatura/evandro_silva/`.
