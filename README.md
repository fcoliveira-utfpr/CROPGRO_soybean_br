# Calibração nacional do CROPGRO-Soybean por grupo de maturidade e grupo de maturidade efetivo da soja no Brasil

Análises, dados e produtos de uma calibração dos coeficientes do CSM-CROPGRO-Soybean (DSSAT 4.8) como **funções do
grupo de maturidade relativa (GMR)** para todo o Brasil, da **inversão da fenologia por satélite (MODIS)** para obter o
GMR efetivo e a data de semeadura da soja entre 2000 e 2024, e da **limitação hídrica da produtividade** simulada com
essas informações. Este repositório é a base de dois artigos científicos (seção [Artigos](#artigos)).

> **English summary.** This repository contains the code, data and products of a nationwide calibration of the
> CROPGRO-Soybean (DSSAT) phenology as functions of the relative maturity group (RMG) for Brazil, using a machine-learning
> emulator and hierarchical Bayesian inference, and of the amortized inversion (neural posterior estimation) of MODIS
> phenology into the effective RMG and sowing date of soybean at 224 grid points from 2000 to 2024. It provides
> ready-to-use DSSAT cultivar files for Brazil, a notebook that runs DSSAT in the browser (Google Colab) without any
> installation, and a static atlas. It supports two papers, one agronomic and one computational. Code comments and
> internal notes are mostly in Portuguese; the products (`calibracao_mg/produtos/`) are in English.

---

## Sumário

- [Artigos](#artigos)
- [Principais resultados](#principais-resultados)
- [Produtos para uso](#produtos-para-uso)
- [Estrutura do repositório](#estrutura-do-repositório)
- [Fluxo de trabalho](#fluxo-de-trabalho)
- [Como reproduzir](#como-reproduzir)
- [Dados e fontes](#dados-e-fontes)
- [O que não está no repositório](#o-que-não-está-no-repositório)
- [Documentação interna](#documentação-interna)
- [Como citar, licença e contato](#como-citar-licença-e-contato)

---

## Artigos

O projeto resultou em dois manuscritos independentes, em preparação para revistas da Springer. Os manuscritos não estão
no repositório enquanto não forem publicados (ficam em `_nao_sobe/artigo/`, fora do controle de versão).

| | Artigo agronômico | Artigo computacional |
|---|---|---|
| Título provisório | *Maturity-group-based phenology of CROPGRO-Soybean for Brazil and the effective maturity group of soybean from 2000 to 2024* | *Emulator-based Bayesian calibration and amortized inference of soybean phenology across Brazil with a process-based crop model* |
| Revista-alvo | Agronomy for Sustainable Development | Machine Learning for Computational Science and Engineering |
| Foco | Coeficientes θ(GMR) nacionais; GMR efetivo e semeadura de 2000 a 2024; limitação hídrica (Yw/Yp) validada com o IBGE; cultivares por ponto para planejamento | Emulador validado por local; verossimilhança hierárquica; identificabilidade ao longo de 13 versões da posterior; restrição fisiológica e consequência para a produtividade; inferência amortizada (NPE) comparada com MCMC por registro |
| Partes do repositório | etapas 1–4, 7–10 | etapas 5–8 |

---

## Principais resultados

**Fenologia calibrada (versão final v13).** Seis coeficientes como funções lineares do GMR (CSDL, PPSEN, EM-FL, FL-SD,
SD-PM e R1PPO), calibrados com 746 termos de verossimilhança: 453 observações de 8 estudos de campo em 13 locais,
médias multilocais de 88 cultivares (Alliprandini et al., 2009) e 117 registros de catálogos de cultivares.
O modelo reproduz as médias multilocais de R1 e R8 com RMSE de 3,2 e 1,9 dias (d de Willmott de 0,93 e 0,99).

| Coeficiente | GMR 5 | GMR 7 | GMR 9 |
|---|---|---|---|
| CSDL (h) | 12,88 | 12,39 | 11,90 |
| PPSEN (1/h) | 0,228 | 0,253 | 0,279 |
| EM-FL (dias fototérmicos) | 21,1 | 28,4 | 35,8 |
| FL-SD (dias fototérmicos) | 13,9 | 10,8 | 10,5 |
| SD-PM (dias fototérmicos) | 25,9 | 25,1 | 25,0 |
| R1PPO (1/h) | 0,95 | 0,86 | 0,77 |

**Emulador.** Gradient boosting treinado com 2,81 milhões de safras simuladas em 73 locais; com locais inteiros fora do
treino, RMSE de 1,2 dia em R1 e 2,0 dias em R7, cerca de 3×10⁴ vezes mais rápido que o DSSAT.

**GMR efetivo e semeadura (MODIS, 4.030 ponto-safra em 224 pontos).** A inversão amortizada recupera o GMR com erro de
0,75 e cobertura de 88% do intervalo de 90%, com posteriores praticamente iguais às da MCMC por registro, a um custo cerca
de 600 vezes menor. O GMR efetivo caiu 0,25 unidade por década entre 2000 e 2024 (de 9,0 para 8,2 em MT e de 8,6 para
7,7 no Centro-Oeste), e a semeadura antecipou 3,9 dias por década. O padrão espacial coincide com o zoneamento de GMR de
Ferreira et al. (2026): 0,136 GMR por grau de latitude, contra 0,14 naquele trabalho.

**Limitação hídrica.** Yw/Yp de 0,74 no Brasil, de 0,92 em MT/RO e Norte a 0,53 no Sul; as anomalias sem tendência do Yw
simulado acompanham as do rendimento do IBGE (r = 0,76 por região e safra).

**Limitações conhecidas.** Os coeficientes são efetivos (a partição entre R1PPO e as durações reprodutivas é fracamente
identificável e foi restringida pela faixa publicada); o crescimento é o da BRS 284 (Battisti e Sentelhas, 2019), sem o
ganho genético recente, então a produtividade simulada serve para comparações relativas; o modelo subestima o ciclo das
semeaduras precoces no Sul; e o GMR efetivo do Sul fica acima do recomendado.

---

## Produtos para uso

Em [`calibracao_mg/produtos/`](calibracao_mg/produtos/) (em inglês; detalhes no [README dos produtos](calibracao_mg/produtos/README.md)):

| Produto | Para quê |
|---|---|
| [`dados/SBGRO048_BR.CUL`](calibracao_mg/produtos/dados/SBGRO048_BR.CUL) e [`.ECO`](calibracao_mg/produtos/dados/SBGRO048_BR.ECO) | Cultivares brasileiras para o DSSAT 4.8: genéricas por GMR (`BRG050`–`BRG095`) e uma por ponto da grade (`BR0001`–`BR0224`). Saídas idênticas às do pipeline (verificado por `src/validar_cul.py`) |
| [`dados/*.csv`](calibracao_mg/produtos/dados/) | Coeficientes por GMR (com intervalo de 90%), por ponto e por ponto-safra; GMR efetivo, semeadura, Yp, Yw e rendimento do IBGE |
| [`notebooks/usage_example.ipynb`](calibracao_mg/produtos/notebooks/usage_example.ipynb) | Roda o DSSAT com essas cultivares para qualquer local. **Feito para quem sabe pouco de programação: roda no navegador pelo Google Colab, sem instalar nada.** Basta mudar local, safras e datas de semeadura numa célula e executar tudo; o notebook baixa o clima da NASA/POWER, monta a cultivar e o solo, roda o DSSAT e faz os gráficos de data de colheita e limitação hídrica por data de semeadura e GMR, como apoio ao planejamento |
| [`pagina/index.html`](calibracao_mg/produtos/pagina/index.html) | Atlas estático (pronto para o GitHub Pages): mapa do GMR efetivo, da tendência, da semeadura e da limitação hídrica; série e cultivar de cada ponto; curvas θ(GMR); downloads |

Para publicar o atlas no GitHub Pages, aponte o Pages para a pasta `calibracao_mg/produtos/pagina/` (ou copie o
conteúdo dela para `docs/`). Para o botão do Colab, use
`https://colab.research.google.com/github/<usuário>/<repositório>/blob/main/calibracao_mg/produtos/notebooks/usage_example.ipynb`
e atualize `REPO_RAW` na célula de dados do notebook.

---

## Estrutura do repositório

```
.
├── README.md                     este arquivo
├── .gitignore                    o que não sobe (ver abaixo)
└── calibracao_mg/
    ├── README.md                 etapas do projeto e notas técnicas (português)
    ├── requirements.txt          ambiente principal (.venv)
    ├── requirements_sbi.txt      ambiente da inferência amortizada (.venv_sbi)
    ├── src/                      todo o código (Python e pipelines em shell)
    ├── notebooks/                01–05 por etapa; 06_figuras_resultados.ipynb com as figuras dos artigos
    ├── dados/
    │   ├── brutos/ibge/          IBGE PAM (tabela 1612) por município, 1980–2024
    │   ├── processados/          pontos da grade, solos do MapBiomas, municípios, contorno dos estados
    │   └── literatura/           tabelas extraídas da literatura (coeficientes, fenologia observada, catálogos)
    ├── resultados/
    │   ├── bs2019/               replicação de Battisti e Sentelhas (2019): comparações e figuras
    │   ├── solo_mapbiomas/       teste dos perfis de solo
    │   ├── morris/               análise de sensibilidade
    │   ├── emulador_fenologia/   modelos do emulador, validação, observações, comparações das cadeias e emuladores
    │   ├── fenologia_satelite/   fenologia MODIS por ponto e safra
    │   ├── sbi/                  NPE treinado, validação e GMR efetivo por ponto-safra
    │   ├── produtividade/        Yp, Yw, IBGE por célula, limitação hídrica e anomalias
    │   ├── planejamento/         demonstração semeadura × GMR em 5 pontos
    │   ├── logs/                 logs das rodadas (calibrações, treinos, pipelines)
    │   ├── figuras/              figuras em português (PNG e PDF)
    │   └── figuras_en/           as mesmas figuras em inglês
    └── produtos/                 cultivares, tabelas, notebook de uso e atlas (inglês)
```

---

## Fluxo de trabalho

As etapas abaixo estão na ordem em que os dados fluem. Todos os comandos partem de `calibracao_mg/`.

| # | Etapa | Scripts principais | Saídas |
|---|---|---|---|
| 1 | Dados: grade de 462 pontos de Battisti e Sentelhas (2019), clima diário BR-DWGD (Xavier) 1980–2025, soja e solo do MapBiomas, IBGE | `digitalizar_pontos_bs2019.py`, `digitalizar_zonas_bs2019.py`, `gee.py`, `baixar_xavier_bs2019.py`, `baixar_xavier_extensao.py`, `solo_mapbiomas.py`, `baixar_ibge_municipios.py` | `dados/processados/`, `dados/brutos/` |
| 2 | Replicação de Battisti e Sentelhas (2019) e escolha da configuração (ET de Priestley-Taylor; 3 solos com peso pela textura) | `bs2019.py`, `rodar_bs2019.py`, `analise_bs2019.py`, `sensibilidade_bs2019.py`, `rodar_solo_mb.py`, `ptf.py` | `resultados/bs2019/`, `resultados/solo_mapbiomas/`; notebooks 02 e 03 |
| 3 | Sensibilidade de Morris (15 coeficientes, 9.600 simulações) | `faixas.py`, `morris.py` | `resultados/morris/`; notebook 04 |
| 4 | Observações de fenologia com GMR conhecido (experimentos, Alliprandini, catálogos) | `catalogos_embrapa.py`, `observacoes_fenologia.py` | `dados/literatura/`, `resultados/emulador_fenologia/observacoes_fenologia.csv` |
| 5 | Emulador de fenologia (2,81 milhões de safras em 73 locais; validação com locais fora do treino) | `treino_fenologia.py`, `emulador_fenologia.py`, `validar_emulador_fenologia.py`, `comparar_emuladores.py`, `sensibilidade_emulador.py` | `resultados/emulador_fenologia/modelo_*.joblib`, `validacao_metricas.csv`, `comparacao_emuladores.csv` |
| 6 | Calibração bayesiana hierárquica (13 versões; final v13) | `calibrar_fenologia.py`, `comparar_cadeias.py`, `decompor_logvero.py`, `pipeline_fenologia_v*.sh` | cadeias `cadeia_v*.h5` (fora do git), `comparacao_cadeias_v*.txt` |
| 7 | Fenologia MODIS e inversão amortizada (NPE) do GMR efetivo e da semeadura; comparação com MCMC por registro | `fenologia_satelite.py`, `rodar_fenologia_satelite.py`, `sbi_simulador.py`, `sbi_treinar.py`, `pipeline_sbi_final.sh`, `mcmc_por_registro.py` | `resultados/fenologia_satelite/`, `resultados/sbi/` |
| 8 | Produtividade potencial e limitada por água por ponto e safra; limitação hídrica e validação com o IBGE | `rodar_produtividade.py`, `lacuna_produtividade.py`, `limitacao_hidrica.py` | `resultados/produtividade/` |
| 9 | Produtos e demonstração para planejamento | `gerar_produtos.py`, `validar_cul.py`, `gerar_notebook_uso.py`, `gerar_pagina.py` (+ `pagina_modelo.html`), `demo_planejamento.py` | `produtos/`, `resultados/planejamento/` |
| 10 | Figuras dos artigos (português e inglês) | `figuras.py`, `gerar_notebook_figuras.py`, `traducao_en.json` | `notebooks/06_figuras_resultados.ipynb`, `resultados/figuras*/` |

Os scripts `emulador.py`, `treino_emulador.py`, `validar_emulador.py`, `validar_emulador_fen.py` e o notebook 05 são do
emulador antigo de produtividade, substituído pelo de fenologia; ficam como registro.

**Variáveis de ambiente da calibração final.** Todo script que lê a posterior v13 deve rodar com
`DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1` (δ8 fixo em 11 dias, prior da literatura e faixa fisiológica de FL-SD e SD-PM).

---

## Como reproduzir

**1. Ambientes (Python 3.11).**

```bash
cd calibracao_mg
python3.11 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python -m ipykernel install --user --name calibracao_mg
python3.11 -m venv .venv_sbi && .venv_sbi/bin/pip install --only-binary=grpcio -r requirements_sbi.txt
```

**2. DSSAT.** O DSSATTools traz executáveis para Linux e Windows. No macOS, compile o DSSAT 4.8.2 (instruções em
[`calibracao_mg/README.md`](calibracao_mg/README.md), seção Ambiente) e substitua `.venv/.../DSSATTools/bin/dscsm048`.
O DSSAT não aceita caminhos longos: as rodadas usam a pasta curta `calibracao_mg/.run/`.

**3. Google Earth Engine** (só para baixar o clima e o MapBiomas de novo): `earthengine authenticate`, projeto
`mapbiomas-brazil` em `src/gee.py`.

**4. Ordem dos comandos** (tempos aproximados num laptop de 2 núcleos):

```bash
# clima, solo e IBGE (etapa 1)
.venv/bin/python src/baixar_xavier_extensao.py && .venv/bin/python src/baixar_ibge_municipios.py
# treino e validação do emulador (etapa 5; ~2 dias de simulação em lotes)
for lote in 1 2 3 4 5 6; do .venv/bin/python src/treino_fenologia.py <amostras_por_local> 4 $lote; done   # ver o cabeçalho de src/treino_fenologia.py
.venv/bin/python src/validar_emulador_fenologia.py
# fenologia MODIS (etapa 7, Earth Engine)
.venv/bin/python src/rodar_fenologia_satelite.py
# calibração final (etapa 6; ~5 h por versão)
src/pipeline_fenologia_v13.sh 3000          # depois estendida a 9.000 passos (pipeline_fenologia_v13_ext*.sh)
# satélite e NPE (etapa 7; ~1 h)
src/pipeline_sbi_final.sh 100000 resultados/emulador_fenologia/cadeia_v13_a.h5
# produtividade, limitação hídrica (etapa 8; ~5 h)
DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1 .venv/bin/python src/rodar_produtividade.py resultados/emulador_fenologia/cadeia_v13_a.h5 final 4
.venv/bin/python src/lacuna_produtividade.py && .venv/bin/python src/limitacao_hidrica.py
# produtos (etapa 9)
DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1 .venv/bin/python src/gerar_produtos.py
.venv/bin/python src/gerar_notebook_uso.py && .venv/bin/python src/gerar_pagina.py
# figuras (etapa 10): português e inglês
.venv/bin/python src/gerar_notebook_figuras.py
cd notebooks && ../.venv/bin/jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.kernel_name=calibracao_mg 06_figuras_resultados.ipynb
FIG_LANG=en ../.venv/bin/jupyter nbconvert --to notebook --execute --output-dir /tmp --ExecutePreprocessor.kernel_name=calibracao_mg 06_figuras_resultados.ipynb
```

Os comandos exatos de cada versão estão nos scripts `src/pipeline_fenologia_v*.sh`, que trazem no cabeçalho o motivo de cada
versão; os logs das rodadas estão em `resultados/logs/`.

---

## Dados e fontes

| Dado | Fonte | Uso |
|---|---|---|
| Clima diário (Tmax, Tmin, chuva, radiação, UR, vento), 0,1°, 1980–2025 | BR-DWGD (Xavier et al., 2022), via Google Earth Engine | simulações e emulador |
| Área de soja e textura do solo | MapBiomas (cobertura e solo) | máscara do satélite e peso dos solos |
| Fenologia por satélite (NDVI, 250 m) | MODIS MOD13Q1 e MYD13Q1 | GMR efetivo e semeadura |
| Produtividade municipal da soja | IBGE, Produção Agrícola Municipal, tabela 1612 do SIDRA | validação das anomalias |
| Perfis de solo e cultivar de referência (BRS 284) | Battisti e Sentelhas (2019), Tabelas S1 e S2 | configuração e crescimento |
| Fenologia observada e coeficientes publicados | 8 estudos de campo, Alliprandini et al. (2009), catálogos da Embrapa e outros, 21 calibrações publicadas | calibração e priors |
| Clima no notebook de uso | NASA/POWER (API diária) | exemplo de uso fora da grade |

As referências completas estão nos artigos. Os PDFs consultados não são redistribuídos (direitos autorais).

---

## O que não está no repositório

**Fica no computador, fora do git, no lugar onde os scripts o leem** (listado no [`.gitignore`](.gitignore)):

| Arquivo ou pasta | Tamanho | Como obter |
|---|---|---|
| `.venv/`, `.venv_sbi/`, `.tools/`, `calibracao_mg/.run/` | 2,8 GB | recriar (seção anterior) |
| `calibracao_mg/codex.md` | — | diário interno dos autores (decisões e comandos de retomada) |
| `calibracao_mg/dados/brutos/xavier/` | 287 MB | `src/baixar_xavier_*.py` (Earth Engine) |
| `resultados/emulador_fenologia/tabela_treino.parquet`, `treino/`, `validacao_predicoes.parquet` | ~380 MB | `src/treino_fenologia.py` e `src/validar_emulador_fenologia.py`; a depositar no Zenodo |
| Cadeias MCMC `resultados/emulador_fenologia/**/*.h5` (todas as versões; a final é `cadeia_v13_a/b.h5`, 94 MB cada) | ~1,1 GB | `src/pipeline_fenologia_v*.sh`; **as cadeias finais devem ir para o Zenodo**, porque reproduzi-las leva dias |
| Simulações do SBI (`resultados/sbi/*.npz`) | 37 MB | `src/sbi_simulador.py` |
| Simulações da replicação e do solo (`resultados/bs2019/sim*/`, `resultados/solo_mapbiomas/sim_*/`) | 160 MB | `src/rodar_bs2019.py`, `src/rodar_solo_mb.py` |

**Foi movido para [`_nao_sobe/`](_nao_sobe/)** (ver o [LEIA-ME](_nao_sobe/LEIA-ME.md) da pasta): os manuscritos dos dois
artigos, os PDFs da literatura e dos catálogos, a exploração inicial (`initial_codes/`), o emulador antigo de
produtividade (669 MB), logs do TensorBoard e restos de rodadas do DSSAT.

---

## Documentação interna

- [`calibracao_mg/README.md`](calibracao_mg/README.md): etapas, decisões técnicas (solos, pedotransferência, ecótipos),
  inventário da literatura e instruções de compilação do DSSAT.
- `src/pipeline_fenologia_v*.sh`: um script por versão da calibração, com o motivo da versão no cabeçalho.
- O diário de laboratório dos autores (`calibracao_mg/codex.md`) é interno e não está no repositório.

---

## Como citar, licença e contato

**Citação.** Até a publicação dos artigos, cite o repositório: Oliveira, F. C. *et al.* (2026). *Calibração nacional do
CROPGRO-Soybean por grupo de maturidade e GMR efetivo da soja no Brasil* [código e dados]. GitHub/Zenodo. DOI a definir.

**Licença.** Os dados e produtos (tabelas, cultivares `.CUL`/`.ECO`, figuras, atlas e resultados) estão sob a licença
[Creative Commons Atribuição 4.0 Internacional (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/deed.pt-br): podem
ser usados, adaptados e redistribuídos, inclusive comercialmente, desde que os autores sejam citados (texto completo em
[`LICENSE-DATA`](LICENSE-DATA)). Dados de terceiros incluídos mantêm as licenças das fontes. A licença do código (`src/`
e notebooks) ainda será definida.

**Contato.** Universidade Tecnológica Federal do Paraná (UTFPR), Campus Santa Helena — fcoliveira@utfpr.edu.br

**Apoio.** CAPES, Fundação Araucária e UTFPR.
