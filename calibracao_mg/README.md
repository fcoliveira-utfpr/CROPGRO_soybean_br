# Calibração do CROPGRO-Soybean por grupo de maturidade (MG) para o Brasil

Os coeficientes do CROPGRO-Soybean (DSSAT 4.8) são estimados como **funções do MG**, θ(MG), e não como uma tabela por cultivar. A calibração é bayesiana hierárquica e combina experimentos da literatura, fenologia por satélite e produtividade municipal (IBGE).

## Etapas

| # | Etapa | Status |
|---|---|---|
| 01 | Faixas dos coeficientes (espaço de busca e priors) | feito: `notebooks/01_faixas_coeficientes.ipynb` |
| 02 | **Replicação de Battisti & Sentelhas (2019)** com o Xavier 0,1° e os 3 solos do artigo | **validada (v2)**: atingível a ±11% e potencial a ±3% do artigo por zona (zona 7: +9%); `notebooks/02_replicacao_bs2019.ipynb` |
| 02b | Solo do MapBiomas nas áreas de soja. v1: pesos dos 3 solos por textura (feito; r com o IBGE 0,40 → 0,44); v2: perfis completos com pedotransferência brasileira | feito (`notebooks/03_solo_mapbiomas.ipynb`). v1 (3 solos com peso pela textura do MapBiomas): r com o IBGE 0,40 → 0,44. v2 (perfis completos MapBiomas + SoilGrids + Medrado & Lima): r cai para 0,24–0,27, porque a função reduz a água disponível nos solos argilosos, o que vale para o Cerrado mas não para o Sul. **Base de solo da Fase B: v1** |
| 02c | Análise de sensibilidade de Morris (15 parâmetros, 20 pontos, 3 semeaduras, 9.600 rodadas) | **feito** (`notebooks/04_sensibilidade_morris.ipynb`). Selecionados 10: CSDL, PPSEN, EM-FL, FL-SD, SD-PM, LFMAX, SFDUR, SIZLF, fator_ad, prof_raiz |
| 03 | Observações de fenologia com GMR conhecido (8 estudos de campo, Alliprandini 2009, catálogos) | **feito**: 453 observações de experimentos, 176 médias multilocais, 117 registros de catálogos (`src/observacoes_fenologia.py`) |
| 04 | **Emulador de fenologia** (gradient boosting; 2,81 milhões de safras em 73 locais) | **feito**: locais fora do treino R1 1,19 d, R5 2,20 d, R7 1,97 d; comparado com MLP e ridge (`comparacao_emuladores.csv`) e sensibilidade da calibração à escolha (`sensibilidade_emulador.txt`) |
| 05 | **Calibração bayesiana** da fenologia (13 versões) | **feito: v13 final** (δ8 = 11 d, prior da literatura, faixa fisiológica FL-SD ≥ 10,5 e SD-PM ≥ 25; 2 conjuntos × 40 caminhantes × 9.000 passos; R-hat ≤ 1,08; ESS ≈ 760). Histórico nos scripts `src/pipeline_fenologia_v*.sh` e em `resultados/logs/` |
| 06 | Fenologia por satélite (MODIS) e **inversão amortizada** (NPE) do GMR efetivo e da semeadura | **feito**: 4.030 ponto-safra em 224 pontos; GMR RMSE 0,75, cobertura 0,88; MCMC por registro com posteriores iguais (`resultados/sbi/mcmc_por_registro.csv`) |
| 07 | Produtividade (opção A: crescimento da BRS 284 de B&S 2019) | **feito**: Yp e Yw por ponto × safra; lacuna absoluta Yw − IBGE descartada (Yw < IBGE em 42%; sem ganho genético); resultado = limitação hídrica Yw/Yp (Brasil 0,74; RS/SC 0,53; MT 0,92) |
| 08 | Validação com o IBGE (variação entre anos) | **feito**: anomalias sem tendência Yw × IBGE r = 0,76 por região × safra (`src/limitacao_hidrica.py`) |
| 09 | Produtos de uso público e demonstração de planejamento | **feito**: cultivares `.CUL/.ECO` por GMR e por ponto, tabelas, notebook para Colab, atlas estático (`produtos/`); demonstração semeadura × GMR em 5 pontos (`resultados/planejamento/`) |
| 10 | Artigos (ASD e MLCSE) | redigidos, em revisão (em `../_nao_sobe/artigo/`, fora do repositório) |

## Dados

- **Clima:** Xavier (BR-DWGD), diário, 0,1°
- **Área de soja:** MapBiomas (classe soja, anual)
- **Produtividade:** IBGE PAM via API do SIDRA, `t/1612 … v/109,216,214,112 … c81/2713`
- **Solos:** BDSolos / SoilGrids (a definir)
- **Literatura:** `dados/literatura/coeficientes_literatura.csv`

## Dados de observação para a calibração

| Arquivo | Conteúdo |
|---|---|
| `dados/brutos/ibge/sidra_1612_soja_municipios_1980_2024.csv` | IBGE PAM, soja por município: área plantada e colhida, produção e rendimento, 1980–2024 (`src/baixar_ibge_municipios.py`) |
| `dados/brutos/xavier/xavier_no_*_1980_2025.csv.gz` | Clima do Xavier nos 462 pontos, 1980–2025 (`src/baixar_xavier_extensao.py`) |
| `dados/literatura/fenologia_observada.csv` | 45 registros de fenologia observada (dias até R1, R5, R7 e R8): Oliveira (2010, Jaboticabal-SP), Bubans et al. (2021, Ijuí-RS), Sampaio et al. (2021, Paragominas-PA, `PARA1801.SBA`) e a média de Battisti et al. (2017, Fig. 1) |
| `dados/literatura/catalogos_embrapa_ciclo.csv` | Ciclo (dias), GMR e hábito por cultivar, região e altitude, de 13 catálogos de cultivares da Embrapa (2010–2019); 88 registros com ciclo (`src/catalogos_embrapa.py`). Os catálogos de 2015 e 2017 têm o ciclo em imagem (sem OCR) |
| `resultados/fenologia_satelite/fenologia_modis_pontos.csv` | Fenologia por satélite (MODIS NDVI 250 m sobre soja do MapBiomas): início do verde, pico e queda de 50% e 20%, mediana e quartis de ~150 pixels por ponto e safra, 2000/01–2023/24 (`src/fenologia_satelite.py`) |

## Literatura

PDFs em `dados/literatura/pdfs/`. Tabelas extraídas:

| Arquivo | Conteúdo |
|---|---|
| `coeficientes_literatura.csv` | 21 conjuntos de coeficientes de 7 trabalhos (colunas `MG`, `mg_origem`, `habito`, `tipo`, `obs`) |
| `fensterseifer2017_calibracoes.csv` | 21 calibrações independentes da BMX Potência RR + agregações (variabilidade entre calibrações) |
| `salmeron2016_genericos.csv` | Coeficientes genéricos por rMG (EUA, MG 3–6) |
| `kaster_farias2011_gmr_ciclo_medio.csv` | Faixa de GMR do ciclo médio por macrorregião sojícola |

Já extraídos: Sampaio et al. (2021), Teixeira et al. (2019), Battisti et al. (2017), Battisti & Sentelhas (2019), Fensterseifer et al. (2017), Mota et al. (2024), Dallacort et al. (2006). Também lidos: Salmerón & Purcell (2016), Kaster & Farias (2011). Tagliapietra et al. (2021) não publica os coeficientes; usa Mercau et al. (2007) e Monzon et al. (2007).

Datas corrigidas em relação ao `ini_infos.md`: Sampaio é de 2021, Dallacort de 2006 (Acta Sci. Agron. 28:583–589) e Kaster & Farias de 2011.

Pedotransferência: Ottoni et al. (2018, HYBRAS; `ottoni2018_hybras.pdf`) indica Medrado & Lima (2014, Geoderma Regional 1:59–66, doi:10.1016/j.geodrs.2014.08.003) como a melhor função para solos brasileiros (RMSE 0,064 cm³/cm³). Tomasella et al. (2000/2003) não é aplicável, porque exige equivalente de umidade. Alternativa só com textura: Tomasella & Hodnett (1998, Soil Sci.).

Faltam: Irene et al. (2024), Rodrigues et al. (2012), Barbieri et al. (2019), Dapper et al. (2026), Battisti (2016, tese ESALQ), Battisti & Sentelhas (2017), Mercau et al. (2007), Monzon et al. (2007).

Cuidados ao usar o CSV:
- **Teixeira (2019) e Battisti & Sentelhas (2019)** reutilizam os coeficientes de crescimento da BRS 284 (Battisti et al., 2017). Não são observações independentes.
- **CSDL/PPSEN** de Dallacort são o MG 5 padrão, não calibrados. Mota e Fensterseifer não informam esses dois coeficientes.
- **MG de Mota (VII) e Fensterseifer (VI)** é a classe; o valor numérico é o centro da classe (`mg_origem`).

## Observações

- As cultivares `RB0002`/`RB0003` do `.CUL` apontam para os ecótipos `SB0777`/`SB0888`, que não existem no `SBGRO048.ECO` do DSSATTools 3.0.2. O artigo (Sampaio et al., 2021) usa os ecótipos padrão SB0701/SB0801.
- Os ecótipos genéricos de MG ≥ 5 são de crescimento determinado (FL-VS = 9). A maioria das cultivares brasileiras atuais é indeterminada, por isso `FL-VS` e `R1PPO` entram na análise.
- Os valores brasileiros de CSDL e FL-SD ficam abaixo do `MINIMA` do DSSAT: as faixas do GENCALC não servem como estão.

## Replicação de Battisti & Sentelhas (2019)

| Arquivo | Função |
|---|---|
| `src/digitalizar_pontos_bs2019.py` | Lê os 462 pontos de clima da Fig. 1 do PDF (o artigo informa 453) |
| `src/digitalizar_zonas_bs2019.py` | Atribui a cada ponto a zona da Fig. 3 pela cor |
| `src/gee.py` | Xavier BR-DWGD (média dos 4 pixels em volta do nó; borda com média focal), soja e solo do MapBiomas, SRTM |
| `src/baixar_xavier_bs2019.py` | Baixa o clima de 1980–2013 dos 462 pontos (cache em `dados/brutos/xavier/`) |
| `src/bs2019.py` | Solos (Tabela S1), cultivares (Tabela S2), manejo e função `simular` |
| `src/rodar_bs2019.py` | Fatorial 462 × 3 solos × 6 semeaduras × (potencial, atingível), em paralelo, com checkpoint |
| `src/analise_bs2019.py` | Estatísticas, agrupamento de Ward, casamento de zonas e comparação com `dados/literatura/bs2019_alvos.csv` |

Suposições (o artigo não informa): ET por Priestley-Taylor, ecótipo do MG inteiro mais próximo, carbono orgânico, albedo, profundidade de semeadura e ponto de orvalho calculado da UR média.

Diferença interna do artigo: pela Tabela S1, o solo argiloso tem 150 mm de água disponível em 0–80 cm; o texto diz 120 mm. A v1 (FAO-56, 150 mm) superestimou a produtividade atingível das zonas secas em 13–45%. A v2 (Priestley-Taylor, 120 mm) ficou dentro de ±12% no teste de sensibilidade e é a configuração adotada.

## Ambiente

```bash
/opt/anaconda3/envs/sklearning/bin/python -m venv .venv    # Python 3.11
.venv/bin/pip install -r requirements.txt
.venv/bin/python -m ipykernel install --user --name calibracao_mg
```

**DSSAT para macOS.** O DSSATTools traz só o executável Linux. O executável nativo foi compilado assim:

```bash
/opt/anaconda3/bin/conda create -y -p ./.tools -c conda-forge gfortran make
git clone https://github.com/DSSAT/dssat-csm-os.git .tools/src_dssat
cd .tools/src_dssat && git checkout v4.8.2.0     # versão cujos arquivos de dados são idênticos aos do DSSATTools 3.0.2
sed -i '' 's/MACOSX_DEPLOYMENT_TARGET} 10.10/MACOSX_DEPLOYMENT_TARGET} 10.15/' CMakeLists.txt
mkdir build && cd build
PATH=../../bin:$PATH SDKROOT=$(xcrun --show-sdk-path) MACOSX_DEPLOYMENT_TARGET=10.15 \
  cmake .. -DCMAKE_Fortran_COMPILER=$(which gfortran) -DCMAKE_BUILD_TYPE=RELEASE && make -j4
```

O executável `build/bin/dscsm048` substitui `.venv/.../DSSATTools/bin/dscsm048` (o Linux original fica em `dscsm048.linux`). O DSSAT não aceita caminhos longos: as rodadas usam `.run/`.

**GEE:** projeto `mapbiomas-brazil`.

Para a MCMC da fenologia, `src/calibrar_fenologia.py` aceita `--output` para salvar uma cadeia em outro arquivo, `--seed` para controlar a inicialização, `--initial-chain` para iniciar pelas posições finais de outra cadeia, `--move` para escolher o movimento e `--resume` para continuar uma cadeia existente. Sem `--resume`, o script recusa sobrescrever o arquivo. Para analisar uma cadeia específica, use `src/analisar_calibracao_fenologia.py --chain caminho.h5`.

## Estrutura

```
src/genotipo.py    leitura do SBGRO048.CUL/.ECO
src/faixas.py      faixas dos coeficientes
notebooks/         uma etapa por notebook (01_, 02_, …)
dados/brutos       downloads (Xavier, MapBiomas, SIDRA)
dados/processados  WTH, SOL, tabelas prontas
dados/literatura   coeficientes e observações compilados de artigos
resultados/        figuras e tabelas
```
