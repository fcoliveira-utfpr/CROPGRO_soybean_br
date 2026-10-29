# Nationwide maturity-group calibration of CROPGRO-Soybean and the effective maturity group of soybean in Brazil

Analyses, data and products of a calibration of the CSM-CROPGRO-Soybean (DSSAT 4.8) coefficients as **functions of the
relative maturity group (RMG)** for the whole of Brazil, of the **inversion of satellite phenology (MODIS)** into the
effective RMG and sowing date of soybean from 2000 to 2024, and of the **water limitation of yield** simulated with this
information. This repository is the basis of two scientific papers (see [Papers](#papers)).

[![Atlas](https://img.shields.io/badge/atlas-GitHub%20Pages-2f6b3c)](https://fcoliveira-utfpr.github.io/CROPGRO_soybean_br/)
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/fcoliveira-utfpr/CROPGRO_soybean_br/blob/main/calibracao_mg/produtos/notebooks/usage_example.ipynb)
[![Data license: CC BY 4.0](https://img.shields.io/badge/data-CC%20BY%204.0-lightgrey)](LICENSE-DATA)

**Quick links**

- Interactive atlas of the effective RMG, sowing dates and cultivars by point: <https://fcoliveira-utfpr.github.io/CROPGRO_soybean_br/>
- Notebook that runs DSSAT with the Brazilian cultivars in the browser, with nothing to install: [open in Google Colab](https://colab.research.google.com/github/fcoliveira-utfpr/CROPGRO_soybean_br/blob/main/calibracao_mg/produtos/notebooks/usage_example.ipynb)
- DSSAT cultivar and ecotype files: [`SBGRO048_BR.CUL`](calibracao_mg/produtos/dados/SBGRO048_BR.CUL), [`SBGRO048_BR.ECO`](calibracao_mg/produtos/dados/SBGRO048_BR.ECO)
- Tables of coefficients, effective RMG and sowing dates: [`calibracao_mg/produtos/dados/`](calibracao_mg/produtos/dados/)

Folder and file names, code comments and the technical notes in `calibracao_mg/README.md` are in Portuguese; this README
and the products in `calibracao_mg/produtos/` are in English.

---

## Contents

- [Papers](#papers)
- [Main results](#main-results)
- [Products](#products)
- [Repository structure](#repository-structure)
- [Workflow](#workflow)
- [How to reproduce](#how-to-reproduce)
- [Data sources](#data-sources)
- [What is not in the repository](#what-is-not-in-the-repository)
- [Internal documentation](#internal-documentation)
- [Citation, license and contact](#citation-license-and-contact)

---

## Papers

The project supports two independent manuscripts in preparation for Springer journals. The manuscripts are not in the
repository until they are published.

| | Agronomic paper | Computational paper |
|---|---|---|
| Working title | *Maturity-group-based phenology of CROPGRO-Soybean for Brazil and the effective maturity group of soybean from 2000 to 2024* | *Emulator-based Bayesian calibration and amortized inference of soybean phenology across Brazil with a process-based crop model* |
| Target journal | Agronomy for Sustainable Development | Machine Learning for Computational Science and Engineering |
| Focus | National θ(RMG) coefficients; effective RMG and sowing date from 2000 to 2024; water limitation of yield (Yw/Yp) validated against IBGE yields; point-specific cultivars for planning | Emulator validated on held-out locations; hierarchical likelihood; identifiability across 13 versions of the posterior; physiological constraint and its consequence for yield; amortized inference (neural posterior estimation) compared with per-record MCMC |
| Parts of the repository | steps 1–4 and 7–10 | steps 5–8 |

---

## Main results

**Calibrated phenology (final version v13).** Six coefficients as linear functions of the RMG (CSDL, PPSEN, EM-FL, FL-SD,
SD-PM and R1PPO), calibrated with 746 likelihood terms: 453 observations from 8 field studies at 13 sites, multi-site
means of 88 cultivars (Alliprandini et al., 2009) and 117 cultivar-catalogue records. The model reproduces the multi-site
days to R1 and R8 with RMSE of 3.2 and 1.9 days (Willmott's d of 0.93 and 0.99).

| Coefficient | RMG 5 | RMG 7 | RMG 9 |
|---|---|---|---|
| CSDL (h) | 12.88 | 12.39 | 11.90 |
| PPSEN (1/h) | 0.228 | 0.253 | 0.279 |
| EM-FL (photothermal days) | 21.1 | 28.4 | 35.8 |
| FL-SD (photothermal days) | 13.9 | 10.8 | 10.5 |
| SD-PM (photothermal days) | 25.9 | 25.1 | 25.0 |
| R1PPO (h) | 0.95 | 0.86 | 0.77 |

**Emulator.** Gradient boosting trained on 2.81 million simulated seasons at 73 locations; on held-out locations, RMSE of
1.2 days for R1 and 2.0 days for R7, about 3×10⁴ times faster than DSSAT.

**Comparison with GLUE.** Reproduced with the same emulator and likelihood, the GLUE procedure of DSSAT degenerated: with
6,000 parameter sets (the DSSAT default) one set carried all the weight (effective sample size 1.0), and with 60,000 sets
the best log-posterior (−1797) remained far from that of the ensemble sampler (−1579).

**Effective RMG and sowing date (MODIS, 4,030 point-seasons at 224 points).** The amortized inversion recovers the RMG
with an error of 0.75 and 88% coverage of the 90% credible interval, with posteriors nearly identical to those of
per-record MCMC at about 600 times lower cost. The effective RMG decreased by 0.25 units per decade between 2000 and 2024
(from 9.0 to 8.2 in Mato Grosso and from 8.6 to 7.7 in the Centre-West), and sowing moved 3.9 days per decade earlier.
The spatial pattern agrees with the RMG zoning of Ferreira et al. (2026): 0.136 RMG units per degree of latitude, against
0.14 in that study.

**Water limitation.** Yw/Yp of 0.74 for Brazil, from 0.92 in Mato Grosso, Rondônia and the North to 0.53 in the South;
the detrended anomalies of the simulated Yw follow those of IBGE yields (r = 0.76 by region and season).

**Default × calibrated cultivars (Paragominas, PA).** In an experiment with time series of leaf area, biomass and grain,
the national cultivars reduced the errors of the default DSSAT cultivars by about half (RMSE of flowering 3.9 against
7.4 days, of maturity 4.4 against 7.8 days, of grain yield 953 against 1,395 kg/ha), without any local adjustment.

**Planning example.** At five representative points, each RMG unit delays harvest by 11–14 days; in Mato Grosso and Goiás
the water risk is low for sowings from October, in MATOPIBA early sowings are risky, and in Paraná and the South no
combination of sowing date and RMG removes the water risk.

**Known limitations.** The coefficients are effective (the partition between R1PPO and the reproductive durations is
weakly identifiable and was constrained by the published range); growth coefficients are those of BRS 284 (Battisti and
Sentelhas, 2019), without the recent genetic gain, so simulated yields are meant for relative comparisons; the model
underestimates the cycle of early sowings in the South; and the effective RMG of the South is higher than that of
recommended cultivars.

---

## Products

In [`calibracao_mg/produtos/`](calibracao_mg/produtos/) (details in the [products README](calibracao_mg/produtos/README.md)):

| Product | Purpose |
|---|---|
| [`dados/SBGRO048_BR.CUL`](calibracao_mg/produtos/dados/SBGRO048_BR.CUL) and [`.ECO`](calibracao_mg/produtos/dados/SBGRO048_BR.ECO) | Brazilian cultivars for DSSAT 4.8: generic by RMG (`BRG050`–`BRG095`) and one per grid point (`BR0001`–`BR0224`). DSSAT outputs are identical to those of the pipeline (verified by `src/validar_cul.py`) |
| [`dados/*.csv`](calibracao_mg/produtos/dados/) | Coefficients by RMG (with 90% interval), by point and by point-season; effective RMG, sowing date, Yp, Yw and IBGE yield |
| [`notebooks/usage_example.ipynb`](calibracao_mg/produtos/notebooks/usage_example.ipynb) | Runs DSSAT with these cultivars for any location. **Made for people with little or no programming experience: it runs in the browser through [Google Colab](https://colab.research.google.com/github/fcoliveira-utfpr/CROPGRO_soybean_br/blob/main/calibracao_mg/produtos/notebooks/usage_example.ipynb), with nothing to install.** Set the location, seasons and sowing dates in one cell and run all; the notebook downloads NASA/POWER weather, builds the cultivar and the soil, runs DSSAT and plots the harvest date and the water limitation for each sowing date and maturity group, to support planning |
| [Atlas](https://fcoliveira-utfpr.github.io/CROPGRO_soybean_br/) ([`pagina/index.html`](calibracao_mg/produtos/pagina/index.html)) | Static atlas: map of the effective RMG, its trend, the sowing date and the water limitation; series and cultivar of each point; θ(RMG) curves; downloads |

The atlas is served by GitHub Pages from the [`docs/`](docs/) folder (branch `main`, folder `/docs`), a copy of
`calibracao_mg/produtos/pagina/` rebuilt by `src/gerar_pagina.py`; the `.nojekyll` file makes GitHub serve the files as
they are. When the notebook runs in Colab, it reads the data products directly from this repository.

---

## Repository structure

```
.
├── README.md                     this file
├── LICENSE-DATA                  CC BY 4.0 license of the data and products
├── .gitignore                    what is kept out of the repository (see below)
├── docs/                         atlas for GitHub Pages (copy of calibracao_mg/produtos/pagina/)
└── calibracao_mg/
    ├── README.md                 project steps and technical notes (Portuguese)
    ├── requirements.txt          main environment (.venv)
    ├── requirements_sbi.txt      environment of the amortized inference (.venv_sbi)
    ├── src/                      all the code (Python scripts and shell pipelines)
    ├── notebooks/                01–05 by step; 06_figuras_resultados.ipynb builds the figures of the papers
    ├── dados/                    data
    │   ├── brutos/               raw data (IBGE municipal yields 1980–2024, field datasets, DSSAT example files)
    │   ├── processados/          grid points, MapBiomas soils, municipalities, state boundaries
    │   └── literatura/           tables extracted from the literature (coefficients, observed phenology, catalogues)
    ├── resultados/               results
    │   ├── bs2019/               replication of Battisti and Sentelhas (2019): comparisons and figures
    │   ├── solo_mapbiomas/       test of soil profiles
    │   ├── morris/               sensitivity analysis
    │   ├── emulador_fenologia/   emulator models, validation, observations, comparison of chains and emulators
    │   ├── fenologia_satelite/   MODIS phenology by point and season
    │   ├── sbi/                  trained neural posterior estimator, validation, effective RMG by point-season
    │   ├── produtividade/        Yp, Yw, IBGE yields by cell, water limitation and anomalies
    │   ├── planejamento/         sowing date × RMG demonstration at 5 points
    │   ├── curvas/               default × calibrated cultivars against observed curves (Paragominas)
    │   ├── glue/                 GLUE reproduced with the emulator, compared with the ensemble sampler
    │   ├── logs/                 logs of the runs (calibrations, training, pipelines)
    │   ├── figuras/              figures in Portuguese (PNG and PDF)
    │   └── figuras_en/           the same figures in English
    └── produtos/                 cultivars, tables, usage notebook and atlas (English)
```

---

## Workflow

The steps follow the flow of the data. All commands are run from `calibracao_mg/`.

| # | Step | Main scripts | Outputs |
|---|---|---|---|
| 1 | Data: grid of 462 points of Battisti and Sentelhas (2019), daily BR-DWGD weather (Xavier) 1980–2025, MapBiomas soybean and soil, IBGE | `digitalizar_pontos_bs2019.py`, `digitalizar_zonas_bs2019.py`, `gee.py`, `baixar_xavier_bs2019.py`, `baixar_xavier_extensao.py`, `solo_mapbiomas.py`, `baixar_ibge_municipios.py` | `dados/processados/`, `dados/brutos/` |
| 2 | Replication of Battisti and Sentelhas (2019) and choice of the setup (Priestley–Taylor evapotranspiration; 3 soils weighted by texture) | `bs2019.py`, `rodar_bs2019.py`, `analise_bs2019.py`, `sensibilidade_bs2019.py`, `rodar_solo_mb.py`, `ptf.py` | `resultados/bs2019/`, `resultados/solo_mapbiomas/`; notebooks 02 and 03 |
| 3 | Morris sensitivity analysis (15 coefficients, 9,600 simulations) | `faixas.py`, `morris.py` | `resultados/morris/`; notebook 04 |
| 4 | Phenological observations with known RMG (experiments, Alliprandini, catalogues) | `catalogos_embrapa.py`, `observacoes_fenologia.py` | `dados/literatura/`, `resultados/emulador_fenologia/observacoes_fenologia.csv` |
| 5 | Phenology emulator (2.81 million seasons at 73 locations; validation on held-out locations) | `treino_fenologia.py`, `emulador_fenologia.py`, `validar_emulador_fenologia.py`, `comparar_emuladores.py`, `sensibilidade_emulador.py` | `resultados/emulador_fenologia/modelo_*.joblib`, `validacao_metricas.csv`, `comparacao_emuladores.csv` |
| 6 | Hierarchical Bayesian calibration (13 versions; final v13) | `calibrar_fenologia.py`, `comparar_cadeias.py`, `decompor_logvero.py`, `pipeline_fenologia_v*.sh`, `glue_comparacao.py` | chains `cadeia_v*.h5` (not in git), `comparacao_cadeias_v*.txt` |
| 7 | MODIS phenology and amortized inversion (neural posterior estimation) of the effective RMG and sowing date; comparison with per-record MCMC | `fenologia_satelite.py`, `rodar_fenologia_satelite.py`, `sbi_simulador.py`, `sbi_treinar.py`, `pipeline_sbi_final.sh`, `mcmc_por_registro.py` | `resultados/fenologia_satelite/`, `resultados/sbi/` |
| 8 | Potential and water-limited yield by point and season; water limitation and validation against IBGE | `rodar_produtividade.py`, `lacuna_produtividade.py`, `limitacao_hidrica.py` | `resultados/produtividade/` |
| 9 | Products, planning demonstration and default × calibrated curves (Paragominas) | `gerar_produtos.py`, `validar_cul.py`, `gerar_notebook_uso.py`, `gerar_pagina.py` (+ `pagina_modelo.html`), `demo_planejamento.py`, `curvas_paragominas.py`, `figuras_curvas.py` | `produtos/`, `docs/`, `resultados/planejamento/`, `resultados/curvas/` |
| 10 | Figures of the papers (Portuguese and English) | `figuras.py`, `gerar_notebook_figuras.py`, `traducao_en.json` | `notebooks/06_figuras_resultados.ipynb`, `resultados/figuras*/` |

The scripts `emulador.py`, `treino_emulador.py`, `validar_emulador.py`, `validar_emulador_fen.py` and notebook 05 belong
to an earlier yield emulator, replaced by the phenology emulator; they are kept as a record.

**Environment variables of the final calibration.** Every script that reads the v13 posterior must run with
`DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1` (R7–R8 interval fixed at 11 days, literature prior, and physiological range of
FL-SD and SD-PM).

---

## How to reproduce

**1. Environments (Python 3.11).**

```bash
cd calibracao_mg
python3.11 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python -m ipykernel install --user --name calibracao_mg
python3.11 -m venv .venv_sbi && .venv_sbi/bin/pip install --only-binary=grpcio -r requirements_sbi.txt
```

**2. DSSAT.** DSSATTools ships DSSAT executables for Linux and Windows. On macOS, compile DSSAT 4.8.2 (instructions in
[`calibracao_mg/README.md`](calibracao_mg/README.md), section *Ambiente*) and replace `.venv/.../DSSATTools/bin/dscsm048`.
DSSAT does not accept long paths, so the runs use the short folder `calibracao_mg/.run/`.

**3. Google Earth Engine** (only to download the weather and MapBiomas data again): `earthengine authenticate`, project
`mapbiomas-brazil` in `src/gee.py`.

**4. Order of the commands** (approximate times on a 2-core laptop):

```bash
# weather, soil and IBGE (step 1)
.venv/bin/python src/baixar_xavier_extensao.py && .venv/bin/python src/baixar_ibge_municipios.py
# emulator training and validation (step 5; ~2 days of simulation in batches)
for lote in 1 2 3 4 5 6; do .venv/bin/python src/treino_fenologia.py <samples_per_location> 4 $lote; done   # see the header of src/treino_fenologia.py
.venv/bin/python src/validar_emulador_fenologia.py
# MODIS phenology (step 7, Earth Engine)
.venv/bin/python src/rodar_fenologia_satelite.py
# final calibration (step 6; ~5 h per version)
src/pipeline_fenologia_v13.sh 3000          # then extended to 9,000 steps (pipeline_fenologia_v13_ext*.sh)
# satellite inversion (step 7; ~1 h)
src/pipeline_sbi_final.sh 100000 resultados/emulador_fenologia/cadeia_v13_a.h5
# yield and water limitation (step 8; ~5 h)
DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1 .venv/bin/python src/rodar_produtividade.py resultados/emulador_fenologia/cadeia_v13_a.h5 final 4
.venv/bin/python src/lacuna_produtividade.py && .venv/bin/python src/limitacao_hidrica.py
# products (step 9)
DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1 .venv/bin/python src/gerar_produtos.py
.venv/bin/python src/gerar_notebook_uso.py && .venv/bin/python src/gerar_pagina.py
# figures (step 10): Portuguese and English
.venv/bin/python src/gerar_notebook_figuras.py
cd notebooks && ../.venv/bin/jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.kernel_name=calibracao_mg 06_figuras_resultados.ipynb
FIG_LANG=en ../.venv/bin/jupyter nbconvert --to notebook --execute --output-dir /tmp --ExecutePreprocessor.kernel_name=calibracao_mg 06_figuras_resultados.ipynb
```

The exact commands of each calibration version are in the scripts `src/pipeline_fenologia_v*.sh`, whose headers explain
why each version was run; the logs of the runs are in `resultados/logs/`.

---

## Data sources

| Data | Source | Use |
|---|---|---|
| Daily weather (Tmax, Tmin, rainfall, radiation, relative humidity, wind), 0.1°, 1980–2025 | BR-DWGD (Xavier et al., 2022), through Google Earth Engine | simulations and emulator |
| Soybean area and soil texture | MapBiomas (land cover and soil) | satellite mask and soil weights |
| Satellite phenology (NDVI, 250 m) | MODIS MOD13Q1 and MYD13Q1 | effective RMG and sowing date |
| Municipal soybean yield | IBGE, Municipal Agricultural Production, SIDRA table 1612 | validation of the anomalies |
| Soil profiles and reference cultivar (BRS 284) | Battisti and Sentelhas (2019), Tables S1 and S2 | setup and growth |
| Observed phenology and published coefficients | 8 field studies, Alliprandini et al. (2009), Embrapa and other cultivar catalogues, 21 published calibrations | calibration and priors |
| Weather in the usage notebook | NASA/POWER (daily API) | example of use outside the grid |

Full references are given in the papers. The PDFs consulted are not redistributed (copyright).

---

## What is not in the repository

**Kept on the authors' computer, outside git, where the scripts read them** (listed in [`.gitignore`](.gitignore)):

| File or folder | Size | How to obtain |
|---|---|---|
| `.venv/`, `.venv_sbi/`, `.tools/`, `calibracao_mg/.run/` | 2.8 GB | recreate (previous section) |
| `calibracao_mg/codex.md` | — | internal lab notebook of the authors |
| `calibracao_mg/dados/brutos/xavier/` | 287 MB | `src/baixar_xavier_*.py` (Earth Engine) |
| `resultados/emulador_fenologia/tabela_treino.parquet`, `treino/`, `validacao_predicoes.parquet` | ~380 MB | `src/treino_fenologia.py` and `src/validar_emulador_fenologia.py`; to be deposited on Zenodo |
| MCMC chains `resultados/emulador_fenologia/**/*.h5` (all versions; the final ones are `cadeia_v13_a/b.h5`, 94 MB each) | ~1.1 GB | `src/pipeline_fenologia_v*.sh`; **the final chains will be deposited on Zenodo**, because reproducing them takes days |
| Simulations of the satellite inversion (`resultados/sbi/*.npz`) | 37 MB | `src/sbi_simulador.py` |
| Simulations of the replication and of the soil test (`resultados/bs2019/sim*/`, `resultados/solo_mapbiomas/sim_*/`) | 160 MB | `src/rodar_bs2019.py`, `src/rodar_solo_mb.py` |

**Also kept out of the repository:** the manuscripts of the two papers (until publication), the PDFs of the literature and
of the cultivar catalogues (copyright), the initial exploratory code, an earlier yield emulator (669 MB), TensorBoard logs
and leftover DSSAT run files.

---

## Internal documentation

- [`calibracao_mg/README.md`](calibracao_mg/README.md) (Portuguese): project steps, technical decisions (soils,
  pedotransfer, ecotypes), inventory of the literature and instructions to compile DSSAT.
- `src/pipeline_fenologia_v*.sh`: one script per calibration version, with the reason for the version in its header.

---

## Citation, license and contact

**Citation.** Until the papers are published, please cite the repository: Oliveira, F. C. *et al.* (2026). *Nationwide
maturity-group calibration of CROPGRO-Soybean and the effective maturity group of soybean in Brazil* [code and data].
GitHub: <https://github.com/fcoliveira-utfpr/CROPGRO_soybean_br>. DOI to be assigned (Zenodo).

**License.** The data and products (tables, `.CUL`/`.ECO` cultivars, figures, atlas and results) are licensed under
[Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/): they may be
used, adapted and redistributed, including commercially, provided the authors are credited (full text in
[`LICENSE-DATA`](LICENSE-DATA)). Third-party data included keep the licenses of their sources. The license of the code
(`src/` and notebooks) is still to be defined.

**Contact.** Federal University of Technology – Paraná (UTFPR), Santa Helena campus — fcoliveira@utfpr.edu.br

**Funding.** CAPES, Fundação Araucária and UTFPR.
