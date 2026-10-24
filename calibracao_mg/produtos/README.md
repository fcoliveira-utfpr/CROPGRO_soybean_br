# Brazilian soybean cultivars for CROPGRO-Soybean (DSSAT)

Data products of the national calibration of CROPGRO-Soybean phenology as functions of the relative maturity group
(RMG), and of the effective RMG and sowing date inferred from MODIS phenology at 224 grid points of Brazil
(harvests 2001–2024).

## Contents

| Path | Content |
|---|---|
| `dados/theta_rmg.csv` | Phenological coefficients (CSDL, PPSEN, EM-FL, FL-SD, SD-PM, R1PPO) from RMG 5.0 to 9.5 (step 0.1): posterior median and 90% interval |
| `dados/coefficients_point.csv` | One row per grid point: effective RMG (mean, median, median of harvests 2021–2024, trend per decade), sowing date and its trend, water-limited/potential yield, soil texture fractions of the soybean area, DSSAT cultivar code and its coefficients |
| `dados/coefficients_point_season.csv` | One row per point and season: effective RMG (mean, SD, 5/50/95% quantiles), sowing date, coefficients at the median RMG, simulated potential and water-limited yield, IBGE yield of the 1° cell |
| `dados/SBGRO048_BR.CUL`, `dados/SBGRO048_BR.ECO` | DSSAT 4.8 cultivar and ecotype files: `BRG050`–`BRG095` (generic cultivars for RMG 5.0–9.5) and `BR0001`–`BR0224` (one cultivar per grid point) |
| `notebooks/usage_example.ipynb` | Runs DSSAT with these cultivars and NASA/POWER weather for any location, including a sowing date × maturity group comparison (runs in Google Colab) |
| `pagina/index.html` | Static atlas: map of the effective RMG, sowing date and water limitation; series and cultivar of each point; coefficient curves |

## The notebook: no programming or installation needed

`notebooks/usage_example.ipynb` was written for agronomists, consultants, students and extension agents with little or
no programming experience. It runs in the browser with Google Colab, so nothing needs to be installed: open it in Colab,
change the location, seasons and sowing dates in the cell *Your settings*, and choose *Runtime → Run all*. The notebook
downloads the weather (NASA/POWER), builds the cultivar and the soil, runs DSSAT and draws the graphs. Its purpose is to
support planning: comparing sowing dates and maturity groups for a location by harvest date (time left for a second
crop) and by water risk.

## How the cultivars are built

* Phenology: six coefficients as linear functions of the RMG, estimated by Bayesian calibration (emulator of DSSAT,
  hierarchical likelihood, two ensembles of 40 walkers and 9,000 steps) against field experiments, multi-site cultivar
  trials and cultivar catalogues. R1PPO is an ecotype coefficient and is written in the `.ECO` file.
* Growth: coefficients of BRS 284 (Battisti and Sentelhas, 2019, Field Crops Research 240:95–105, Table S2); FL-SH and
  FL-LF of the nearest maturity group of that study; ecotype FL-VS = 26.
* Base generic cultivar and ecotype: nearest integer maturity group (990005–990009, SB0501–SB0901).
* Point cultivars use the median effective RMG of the harvests 2021–2024 (the median is robust to seasons with an
  anomalous NDVI curve, such as the drought of 2021/22 in Paraná).

The DSSAT outputs obtained with `SBGRO048_BR.CUL`/`.ECO` were verified to be identical to those of the cultivars built
in the calibration pipeline (`src/validar_cul.py`).

## Using the DSSAT files

Copy both files to the folder of the experiment, or to `DSSAT48/Genotype` renamed as `SBGRO048.CUL` and `SBGRO048.ECO`
(back up the originals first, or append the lines to them), and use the cultivar code in the FileX.

## Limitations

* The cultivars are effective cultivars: each represents the mix of cultivars grown in a 1° cell.
* Growth coefficients are those of a single reference cultivar and do not include the genetic gain of recent cultivars;
  use simulated yields for relative comparisons (sowing dates, maturity groups, seasons, regions).
* In southern Brazil the effective RMG (about 7.0) is higher than that of recommended cultivars (5.5–6.4), because the
  model underestimates the cycle of early sowings there and the satellite lag between maturity and NDVI decline is
  partly confounded with the RMG.
* Sowing dates inferred from the satellite have not yet been validated against field records.

## Regenerating

From `calibracao_mg/`: `DELTA8_FIXO=11 PRIOR_LIT=1 FAIXA_FISIO=1 .venv/bin/python src/gerar_produtos.py`, then
`src/validar_cul.py`, `src/gerar_notebook_uso.py` and `src/gerar_pagina.py`.

## License

Data and products: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) (see `LICENSE-DATA` at the root of the repository).
Please cite the authors when using the cultivars, tables or the atlas.
