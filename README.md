# Seasonality of tuberculosis notification in India: disease signal versus public-sector reporting artefact

Data and analysis code for the paper:

> **Disease or calendar? Seasonality of tuberculosis notification in India and a public-sector reporting artefact, Ni-kshay 2022–2025.**
> Arun Mitra et al. (under review).

A longitudinal ecological analysis of monthly TB notification at the district level (747 Local
Government Directory districts), separating the genuine seasonal signal from a reporting-calendar
and campaign artefact, using the public-versus-private sector contrast as the discriminating test.

## What the analysis does

1. **Objective 1 — describe seasonality.** Seasonal-trend decomposition (STL) of monthly
   notification, nationally and by sector, with moving-block bootstrap confidence intervals.
2. **Objective 2 — signal vs artefact.** Sector-specific negative-binomial models with population
   and working-day offsets, harmonic and lagged-climate terms, and reporting-calendar/campaign
   terms; the public-vs-private contrast identifies the artefact. The modelled seasonal variance is
   partitioned into harmonic, climate, and calendar-campaign blocks.
3. **Objective 3 — monitoring consequence and heterogeneity.** The share of month-on-month
   movements that reverse after seasonal and working-day adjustment, and a beta regression of where
   the artefact concentrates (with a Census-2011 socio-demographic robustness check and a spatial
   autocorrelation test).

## Repository structure

```
.
├── code/     R analysis scripts (sequential; see the map below)
├── data/     district-level aggregates, covariates, and cached model objects (RDS)
├── files/    district geometry (LGD .gpkg) and lookups
├── out/      result tables (CSV) and figures (PNG) behind the paper
├── paper/    manuscript and supplement (Quarto) + bibliographies
└── README.md
```

### Code → output map

| Script | Does | Reads | Writes | Runs from quick-start data? |
|---|---|---|---|---|
| `00_setup.R` | packages, `here()` paths, seed | — | — | always sourced first |
| `_install_packages.R` | one-time dependency install | — | — | run once |
| `01_build_district_panel.R` | TU→district crosswalk, fold newly-created districts to parents, build the district-month panel | TU-level Ni-kshay (**non-public**) | `data/panel_dm.rds`, `files/district_anl.gpkg`, lookups | **No — provenance only** |
| `02_covariates.R` | NASA POWER climate, working-day/holiday/quarter-end/campaign calendar, district modifiers | panel, geometry, POWER API, TU-level (**non-public** for the modifier part) | `data/climate_dm.rds`, `data/calendar.rds`, `data/campaign_dm.rds`, `data/modifiers*.rds` | **Partly — outputs shipped** |
| `03_seasonality_stl.R` | Objective 1 STL by sector + bootstrap | `data/panel_dm.rds` | `out/o1_seasonal_national.csv`, `out/o1_seasonal_plot.png` | **Yes** |
| `04_sector_models.R` | assemble model frame, fit sector NB models, variance partition | `data/panel_dm.rds` + covariates | `data/model_frame.rds`, `data/fit_*.rds`, `out/o2_*` | **Yes** (cached fits shipped) |
| `05_monitoring_and_heterogeneity.R` | payoff (sign-flip / surge artefact), heterogeneity beta regression, Census-adjusted check, Moran's I | panel, calendar, modifiers, geometry | `out/o3_payoff*`, `out/o3_heterogeneity*` | **Yes** |
| `06_figures_and_tables.R` | assemble manuscript figures/tables | panel, fits | `out/` figures + coefficient table | **Yes** |
| `07_sensitivity.R` | STL spec, climate lag, working-day handling, payoff method | panel, model frame, fits | `out/sensitivity.csv` | **Yes** (cached sensitivity fits shipped) |

## Reproducing the analysis

Requires **R ≥ 4.4.1**. Paths resolve relative to the repository root via `here()`, so the project
runs unchanged after `git clone`.

```r
# 1. install dependencies (once)
source("code/_install_packages.R")

# 2. quick start — reproduce all models and figures from the shipped district data
source("code/03_seasonality_stl.R")
source("code/04_sector_models.R")
source("code/05_monitoring_and_heterogeneity.R")
source("code/06_figures_and_tables.R")
source("code/07_sensitivity.R")
```

Each script sources `code/00_setup.R` itself. The model fits and the panel are shipped in `data/`,
so scripts `03`–`07` reproduce every result table in `out/` (and every figure in the paper) without
re-running the multi-hour model fits or the climate download.

### Full rebuild from source (optional)

Scripts `01` and the climate-API / TU-level parts of `02` rebuild the inputs from TU-level Ni-kshay
data, which is **not redistributable** (see below). They are provided for transparency. To run them,
set the environment variable `TU_SOURCE_ROOT` to your local source tree before sourcing
`00_setup.R`; otherwise script `01` stops with a message directing you to the shipped data.

### Rendering the paper

```bash
cd paper
QUARTO_R="/path/to/R-4.4.1" quarto render manuscript_seasonality.qmd
QUARTO_R="/path/to/R-4.4.1" quarto render supplement_seasonality.qmd
```

Both read their figures from `../out/` and their tables from `../out/` and `../files/`.

## Dependencies

Core packages (see `code/_install_packages.R` for the full install): `here`, `dplyr`, `tidyr`,
`readr`, `stringr`, `tibble`, `readxl`, `scales`; `sf`, `spdep`, `spatialreg`; `tsibble`,
`lubridate`, `feasts`, `fabletools`, `fable`; `glmmTMB`, `DHARMa`, `MuMIn`, `tscount`, `betareg`;
`duckdb`, `DBI`; `nasapower`; `ggplot2`. Rendering the paper additionally needs Quarto and `knitr`.

## Data availability

- **TB notification (Ni-kshay).** Facility-level (TB Unit) notification counts are governed by the
  National TB Elimination Programme and are **not publicly redistributable**. This repository ships
  the **district-level monthly aggregates** (`data/panel_dm.rds`) that underlie every result; no
  facility-level rows are included.
- **Climate.** NASA POWER monthly data, publicly available (https://power.larc.nasa.gov/); the
  district-month extract is shipped in `data/climate_dm.rds`.
- **Population.** WorldPop, publicly available (https://www.worldpop.org/); district denominators in
  `data/district_pop_anl.rds`.
- **Geometry.** `files/district_anl.gpkg` (747 districts) is a derived layer: the 2011 Local
  Government Directory (LGD) sub-district boundaries dissolved to district level, with 18
  newly-created districts folded to their pre-split parents (`01_build_district_panel.R`). The
  shipped copy is **topology-preserving simplified** (mapshaper, ~5% of vertices, ~4.6 MB instead of
  73 MB); this is for display and the queen-contiguity graph only and leaves the contiguity
  structure used for Moran's I unchanged (4,078 neighbour links, identical to the full-resolution
  layer). The raw LGD boundaries are from the Local Government Directory
  (https://lgdirectory.gov.in/); the full-resolution district layer is regenerated by
  `01_build_district_panel.R`.
- **Socio-demography.** 2011 Census district indicators (`data/pc11_pca_district.tab`).

## Ethics

The study used aggregate programmatic data with no individual identifiers. Institutional Ethics
Committee approval number and date: *to be inserted*.

## Citation

If you use this code or the district-level data, please cite the paper (citation and archived DOI to
be added on publication).

```
Mitra A, et al. Disease or calendar? Seasonality of tuberculosis notification in India and a
public-sector reporting artefact, Ni-kshay 2022–2025. [Journal], [Year]. DOI: [to be added]
```

## License

Code released under the MIT License (see `LICENSE`). The district-level data are released for
research and replication of this analysis.
