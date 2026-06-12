# 00_setup.R — packages, paths, and global settings.
# Source this first in every script: source(here::here("code", "00_setup.R"))
#
# Paths are resolved relative to the repository root with here(), so the project
# runs unchanged on any machine after `git clone`. Requires R >= 4.4.1.

stopifnot(getRversion() >= "4.4.1")

## --- packages -------------------------------------------------------------
## One-time install: run code/_install_packages.R (see README for the version list).
suppressMessages({
  library(here)
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(tsibble); library(lubridate)
})

## --- paths ----------------------------------------------------------------
DER   <- here("data")    # district-level inputs + model intermediates (the quick-start data)
OUT   <- here("out")     # analysis result tables and figures
FILES <- here("files")   # district geometry and lookups
GEO   <- file.path(FILES, "district_anl.gpkg")   # 747-district LGD polygons
for (d in c(DER, OUT, FILES)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

## Provenance roots — only needed to REBUILD the panel from the non-public TU-level
## Ni-kshay data (script 01, and the climate-API / TU parts of 02). Leave unset to run
## from the shipped district aggregates. Point TU_SOURCE_ROOT at your local source tree.
ROOT <- Sys.getenv("TU_SOURCE_ROOT", unset = NA_character_)
SQ   <- if (!is.na(ROOT)) file.path(ROOT, "side_quest_subdistrict") else NA_character_
P2   <- here()

## --- global settings ------------------------------------------------------
SEED  <- 1L
B     <- 2000L          # bootstrap resamples
YEARS <- 2022:2026
set.seed(SEED)
