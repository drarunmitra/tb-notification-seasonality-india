# _install_packages.R — one-time dependency install. Run once before the pipeline.
# R >= 4.4.1. MASS and splines ship with R and are not installed here.

pkgs <- c(
  # data wrangling / IO
  "here", "dplyr", "tidyr", "readr", "stringr", "tibble", "readxl", "scales",
  # spatial
  "sf", "spdep", "spatialreg",
  # time series / seasonality
  "tsibble", "lubridate", "feasts", "fabletools", "fable",
  # count / mixed / beta models and diagnostics
  "glmmTMB", "DHARMa", "MuMIn", "tscount", "betareg",
  # database backend for panel aggregation (provenance step only)
  "duckdb", "DBI",
  # climate extraction (provenance step only)
  "nasapower",
  # plotting
  "ggplot2"
)

to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
} else {
  message("All required packages are already installed.")
}
