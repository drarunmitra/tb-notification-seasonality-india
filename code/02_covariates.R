# ----------------------------------------------------------------------------
# PROVENANCE STEP. This script rebuilds inputs from TU-level Ni-kshay data that
# is NOT redistributable. The shipped files in data/ and files/ are the
# quick-start entry point; you only need this script to rebuild from source.
# To run it, set the environment variable TU_SOURCE_ROOT to your local source
# tree (the parent of side_quest_subdistrict/) before sourcing 00_setup.R.
# ----------------------------------------------------------------------------

# 02_covariates.R — Climate, calendar/campaign, and district modifier covariates
#
# Renumbered, path-normalized reorganization of the working pipeline scripts: 02_climate.R, 05_calendar.R, 06_modifiers.R.
# Computation is unchanged from the scripts that produced the shipped results in out/.

source(here::here("code", "00_setup.R"))


## ===========================================================================
## from 02_climate.R
## ===========================================================================

# 02_climate.R
# Phase 0, Task 0.3 - Build biological-side climate covariates per district per month
# from NASA POWER (free, no auth). Covers 2022-01 to 2025-12.
# 2026 monthly data is NOT available from NASA POWER (API confirmed: data ends 2025-12-31).
# Output: derived/climate_dm.rds  cols: dist_lgd, month, t2m, t2m_min, rh, ah, precip, solar


library(sf)
library(nasapower)

# ── Paths ─────────────────────────────────────────────────────────────────────
CACHE_DIR  <- file.path(DER, "power_cache")
OUT_FILE   <- file.path(DER, "climate_dm.rds")
GPKG_FILE  <- GEO
PANEL_FILE <- file.path(DER, "panel_dm.rds")

dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. District centroids (EPSG:4326) ─────────────────────────────────────────
cat("[1/4] Computing district centroids...\n")
districts <- sf::st_read(GPKG_FILE, quiet = TRUE)

# Ensure CRS is geographic before point_on_surface; transform if needed
if (!sf::st_is_longlat(districts)) {
  districts <- sf::st_transform(districts, 4326)
} else {
  districts <- sf::st_transform(districts, 4326)  # force canonical
}

# Suppress s2 warnings for point_on_surface on geographic CRS
suppressMessages(
  centroids <- sf::st_point_on_surface(districts)
)
coords <- sf::st_coordinates(centroids)
dist_coords <- data.frame(
  dist_lgd = as.integer(districts$dist_lgd_anl),
  lon      = round(coords[, "X"], 6),
  lat      = round(coords[, "Y"], 6),
  stringsAsFactors = FALSE
)

cat(sprintf("  %d district centroids computed\n", nrow(dist_coords)))

# ── 2. Pull NASA POWER monthly data (with per-district cache) ─────────────────
cat("[2/4] Fetching NASA POWER data (2022-2025)...\n")
cat("  Cache dir:", CACHE_DIR, "\n")

PARS      <- c("T2M", "T2M_MIN", "RH2M", "PRECTOTCORR", "ALLSKY_SFC_SW_DWN")
YEAR_FROM <- 2022L
YEAR_TO   <- 2025L   # 2026 confirmed unavailable; data ends 2025-12-31

failures   <- character(0)
n_cached   <- 0L
n_fetched  <- 0L
n_total    <- nrow(dist_coords)

for (i in seq_len(n_total)) {
  did  <- dist_coords$dist_lgd[i]
  lon  <- dist_coords$lon[i]
  lat  <- dist_coords$lat[i]
  cache_file <- file.path(CACHE_DIR, paste0(did, ".rds"))

  if (file.exists(cache_file)) {
    n_cached <- n_cached + 1L
    if (n_cached %% 100 == 0) cat(sprintf("  ... %d/%d cached\n", i, n_total))
    next
  }

  result <- tryCatch({
    raw <- nasapower::get_power(
      community    = "ag",
      temporal_api = "monthly",
      pars         = PARS,
      lonlat       = c(lon, lat),
      dates        = c(YEAR_FROM, YEAR_TO)
    )
    raw
  }, error = function(e) {
    cat(sprintf("  FAIL district %d: %s\n", did, conditionMessage(e)))
    NULL
  })

  if (is.null(result)) {
    failures <- c(failures, as.character(did))
    Sys.sleep(1.0)   # longer back-off on failure
    next
  }

  saveRDS(result, cache_file)
  n_fetched <- n_fetched + 1L

  if (n_fetched %% 50 == 0) {
    cat(sprintf("  ... fetched %d new (total %d/%d processed)\n", n_fetched, i, n_total))
  }

  # Polite delay: 0.2-0.5 s
  Sys.sleep(runif(1, 0.2, 0.5))
}

cat(sprintf("  Done: %d cached, %d newly fetched, %d failed\n",
            n_cached, n_fetched, length(failures)))

# ── 3. Reshape to long: dist_lgd x year x month ───────────────────────────────
cat("[3/4] Reshaping to long format...\n")

MONTH_COLS <- c("JAN","FEB","MAR","APR","MAY","JUN",
                "JUL","AUG","SEP","OCT","NOV","DEC")

all_rows <- vector("list", length = n_total)

for (i in seq_len(n_total)) {
  did        <- dist_coords$dist_lgd[i]
  cache_file <- file.path(CACHE_DIR, paste0(did, ".rds"))

  if (!file.exists(cache_file)) next  # failed district - skip

  raw <- readRDS(cache_file)

  # Wide format: PARAMETER x YEAR, months as columns (JAN..DEC, ANN)
  # Pivot to long: one row per dist_lgd x year x month x parameter
  month_present <- intersect(MONTH_COLS, colnames(raw))

  long <- raw[, c("PARAMETER", "YEAR", month_present), drop = FALSE]
  long <- tidyr::pivot_longer(long,
    cols      = dplyr::all_of(month_present),
    names_to  = "month_abbr",
    values_to = "value"
  )

  # Add district id
  long$dist_lgd <- did
  all_rows[[i]] <- long
}

combined <- dplyr::bind_rows(all_rows)
cat(sprintf("  Combined long rows: %d\n", nrow(combined)))

# Map month abbreviation to integer
month_map <- setNames(1:12, MONTH_COLS)
combined$month_num <- month_map[combined$month_abbr]

# Pivot wider: one row per dist_lgd x year x month
wide <- tidyr::pivot_wider(
  combined[, c("dist_lgd", "YEAR", "month_num", "PARAMETER", "value")],
  names_from  = "PARAMETER",
  values_from = "value"
)

# Rename to lower-case clean names
wide <- dplyr::rename(wide,
  year     = YEAR,
  t2m      = T2M,
  t2m_min  = T2M_MIN,
  rh       = RH2M,
  precip   = PRECTOTCORR,
  solar    = ALLSKY_SFC_SW_DWN
)

# ── Build yearmonth index matching panel ──────────────────────────────────────
wide$month <- tsibble::yearmonth(
  paste0(wide$year, "-", formatC(wide$month_num, width = 2, flag = "0"))
)

# ── 3b. Derive absolute humidity ─────────────────────────────────────────────
# AH (g/m3) = 6.112 * exp(17.67*T/(T+243.5)) * RH * 2.1674 / (273.15 + T)
wide$ah <- with(wide,
  6.112 * exp(17.67 * t2m / (t2m + 243.5)) * rh * 2.1674 / (273.15 + t2m)
)

# ── Final output columns ──────────────────────────────────────────────────────
climate_dm <- dplyr::select(wide,
  dist_lgd, month, t2m, t2m_min, rh, ah, precip, solar
)

# Sort
climate_dm <- dplyr::arrange(climate_dm, dist_lgd, month)

saveRDS(climate_dm, OUT_FILE)
cat(sprintf("  Saved: %s\n", OUT_FILE))
cat(sprintf("  Rows: %d, Districts: %d, Months: %d\n",
            nrow(climate_dm),
            dplyr::n_distinct(climate_dm$dist_lgd),
            dplyr::n_distinct(climate_dm$month)))

# ── 4. VERIFICATION ───────────────────────────────────────────────────────────
cat("\n[4/4] VERIFICATION\n")
cat(rep("-", 60), "\n", sep = "")

panel <- readRDS(PANEL_FILE)

# Coverage
n_districts_climate <- dplyr::n_distinct(climate_dm$dist_lgd)
n_distmonths        <- nrow(climate_dm)
cat(sprintf("Districts covered: %d / 747\n", n_districts_climate))
cat(sprintf("Total district-months in climate_dm: %d\n", n_distmonths))

# Panel keys for 2022-01..2025-12 (complete window, excluding 2026)
panel_complete <- dplyr::filter(panel,
  month >= tsibble::yearmonth("2022 Jan"),
  month <= tsibble::yearmonth("2025 Dec")
)
panel_keys <- dplyr::distinct(as.data.frame(panel_complete)[, c("dist_lgd","month")])

climate_keys <- dplyr::distinct(climate_dm[, c("dist_lgd","month")])

# Anti-join: panel keys missing from climate
missing_in_climate <- dplyr::anti_join(panel_keys, climate_keys,
                                        by = c("dist_lgd","month"))
n_missing_complete <- nrow(missing_in_climate)
cat(sprintf("\nComplete window (2022-01..2025-12):\n"))
cat(sprintf("  Panel district-months: %d\n", nrow(panel_keys)))
cat(sprintf("  Missing in climate_dm (must be 0 for OK): %d\n", n_missing_complete))
if (n_missing_complete == 0) {
  cat("  [OK] All panel district-months covered in complete window\n")
} else {
  cat("  [WARN] Missing district-months:\n")
  print(head(missing_in_climate, 20))
}

# 2026 coverage
panel_2026 <- dplyr::filter(panel,
  month >= tsibble::yearmonth("2026 Jan")
)
climate_2026 <- dplyr::filter(climate_dm,
  month >= tsibble::yearmonth("2026 Jan")
)
cat(sprintf("\n2026 coverage:\n"))
cat(sprintf("  Panel district-months in 2026: %d\n", nrow(panel_2026)))
cat(sprintf("  Climate district-months in 2026: %d (expected: 0 - API unavailable)\n",
            nrow(climate_2026)))

# Failures
cat(sprintf("\nFailed fetches: %d (must be 0 for OK)\n", length(failures)))
if (length(failures) > 0) {
  cat("  Failed districts:", paste(failures, collapse = ", "), "\n")
}

# Sanity check ranges
cat("\nSanity ranges (complete window 2022-2025):\n")
clim_complete <- dplyr::filter(climate_dm,
  month >= tsibble::yearmonth("2022 Jan"),
  month <= tsibble::yearmonth("2025 Dec")
)

sanity_check <- function(x, label, lo, hi) {
  mn  <- min(x, na.rm = TRUE)
  med <- median(x, na.rm = TRUE)
  mx  <- max(x, na.rm = TRUE)
  ok  <- mn >= lo & mx <= hi
  flag <- if (ok) "[OK]" else "[WARN]"
  cat(sprintf("  %s %s: min=%.2f, median=%.2f, max=%.2f (expected %g..%g)\n",
              flag, label, mn, med, mx, lo, hi))
}

sanity_check(clim_complete$t2m,     "t2m (C)",    -5,  40)
sanity_check(clim_complete$t2m_min, "t2m_min (C)",-10, 40)
sanity_check(clim_complete$rh,      "rh (%)",      0,  100)
sanity_check(clim_complete$ah,      "ah (g/m3)",   0,   60)
sanity_check(clim_complete$precip,  "precip (mm/d)",0, 30)
sanity_check(clim_complete$solar,   "solar (MJ/m2/d)", 0, 35)

# Overall verdict
verdict_ok <- (n_missing_complete == 0) && (length(failures) == 0) &&
              (n_districts_climate == 747)
cat(sprintf("\n%s\n", rep("-", 60)))
if (verdict_ok) {
  cat("VERIFICATION: PASS\n")
} else {
  cat("VERIFICATION: FAIL - see warnings above\n")
}

cat("\n2026 handling note: NASA POWER monthly data confirmed unavailable beyond\n")
cat("2025-12-31. Script fetched 2022-2025 only. 2026 rows absent from climate_dm\n")
cat("as expected - the 2026 tail is largely excluded from seasonal analysis.\n")

## ===========================================================================
## from 05_calendar.R
## ===========================================================================

# 05_calendar.R
# Phase 0, Task 0.6 - Reporting-calendar / programme covariates
#
# Produces:
#   derived/calendar.rds     - one row per month 2022-01..2026-12
#   derived/campaign_dm.rds  - district x month campaign flag
#
# Sources:
#   - DoPT OM F.No.12/5/2021-JCA-2 (2022 holidays) via igecorner.com
#   - DoPT OM for 2023 via 7thpaycommissionnews.in
#   - DoPT OM for 2024 via 7thpaycommissionnews.in
#   - DoPT OM for 2025 via gconnect.in (DOPT OM dated 9 Jul 2024)
#   - 2026 via cleartax.in / bankbazaar.com (official DoPT circular)
#   - officeholidays.com cross-checked for all years


library(tsibble)   # yearmonth
library(splines)   # ns()

# ── 1. National gazetted holiday table (weekday-only closures) ──────────────
#
# FIXED holidays (same date every year):
#   Republic Day: 26 Jan
#   Independence Day: 15 Aug
#   Gandhi Jayanti: 2 Oct
#   Christmas: 25 Dec
#
# MOVABLE holidays - actual gazetted dates sourced from DoPT OMs per year:
#   (State-specific holidays are OUT OF SCOPE - documented refinement)
#
# NOTE: Only weekday occurrences are subtracted from working_days.
#       Holidays that fall on Sat/Sun are noted but not subtracted.
#
# Muharram: Islamic New Year (Ashura/Muharram 10th)
# Id-e-Milad / Milad-un-Nabi: Prophet's Birthday
# Janmashtami: Krishna's birthday (Vaishnava date used, as per DoPT)
# Dussehra / Vijaya Dashami
# Guru Nanak Jayanti
# Mahavir Jayanti
# Buddha Purnima
# Holi
# Good Friday (Easter-linked)
# Id-ul-Fitr (Eid)
# Id-ul-Zuha / Bakrid
#
# 2022 -----------------------------------------------------------------------
# Source: DoPT OM F.No.12/5/2021-JCA-2; confirmed via igecorner.com
hols_2022 <- as.Date(c(
  "2022-01-26",  # Republic Day (fixed)
  "2022-03-18",  # Holi
  "2022-04-14",  # Mahavir Jayanti
  "2022-04-15",  # Good Friday
  "2022-05-03",  # Id-ul-Fitr
  "2022-05-16",  # Buddha Purnima
  "2022-07-10",  # Id-ul-Zuha (Bakrid)  -- falls Sunday; no weekday deduction
  "2022-08-09",  # Muharram
  "2022-08-15",  # Independence Day (fixed)
  "2022-08-19",  # Janmashtami (Vaishnava)
  "2022-10-02",  # Gandhi Jayanti (fixed)  -- also Dussehra sub; see below
  "2022-10-05",  # Dussehra (Vijaya Dashami)
  "2022-10-09",  # Milad-un-Nabi (Id-e-Milad)  -- falls Sunday; no wd deduction
  "2022-10-24",  # Diwali (Deepavali)
  "2022-11-08",  # Guru Nanak Jayanti
  "2022-12-25"   # Christmas (fixed)
))

# 2023 -----------------------------------------------------------------------
# Source: DoPT OM issued June 2022; confirmed via 7thpaycommissionnews.in
hols_2023 <- as.Date(c(
  "2023-01-26",  # Republic Day (fixed)
  "2023-03-08",  # Holi
  "2023-03-30",  # Ram Navami (gazetted 2023)
  "2023-04-04",  # Mahavir Jayanti
  "2023-04-07",  # Good Friday
  "2023-04-22",  # Id-ul-Fitr  -- falls Saturday; no weekday deduction
  "2023-05-05",  # Buddha Purnima
  "2023-06-29",  # Id-ul-Zuha (Bakrid)
  "2023-07-29",  # Muharram  -- falls Saturday; no weekday deduction
  "2023-08-15",  # Independence Day (fixed)
  "2023-09-07",  # Janmashtami (Vaishnava)
  "2023-09-28",  # Milad-un-Nabi (Id-e-Milad)
  "2023-10-02",  # Gandhi Jayanti (fixed)
  "2023-10-24",  # Dussehra
  "2023-11-12",  # Diwali  -- falls Sunday; no weekday deduction
  "2023-11-27",  # Guru Nanak Jayanti
  "2023-12-25"   # Christmas (fixed)
))

# 2024 -----------------------------------------------------------------------
# Source: DoPT OM; confirmed via 7thpaycommissionnews.in
hols_2024 <- as.Date(c(
  "2024-01-26",  # Republic Day (fixed)
  "2024-03-25",  # Holi
  "2024-03-29",  # Good Friday
  "2024-04-11",  # Id-ul-Fitr
  "2024-04-17",  # Ram Navami (gazetted 2024)
  "2024-04-21",  # Mahavir Jayanti  -- falls Sunday; no weekday deduction
  "2024-05-23",  # Buddha Purnima
  "2024-06-17",  # Id-ul-Zuha (Bakrid)
  "2024-07-17",  # Muharram
  "2024-08-15",  # Independence Day (fixed)
  "2024-08-26",  # Janmashtami (Vaishnava)
  "2024-09-16",  # Milad-un-Nabi (Id-e-Milad)
  "2024-10-02",  # Gandhi Jayanti (fixed)
  "2024-10-12",  # Dussehra  -- falls Saturday; no weekday deduction
  "2024-10-31",  # Diwali (Deepavali)
  "2024-11-15",  # Guru Nanak Jayanti
  "2024-12-25"   # Christmas (fixed)
))

# 2025 -----------------------------------------------------------------------
# Source: DoPT OM dated 9 Jul 2024; confirmed via gconnect.in
# Note: Dussehra (Oct 2) coincides with Gandhi Jayanti; listed once.
hols_2025 <- as.Date(c(
  "2025-01-26",  # Republic Day (fixed)
  "2025-03-14",  # Holi
  "2025-03-31",  # Id-ul-Fitr
  "2025-04-10",  # Mahavir Jayanti
  "2025-04-18",  # Good Friday
  "2025-05-12",  # Buddha Purnima
  "2025-06-07",  # Id-ul-Zuha (Bakrid)  -- falls Saturday; no weekday deduction
  "2025-07-06",  # Muharram
  "2025-08-15",  # Independence Day (fixed)
  "2025-08-16",  # Janmashtami (Vaishnava)
  "2025-09-05",  # Milad-un-Nabi (Id-e-Milad)
  "2025-10-02",  # Gandhi Jayanti (fixed) / Dussehra (same date 2025)
  "2025-10-20",  # Diwali (Deepavali)
  "2025-11-05",  # Guru Nanak Jayanti
  "2025-12-25"   # Christmas (fixed)
))

# 2026 -----------------------------------------------------------------------
# Source: bankbazaar.com / cleartax.in referencing DoPT circular 2026
# Note: Diwali 2026 falls on Sunday (Nov 8); no weekday deduction.
hols_2026 <- as.Date(c(
  "2026-01-26",  # Republic Day (fixed)
  "2026-03-04",  # Holi
  "2026-03-21",  # Id-ul-Fitr  -- falls Saturday; no weekday deduction
  "2026-03-26",  # Ram Navami (gazetted 2026)
  "2026-03-31",  # Mahavir Jayanti
  "2026-04-03",  # Good Friday
  "2026-05-01",  # Buddha Purnima
  "2026-05-27",  # Id-ul-Zuha (Bakrid)
  "2026-06-26",  # Muharram
  "2026-08-15",  # Independence Day (fixed)
  "2026-08-25",  # Id-e-Milad (Milad-un-Nabi)  -- some sources say 26 Aug
  "2026-09-04",  # Janmashtami
  "2026-10-02",  # Gandhi Jayanti (fixed)
  "2026-10-20",  # Dussehra
  "2026-11-08",  # Diwali  -- falls Sunday; no weekday deduction
  "2026-11-24",  # Guru Nanak Jayanti
  "2026-12-25"   # Christmas (fixed)
))

all_holidays <- c(hols_2022, hols_2023, hols_2024, hols_2025, hols_2026)

# ── 2. Working-days calculation ─────────────────────────────────────────────

# Weekday holidays only (Mon=2 .. Fri=6 in R's weekdays() / as.POSIXlt$wday)
weekday_holidays <- all_holidays[
  !weekdays(all_holidays) %in% c("Saturday", "Sunday")
]

# Count weekdays (Mon-Fri) in a month
count_weekdays <- function(yr, mo) {
  first_day <- as.Date(sprintf("%d-%02d-01", yr, mo))
  # Last day of month: first day of next month minus 1
  next_first <- if (mo < 12L) {
    as.Date(sprintf("%d-%02d-01", yr, mo + 1L))
  } else {
    as.Date(sprintf("%d-01-01", yr + 1L))
  }
  days <- seq(first_day, next_first - 1L, by = "day")
  sum(!weekdays(days) %in% c("Saturday", "Sunday"))
}

# ── 3. Build calendar tibble 2022-01 .. 2026-12 ─────────────────────────────

months_seq <- seq(
  as.Date("2022-01-01"),
  as.Date("2026-12-01"),
  by = "month"
)

# Month index for spline (1 = Jan 2022)
month_index <- seq_along(months_seq)

# Spline basis for COVID recovery (natural spline, df=3)
covid_basis <- splines::ns(month_index, df = 3)

calendar <- tibble::tibble(
  date       = months_seq,
  month      = tsibble::yearmonth(months_seq),
  yr         = as.integer(format(months_seq, "%Y")),
  mo         = as.integer(format(months_seq, "%m")),
  month_idx  = month_index
) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    # Total weekdays in month
    total_wd = count_weekdays(yr, mo),
    # Gazetted weekday holidays falling in this month
    hols_in_month = sum(
      format(weekday_holidays, "%Y-%m") == sprintf("%d-%02d", yr, mo)
    ),
    working_days = total_wd - hols_in_month,
    fy_end = (mo == 3L),
    q_end  = (mo %in% c(3L, 6L, 9L, 12L))
  ) |>
  dplyr::ungroup() |>
  dplyr::bind_cols(
    tibble::as_tibble(covid_basis, .name_repair = "minimal") |>
      stats::setNames(c("covid_s1", "covid_s2", "covid_s3"))
  ) |>
  dplyr::select(month, working_days, fy_end, q_end, covid_s1, covid_s2, covid_s3)

# ── 4. Verification ─────────────────────────────────────────────────────────

panel_dm       <- readRDS(file.path(DER, "panel_dm.rds"))
# index() returns the index name; pull the actual column values
panel_months   <- unique(panel_dm[[as.character(tsibble::index(panel_dm))]])
district_frame <- readRDS(file.path(DER, "district_frame.rds"))

cat("\n=== VERIFICATION ===\n")

# V1: All panel months covered in calendar
missing_months <- setdiff(as.character(panel_months), as.character(calendar$month))
cat(sprintf("Calendar rows: %d\n", nrow(calendar)))
cat(sprintf("Panel months missing from calendar: %d  [OK if 0]\n",
            length(missing_months)))
if (length(missing_months) > 0) cat("  Missing:", paste(missing_months, collapse=", "), "\n")

# V2: working_days range
cat(sprintf("working_days min=%d  max=%d  [OK if 18..27]\n",
            min(calendar$working_days), max(calendar$working_days)))

# Per-year holiday count
hol_by_year <- tapply(
  weekday_holidays,
  format(weekday_holidays, "%Y"),
  length
)
cat("Weekday holidays used per year:\n")
for (yr in names(hol_by_year)) cat(sprintf("  %s: %d\n", yr, hol_by_year[[yr]]))

# V3: fy_end / q_end flags
n_fy_end <- sum(calendar$fy_end)
cat(sprintf("fy_end TRUE count=%d (expect 5, one March per year 2022-2026)\n", n_fy_end))
# Extract month number from yearmonth for display
fy_mo_nums <- sort(unique(as.integer(format(as.Date(calendar$month[calendar$fy_end]), "%m"))))
cat(sprintf("fy_end months (month numbers): %s  [OK if only 3]\n", paste(fy_mo_nums, collapse=",")))

n_q_end <- sum(calendar$q_end)
q_mo_nums <- sort(unique(as.integer(format(as.Date(calendar$month[calendar$q_end]), "%m"))))
cat(sprintf("q_end TRUE count=%d (expect 20 rows)\n", n_q_end))
cat(sprintf("q_end months (month numbers): %s  [OK if 3,6,9,12]\n", paste(q_mo_nums, collapse=",")))

# ── 5. Save calendar ────────────────────────────────────────────────────────
saveRDS(calendar, file.path(DER, "calendar.rds"))
cat(sprintf("\nSaved: %s\n", file.path(DER, "calendar.rds")))

# ── 6. Build campaign_dm ────────────────────────────────────────────────────
#
# NTEP 100-Day Intensified Active Case-Finding TB Campaign
# Verified window: ~7 Dec 2024 to 17 Mar 2025
# Flagged months: 2024-12, 2025-01, 2025-02, 2025-03
#
# NOTE: The campaign initially targeted 347 specific districts.
# The curated 347-district list has NOT yet been compiled.
# This implementation uses a NATIONAL APPROXIMATION: all 747 districts
# are flagged for those 4 months.
# TODO: Replace with 347-district curated list when available.

campaign_months <- tsibble::yearmonth(c("2024 Dec", "2025 Jan", "2025 Feb", "2025 Mar"))

# district_frame may be a tibble or tsibble; extract the dist_lgd column
if ("dist_lgd" %in% names(district_frame)) {
  dist_ids <- district_frame$dist_lgd
} else {
  stop("district_frame does not have a 'dist_lgd' column")
}

campaign_dm <- tidyr::expand_grid(
  dist_lgd = dist_ids,
  month    = tsibble::yearmonth(
    seq(as.Date("2022-01-01"), as.Date("2026-05-01"), by = "month")
  )
) |>
  dplyr::mutate(
    campaign = as.integer(month %in% campaign_months)
  )

# V4: Campaign verification
cat("\n--- Campaign verification ---\n")
n_campaign_1 <- sum(campaign_dm$campaign == 1L)
campaign_months_found <- sort(unique(campaign_dm$month[campaign_dm$campaign == 1L]))
cat(sprintf("campaign==1 rows: %d  [expect 2988 = 747 * 4]\n", n_campaign_1))
cat(sprintf("Campaign months: %s  [expect 2024-Dec..2025-Mar]\n",
            paste(as.character(campaign_months_found), collapse=", ")))
cat(sprintf("Districts flagged: %d  [NATIONAL APPROX - replace with 347-district list]\n",
            length(unique(campaign_dm$dist_lgd[campaign_dm$campaign == 1L]))))

saveRDS(campaign_dm, file.path(DER, "campaign_dm.rds"))
cat(sprintf("Saved: %s\n", file.path(DER, "campaign_dm.rds")))

cat("\n=== DONE ===\n")

## ===========================================================================
## from 06_modifiers.R
## ===========================================================================

# 06_modifiers.R
# Build derived/modifiers.rds: one row per analysis district (747-frame)
# with district covariates that predict WHERE the reporting artefact bites.
# Keys: dist_lgd (== dist_lgd_anl in the analysis frame).
#
# PART A: In-hand predictors (~100% coverage expected)
#   1. latitude (centroid, EPSG:4326)
#   2. area_km2 (st_area), density_2024
#   3. public_share (sum public / sum total over complete months 2022-2025)
#   4. pct_urban_notif (notification-weighted urban share via DuckDB + tu_panel.parquet)
#   5. region, stname, total_notif_2022_2025
#
# PART B: Best-effort downloadable structural modifiers
#   6. dist_lgd <-> dtcode11 crosswalk from subdistrict_analysis.csv
#   7. Census-2011 district indicators (data.gov.in PCA / household amenities)
#   8. NFHS-5 (2019-21) district indicators


library(sf)
library(duckdb); library(DBI)

cat("===== 06_modifiers.R =====\n")
cat("Working dir:", getwd(), "\n\n")

# =========================================================
# 0. Load lookup: dist_lgd -> dist_lgd_anl
# =========================================================
lookup <- readRDS(file.path(FILES, "dist_lgd_anl_lookup.rds"))
# PoK / unresolvable districts map to NA dist_lgd_anl -> drop them
lookup <- lookup[!is.na(lookup$dist_lgd_anl), ]
cat("Lookup rows (non-NA anl):", nrow(lookup), "\n")
cat("Unique dist_lgd_anl:", length(unique(lookup$dist_lgd_anl)), "\n\n")

# Helper: recode dist_lgd -> dist_lgd_anl
recode_anl <- function(df, col = "dist_lgd") {
  df <- merge(df, lookup[, c("dist_lgd", "dist_lgd_anl")],
              by.x = col, by.y = "dist_lgd", all.x = FALSE)
  df
}

# =========================================================
# PART A.1-2  Geometry: latitude & area from district_anl.gpkg
# =========================================================
cat("--- PART A.1-2: Geometry ---\n")
gpkg_path <- GEO
gpkg <- sf::st_read(gpkg_path, quiet = TRUE)
# gpkg key: dist_lgd_anl (already the analysis code)
cat("gpkg dims:", dim(gpkg), "\n")
cat("CRS:", sf::st_crs(gpkg)$Name, "\n")

# Ensure WGS84 for latitude
if (is.na(sf::st_crs(gpkg)) || sf::st_crs(gpkg)$epsg != 4326) {
  gpkg <- sf::st_transform(gpkg, 4326)
}

# Point on surface for stable centroid (handles non-convex polygons)
pts <- sf::st_point_on_surface(gpkg)
coords <- sf::st_coordinates(pts)
lat_df <- data.frame(
  dist_lgd_anl = gpkg$dist_lgd_anl,
  latitude     = coords[, "Y"]
)

# Area in km2: use equal-area projection (EPSG:32644 UTM zone 44N covers India)
gpkg_ea <- sf::st_transform(gpkg, 32644)
area_m2 <- as.numeric(sf::st_area(gpkg_ea))
area_df <- data.frame(
  dist_lgd_anl = gpkg$dist_lgd_anl,
  area_km2     = area_m2 / 1e6
)

geom_df <- merge(lat_df, area_df, by = "dist_lgd_anl")
cat("Geometry frame n:", nrow(geom_df), "\n")
cat("Latitude range:", round(range(geom_df$latitude), 2), "\n")
cat("Area_km2 range:", round(range(geom_df$area_km2), 1), "\n\n")

# =========================================================
# PART A.2 (cont.)  Population 2024 -> density
# =========================================================
cat("--- PART A.2 (cont.): Population & density ---\n")
pop_long <- readRDS(file.path(DER, "district_pop_anl.rds"))
pop24 <- pop_long[pop_long$yr == 2024, c("dist_lgd_anl", "pop")]
names(pop24)[2] <- "pop_2024"
cat("Pop 2024 n districts:", nrow(pop24), "\n")

# Merge geom + pop -> density
geom_pop <- merge(geom_df, pop24, by = "dist_lgd_anl", all.x = TRUE)
geom_pop$density_2024 <- geom_pop$pop_2024 / geom_pop$area_km2
cat("density_2024 range (pop/km2):", round(range(geom_pop$density_2024, na.rm=TRUE), 1), "\n\n")

# =========================================================
# PART A.3  public_share from panel_dm (complete months 2022-2025)
# =========================================================
cat("--- PART A.3: public_share ---\n")
panel <- readRDS(file.path(DER, "panel_dm.rds"))
# panel key: dist_lgd (raw, pre-recode), cols: public, private, total, month, complete
panel_df <- as.data.frame(panel)

# Filter complete months, years 2022-2025
panel_df$year <- as.integer(format(panel_df$month, "%Y"))
panel_cmp <- panel_df[panel_df$complete & panel_df$year >= 2022 & panel_df$year <= 2025, ]
cat("Complete obs 2022-2025:", nrow(panel_cmp), "\n")

pub_share <- aggregate(
  cbind(public = panel_cmp$public, total = panel_cmp$total),
  by = list(dist_lgd = panel_cmp$dist_lgd),
  FUN = sum, na.rm = TRUE
)
pub_share$public_share <- pub_share$public / pub_share$total

# Also total notifications 2022-2025
total_notif <- data.frame(
  dist_lgd           = pub_share$dist_lgd,
  total_notif_2022_2025 = pub_share$total
)

# Recode to anl frame
pub_share_anl <- recode_anl(pub_share[, c("dist_lgd","public","total","public_share")])
# If multiple old districts -> one analysis district: aggregate public+total, recalc share
pub_share_anl2 <- aggregate(
  cbind(public = pub_share_anl$public, total = pub_share_anl$total),
  by = list(dist_lgd_anl = pub_share_anl$dist_lgd_anl),
  FUN = sum, na.rm = TRUE
)
pub_share_anl2$public_share <- pub_share_anl2$public / pub_share_anl2$total
pub_share_anl2$total_notif_2022_2025 <- pub_share_anl2$total
pub_share_anl2 <- pub_share_anl2[, c("dist_lgd_anl","public_share","total_notif_2022_2025")]
cat("public_share n:", nrow(pub_share_anl2), "\n")
cat("public_share range:", round(range(pub_share_anl2$public_share, na.rm=TRUE), 3), "\n\n")

# =========================================================
# PART A.4  pct_urban_notif via DuckDB + tu_panel.parquet
# =========================================================
cat("--- PART A.4: pct_urban_notif (DuckDB) ---\n")

tu2dist <- readRDS(file.path(DER, "tu2district.rds"))
cat("tu2district urban_rural values:\n")
print(table(tu2dist$urban_rural, useNA = "always"))
# urban_rural: Rural / Urban / Water
# Treat "Urban" as urban; "Rural" and "Water" as non-urban
tu2dist$is_urban <- tu2dist$urban_rural == "Urban"
cat("Urban mapping: Urban=TRUE, Rural/Water=FALSE\n")
tu2dist_slim <- tu2dist[, c("tu_unique", "dist_lgd", "is_urban")]

# Use DuckDB to aggregate notifications by district x urban_rural
pq_path <- file.path(ROOT, "data/derived/tu_panel.parquet")
con <- dbConnect(duckdb(), dbdir = ":memory:")

# Write tu2dist_slim as a DuckDB table
dbWriteTable(con, "tu2dist", tu2dist_slim, overwrite = TRUE)

# Read parquet; join to tu2dist; aggregate 2022-2025
# Note: tu_panel.parquet has tu_unique+period but NOT dist_lgd;
# dist_lgd comes from tu2dist (d.dist_lgd).
q_urban <- sprintf(
  "SELECT d.dist_lgd,
          SUM(CASE WHEN d.is_urban THEN t.total ELSE 0 END) AS urban_notif,
          SUM(t.total) AS all_notif
   FROM read_parquet('%s') AS t
   JOIN tu2dist AS d ON t.tu_unique = d.tu_unique
   WHERE t.period >= '2022-01' AND t.period <= '2025-12'
   GROUP BY d.dist_lgd",
  pq_path
)
urban_agg <- dbGetQuery(con, q_urban)
dbDisconnect(con)

urban_agg$pct_urban_notif <- urban_agg$urban_notif / urban_agg$all_notif
cat("DuckDB urban agg n districts:", nrow(urban_agg), "\n")

# Check how many tu_unique in panel are matched to dist_lgd via tu2dist
cat("pct_urban_notif range:", round(range(urban_agg$pct_urban_notif, na.rm=TRUE), 3), "\n")
cat("NAs in pct_urban_notif:", sum(is.na(urban_agg$pct_urban_notif)), "\n")

# Recode dist_lgd -> dist_lgd_anl; pop-weighted mean for merged districts
urban_anl <- recode_anl(urban_agg[, c("dist_lgd","urban_notif","all_notif")])
urban_anl2 <- aggregate(
  cbind(urban_notif = urban_anl$urban_notif, all_notif = urban_anl$all_notif),
  by = list(dist_lgd_anl = urban_anl$dist_lgd_anl),
  FUN = sum, na.rm = TRUE
)
urban_anl2$pct_urban_notif <- urban_anl2$urban_notif / urban_anl2$all_notif
urban_anl2 <- urban_anl2[, c("dist_lgd_anl","pct_urban_notif")]
cat("After recode, n:", nrow(urban_anl2), "\n")
cat("NOTE: pct_urban_notif reflects the TU-level urban_rural classification in tu2district.\n")
cat("  ~93% of TUs are classified Urban; districts with only Urban TUs get pct=1.0\n")
cat("  (", sum(urban_anl2$pct_urban_notif == 1.0, na.rm=TRUE), "of", nrow(urban_anl2),
    "have pct_urban_notif=1.0; variable still discriminates",
    sum(urban_anl2$pct_urban_notif > 0 & urban_anl2$pct_urban_notif < 1, na.rm=TRUE),
    "mixed-coverage districts)\n\n")

# =========================================================
# PART A.5  region, stname from tu2district
# =========================================================
cat("--- PART A.5: region, stname ---\n")
# Get district-level region + stname (mode per district if multiple TUs)
# Most districts have one stname/region; take the first unique per dist_lgd
reg_st <- unique(tu2dist[, c("dist_lgd","stname","region")])
# In case of duplicates (should not happen for stname/region at district level)
reg_st <- reg_st[!duplicated(reg_st$dist_lgd), ]
cat("reg_st n:", nrow(reg_st), "\n")

reg_st_anl <- recode_anl(reg_st)
# If a merge produced duplicates (very unlikely), keep first
reg_st_anl <- reg_st_anl[!duplicated(reg_st_anl$dist_lgd_anl), ]
reg_st_anl <- reg_st_anl[, c("dist_lgd_anl","stname","region")]
cat("reg_st_anl n:", nrow(reg_st_anl), "\n\n")

# =========================================================
# Assemble Part A base frame on the 747-analysis frame
# =========================================================
cat("--- Assembling Part A ---\n")
base <- geom_pop[, c("dist_lgd_anl","latitude","area_km2","pop_2024","density_2024")]
base <- merge(base, pub_share_anl2, by = "dist_lgd_anl", all.x = TRUE)
base <- merge(base, urban_anl2,    by = "dist_lgd_anl", all.x = TRUE)
base <- merge(base, reg_st_anl,    by = "dist_lgd_anl", all.x = TRUE)
cat("Part A frame rows:", nrow(base), "\n")
cat("Part A frame cols:", ncol(base), "\n")
# Rename dist_lgd_anl to dist_lgd (== dist_lgd_anl in the 747-frame)
names(base)[names(base) == "dist_lgd_anl"] <- "dist_lgd"

# =========================================================
# PART B.6  dist_lgd <-> dtcode11 crosswalk from subdistrict_analysis.csv
# =========================================================
cat("\n--- PART B.6: Census crosswalk ---\n")
sq_csv <- read.csv(file.path(SQ, "subdistrict_analysis.csv"), stringsAsFactors = FALSE)
# Build distinct crosswalk: one row per dist_lgd
census_xwalk <- unique(sq_csv[, c("dist_lgd","dtcode11","stcode11","dtname","stname")])
census_xwalk <- census_xwalk[!duplicated(census_xwalk$dist_lgd), ]
cat("Census xwalk n:", nrow(census_xwalk), "\n")
# Recode to anl frame
census_xwalk_anl <- merge(census_xwalk, lookup, by = "dist_lgd", all.x = FALSE)
census_xwalk_anl <- census_xwalk_anl[!is.na(census_xwalk_anl$dist_lgd_anl), ]
census_xwalk_anl <- census_xwalk_anl[!duplicated(census_xwalk_anl$dist_lgd_anl), ]
census_xwalk_anl <- census_xwalk_anl[, c("dist_lgd_anl","dtcode11","stcode11")]
names(census_xwalk_anl)[names(census_xwalk_anl) == "dist_lgd_anl"] <- "dist_lgd"
cat("Census xwalk after recode n:", nrow(census_xwalk_anl), "\n")

# Merge into base
base <- merge(base, census_xwalk_anl, by = "dist_lgd", all.x = TRUE)
cat("dtcode11 coverage:", sum(!is.na(base$dtcode11)), "/", nrow(base), "\n\n")

# =========================================================
# PART B.7  Census-2011 district indicators from data.gov.in
# =========================================================
cat("--- PART B.7: Census-2011 district indicators (best-effort download) ---\n")
b7_fetched <- FALSE
b7_joined  <- FALSE
b7_skip_reason <- NA_character_

tryCatch({
  # Try SHRUG (Socioeconomic High-resolution Rural-Urban Geographic) dataset:
  # a well-known research-grade district-level Census-2011 CSV published by IDFC.
  # Also try datameet compiled PCA CSV.
  census_urls <- list(
    # SHRUG district-level CSV (literacy, SC, ST, etc.) - GitHub mirror
    list(url  = "https://raw.githubusercontent.com/devdatalab/shrug/main/data/shrug-v1.5/pd01_district_level.csv",
         ext  = ".csv",
         desc = "SHRUG pd01_district_level"),
    # datameet India compiled Census 2011 district socio data
    list(url  = "https://raw.githubusercontent.com/datameet/india-district-boundaries/master/census2011/primary-census-abstract-data-tables/DDW_PCA0000_2011_Indiastatedist.xlsx",
         ext  = ".xlsx",
         desc = "datameet DDW_PCA0000_2011")
  )
  census_df <- NULL
  for (src in census_urls) {
    tmp_f <- tempfile(fileext = src$ext)
    dl_res <- tryCatch(
      download.file(src$url, destfile = tmp_f,
                    quiet = TRUE, method = "curl", timeout = 40),
      error = function(e) -1L
    )
    if (!identical(dl_res, 0L)) {
      cat("  Failed:", src$desc, "\n"); next
    }
    # Validate: file must be > 5KB and parseable as CSV (or xlsx)
    fsz <- file.info(tmp_f)$size
    if (is.na(fsz) || fsz < 5000) {
      cat("  Too small (likely error page):", src$desc, "size=", fsz, "\n"); next
    }
    if (src$ext == ".csv") {
      df_try <- tryCatch(read.csv(tmp_f, stringsAsFactors=FALSE, nrows=5),
                         error = function(e) NULL)
      if (is.null(df_try) || nrow(df_try) == 0) {
        cat("  Not valid CSV:", src$desc, "\n"); next
      }
      census_df <- read.csv(tmp_f, stringsAsFactors=FALSE)
      cat("  SUCCESS:", src$desc, "dims:", paste(dim(census_df), collapse="x"), "\n")
      b7_fetched <- TRUE
      break
    } else if (src$ext == ".xlsx") {
      if (!requireNamespace("readxl", quietly=TRUE)) {
        cat("  readxl not available, skipping xlsx:", src$desc, "\n"); next
      }
      df_try <- tryCatch(readxl::read_excel(tmp_f, n_max=5), error=function(e) NULL)
      if (is.null(df_try) || nrow(df_try) == 0) {
        cat("  Not valid xlsx:", src$desc, "\n"); next
      }
      census_df <- as.data.frame(readxl::read_excel(tmp_f))
      cat("  SUCCESS:", src$desc, "dims:", paste(dim(census_df), collapse="x"), "\n")
      b7_fetched <- TRUE
      break
    }
  }
  if (!b7_fetched) stop("All Census-2011 URLs failed or returned unusable files")

  # If fetched: try to identify and join dtcode11-keyed columns
  if (b7_fetched && !is.null(census_df)) {
    cat("Census-2011 cols (first 20):", paste(names(census_df)[1:min(20,ncol(census_df))], collapse=", "), "\n")
    # Look for dtcode / district code column
    code_col <- grep("^(district_code|dtcode|dist_code|D_CODE|PC11_D_ID|distcode)",
                     names(census_df), ignore.case=TRUE, value=TRUE)[1]
    cat("District code column found:", ifelse(is.na(code_col), "NONE", code_col), "\n")
    # Further processing skipped if code column not found - document only
    if (is.na(code_col)) {
      b7_skip_reason <- "Fetched Census-2011 file but no dtcode11-compatible column identified for join"
      cat("Census-2011: FETCHED but not joined -", b7_skip_reason, "\n")
    }
  }
}, error = function(e) {
  b7_skip_reason <<- paste("Census-2011:", conditionMessage(e))
  cat("SKIPPED - Census-2011:", b7_skip_reason, "\n")
})
if (!b7_fetched && is.na(b7_skip_reason)) {
  b7_skip_reason <- "Census-2011 download failed or produced no usable file"
}
if (!b7_fetched) {
  cat("Census-2011 indicators: SKIPPED (", b7_skip_reason, ")\n")
} else if (!b7_joined) {
  cat("Census-2011 indicators: FETCHED but not joined (",
      ifelse(is.na(b7_skip_reason), "no join key match", b7_skip_reason), ")\n")
}

# =========================================================
# PART B.8  NFHS-5 (2019-21) district indicators (best-effort)
# =========================================================
cat("\n--- PART B.8: NFHS-5 district indicators (best-effort download) ---\n")
b8_fetched  <- FALSE
b8_joined   <- FALSE
b8_skip_reason <- NA_character_
b8_match_rate <- NA_real_

tryCatch({
  # Try multiple compiled NFHS-5 district CSV sources
  nfhs5_urls <- c(
    "https://raw.githubusercontent.com/pratapvardhan/NFHS-5/main/data/nfhs5_district.csv",
    "https://raw.githubusercontent.com/pratapvardhan/NFHS-5/master/data/nfhs5_district.csv",
    "https://raw.githubusercontent.com/cegis-org/nfhs-5/main/data/nfhs5_district_indicators.csv"
  )
  tmp_nfhs <- tempfile(fileext = ".csv")
  dl_res <- -1L
  for (nfhs5_url in nfhs5_urls) {
    dl_res <- tryCatch(
      download.file(nfhs5_url, destfile = tmp_nfhs,
                    quiet = TRUE, method = "curl", timeout = 45),
      error = function(e) -1L
    )
    if (identical(dl_res, 0L)) {
      fsz <- file.info(tmp_nfhs)$size
      if (!is.na(fsz) && fsz > 50000) { cat("  Fetched:", nfhs5_url, "size=", fsz, "\n"); break }
      cat("  Too small (error page?), size=", fsz, "from:", nfhs5_url, "\n")
      dl_res <- -1L
    }
  }
  if (!identical(dl_res, 0L)) stop("NFHS-5 download failed from all tried URLs")

  nfhs5_raw <- read.csv(tmp_nfhs, stringsAsFactors = FALSE)
  b8_fetched <- TRUE
  cat("NFHS-5 raw dims:", dim(nfhs5_raw), "\n")
  cat("NFHS-5 names (first 20):", paste(names(nfhs5_raw)[1:min(20,ncol(nfhs5_raw))], collapse=", "), "\n")

  # Find district name + state name columns for matching
  name_cols <- grep("district|state", names(nfhs5_raw), ignore.case=TRUE, value=TRUE)
  cat("Name-like cols:", paste(name_cols, collapse=", "), "\n")

  # Normalise function: lower, strip punctuation/extra spaces
  norm_str <- function(x) {
    x <- tolower(trimws(x))
    x <- gsub("[^a-z0-9 ]", " ", x)
    x <- gsub("\\s+", " ", x)
    trimws(x)
  }

  # Build a matching key from subdistrict_analysis district name + stname at district level
  dist_names <- unique(sq_csv[, c("dist_lgd","dtname","stname")])
  dist_names <- dist_names[!duplicated(dist_names$dist_lgd), ]
  dist_names$key_dtname <- norm_str(dist_names$dtname)
  dist_names$key_stname <- norm_str(dist_names$stname)

  # Detect district and state columns in NFHS-5
  dcol <- name_cols[grepl("district", name_cols, ignore.case=TRUE)][1]
  scol <- name_cols[grepl("state",    name_cols, ignore.case=TRUE)][1]
  if (is.na(dcol) || is.na(scol)) stop("Cannot identify district/state columns in NFHS-5 file")

  nfhs5_raw$key_dtname <- norm_str(nfhs5_raw[[dcol]])
  nfhs5_raw$key_stname <- norm_str(nfhs5_raw[[scol]])

  # Target cols: undernutrition/BMI, tobacco, clean cooking fuel, diabetes proxy
  # Look for these by pattern in column names
  pick_cols <- function(df, patterns) {
    unlist(lapply(patterns, function(p) {
      grep(p, names(df), ignore.case=TRUE, value=TRUE)[1]
    }))
  }
  want_patterns <- c("bmi|thin|undernut|wasting",
                     "tobacco|smoke|smok",
                     "clean.*fuel|lpg|cooking",
                     "blood.*sugar|diabetes|hyperglycaemia")
  sel_cols <- pick_cols(nfhs5_raw, want_patterns)
  sel_cols <- sel_cols[!is.na(sel_cols)]
  cat("Selected NFHS-5 cols:", paste(sel_cols, collapse=", "), "\n")

  # Join by district+state key
  nfhs5_slim <- nfhs5_raw[, c("key_dtname","key_stname", sel_cols), drop=FALSE]
  joined <- merge(dist_names, nfhs5_slim,
                  by = c("key_dtname","key_stname"), all.x = FALSE)
  b8_match_rate <- nrow(joined) / nrow(dist_names)
  cat("NFHS-5 match rate:", round(b8_match_rate*100, 1), "% (", nrow(joined), "/", nrow(dist_names), ")\n")

  if (nrow(joined) > 0 && b8_match_rate > 0.3) {
    # Recode dist_lgd -> dist_lgd_anl; pop-weighted mean for merged districts
    nfhs5_coded <- merge(joined, lookup, by = "dist_lgd", all.x = FALSE)
    nfhs5_coded <- nfhs5_coded[!is.na(nfhs5_coded$dist_lgd_anl), ]
    # Aggregate: for numeric cols, take mean (pop-weighted ideally, but pop available)
    # Merge in pop2024 for weighting
    nfhs5_coded <- merge(nfhs5_coded, pop24, by.x = "dist_lgd_anl", by.y = "dist_lgd_anl", all.x = TRUE)
    num_sel <- intersect(sel_cols, names(nfhs5_coded))
    if (length(num_sel) > 0) {
      # Weighted mean per dist_lgd_anl
      nfhs5_agg <- do.call(rbind, lapply(split(nfhs5_coded, nfhs5_coded$dist_lgd_anl), function(d) {
        w <- d$pop_2024; if (all(is.na(w))) w <- rep(1, nrow(d))
        w[is.na(w)] <- median(w, na.rm=TRUE)
        out <- data.frame(dist_lgd_anl = d$dist_lgd_anl[1])
        for (col in num_sel) {
          v <- suppressWarnings(as.numeric(d[[col]]))
          out[[col]] <- if (all(is.na(v))) NA_real_ else weighted.mean(v, w, na.rm=TRUE)
        }
        out
      }))
      nfhs5_agg$dist_lgd <- nfhs5_agg$dist_lgd_anl
      nfhs5_agg$dist_lgd_anl <- NULL
      base <- merge(base, nfhs5_agg, by = "dist_lgd", all.x = TRUE)
      b8_joined <- TRUE
      cat("NFHS-5 indicators joined to base frame\n")
    }
  } else {
    b8_skip_reason <- paste0("Match rate too low (", round(b8_match_rate*100,1), "%) or no matches")
    cat("NFHS-5: not joined -", b8_skip_reason, "\n")
  }
}, error = function(e) {
  b8_skip_reason <<- paste("NFHS-5:", conditionMessage(e))
  cat("SKIPPED - NFHS-5:", b8_skip_reason, "\n")
})

if (!b8_fetched) {
  b8_skip_reason <- ifelse(is.na(b8_skip_reason),
    "NFHS-5 download not attempted or failed",
    b8_skip_reason)
}

# =========================================================
# Finalise and save
# =========================================================
cat("\n===== FINAL ASSEMBLY =====\n")
# Ensure one row per dist_lgd (should already be)
base <- base[!duplicated(base$dist_lgd), ]
# Sort by dist_lgd
base <- base[order(base$dist_lgd), ]
cat("Final rows:", nrow(base), "\n")
cat("Final distinct dist_lgd:", length(unique(base$dist_lgd)), "\n")
cat("Columns:", paste(names(base), collapse=", "), "\n\n")

# Save
out_path <- file.path(DER, "modifiers.rds")
saveRDS(base, out_path)
cat("Saved to:", out_path, "\n\n")

# =========================================================
# VERIFICATION & COVERAGE TABLE
# =========================================================
cat("===== COVERAGE TABLE =====\n")
n_total <- nrow(base)
cov <- sapply(names(base), function(col) {
  n_present <- sum(!is.na(base[[col]]))
  pct <- round(100 * n_present / n_total, 1)
  c(n_present = n_present, n_total = n_total, pct = pct)
})
cov_df <- as.data.frame(t(cov))
cov_df$column <- rownames(cov_df)
cov_df <- cov_df[, c("column","n_present","n_total","pct")]
rownames(cov_df) <- NULL
print(cov_df)

cat("\n===== PART A SANITY CHECKS =====\n")
part_a_cols <- c("latitude","area_km2","density_2024","public_share","pct_urban_notif",
                 "region","stname","total_notif_2022_2025")
cat("Part A coverage and summary statistics:\n\n")
for (col in part_a_cols) {
  if (!col %in% names(base)) {
    cat(col, ": MISSING FROM FRAME\n")
    next
  }
  n_ok  <- sum(!is.na(base[[col]]))
  pct   <- round(100 * n_ok / n_total, 1)
  cat(sprintf("%-28s  n=%d/%d (%.1f%%)", col, n_ok, n_total, pct))
  vals <- base[[col]]
  if (is.numeric(vals)) {
    cat(sprintf("  min=%.3g  median=%.3g  max=%.3g",
                min(vals, na.rm=TRUE), median(vals, na.rm=TRUE), max(vals, na.rm=TRUE)))
  }
  cat("\n")
}

cat("\n===== PART B SOURCE STATUS =====\n")
cat("Part B.6 dtcode11 crosswalk:    DONE (built from subdistrict_analysis.csv in-hand)\n")
cat(sprintf("  Coverage: %d / %d (%.1f%%)\n",
    sum(!is.na(base$dtcode11)), nrow(base),
    100*sum(!is.na(base$dtcode11))/nrow(base)))
if (b7_fetched) {
  cat("Part B.7 Census-2011 indicators: FETCHED (but not joined - data validation needed)\n")
} else {
  cat("Part B.7 Census-2011 indicators: SKIPPED\n")
  cat("  Reason:", b7_skip_reason, "\n")
}
if (b8_joined) {
  cat(sprintf("Part B.8 NFHS-5 indicators:       JOINED (match rate: %.1f%%)\n",
              100 * b8_match_rate))
} else if (b8_fetched) {
  cat(sprintf("Part B.8 NFHS-5 indicators:       FETCHED but not joined (%.1f%% match)\n",
              100 * b8_match_rate))
  cat("  Reason:", b8_skip_reason, "\n")
} else {
  cat("Part B.8 NFHS-5 indicators:       SKIPPED\n")
  cat("  Reason:", b8_skip_reason, "\n")
}

cat("\n===== DONE =====\n")
