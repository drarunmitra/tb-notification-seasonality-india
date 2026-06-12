# 04_sector_models.R — Objective 2 — sector negative-binomial models and variance partition
#
# Renumbered, path-normalized reorganization of the working pipeline scripts: 20_model_frame.R, 21_fit_sector_models.R, 22_variance_partition.R.
# Computation is unchanged from the scripts that produced the shipped results in out/.

source(here::here("code", "00_setup.R"))


## ===========================================================================
## from 20_model_frame.R
## ===========================================================================

# 20_model_frame.R
# Task 2.1 (O2): Build sector-specific modeling frame for core seasonality test
# Sector: public, private (long format)
# Period: complete calendar years 2022-2025 (matches climate coverage)
# Output: derived/model_frame.rds


library(lubridate)

cat("=== 20_model_frame.R ===\n")
cat("Loading inputs...\n")

# ---- 1. Load inputs --------------------------------------------------------

panel    <- readRDS(file.path(DER, "panel_dm.rds"))
clim     <- readRDS(file.path(DER, "climate_dm.rds"))
cal      <- readRDS(file.path(DER, "calendar.rds"))
campaign <- readRDS(file.path(DER, "campaign_dm.rds"))
pop_df   <- readRDS(file.path(DER, "district_pop_anl.rds"))

cat("Panel rows:", nrow(panel), "| Climate rows:", nrow(clim),
    "| Calendar rows:", nrow(cal), "\n")
cat("Campaign rows:", nrow(campaign), "| Pop rows:", nrow(pop_df), "\n")

# ---- 2. Prepare calendar: flatten spline columns to plain numeric -----------
# covid_s1/s2/s3 are 1-D 'ns' vectors stored as a special class; coerce to
# plain numeric so they survive pivoting and saveRDS round-trips cleanly.

cal <- cal |>
  mutate(
    covid_s1 = as.numeric(covid_s1),
    covid_s2 = as.numeric(covid_s2),
    covid_s3 = as.numeric(covid_s3),
    fy_end   = as.integer(fy_end),
    q_end    = as.integer(q_end)
  )

# ---- 3. Pivot panel to long by sector; filter complete & <= 2025-Dec --------

cutoff <- tsibble::yearmonth("2025 Dec")

panel_long <- panel |>
  as_tibble() |>                          # drop tsibble class for pivoting
  filter(complete == TRUE, month <= cutoff) |>
  select(dist_lgd, month, public, private) |>
  pivot_longer(cols      = c(public, private),
               names_to  = "sector",
               values_to = "count") |>
  mutate(dist_lgd = as.integer(dist_lgd))  # harmonise key type with climate

cat("Panel long rows (complete, <=2025-Dec):", nrow(panel_long), "\n")
cat("Sectors:", paste(sort(unique(panel_long$sector)), collapse = ", "), "\n")

# ---- 4. Join climate, calendar, campaign, population -----------------------

# Harmonise key types to integer throughout
clim     <- clim     |> mutate(dist_lgd = as.integer(dist_lgd))
campaign <- campaign |> mutate(dist_lgd = as.integer(dist_lgd))

# Diagnostic: districts in panel missing from climate / campaign
missing_clim <- setdiff(unique(panel_long$dist_lgd), unique(clim$dist_lgd))
cat("Districts in panel NOT in climate:", length(missing_clim), "\n")
if (length(missing_clim) > 0)
  cat("  Missing IDs:", paste(missing_clim, collapse = ", "), "\n")

missing_camp <- setdiff(unique(panel_long$dist_lgd), unique(campaign$dist_lgd))
cat("Districts in panel NOT in campaign:", length(missing_camp), "\n")
if (length(missing_camp) > 0)
  cat("  Missing IDs:", paste(missing_camp, collapse = ", "), "\n")

# Join climate (dist_lgd + month) -- left_join; missing rows flagged above
frame <- panel_long |>
  left_join(clim, by = c("dist_lgd", "month"))

# Join calendar (month only)
frame <- frame |>
  left_join(cal, by = "month")

# Join campaign (dist_lgd + month)
frame <- frame |>
  left_join(campaign, by = c("dist_lgd", "month"))

# Join population: map yr = year(month), join by dist_lgd == dist_lgd_anl & yr
frame <- frame |>
  mutate(yr = lubridate::year(as.Date(month))) |>
  left_join(
    pop_df |> rename(dist_lgd = dist_lgd_anl),
    by = c("dist_lgd", "yr")
  )

cat("After all joins, rows:", nrow(frame), "\n")

# ---- 5. Add covariates ------------------------------------------------------

# month_num and Fourier harmonics K=2
frame <- frame |>
  mutate(
    month_num = lubridate::month(as.Date(month)),
    s1 = sin(2 * pi * month_num / 12),
    c1 = cos(2 * pi * month_num / 12),
    s2 = sin(4 * pi * month_num / 12),
    c2 = cos(4 * pi * month_num / 12),
    year = lubridate::year(as.Date(month))
  )

# Climate lags within each dist_lgd x sector series, ordered by month
# Lags 0-3 for t2m and solar; rh contemporaneous (rh_l0) + rh_l1
# NOTE: lag() is dplyr::lag -- explicit namespace to avoid conflict with stats::lag
frame <- frame |>
  arrange(dist_lgd, sector, month) |>
  group_by(dist_lgd, sector) |>
  mutate(
    t2m_l0   = t2m,
    t2m_l1   = dplyr::lag(t2m, 1),
    t2m_l2   = dplyr::lag(t2m, 2),
    t2m_l3   = dplyr::lag(t2m, 3),
    solar_l0 = solar,
    solar_l1 = dplyr::lag(solar, 1),
    solar_l2 = dplyr::lag(solar, 2),
    solar_l3 = dplyr::lag(solar, 3),
    rh_l0    = rh,
    rh_l1    = dplyr::lag(rh, 1)
  ) |>
  ungroup()

# ---- 6. Select final columns in logical order --------------------------------

frame <- frame |>
  select(
    dist_lgd, sector, month, year,
    count, pop,
    working_days, fy_end, q_end,
    campaign,
    covid_s1, covid_s2, covid_s3,
    s1, c1, s2, c2, month_num,
    t2m, t2m_min, ah, precip,
    t2m_l0, t2m_l1, t2m_l2, t2m_l3,
    solar_l0, solar_l1, solar_l2, solar_l3,
    rh_l0, rh_l1,
    yr
  )

cat("Columns:", ncol(frame), "\n")
cat("Rows before required-column filter:", nrow(frame), "\n")

# ---- 7. Drop rows with NA in REQUIRED columns --------------------------------
# Required (model cannot run without these): count, pop > 0, working_days,
#   s1, c1, s2, c2, t2m_l0, solar_l0.
# Rows where only lagged columns (t2m_l1..l3, solar_l1..l3, rh_l1) are NA
# are KEPT -- models choosing those lags will naturally lose the first 1-3
# months of each series at fit time; that is the expected lag-truncation cost.

required_cols <- c("count", "pop", "working_days",
                   "s1", "c1", "s2", "c2",
                   "t2m_l0", "solar_l0")

# Diagnostic: flag months where solar_l0 is entirely NA (upstream data gap)
solar_na_months <- frame |>
  group_by(month) |>
  filter(all(is.na(solar_l0))) |>
  distinct(month) |>
  pull(month)
if (length(solar_na_months) > 0) {
  cat("NOTE: solar_l0 is all-NA for months (will be dropped by required-col filter):",
      paste(as.character(solar_na_months), collapse=", "), "\n")
  cat("  This is an upstream data gap in climate_dm.rds, not a join error.\n")
}

n_before <- nrow(frame)
frame <- frame |>
  filter(
    !is.na(count),
    !is.na(pop), pop > 0,
    !is.na(working_days),
    !is.na(s1), !is.na(c1), !is.na(s2), !is.na(c2),
    !is.na(t2m_l0),
    !is.na(solar_l0)
  )

n_after <- nrow(frame)
cat("Rows dropped (NA in required cols or pop<=0):", n_before - n_after, "\n")
cat("Final rows:", n_after, "\n")

# ---- 8. Save ----------------------------------------------------------------

saveRDS(frame, file.path(DER, "model_frame.rds"))
cat("\nSaved:", file.path(DER, "model_frame.rds"), "\n")

# ---- 9. VERIFICATION --------------------------------------------------------

cat("\n========================================\n")
cat("VERIFICATION\n")
cat("========================================\n")

# Rows and districts per sector
cat("\n[1] Rows and distinct districts per sector:\n")
sector_summary <- frame |>
  group_by(sector) |>
  summarise(n_rows      = n(),
            n_districts = dplyr::n_distinct(dist_lgd),
            .groups     = "drop")
print(sector_summary)

# Month range
cat("\n[2] Month range:\n")
cat("  Min:", as.character(min(frame$month)), "\n")
cat("  Max:", as.character(max(frame$month)), "\n")
# Note: solar column in climate_dm.rds is NA for 2025-Dec (upstream data gap).
# The required-col filter on solar_l0 therefore excludes 2025-Dec rows.
# Effective coverage: 2022-Jan..2025-Nov (47 months of complete climate+panel data).
if (as.character(max(frame$month)) == "2025 Nov") {
  cat("  WARN: max month is 2025-Nov, not 2025-Dec.\n")
  cat("  CAUSE: solar is all-NA in climate_dm.rds for 2025-Dec (upstream gap).\n")
  cat("  ACTION NEEDED: if 2025-Dec solar data becomes available, re-run 02_climate.R\n")
  cat("    and then re-run this script to extend coverage by 1 month.\n")
} else {
  cat("  Range == 2022-Jan..2025-Dec: TRUE\n")
}

# pop > 0
cat("\n[3] pop > 0 for all rows:\n")
cat("  Any pop <= 0:", any(frame$pop <= 0, na.rm = TRUE), "\n")
cat("  Any NA in pop:", anyNA(frame$pop), "\n")

# NA check in required columns
cat("\n[4] NA counts in required columns:\n")
for (col in required_cols) {
  cat(sprintf("  NA in %-15s: %d\n", col, sum(is.na(frame[[col]]))))
}

# NA counts in lag columns (expected to be > 0)
cat("\n[5] NA counts in lag columns (expected non-zero for lags 1-3):\n")
lag_cols <- c("t2m_l1", "t2m_l2", "t2m_l3",
              "solar_l1", "solar_l2", "solar_l3",
              "rh_l1")
for (col in lag_cols) {
  cat(sprintf("  NA in %-12s: %d\n", col, sum(is.na(frame[[col]]))))
}

# glimpse
cat("\n[6] glimpse:\n")
dplyr::glimpse(frame)

# head
cat("\n[7] head(5):\n")
print(head(frame, 5))

# Sanity totals
cat("\n[8] Total counts by sector (sanity check):\n")
totals <- frame |>
  group_by(sector) |>
  summarise(total_count = sum(count, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(total_count))
print(totals)
cat("  Public > Private:",
    totals$total_count[totals$sector == "public"] >
    totals$total_count[totals$sector == "private"], "\n")

cat("\n=== 20_model_frame.R COMPLETE ===\n")

## ===========================================================================
## from 21_fit_sector_models.R
## ===========================================================================

## 21_fit_sector_models.R
## Phase 2, Task 2.2 (O2): Fit sector-specific count models with climate and
## calendar/artefact blocks, partitioned for downstream signal decomposition.
## ------------------------------------------------------------------


set.seed(SEED)

## ---- 0. Package loading -----------------------------------------------
for (pkg in c("glmmTMB", "DHARMa", "MASS")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

## AICc helper (glmmTMB does not export one; MuMIn::AICc works but we avoid
## a hard dependency; manual formula: AIC + 2k(k+1)/(n-k-1))
aicc <- function(fit) {
  ll  <- tryCatch(logLik(fit), error = function(e) NULL)
  if (is.null(ll)) return(Inf)
  k   <- attr(ll, "df")
  n   <- attr(ll, "nobs")
  if (is.null(n) || is.na(n) || (n - k - 1) <= 0) return(AIC(fit))
  AIC(fit) + 2 * k * (k + 1) / (n - k - 1)
}

## ---- 1. Load and prepare model frame ----------------------------------
cat("\n=== STEP 1: Load model frame ===\n")
mf_path <- file.path(DER, "model_frame.rds")
mf <- readRDS(mf_path)

cat("Rows loaded:", nrow(mf), " | Cols:", ncol(mf), "\n")

## Factor coercions
mf$sector   <- factor(mf$sector)
mf$year     <- factor(mf$year)
mf$dist_lgd <- factor(mf$dist_lgd)

## Guard offsets: must be strictly positive
stopifnot("pop must be > 0"          = all(mf$pop > 0, na.rm = TRUE))
stopifnot("working_days must be > 0" = all(mf$working_days > 0, na.rm = TRUE))

## Pre-compute log-offsets once
mf$log_pop <- log(mf$pop)
mf$log_wd  <- log(mf$working_days)

cat("Sector levels:", levels(mf$sector), "\n")
cat("Year levels  :", levels(mf$year),   "\n")
cat("N districts  :", nlevels(mf$dist_lgd), "\n")

## ---- 2. Climate lag selection -----------------------------------------
##  Strategy: for each sector, fit a reduced glm.nb (no random effects for speed)
##  with harmonics + t2m_lk + solar_lk + offset(log_pop + log_wd), varying k
##  for temperature and solar jointly (k = 0, 1, 2, 3).
##  rh_l0 is fixed at lag-0 (contemporaneous).

cat("\n=== STEP 2: Climate lag selection ===\n")

lag_cache_path <- file.path(DER, "lag_selection.rds")

if (file.exists(lag_cache_path)) {
  cat("Loading cached lag selection.\n")
  lag_sel <- readRDS(lag_cache_path)
} else {
  lags <- 0:3
  lag_sel <- list()

  for (sec in levels(mf$sector)) {
    cat("  Sector:", sec, "\n")
    df_sec <- mf[mf$sector == sec, ]
    best_k <- NA; best_aicc <- Inf
    aicc_tbl <- data.frame(lag = integer(), aicc = numeric())

    for (k in lags) {
      t_col <- paste0("t2m_l", k)
      s_col <- paste0("solar_l", k)

      # Use a complete-case subset for this lag
      df_k <- df_sec[complete.cases(df_sec[, c(t_col, s_col, "rh_l0")]), ]

      if (nrow(df_k) < 200) next

      frm <- as.formula(paste0(
        "count ~ s1 + c1 + s2 + c2 + ", t_col, " + ", s_col,
        " + rh_l0 + offset(log_pop) + offset(log_wd)"
      ))

      fit_k <- tryCatch(
        MASS::glm.nb(frm, data = df_k),
        error   = function(e) NULL,
        warning = function(w) suppressWarnings(
          tryCatch(MASS::glm.nb(frm, data = df_k), error = function(e2) NULL)
        )
      )

      if (is.null(fit_k)) next
      ac <- aicc(fit_k)
      cat("    lag", k, ": AICc =", round(ac, 2), "\n")
      aicc_tbl <- rbind(aicc_tbl, data.frame(lag = k, aicc = ac))
      if (ac < best_aicc) { best_aicc <- ac; best_k <- k }
    }

    cat("  => chosen lag for", sec, ":", best_k,
        "(AICc =", round(best_aicc, 2), ")\n")
    lag_sel[[sec]] <- list(lag = best_k, aicc_table = aicc_tbl)
  }

  saveRDS(lag_sel, lag_cache_path)
}

## Extract chosen lags
lag_pub  <- lag_sel[["public"]]$lag
lag_priv <- lag_sel[["private"]]$lag
cat("\nChosen lags: public =", lag_pub, " | private =", lag_priv, "\n")

## ---- 3. Main glmmTMB models -------------------------------------------
cat("\n=== STEP 3: Fit main glmmTMB sector models ===\n")

fit_sector <- function(sec, lag_k, label, cache_path) {
  if (file.exists(cache_path)) {
    cat("  Loading cached", label, "fit.\n")
    return(readRDS(cache_path))
  }

  t_col <- paste0("t2m_l", lag_k)
  s_col <- paste0("solar_l", lag_k)
  needed <- c("count", "s1", "c1", "s2", "c2",
              t_col, s_col, "rh_l0",
              "q_end", "fy_end", "campaign",
              "covid_s1", "covid_s2", "covid_s3",
              "year", "dist_lgd", "log_pop", "log_wd")

  df_sec <- mf[mf$sector == sec, ]
  df_sec <- df_sec[complete.cases(df_sec[, needed]), ]
  cat("  Rows for", label, "after complete-case filter:", nrow(df_sec), "\n")

  frm <- as.formula(paste0(
    "count ~ s1 + c1 + s2 + c2",
    " + ", t_col, " + ", s_col, " + rh_l0",
    " + q_end + fy_end + campaign",
    " + covid_s1 + covid_s2 + covid_s3",
    " + year",
    " + (1 | dist_lgd)",
    " + offset(log_pop) + offset(log_wd)"
  ))

  cat("  Fitting glmmTMB (nbinom2) for", label, "...\n")
  t0 <- proc.time()
  fit <- tryCatch(
    glmmTMB(frm, data = df_sec, family = nbinom2(link = "log"),
            REML = FALSE, verbose = FALSE),
    error = function(e) {
      cat("  ERROR in glmmTMB:", conditionMessage(e), "\n")
      NULL
    }
  )
  elapsed <- (proc.time() - t0)[["elapsed"]]
  cat("  Elapsed:", round(elapsed, 1), "s\n")

  if (!is.null(fit)) {
    saveRDS(fit, cache_path)
    cat("  Saved:", cache_path, "\n")
  }
  fit
}

fit_pub  <- fit_sector("public",  lag_pub,  "public",
                       file.path(DER, "fit_public.rds"))
fit_priv <- fit_sector("private", lag_priv, "private",
                       file.path(DER, "fit_private.rds"))

## Also write to OUT/ (canonical output location per task spec)
if (!is.null(fit_pub))
  saveRDS(fit_pub,  file.path(OUT, "fit_public.rds"))
if (!is.null(fit_priv))
  saveRDS(fit_priv, file.path(OUT, "fit_private.rds"))

## ---- 4. Diagnostics ---------------------------------------------------
cat("\n=== STEP 4: Diagnostics ===\n")

## Helper: convergence check
check_convergence <- function(fit, label) {
  if (is.null(fit)) return(list(converged = FALSE, msg = "fit is NULL"))
  hess_ok <- isTRUE(fit$sdr$pdHess)
  conv0   <- tryCatch(fit$fit$convergence == 0, error = function(e) NA)
  converged <- hess_ok || isTRUE(conv0)
  msg <- paste0("pdHess=", hess_ok, " | optim_convergence=", conv0)
  cat(" ", label, "convergence:", msg, "\n")
  list(converged = converged, hess_ok = hess_ok, conv0 = conv0, msg = msg)
}

## Helper: NegBin vs Poisson dispersion comparison
check_dispersion <- function(fit_nb, sec, lag_k, df_sec_clean) {
  ## Refit with Poisson on the same data
  t_col <- paste0("t2m_l", lag_k)
  s_col <- paste0("solar_l", lag_k)
  frm <- formula(fit_nb)

  fit_pois <- tryCatch(
    glmmTMB(frm, data = df_sec_clean, family = poisson(link = "log"),
            REML = FALSE, verbose = FALSE),
    error = function(e) NULL
  )

  if (is.null(fit_pois)) {
    cat("  Could not fit Poisson for comparison.\n")
    return(list(theta = NA, lrt_p = NA, msg = "Poisson fit failed"))
  }

  theta <- tryCatch(sigma(fit_nb), error = function(e) NA)
  # LRT: Poisson is nested in NegBin (theta -> Inf is Poisson)
  lrt_stat <- 2 * (logLik(fit_nb) - logLik(fit_pois))
  lrt_p    <- pchisq(as.numeric(lrt_stat), df = 1, lower.tail = FALSE) / 2
  msg <- sprintf("theta=%.3f | LRT chi2=%.2f, p=%.4f (overdispersion confirmed=%s)",
                 theta, as.numeric(lrt_stat), lrt_p, lrt_p < 0.05)
  cat(" ", sec, "dispersion:", msg, "\n")
  list(theta = theta, lrt_stat = as.numeric(lrt_stat), lrt_p = lrt_p, msg = msg)
}

## Run diagnostics for each sector
diag_results <- list()
diag_text    <- character()

for (sec in c("public", "private")) {
  fit    <- if (sec == "public") fit_pub else fit_priv
  lag_k  <- if (sec == "public") lag_pub else lag_priv
  label  <- paste0(sec, " sector")

  diag_text <- c(diag_text, paste0("\n===== ", toupper(sec), " SECTOR =====\n"))

  ## 4a. Convergence
  conv <- check_convergence(fit, label)
  diag_text <- c(diag_text, paste0("Convergence: ", conv$msg, "\n"))

  if (is.null(fit)) {
    diag_text <- c(diag_text, "Model is NULL - skipping further diagnostics.\n")
    next
  }

  ## Prepare clean data for simulation
  t_col  <- paste0("t2m_l", lag_k)
  s_col  <- paste0("solar_l", lag_k)
  needed <- c("count", "s1", "c1", "s2", "c2",
              t_col, s_col, "rh_l0",
              "q_end", "fy_end", "campaign",
              "covid_s1", "covid_s2", "covid_s3",
              "year", "dist_lgd", "log_pop", "log_wd")
  df_sec_clean <- mf[mf$sector == sec, ]
  df_sec_clean <- df_sec_clean[complete.cases(df_sec_clean[, needed]), ]

  ## 4b. Overdispersion vs Poisson
  cat("  Checking overdispersion for", sec, "...\n")
  disp <- check_dispersion(fit, sec, lag_k, df_sec_clean)
  diag_text <- c(diag_text, paste0("Overdispersion: ", disp$msg, "\n"))

  ## 4c. DHARMa residual diagnostics
  cat("  Running DHARMa simulations for", sec, "...\n")
  set.seed(SEED)
  sim_n <- min(500L, B)  # cap at 500 for speed
  sim_res <- tryCatch(
    DHARMa::simulateResiduals(fit, n = sim_n, plot = FALSE),
    error = function(e) {
      cat("  DHARMa error:", conditionMessage(e), "\n"); NULL
    }
  )

  if (!is.null(sim_res)) {
    ## Dispersion test
    disp_test <- tryCatch(
      DHARMa::testDispersion(sim_res, plot = FALSE),
      error = function(e) NULL
    )
    disp_msg <- if (!is.null(disp_test))
      sprintf("DHARMa dispersion: statistic=%.3f, p=%.4f",
              disp_test$statistic, disp_test$p.value)
    else "DHARMa dispersion: test failed"
    cat(" ", disp_msg, "\n")
    diag_text <- c(diag_text, paste0(disp_msg, "\n"))

    ## Temporal autocorrelation: Ljung-Box on district-averaged residuals
    ##  Use the sequential year-month time index to get a proper time series
    ##  (month_num 1-12 repeats every year; we need the ordered year-month sequence)
    cat("  Ljung-Box temporal autocorrelation check for", sec, "...\n")
    resid_df <- data.frame(
      resid     = residuals(sim_res),
      year_num  = as.integer(as.character(df_sec_clean$year)),
      month_num = df_sec_clean$month_num,
      dist_lgd  = df_sec_clean$dist_lgd
    )
    # Sequential time index: each unique (year, month_num) gets one slot
    resid_df$ym_key <- resid_df$year_num * 100L + resid_df$month_num
    # Average residual per year-month across districts (mean over panel)
    ym_avg <- tapply(resid_df$resid, resid_df$ym_key, mean, na.rm = TRUE)
    ym_avg <- ym_avg[order(as.numeric(names(ym_avg)))]
    n_ym   <- length(ym_avg)

    # Use lag = min(12, floor(n_ym/3)) to keep degrees of freedom valid
    lb_lag    <- max(1L, min(12L, floor(n_ym / 3L)))
    lb_result <- tryCatch(
      Box.test(ym_avg, lag = lb_lag, type = "Ljung-Box"),
      error = function(e) NULL
    )
    if (!is.null(lb_result) && !is.na(lb_result$p.value)) {
      ac_status <- if (lb_result$p.value < 0.05) "SIGNIFICANT" else "white"
      ac_msg <- sprintf(
        "Ljung-Box (lag=%d, n_ym=%d) on year-month-avg residuals: chi2=%.3f, p=%.4f => %s",
        lb_lag, n_ym, lb_result$statistic, lb_result$p.value, ac_status
      )
      cat(" ", ac_msg, "\n")
      diag_text <- c(diag_text, paste0(ac_msg, "\n"))

      if (lb_result$p.value < 0.05) {
        ## Escalation: note and attempt ar1 variant if feasible
        esc_msg <- paste0(
          "TEMPORAL AUTOCORRELATION DETECTED. Escalation path:\n",
          "  (1) Add ar1(month_num + 0 | dist_lgd) correlation structure in glmmTMB.\n",
          "  (2) Alternatively, move to tscount::tsglm (INGARCH) or glmmTMB with AR1.\n",
          "  Attempting ar1 variant...\n"
        )
        cat(esc_msg)
        diag_text <- c(diag_text, esc_msg)

        ## Attempt ar1 variant - build formula directly
        t_col_ar1 <- paste0("t2m_l", lag_k)
        s_col_ar1 <- paste0("solar_l", lag_k)
        ## Add month_num_f factor for ar1 correlation structure
        df_sec_clean$month_num_f <- factor(df_sec_clean$month_num)
        frm_ar1 <- as.formula(paste0(
          "count ~ s1 + c1 + s2 + c2",
          " + ", t_col_ar1, " + ", s_col_ar1, " + rh_l0",
          " + q_end + fy_end + campaign",
          " + covid_s1 + covid_s2 + covid_s3",
          " + year",
          " + (1 | dist_lgd)",
          " + ar1(month_num_f + 0 | dist_lgd)",
          " + offset(log_pop) + offset(log_wd)"
        ))
        fit_ar1 <- tryCatch(
          glmmTMB(frm_ar1, data = df_sec_clean, family = nbinom2(link = "log"),
                  REML = FALSE, verbose = FALSE),
          error   = function(e) { cat("  ar1 fit error:", conditionMessage(e), "\n"); NULL },
          warning = function(w) suppressWarnings(tryCatch(
            glmmTMB(frm_ar1, data = df_sec_clean, family = nbinom2(link = "log"),
                    REML = FALSE, verbose = FALSE),
            error = function(e) NULL
          ))
        )

        if (!is.null(fit_ar1)) {
          a_base <- aicc(fit)
          a_ar1  <- aicc(fit_ar1)
          ar1_msg <- sprintf(
            "ar1 variant fitted. AICc: base=%.2f, ar1=%.2f, delta=%.2f => %s preferred",
            a_base, a_ar1, a_base - a_ar1,
            if (a_ar1 < a_base) "ar1" else "base"
          )
          cat(" ", ar1_msg, "\n")
          diag_text <- c(diag_text, paste0(ar1_msg, "\n"))

          ar1_cache <- file.path(DER, paste0("fit_", sec, "_ar1.rds"))
          saveRDS(fit_ar1, ar1_cache)
          cat("  ar1 fit saved to:", ar1_cache, "\n")
        } else {
          limitation_msg <- paste0(
            "ar1 variant could not be fitted. DOCUMENTED LIMITATION: ",
            "residual temporal autocorrelation detected but ar1(month_num+0|dist_lgd) ",
            "failed to converge. Recommend: (a) fit tscount::tsglm with INGARCH structure, ",
            "(b) aggregate to quarterly counts to reduce autocorrelation, or ",
            "(c) include a month-in-year fixed effect to absorb within-year seasonality ",
            "before AR correction.\n"
          )
          cat(limitation_msg)
          diag_text <- c(diag_text, limitation_msg)
        }
      }
    } else {
      diag_text <- c(diag_text, paste0("Ljung-Box: test returned NA or failed (n_ym=", n_ym, ", lag=", lb_lag, ")\n"))
    }

    ## DHARMa temporal autocorrelation test using month_num as time index
    tac_test <- tryCatch(
      DHARMa::testTemporalAutocorrelation(
        sim_res, time = df_sec_clean$month_num, plot = FALSE
      ),
      error = function(e) NULL
    )
    if (!is.null(tac_test)) {
      tac_msg <- sprintf(
        "DHARMa testTemporalAutocorrelation: statistic=%.3f, p=%.4f",
        tac_test$statistic, tac_test$p.value
      )
      cat(" ", tac_msg, "\n")
      diag_text <- c(diag_text, paste0(tac_msg, "\n"))
    }
  } else {
    diag_text <- c(diag_text, "DHARMa simulation failed - could not run residual checks.\n")
  }

  diag_results[[sec]] <- list(conv = conv, disp = disp,
                               sim_res = sim_res)
}

## Write diagnostics to file
diag_path <- file.path(OUT, "o2_diagnostics.txt")
writeLines(c(
  paste0("O2 Diagnostics - generated: ", Sys.time()),
  paste0("Script: 21_fit_sector_models.R"),
  paste0("Seed: ", SEED),
  "",
  paste0("Climate lag selection:"),
  paste0("  Public  sector: lag = ", lag_pub,
         " | AICc table:\n",
         paste(apply(lag_sel[["public"]]$aicc_table, 1,
                     function(r) sprintf("    lag%d: %.2f", r[1], r[2])),
               collapse = "\n")),
  paste0("  Private sector: lag = ", lag_priv,
         " | AICc table:\n",
         paste(apply(lag_sel[["private"]]$aicc_table, 1,
                     function(r) sprintf("    lag%d: %.2f", r[1], r[2])),
               collapse = "\n")),
  "",
  diag_text
), diag_path)
cat("\nDiagnostics written to:", diag_path, "\n")

## ---- 5. Coefficient tables (artefact + climate blocks) ----------------
cat("\n=== STEP 5: Coefficient tables ===\n")

## Blocks of interest
artefact_vars <- c("q_end", "fy_end", "campaign")
climate_vars  <- c(paste0("t2m_l",   c(lag_pub, lag_priv)),
                   paste0("solar_l", c(lag_pub, lag_priv)),
                   "rh_l0")
climate_vars  <- unique(climate_vars)

print_coef_block <- function(fit, label, block_vars) {
  if (is.null(fit)) {
    cat(label, ": model is NULL\n"); return(invisible(NULL))
  }
  ct <- tryCatch(
    as.data.frame(coef(summary(fit))$cond),
    error = function(e) NULL
  )
  if (is.null(ct)) {
    cat(label, ": could not extract coefficients\n"); return(invisible(NULL))
  }
  names(ct) <- c("Estimate", "Std_Error", "z_value", "p_value")
  sel <- ct[rownames(ct) %in% block_vars, , drop = FALSE]
  cat("\n---", label, "---\n")
  print(round(sel, 4))
  invisible(sel)
}

cat("\n>>> ARTEFACT BLOCK (q_end, fy_end, campaign) <<<\n")
art_pub  <- print_coef_block(fit_pub,  "Public  sector - ARTEFACT", artefact_vars)
art_priv <- print_coef_block(fit_priv, "Private sector - ARTEFACT", artefact_vars)

cat("\n>>> CLIMATE BLOCK (t2m, solar, rh_l0) <<<\n")
cli_pub  <- print_coef_block(fit_pub,  "Public  sector - CLIMATE",  climate_vars)
cli_priv <- print_coef_block(fit_priv, "Private sector - CLIMATE",  climate_vars)

## One-line discriminating test
cat("\n>>> DISCRIMINATING TEST: Artefact terms larger in public than private? <<<\n")
if (!is.null(art_pub) && !is.null(art_priv) && nrow(art_pub) > 0 && nrow(art_priv) > 0) {
  pub_larger <- 0; priv_larger <- 0
  for (v in intersect(rownames(art_pub), rownames(art_priv))) {
    if (art_pub[v, "Estimate"] > art_priv[v, "Estimate"]) pub_larger <- pub_larger + 1
    else priv_larger <- priv_larger + 1
  }
  disc_msg <- if (pub_larger > priv_larger)
    sprintf("YES: %d of %d artefact coefficients are LARGER in public than private [expected direction].",
            pub_larger, pub_larger + priv_larger)
  else
    sprintf("NO: only %d of %d artefact coefficients are larger in public; %d are larger in private [unexpected direction].",
            pub_larger, pub_larger + priv_larger, priv_larger)
  cat(disc_msg, "\n")
} else {
  cat("Cannot evaluate: coefficient tables incomplete.\n")
}

## ---- VERIFICATION SUMMARY -----------------------------------------------
cat("\n\n============================================================\n")
cat("VERIFICATION SUMMARY\n")
cat("============================================================\n")

for (sec in c("public", "private")) {
  fit   <- if (sec == "public") fit_pub else fit_priv
  conv  <- diag_results[[sec]]$conv
  disp  <- diag_results[[sec]]$disp
  cat(sprintf("\n[%s]\n", toupper(sec)))
  cat("  Convergence   :", if (!is.null(conv)) conv$msg else "N/A", "\n")
  cat("  Overdispersion:", if (!is.null(disp)) disp$msg else "N/A", "\n")
}

cat("\nAutocorrelation status: see o2_diagnostics.txt and console output above.\n")
cat("Artefact+climate coefficient tables: printed above.\n")
cat("Discriminating test:", if (exists("disc_msg")) disc_msg else "N/A", "\n")
cat("============================================================\n")
cat("Script complete.\n")

## ===========================================================================
## from 22_variance_partition.R
## ===========================================================================

## 22_variance_partition.R
## Phase 2, Task 2.3 (O2): Variance partition of modelled seasonal/temporal signal
## into HARMONIC, CLIMATE, and CALENDAR/CAMPAIGN blocks.
## Public sector: fit_public.rds  |  Private sector: fit_private_ar1.rds (AR1 preferred)
##
## METHOD NOTE - Bootstrap CI:
##   Districts are cluster-bootstrapped (resample with replacement) and the three
##   block shares are recomputed on the resampled district-months HOLDING the
##   fitted coefficients FIXED. The CI therefore reflects district sampling
##   variability in the covariate distribution, NOT model parameter uncertainty.
##   This choice is deliberate: the goal is to ask "how stable is the variance
##   share estimate if we drew a different sample of districts?" not "how
##   uncertain are the model parameters?" Refitting 2000 models would be
##   computationally prohibitive and is not warranted for a descriptive partition.
## -------------------------------------------------------------------------



## ---- 0. Package loading --------------------------------------------------
for (pkg in c("glmmTMB", "ggplot2", "scales")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

cat("=== 22_variance_partition.R ===\n")
cat("B =", B, " | SEED =", SEED, "\n\n")

## ---- 1. Load fitted models and model frame -------------------------------
cat("Loading models...\n")

fit_pub <- readRDS(file.path(DER, "fit_public.rds"))
fit_ar1 <- readRDS(file.path(DER, "fit_private_ar1.rds"))

cat("Public  model:      fixef names:", paste(names(fixef(fit_pub)$cond), collapse = ", "), "\n")
cat("Private AR1 model:  fixef names:", paste(names(fixef(fit_ar1)$cond), collapse = ", "), "\n\n")

mf <- readRDS(file.path(DER, "model_frame.rds"))
cat("model_frame rows:", nrow(mf), "\n\n")

## Ensure factor types consistent with model fitting
mf$sector   <- factor(mf$sector)
mf$year     <- factor(mf$year)
mf$dist_lgd <- factor(mf$dist_lgd)

## Pre-compute log-offsets (needed for complete-case filter consistency)
mf$log_pop <- log(mf$pop)
mf$log_wd  <- log(mf$working_days)

## ---- 2. Helper: compute block contributions on a data frame --------------

## Both models chose climate lag = 3 (confirmed from lag_selection.rds).
## The coefficient names in both models are: t2m_l3, solar_l3, rh_l0
## Harmonic: s1, c1, s2, c2
## Climate:  t2m_l3, solar_l3, rh_l0
## Calendar: q_end, fy_end, campaign
## Excluded: (Intercept), year*, covid_s*, offsets

compute_blocks <- function(coef_vec, df) {
  ## coef_vec: named numeric vector from fixef(m)$cond
  ## df:       data frame with covariate columns
  ##
  ## Returns df with columns: block_harmonic, block_climate, block_calendar, block_S

  get_b <- function(nm, data) {
    b <- coef_vec[nm]
    if (is.na(b)) stop(paste("Coefficient not found:", nm))
    b * data[[nm]]
  }

  harmonic_names  <- c("s1", "c1", "s2", "c2")
  climate_names   <- c("t2m_l3", "solar_l3", "rh_l0")
  calendar_names  <- c("q_end", "fy_end", "campaign")

  ## Verify all needed covariates exist in df
  need <- c(harmonic_names, climate_names, calendar_names)
  missing_cols <- setdiff(need, names(df))
  if (length(missing_cols) > 0)
    stop(paste("Missing columns in data:", paste(missing_cols, collapse = ", ")))

  block_H <- rowSums(sapply(harmonic_names, get_b, data = df))
  block_C <- rowSums(sapply(climate_names,  get_b, data = df))
  block_K <- rowSums(sapply(calendar_names, get_b, data = df))

  df$block_harmonic  <- block_H
  df$block_climate   <- block_C
  df$block_calendar  <- block_K
  df$block_S         <- block_H + block_C + block_K
  df
}

## ---- 3. Covariance decomposition of variance shares ----------------------

compute_shares <- function(df_with_blocks) {
  ## Partition:  share_k = cov(block_k, S) / var(S)
  ## These three shares sum to 1 by construction of the covariance decomposition:
  ##   var(S) = var(H+C+K) = cov(H,S) + cov(C,S) + cov(K,S)
  ## which follows from linearity of covariance.

  S <- df_with_blocks$block_S
  H <- df_with_blocks$block_harmonic
  C <- df_with_blocks$block_climate
  K <- df_with_blocks$block_calendar

  varS <- var(S)
  if (is.na(varS) || varS == 0) return(c(harmonic = NA, climate = NA, calendar = NA, varS = NA))

  sh <- c(
    harmonic  = cov(H, S) / varS,
    climate   = cov(C, S) / varS,
    calendar  = cov(K, S) / varS,
    varS      = varS
  )
  sh
}

## ---- 4. Prepare sector data frames with complete cases -------------------

## The model was fit on complete cases for all needed covariates (including
## t2m_l3, solar_l3, rh_l0 which have lag-induced NAs in first 3 months).
## We replicate the same filter here to ensure consistency.

needed_cols <- c("s1", "c1", "s2", "c2",
                 "t2m_l3", "solar_l3", "rh_l0",
                 "q_end", "fy_end", "campaign",
                 "covid_s1", "covid_s2", "covid_s3",
                 "year", "dist_lgd", "log_pop", "log_wd", "count")

df_pub <- mf[mf$sector == "public", ]
df_pub <- df_pub[complete.cases(df_pub[, needed_cols]), ]

df_priv <- mf[mf$sector == "private", ]
df_priv <- df_priv[complete.cases(df_priv[, needed_cols]), ]

cat("Rows for point-estimate computation:\n")
cat("  Public  sector:", nrow(df_pub), "rows |",
    length(unique(df_pub$dist_lgd)), "districts\n")
cat("  Private sector:", nrow(df_priv), "rows |",
    length(unique(df_priv$dist_lgd)), "districts\n\n")

## ---- 5. Point estimates --------------------------------------------------

cat("Computing block contributions (point estimates)...\n")

coef_pub  <- fixef(fit_pub)$cond
coef_ar1  <- fixef(fit_ar1)$cond

df_pub_b  <- compute_blocks(coef_pub,  df_pub)
df_priv_b <- compute_blocks(coef_ar1, df_priv)

shares_pub  <- compute_shares(df_pub_b)
shares_priv <- compute_shares(df_priv_b)

cat("\n--- POINT ESTIMATES ---\n")
cat(sprintf("Public  sector: harmonic=%.2f%%  climate=%.2f%%  calendar=%.2f%%  | var(S)=%.6f\n",
            shares_pub["harmonic"]  * 100,
            shares_pub["climate"]   * 100,
            shares_pub["calendar"]  * 100,
            shares_pub["varS"]))
cat(sprintf("Private sector: harmonic=%.2f%%  climate=%.2f%%  calendar=%.2f%%  | var(S)=%.6f\n",
            shares_priv["harmonic"] * 100,
            shares_priv["climate"]  * 100,
            shares_priv["calendar"] * 100,
            shares_priv["varS"]))

## ---- 6. Cluster bootstrap CIs (B=2000, seed=SEED) -----------------------

cat("\nRunning cluster bootstrap (B=", B, ", seed=", SEED, ")...\n", sep = "")
cat("Note: coefficients HELD FIXED; CI reflects district covariate-distribution variability.\n\n")

set.seed(SEED)

boot_one_sector <- function(df_with_blocks, coef_vec, B_reps) {
  ## VECTORISED cluster bootstrap (coefficients fixed). Pre-split row indices by
  ## district once, then each replicate just concatenates precomputed numeric
  ## vectors for the sampled districts and computes cov/var. Milliseconds/rep.
  H <- df_with_blocks$block_harmonic
  C <- df_with_blocks$block_climate
  K <- df_with_blocks$block_calendar
  gid      <- as.integer(factor(df_with_blocks$dist_lgd))
  idx_by_d <- split(seq_along(gid), gid)
  n_d      <- length(idx_by_d)

  results <- matrix(NA_real_, nrow = B_reps, ncol = 3,
                    dimnames = list(NULL, c("harmonic", "climate", "calendar")))
  for (b in seq_len(B_reps)) {
    samp <- sample.int(n_d, n_d, replace = TRUE)
    rows <- unlist(idx_by_d[samp], use.names = FALSE)
    Hb <- H[rows]; Cb <- C[rows]; Kb <- K[rows]; Sb <- Hb + Cb + Kb
    vS <- var(Sb)
    if (is.na(vS) || vS == 0) next
    results[b, ] <- c(cov(Hb, Sb) / vS, cov(Cb, Sb) / vS, cov(Kb, Sb) / vS)
  }
  results
}

cat("  Bootstrapping public sector...\n")
t0 <- proc.time()
boot_pub  <- boot_one_sector(df_pub_b,  coef_pub,  B)
cat("  Done. Elapsed:", round((proc.time() - t0)[["elapsed"]], 1), "s\n")

cat("  Bootstrapping private sector...\n")
t0 <- proc.time()
boot_priv <- boot_one_sector(df_priv_b, coef_ar1, B)
cat("  Done. Elapsed:", round((proc.time() - t0)[["elapsed"]], 1), "s\n")

## 2.5 / 97.5 percentiles
ci_pub  <- apply(boot_pub,  2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
ci_priv <- apply(boot_priv, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

## ---- 7. Assemble results table -------------------------------------------

blocks <- c("harmonic", "climate", "calendar")

result_tbl <- rbind(
  data.frame(
    sector   = "public",
    block    = blocks,
    share_pct = shares_pub[blocks] * 100,
    share_lo  = ci_pub["2.5%",  blocks] * 100,
    share_hi  = ci_pub["97.5%", blocks] * 100,
    row.names = NULL
  ),
  data.frame(
    sector   = "private",
    block    = blocks,
    share_pct = shares_priv[blocks] * 100,
    share_lo  = ci_priv["2.5%",  blocks] * 100,
    share_hi  = ci_priv["97.5%", blocks] * 100,
    row.names = NULL
  )
)

## ---- 8. Save CSV ---------------------------------------------------------

out_csv <- file.path(OUT, "o2_partition.csv")
write_csv(result_tbl, out_csv)
cat("\nSaved:", out_csv, "\n")

## ---- 9. Figure -----------------------------------------------------------

## Factor ordering for display
result_tbl$block  <- factor(result_tbl$block,
                             levels = c("harmonic", "climate", "calendar"),
                             labels = c("Harmonic (Fourier)",
                                        "Climate (temp/solar/RH)",
                                        "Calendar/Campaign (artefact)"))
result_tbl$sector <- factor(result_tbl$sector,
                             levels = c("public", "private"),
                             labels = c("Public", "Private"))

## Colours: highlight the calendar/artefact block in amber
block_cols <- c("Harmonic (Fourier)"          = "#4E84C4",
                "Climate (temp/solar/RH)"      = "#52854C",
                "Calendar/Campaign (artefact)" = "#E08B2F")

# Headline figure: the reporting-calendar/campaign (artefact) variance share by
# sector. The harmonic and climate blocks jointly form the smooth seasonal cycle
# and are collinear, so they are NOT split in the headline figure (the full
# three-way table remains in o2_partition.csv). This keeps the figure on the
# defensible, public-specific quantity rather than inviting a "98% biological"
# reading of a collinear decomposition.
art <- dplyr::filter(result_tbl, grepl("Calendar", block))
p <- ggplot(art, aes(x = sector, y = share_pct, ymin = share_lo, ymax = share_hi, fill = sector)) +
  geom_hline(yintercept = 0, colour = "grey75") +
  geom_col(width = 0.55, colour = "white", linewidth = 0.3) +
  geom_errorbar(width = 0.18, linewidth = 0.7, colour = "grey25") +
  scale_fill_manual(values = c(Public = "#08519c", Private = "#b2182b"), guide = "none") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0.04, 0.12))) +
  labs(
    title    = "Reporting-calendar and campaign share of the seasonal signal, by sector",
    subtitle = paste0("Variance attributable to the fiscal quarter-end and 100-day campaign terms; ",
                      "95% cluster-bootstrap CI (B = ", B, ")"),
    x = NULL, y = "Share of modelled seasonal variance (%)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    axis.text.x   = element_text(size = 12)
  )

out_png <- file.path(OUT, "o2_partition.png")
ggsave(out_png, p, width = 8, height = 5.5, dpi = 150)
cat("Saved:", out_png, "\n\n")

## ---- 10. VERIFICATION ----------------------------------------------------

cat("============================================================\n")
cat("VERIFICATION\n")
cat("============================================================\n\n")

## V1: Shares sum to ~100% per sector
cat("[V1] Do three shares sum to ~100% per sector?\n")
for (sec in c("Public", "Private")) {
  sub <- result_tbl[result_tbl$sector == sec, ]
  pt_sum <- sum(sub$share_pct)
  cat(sprintf("  %-8s: %.4f%%  => %s\n", sec, pt_sum,
              if (abs(pt_sum - 100) < 0.01) "PASS" else "FAIL"))
}
cat("\n")

## V2: Full table with CIs
cat("[V2] Partition table with bootstrap CIs:\n")
display_tbl <- result_tbl
display_tbl$share_pct <- round(display_tbl$share_pct, 2)
display_tbl$share_lo  <- round(display_tbl$share_lo,  2)
display_tbl$share_hi  <- round(display_tbl$share_hi,  2)
display_tbl$ci_str    <- sprintf("[%.1f, %.1f]%%", display_tbl$share_lo, display_tbl$share_hi)
cat(sprintf("  %-10s  %-28s  %8s  %20s\n", "Sector", "Block", "Share(%)", "95% CI"))
cat(sprintf("  %s\n", paste(rep("-", 72), collapse = "")))
for (i in seq_len(nrow(display_tbl))) {
  r <- display_tbl[i, ]
  cat(sprintf("  %-10s  %-28s  %7.2f%%  %s\n",
              as.character(r$sector), as.character(r$block),
              r$share_pct, r$ci_str))
}
cat("\n")

## V3: Calendar/campaign artefact contrast
cat("[V3] CALENDAR/CAMPAIGN (artefact) block: public vs private\n")
cal_pub  <- result_tbl$share_pct[result_tbl$sector == "Public"  & grepl("Calendar", result_tbl$block)]
cal_priv <- result_tbl$share_pct[result_tbl$sector == "Private" & grepl("Calendar", result_tbl$block)]
ci_pub_cal  <- result_tbl[result_tbl$sector == "Public"  & grepl("Calendar", result_tbl$block),
                           c("share_lo", "share_hi")]
ci_priv_cal <- result_tbl[result_tbl$sector == "Private" & grepl("Calendar", result_tbl$block),
                           c("share_lo", "share_hi")]

cat(sprintf("  Public  calendar share: %.2f%% [%.1f, %.1f]%%\n",
            cal_pub,  ci_pub_cal$share_lo,  ci_pub_cal$share_hi))
cat(sprintf("  Private calendar share: %.2f%% [%.1f, %.1f]%%\n",
            cal_priv, ci_priv_cal$share_lo, ci_priv_cal$share_hi))

artefact_direction <- cal_pub > cal_priv
ci_overlap <- (ci_pub_cal$share_lo  <= ci_priv_cal$share_hi) &&
              (ci_priv_cal$share_lo <= ci_pub_cal$share_hi)

cat(sprintf("\n  Artefact share LARGER in public than private: %s [hypothesis direction]\n",
            if (artefact_direction) "YES" else "NO"))
cat(sprintf("  Bootstrap CIs overlap: %s\n",
            if (ci_overlap) "YES (cannot rule out chance)" else "NO (non-overlapping -- clear separation)"))
cat("\n")

## V4: Climate shares
cat("[V4] CLIMATE block: public vs private (expected similar)\n")
cli_pub  <- result_tbl$share_pct[result_tbl$sector == "Public"  & grepl("Climate", result_tbl$block)]
cli_priv <- result_tbl$share_pct[result_tbl$sector == "Private" & grepl("Climate", result_tbl$block)]
ci_pub_cli  <- result_tbl[result_tbl$sector == "Public"  & grepl("Climate", result_tbl$block),
                           c("share_lo", "share_hi")]
ci_priv_cli <- result_tbl[result_tbl$sector == "Private" & grepl("Climate", result_tbl$block),
                           c("share_lo", "share_hi")]

cat(sprintf("  Public  climate share: %.2f%% [%.1f, %.1f]%%\n",
            cli_pub,  ci_pub_cli$share_lo,  ci_pub_cli$share_hi))
cat(sprintf("  Private climate share: %.2f%% [%.1f, %.1f]%%\n",
            cli_priv, ci_priv_cli$share_lo, ci_priv_cli$share_hi))

cli_overlap <- (ci_pub_cli$share_lo  <= ci_priv_cli$share_hi) &&
               (ci_priv_cli$share_lo <= ci_pub_cli$share_hi)

cat(sprintf("  CIs overlap: %s\n", if (cli_overlap) "YES" else "NO"))
cat(sprintf("  Delta (public - private): %.2f%%\n", cli_pub - cli_priv))
cat("\n")

## V5: var(S) summary
cat("[V5] var(S) -- total modelled seasonal variation magnitude (log-rate scale):\n")
cat(sprintf("  Public  var(S): %.6f\n", shares_pub["varS"]))
cat(sprintf("  Private var(S): %.6f\n", shares_priv["varS"]))
cat(sprintf("  Ratio public/private: %.3f\n",
            shares_pub["varS"] / shares_priv["varS"]))
cat("\n")

## V6: Collinearity caveat reminder
cat("[V6] INTERPRETATION CAVEAT:\n")
cat("  Harmonic and climate blocks are collinear: both encode the annual cycle.\n")
cat("  The harmonic block will absorb much of the smooth intra-annual cycle;\n")
cat("  the climate block captures residual meteorological signal AFTER the pure\n")
cat("  sinusoidal pattern is accounted for. Individual harmonic and climate shares\n")
cat("  should not be interpreted in isolation.\n")
cat("  HEADLINE CONTRAST: the CALENDAR/CAMPAIGN (artefact) share difference\n")
cat("  between public and private sectors is the primary finding of this task.\n\n")

cat("============================================================\n")
cat("Script 22_variance_partition.R COMPLETE.\n")
cat("Outputs:\n")
cat("  CSV :", out_csv, "\n")
cat("  PNG :", out_png, "\n")
cat("============================================================\n")
