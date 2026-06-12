# 03_seasonality_stl.R — Objective 1 — national and sector STL seasonality
#
# Renumbered, path-normalized reorganization of the working pipeline scripts: 10_stl_national.R.
# Computation is unchanged from the scripts that produced the shipped results in out/.

source(here::here("code", "00_setup.R"))


## ===========================================================================
## from 10_stl_national.R
## ===========================================================================

# 10_stl_national.R
# Phase 1 Task 1.1 (O1): National STL seasonal decomposition for TB notifications
# Quantify seasonal pattern by sector (total / public / private) with bootstrap CIs.
# Method: STL on log(VAR+1), multiplicative interpretation via back-transform.
# Bootstrap: moving-block bootstrap (block length = 12 months) on the STL
#            remainder-adjusted series, B=2000 reps.

# paths and base packages come from 00_setup.R (sourced above)

## ---- install / load feasts + fabletools ------------------------------------
needed <- c("feasts", "fabletools", "fable")
to_inst <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_inst) > 0) {
  message("Installing: ", paste(to_inst, collapse = ", "))
  install.packages(to_inst, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(feasts)
  library(fabletools)
  library(fable)
})

set.seed(SEED)

## ---- 1. Load & aggregate to national monthly series -----------------------
panel <- readRDS(file.path(DER, "panel_dm.rds"))
# panel is a tsibble keyed by dist_lgd, index = month (yearmonth), 747 districts x 53 months
# Drop incomplete months (lag-affected 2026-05; complete==FALSE)
panel_complete <- panel |> dplyr::filter(complete)

# Aggregate to national: one row per month, sums across all districts
nat <- panel_complete |>
  index_by(month) |>
  summarise(
    public  = sum(public,  na.rm = TRUE),
    private = sum(private, na.rm = TRUE),
    total   = sum(total,   na.rm = TRUE),
    .groups = "drop"
  )
# nat is a 1-key-free tsibble (regular monthly, ~52 months)
stopifnot(nrow(nat) >= 24, !any(is.na(nat$total)))
message("National series: ", nrow(nat), " complete months")

## ---- 2. Pivot to long & fit STL on log(VAR+1) ----------------------------
# Pivot to a keyed tsibble: key = series (total/public/private)
nat_long <- nat |>
  as_tibble() |>
  tidyr::pivot_longer(cols = c(total, public, private),
                      names_to = "series", values_to = "count") |>
  mutate(log_count = log(count + 1)) |>
  as_tsibble(index = month, key = series)

# Fit STL: log(count+1) ~ trend() + season(window=11), robust=TRUE
stl_fits <- nat_long |>
  model(
    stl = STL(log_count ~ trend() + season(window = 11), robust = TRUE)
  )

# Extract components from all series
stl_components <- components(stl_fits)
# columns: month, series, log_count, trend, season_year, remainder, season_adjust

# Save STL components to derived/
saveRDS(stl_components, file.path(DER, "stl_national.rds"))
message("STL components saved to ", file.path(DER, "stl_national.rds"))

## ---- 3. Point estimates: peak/trough, amplitude, seasonal strength --------
# Add calendar month number and name
comp_df <- stl_components |>
  as_tibble() |>
  mutate(cal_month = lubridate::month(tsibble::yearmonth(month)))

# Mean seasonal factor by calendar month (on log scale)
season_means <- comp_df |>
  group_by(series, cal_month) |>
  summarise(mean_season = mean(season_year, na.rm = TRUE), .groups = "drop")

# Verification check 1: seasonal component means ~0 on log scale
cat("\n=== VERIFICATION: Mean seasonal component per series (should be ~0) ===\n")
overall_means <- comp_df |>
  group_by(series) |>
  summarise(mean_season_overall = mean(season_year, na.rm = TRUE))
print(overall_means)

# Peak/trough by calendar month
peak_trough <- season_means |>
  group_by(series) |>
  summarise(
    peak_month   = month.name[cal_month[which.max(mean_season)]],
    trough_month = month.name[cal_month[which.min(mean_season)]],
    max_mean_season = max(mean_season),
    min_mean_season = min(mean_season),
    amp_pct = 100 * (exp(max_mean_season) - exp(min_mean_season)),
    .groups = "drop"
  )

# Seasonal strength via feasts::feat_stl on the log series
strength_df <- nat_long |>
  features(log_count, feat_stl) |>
  select(series, seasonal_strength_year)

## ---- 4. Moving-block bootstrap for CI on amp_pct -------------------------
# Strategy: moving-block bootstrap (block length L=12) on the REMAINDER of each
# STL decomposition (the detrended, de-seasonalised residuals). For each resample,
# reconstruct the series as trend + season + bootstrap_remainder, then re-estimate
# seasonal month means from that reconstructed series (via another STL fit).
# This propagates uncertainty in the seasonal estimates driven by within-series
# autocorrelation, without re-running full STL 2000 times on observed data.
#
# NOTE: We avoid a full STL re-fit per bootstrap rep for speed. Instead we
# bootstrap via: construct resample of log(count+1) = trend + season + boot_remainder,
# then fit a new STL on that resample. This IS a full STL refit — but on a
# short (52-month) univariate series STL is very fast (< 1s per rep).

do_bootstrap_amp <- function(series_name, comp_df, B, L = 12) {
  sub <- comp_df |>
    dplyr::filter(series == series_name) |>
    dplyr::arrange(month)

  trend_v   <- sub$trend
  season_v  <- sub$season_year
  remainder <- sub$remainder
  n         <- length(remainder)
  months_v  <- sub$month

  # Number of blocks that fit; blocks can overlap (moving block)
  n_blocks_pool <- n - L + 1
  blocks_needed <- ceiling(n / L)

  amp_boot <- numeric(B)

  for (b in seq_len(B)) {
    # Sample block start indices with replacement
    starts <- sample(seq_len(n_blocks_pool), size = blocks_needed, replace = TRUE)
    boot_rem <- unlist(lapply(starts, function(s) remainder[s:(s + L - 1)]))[seq_len(n)]

    # Reconstruct log-scale series
    boot_log <- trend_v + season_v + boot_rem

    # Build a temporary tsibble for this replicate
    tmp <- tsibble::tsibble(
      month     = months_v,
      log_count = boot_log,
      index     = month
    )

    # Fit STL
    fit_b <- tmp |>
      fabletools::model(
        stl = feasts::STL(log_count ~ trend() + season(window = 11), robust = TRUE)
      )
    comp_b <- fabletools::components(fit_b) |>
      dplyr::as_tibble() |>
      dplyr::mutate(cal_month = lubridate::month(tsibble::yearmonth(month)))

    # Seasonal month means
    sm_b <- comp_b |>
      dplyr::group_by(cal_month) |>
      dplyr::summarise(mean_s = mean(season_year, na.rm = TRUE), .groups = "drop")

    amp_boot[b] <- 100 * (exp(max(sm_b$mean_s)) - exp(min(sm_b$mean_s)))
  }
  amp_boot
}

message("Running moving-block bootstrap (B=", B, ", L=12) for 3 series...")
message("  series: total ...")
boot_total   <- do_bootstrap_amp("total",   comp_df, B)
message("  series: public ...")
boot_public  <- do_bootstrap_amp("public",  comp_df, B)
message("  series: private ...")
boot_private <- do_bootstrap_amp("private", comp_df, B)
message("Bootstrap complete.")

ci_df <- tibble::tibble(
  series  = c("total", "public", "private"),
  amp_lo  = c(quantile(boot_total,   0.025),
               quantile(boot_public,  0.025),
               quantile(boot_private, 0.025)),
  amp_hi  = c(quantile(boot_total,   0.975),
               quantile(boot_public,  0.975),
               quantile(boot_private, 0.975))
)

## ---- 5. Assemble output table --------------------------------------------
results <- peak_trough |>
  dplyr::left_join(strength_df, by = "series") |>
  dplyr::left_join(ci_df,       by = "series") |>
  dplyr::select(series, peak_month, trough_month, amp_pct, amp_lo, amp_hi,
                seasonal_strength = seasonal_strength_year) |>
  dplyr::arrange(factor(series, levels = c("total", "public", "private")))

## ---- Verification 2: finite amp, amp_lo <= amp_pct <= amp_hi -------------
cat("\n=== VERIFICATION: Seasonal summary table ===\n")
print(results, digits = 4)

checks_pass <- all(
  is.finite(results$amp_pct) &
  results$amp_lo <= results$amp_pct &
  results$amp_pct <= results$amp_hi
)
cat("\namp_lo <= amp_pct <= amp_hi for all series:", checks_pass, "\n")

## ---- Verification 3: public vs private comparison ------------------------
pub_row  <- results |> dplyr::filter(series == "public")
priv_row <- results |> dplyr::filter(series == "private")

cat("\n=== PUBLIC vs PRIVATE AMPLITUDE COMPARISON ===\n")
cat(sprintf("  Public  amp_pct = %.2f%% [%.2f, %.2f]\n",
            pub_row$amp_pct, pub_row$amp_lo, pub_row$amp_hi))
cat(sprintf("  Private amp_pct = %.2f%% [%.2f, %.2f]\n",
            priv_row$amp_pct, priv_row$amp_lo, priv_row$amp_hi))

pub_gt_priv <- pub_row$amp_pct > priv_row$amp_pct
ci_overlap  <- pub_row$amp_lo  <= priv_row$amp_hi & priv_row$amp_lo <= pub_row$amp_hi

cat(sprintf(
  "  Public-sector amplitude IS %s than private-sector (artefact hypothesis direction: %s).\n",
  ifelse(pub_gt_priv, "LARGER", "smaller"),
  ifelse(pub_gt_priv, "SUPPORTED", "NOT supported")
))
cat(sprintf(
  "  Bootstrap 95%% CIs %s.\n",
  ifelse(ci_overlap, "OVERLAP (evidence is not conclusive)", "DO NOT overlap (clear separation)")
))

## ---- Write CSV output ----------------------------------------------------
readr::write_csv(results, file.path(OUT, "o1_seasonal_national.csv"))
message("Output written to ", file.path(OUT, "o1_seasonal_national.csv"))
message("\nDone. Script completed successfully.")
