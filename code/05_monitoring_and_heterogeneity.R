# 05_monitoring_and_heterogeneity.R — Objective 3 — monitoring payoff and district heterogeneity
#
# Renumbered, path-normalized reorganization of the working pipeline scripts: 31_payoff.R, 30_heterogeneity.R, 35_heterogeneity_adjusted.R.
# Computation is unchanged from the scripts that produced the shipped results in out/.

source(here::here("code", "00_setup.R"))


## ===========================================================================
## from 31_payoff.R
## ===========================================================================

# 31_payoff.R
# Phase 3, Task 3.2 (O3) - Surveillance Payoff Analysis
#
# Quantifies the practical cost of NOT adjusting TB notifications for
# seasonal and working-day artefacts: how often a district's month-on-month
# change that a programme would react to (a "surge" or "dip") disappears
# after seasonal + working-day adjustment.
#
# Outputs:
#   out/o3_payoff.csv        - per-district artefact rates
#   out/o3_payoff_map.png    - choropleth of sign-flip rate


suppressMessages({
  library(feasts)
  library(fabletools)
  library(sf)
  library(ggplot2)
})

sf::sf_use_s2(FALSE)

cat("=== 31_payoff.R: Surveillance Payoff Analysis ===\n\n")

# ── 1. Load data ─────────────────────────────────────────────────────────────

cat("Loading inputs...\n")
panel_dm <- readRDS(file.path(DER, "panel_dm.rds"))
calendar  <- readRDS(file.path(DER, "calendar.rds"))

# Build per-district TOTAL monthly series from complete months only
# We keep all rows but will use complete flag for eligibility, not for filtering
# the series (we need consecutive series for STL)
dm_tbl <- as_tibble(panel_dm)

# Join working_days from calendar
dm_tbl <- dm_tbl |>
  dplyr::left_join(
    dplyr::select(as_tibble(calendar), month, working_days),
    by = "month"
  )

# ── 2. Working-day normalisation ─────────────────────────────────────────────
# total_wd = total / working_days * mean(working_days)
mean_wd <- mean(calendar$working_days)
cat(sprintf("Mean working days: %.4f\n", mean_wd))

dm_tbl <- dm_tbl |>
  dplyr::mutate(
    total_wd = total / working_days * mean_wd
  )

# ── 3. Determine which districts use STL vs simple index vs excluded ──────────

# Eligibility: STL needs >=24 complete months (two full seasonal cycles); below
# that a seasonal pattern cannot be estimated, so the district is flagged as
# insufficient data rather than excluded silently. (Threshold lowered from 36 to
# 24 so a 35-month district is no longer dropped by a single month.)
MIN_MONTHS <- 24L
elig <- dm_tbl |>
  dplyr::filter(complete == TRUE) |>
  dplyr::group_by(dist_lgd) |>
  dplyr::summarise(
    n_complete  = dplyr::n(),
    med_total   = median(total),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    method = dplyr::case_when(
      n_complete >= MIN_MONTHS & med_total >= 5 ~ "STL",
      n_complete >= MIN_MONTHS & med_total <  5 ~ "simple",
      TRUE                                      ~ "insufficient"
    )
  )

n_stl      <- sum(elig$method == "STL")
n_simple   <- sum(elig$method == "simple")
n_insuff   <- sum(elig$method == "insufficient")
cat(sprintf("\nDistrict method breakdown:\n"))
cat(sprintf("  STL (>=24 complete months AND med>=5): %d\n", n_stl))
cat(sprintf("  Simple index (>=24 complete, med<5):   %d\n", n_simple))
cat(sprintf("  Insufficient data (<24 months):         %d\n", n_insuff))
cat(sprintf("  Total districts in panel:              %d\n", nrow(elig)))

stl_dists    <- elig$dist_lgd[elig$method == "STL"]
simple_dists <- elig$dist_lgd[elig$method == "simple"]
insuff_dists <- elig$dist_lgd[elig$method == "insufficient"]

# ── 4. STL seasonal adjustment ───────────────────────────────────────────────

cat("\nFitting STL for", n_stl, "districts (this may take a minute)...\n")

# Build tsibble for STL districts using all rows (not just complete),
# log-transform total_wd + 1
dm_stl <- dm_tbl |>
  dplyr::filter(dist_lgd %in% stl_dists) |>
  dplyr::mutate(log_total_wd = log(total_wd + 1)) |>
  tsibble::as_tsibble(key = dist_lgd, index = month)

# Fit STL: decompose log_total_wd
# feasts::STL with robust=TRUE; period=12 for monthly data
stl_fits <- dm_stl |>
  fabletools::model(
    stl = feasts::STL(log_total_wd ~ trend(window = 13) + season(window = "periodic"), robust = TRUE)
  )

# Extract components
stl_comps <- fabletools::components(stl_fits)

# Seasonally-adjusted: remove season_year component from log space
# sa_log = trend + remainder = log_total_wd - season_year
stl_sa <- stl_comps |>
  as_tibble() |>
  dplyr::mutate(
    sa = exp(trend + remainder)   # back-transform: removes seasonal component
  ) |>
  dplyr::select(dist_lgd, month, sa)

cat(sprintf("  STL fits complete. SA series rows: %d\n", nrow(stl_sa)))

# ── 5. Simple multiplicative seasonal index ───────────────────────────────────

cat("Computing simple seasonal index for", n_simple, "districts...\n")

# For simple: ratio-to-12-month-moving-average by calendar month
# Only for districts with n_complete >= 36 (simple group)
# For excluded: will not produce SA series

compute_simple_sa <- function(sub) {
  # sub is a data frame for a single district, sorted by month
  sub <- sub |> dplyr::arrange(month)
  n <- nrow(sub)
  # 12-month centred moving average
  ma <- stats::filter(sub$total_wd, rep(1/12, 12), sides = 2)
  # ratio series
  ratio <- sub$total_wd / ma
  # seasonal index: average ratio by calendar month (integer 1-12)
  sub$mo_num <- as.integer(format(as.Date(sub$month), "%m"))
  si <- tapply(ratio, sub$mo_num, mean, na.rm = TRUE)
  # normalize so indices sum to 12
  si <- si / mean(si, na.rm = TRUE)
  # deseasonalise: sa = total_wd / seasonal_index_for_that_month
  sub$sa <- sub$total_wd / si[sub$mo_num]
  sub[, c("dist_lgd", "month", "sa")]
}

if (n_simple > 0) {
  simple_list <- lapply(simple_dists, function(d) {
    sub <- dm_tbl |>
      dplyr::filter(dist_lgd == d) |>
      dplyr::arrange(month)
    compute_simple_sa(sub)
  })
  simple_sa <- dplyr::bind_rows(simple_list)
  cat(sprintf("  Simple SA series rows: %d\n", nrow(simple_sa)))
} else {
  simple_sa <- tibble::tibble(dist_lgd = character(), month = tsibble::yearmonth(character()), sa = numeric())
}

# Combine SA series from both methods
sa_all <- dplyr::bind_rows(stl_sa, simple_sa)
cat(sprintf("Combined SA series: %d rows across %d districts\n", nrow(sa_all), dplyr::n_distinct(sa_all$dist_lgd)))

# ── 6. Compute month-on-month percentage changes ──────────────────────────────

cat("\nComputing month-on-month changes...\n")

# Merge raw total and SA series
mom_data <- dm_tbl |>
  dplyr::filter(dist_lgd %in% c(stl_dists, simple_dists)) |>
  dplyr::select(dist_lgd, month, total) |>
  dplyr::left_join(sa_all, by = c("dist_lgd", "month")) |>
  dplyr::arrange(dist_lgd, month) |>
  dplyr::group_by(dist_lgd) |>
  dplyr::mutate(
    mom_raw = (total - dplyr::lag(total)) / dplyr::lag(total),
    mom_adj = (sa    - dplyr::lag(sa))    / dplyr::lag(sa)
  ) |>
  dplyr::ungroup()

# ── 7. Define artefact indicators per district-month ─────────────────────────

mom_flags <- mom_data |>
  dplyr::filter(!is.na(mom_raw), !is.na(mom_adj), is.finite(mom_raw), is.finite(mom_adj)) |>
  dplyr::mutate(
    sign_flip      = sign(mom_raw) != sign(mom_adj),
    surge_artefact = (mom_raw >  0.10) & !(mom_adj >  0.10),
    dip_artefact   = (mom_raw < -0.10) & !(mom_adj < -0.10)
  )

cat(sprintf("District-months available for analysis: %d\n", nrow(mom_flags)))

# ── 8. Aggregate per district ─────────────────────────────────────────────────

per_district <- mom_flags |>
  dplyr::group_by(dist_lgd) |>
  dplyr::summarise(
    n_months          = dplyr::n(),
    n_raw_surges      = sum(mom_raw >  0.10),
    n_raw_dips        = sum(mom_raw < -0.10),
    sign_flip_rate    = mean(sign_flip),
    surge_artefact_rate = mean(surge_artefact),
    dip_artefact_rate   = mean(dip_artefact),
    # among raw surges: fraction that are artefactual
    surge_artefact_among_surges = ifelse(n_raw_surges > 0,
                                         sum(surge_artefact) / n_raw_surges,
                                         NA_real_),
    dip_artefact_among_dips     = ifelse(n_raw_dips > 0,
                                         sum(dip_artefact) / n_raw_dips,
                                         NA_real_),
    .groups = "drop"
  ) |>
  dplyr::left_join(
    dplyr::select(elig, dist_lgd, method, n_complete, med_total),
    by = "dist_lgd"
  )

cat(sprintf("Rows in per-district output: %d\n", nrow(per_district)))

# ── 9. National headline numbers ──────────────────────────────────────────────

# Overall district-months: sign flip rate
n_dm_total      <- nrow(mom_flags)
n_sign_flip     <- sum(mom_flags$sign_flip)
pct_sign_flip   <- 100 * n_sign_flip / n_dm_total

# Raw >10% surges: fraction artefactual
n_raw_surges_all    <- sum(mom_flags$mom_raw >  0.10)
n_surge_artefact_all <- sum(mom_flags$surge_artefact)
pct_surge_artefact  <- 100 * n_surge_artefact_all / n_raw_surges_all

# Raw >10% dips: fraction artefactual
n_raw_dips_all      <- sum(mom_flags$mom_raw < -0.10)
n_dip_artefact_all  <- sum(mom_flags$dip_artefact)
pct_dip_artefact    <- 100 * n_dip_artefact_all / n_raw_dips_all

cat("\n=== NATIONAL HEADLINE NUMBERS ===\n")
cat(sprintf("Total district-months analysed: %d\n", n_dm_total))
cat(sprintf("  (1) Sign-flip rate (MoM direction reverses after adjustment): %.1f%% of district-months\n", pct_sign_flip))
cat(sprintf("  (2) Raw >10%% surges: %d total | %d artefactual (%.1f%%)\n",
            n_raw_surges_all, n_surge_artefact_all, pct_surge_artefact))
cat(sprintf("  (3) Raw >10%% dips:   %d total | %d artefactual (%.1f%%)\n",
            n_raw_dips_all, n_dip_artefact_all, pct_dip_artefact))

# Per-district summary of sign-flip rate
sfr <- per_district$sign_flip_rate
cat(sprintf("\nPer-district sign-flip rate — median: %.3f  IQR: [%.3f, %.3f]  range: [%.3f, %.3f]\n",
            median(sfr), quantile(sfr, 0.25), quantile(sfr, 0.75),
            min(sfr), max(sfr)))

# ── 10. Save per-district CSV ──────────────────────────────────────────────────

out_csv <- file.path(OUT, "o3_payoff.csv")
readr::write_csv(per_district, out_csv)
cat(sprintf("\nSaved: %s (%d rows)\n", out_csv, nrow(per_district)))

# ── 11. Map: choropleth of sign_flip_rate ────────────────────────────────────

cat("\nGenerating map...\n")
gpkg_path <- GEO
dist_sf <- sf::st_read(gpkg_path, quiet = TRUE)

cat(sprintf("Spatial file: %d features, key column: dist_lgd_anl\n", nrow(dist_sf)))

# Join sign_flip_rate to spatial data
dist_map <- dist_sf |>
  dplyr::left_join(
    dplyr::select(per_district, dist_lgd, sign_flip_rate),
    by = c("dist_lgd_anl" = "dist_lgd")
  )

matched    <- sum(!is.na(dist_map$sign_flip_rate))
unmatched  <- sum( is.na(dist_map$sign_flip_rate))
cat(sprintf("Spatial join: %d matched, %d unmatched\n", matched, unmatched))

map_plot <- ggplot2::ggplot(dist_map) +
  ggplot2::geom_sf(
    ggplot2::aes(fill = sign_flip_rate),
    colour = NA
  ) +
  ggplot2::scale_fill_viridis_c(
    name   = "Sign-flip\nrate",
    option = "plasma",
    na.value = "grey80",
    labels = scales::percent_format(accuracy = 1)
  ) +
  ggplot2::labs(
    title    = "TB notification MoM direction artefact rate by district",
    subtitle = "Fraction of district-months where MoM sign reverses after seasonal + working-day adjustment",
    caption  = paste0(nrow(per_district), " districts mapped; grey = insufficient data for stable seasonal adjustment: ",
                      unmatched, " districts (Jammu & Kashmir, Arunachal Pradesh, Sikkim)")
  ) +
  ggplot2::theme_void(base_size = 11) +
  ggplot2::theme(
    plot.title    = ggplot2::element_text(size = 13, face = "bold", margin = ggplot2::margin(b = 4)),
    plot.subtitle = ggplot2::element_text(size = 9,  colour = "grey30", margin = ggplot2::margin(b = 8)),
    plot.caption  = ggplot2::element_text(size = 7,  colour = "grey50"),
    legend.position = "right"
  )

map_path <- file.path(OUT, "o3_payoff_map.png")
ggplot2::ggsave(map_path, map_plot, width = 10, height = 7, dpi = 150, bg = "white")
cat(sprintf("Map saved: %s\n", map_path))

# ── 12. VERIFICATION ──────────────────────────────────────────────────────────

cat("\n=== VERIFICATION ===\n")

# V1: o3_payoff.csv correctness
v1_n_rows     <- nrow(per_district)
v1_rates_ok   <- all(
  is.finite(per_district$sign_flip_rate)    & per_district$sign_flip_rate    >= 0 & per_district$sign_flip_rate    <= 1 &
  is.finite(per_district$surge_artefact_rate) & per_district$surge_artefact_rate >= 0 & per_district$surge_artefact_rate <= 1 &
  is.finite(per_district$dip_artefact_rate)   & per_district$dip_artefact_rate   >= 0 & per_district$dip_artefact_rate   <= 1
)
cat(sprintf("V1: o3_payoff.csv rows = %d | All rates finite in [0,1]: %s\n",
            v1_n_rows, if (v1_rates_ok) "OK" else "FAIL"))
cat(sprintf("    STL=%d, simple=%d, insufficient=%d\n", n_stl, n_simple, n_insuff))

# V2: National headline
cat(sprintf("V2: NATIONAL HEADLINE\n"))
cat(sprintf("    Sign-flip: %.1f%% of district-months\n",    pct_sign_flip))
cat(sprintf("    Surge artefact: %.1f%% of raw >10%% surges\n", pct_surge_artefact))
cat(sprintf("    Dip artefact:   %.1f%% of raw >10%% dips\n",   pct_dip_artefact))

# V3: Map file written and non-empty
map_exists   <- file.exists(map_path)
map_size     <- if (map_exists) file.info(map_path)$size else 0L
cat(sprintf("V3: Map file exists: %s | Size: %d bytes | %s\n",
            map_exists, map_size,
            if (map_exists && map_size > 10000) "OK" else "FAIL"))

cat("\n=== DONE ===\n")

## ===========================================================================
## from 30_heterogeneity.R
## ===========================================================================

# 30_heterogeneity.R
# Phase 3, Task 3.1 (O3): What district features predict where the reporting
# artefact bites?
#
# PRIMARY HYPOTHESIS: artefact is larger where public_share is higher
# (the artefact is public-sector-driven).
#
# Outputs:
#   out/o3_heterogeneity.csv  - coefficient table (term, estimate, ci_lo, ci_hi, p, model)
#   out/o3_heterogeneity.png  - forest plot of primary model + partial-effect plot



# ── Packages ─────────────────────────────────────────────────────────────────
for (pkg in c("betareg", "sf", "spdep", "spatialreg", "ggplot2", "patchwork")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# ── 1. Load and join inputs ───────────────────────────────────────────────────
payoff   <- read_csv(file.path(OUT, "o3_payoff.csv"), show_col_types = FALSE)
modifiers <- readRDS(file.path(DER, "modifiers.rds"))

cat("Payoff rows:", nrow(payoff), "\n")
cat("Modifiers rows:", nrow(modifiers), "\n")

# Inner join on dist_lgd
dat <- inner_join(payoff, modifiers, by = "dist_lgd")
cat("After inner join:", nrow(dat), "districts\n")

# Drop districts missing key predictor fields
key_vars <- c("public_share", "latitude", "density_2024",
              "pct_urban_notif", "region", "total_notif_2022_2025",
              "surge_artefact_rate", "sign_flip_rate")
n_before <- nrow(dat)
dat <- dat[complete.cases(dat[, key_vars]), ]
cat("After dropping missing key fields:", nrow(dat), "districts",
    "(dropped", n_before - nrow(dat), ")\n")

# ── 2. Construct predictors ───────────────────────────────────────────────────
# NOTE on pct_urban_notif: 434 of 747 districts have pct_urban_notif == 1
# (right-censored – all notifications from urban areas). This limits its
# discrimination power; it is included but its coefficient should be
# interpreted cautiously.

dat <- dat |>
  mutate(
    region          = factor(region),          # reference = "Central"
    lat_sc          = as.numeric(scale(latitude)),
    logdens_sc      = as.numeric(scale(log(density_2024 + 1))),
    lognotif_sc     = as.numeric(scale(log(total_notif_2022_2025 + 1)))
  )

cat("\npct_urban_notif == 1:", sum(dat$pct_urban_notif == 1),
    "of", nrow(dat), "districts (right-censored; flagged)\n")
cat("public_share range:", round(range(dat$public_share), 3), "\n")

# ── 3. PRIMARY MODEL: surge_artefact_rate ~ predictors ───────────────────────
# Response is in [0, 1). Has some zeros (6 districts) but no ones.
# Strategy: beta regression with boundary-squeeze transform y' = (y*(n-1)+0.5)/n
# where n = number of observations per district. We use a fixed squeeze to
# maintain comparability across districts: y' = (y*(N-1) + 0.5) / N where
# N = nrow(dat). This is the Smithson & Verkuilen (2006) transformation.
# We prefer betareg over quasibinomial because the response is a continuous
# rate (not a count proportion) and beta regression is the natural model.

N <- nrow(dat)
dat$y_surge <- (dat$surge_artefact_rate * (N - 1) + 0.5) / N
dat$y_flip  <- (dat$sign_flip_rate      * (N - 1) + 0.5) / N

cat("\ny_surge range after squeeze:", round(range(dat$y_surge), 5), "\n")
cat("y_flip  range after squeeze:", round(range(dat$y_flip), 5), "\n")

# Formula
f_primary <- y_surge ~ public_share + lat_sc + logdens_sc +
             pct_urban_notif + region + lognotif_sc

cat("\n=== PRIMARY MODEL: betareg (surge_artefact_rate) ===\n")
fit_primary <- betareg(f_primary, data = dat, link = "logit")
summary(fit_primary)

# ── 4. SECONDARY MODEL: sign_flip_rate ───────────────────────────────────────
f_secondary <- y_flip ~ public_share + lat_sc + logdens_sc +
               pct_urban_notif + region + lognotif_sc

cat("\n=== SECONDARY MODEL: betareg (sign_flip_rate) ===\n")
fit_secondary <- betareg(f_secondary, data = dat, link = "logit")
summary(fit_secondary)

# ── 5. Extract coefficient tables ─────────────────────────────────────────────
extract_betareg <- function(fit, model_name) {
  cf <- coef(summary(fit))$mean
  ci <- confint(fit)
  # confint for betareg returns full matrix incl phi; keep only mean params
  mean_terms <- rownames(cf)
  ci_mean <- ci[mean_terms, , drop = FALSE]
  tibble(
    term    = mean_terms,
    estimate = cf[, "Estimate"],
    se       = cf[, "Std. Error"],
    z        = cf[, "z value"],
    p        = cf[, "Pr(>|z|)"],
    ci_lo    = ci_mean[, 1],
    ci_hi    = ci_mean[, 2],
    model    = model_name
  )
}

tbl_primary   <- extract_betareg(fit_primary,   "surge_artefact_betareg")
tbl_secondary <- extract_betareg(fit_secondary, "sign_flip_betareg")

cat("\n--- PRIMARY model coefficient table ---\n")
print(as.data.frame(tbl_primary[, c("term","estimate","ci_lo","ci_hi","p")]),
      digits = 4)

cat("\n--- SECONDARY model coefficient table ---\n")
print(as.data.frame(tbl_secondary[, c("term","estimate","ci_lo","ci_hi","p")]),
      digits = 4)

# Effect size: public_share 0.5 -> 0.9 (mean and link scale)
ps_coef <- tbl_primary$estimate[tbl_primary$term == "public_share"]
# In logit-link betareg, the linear predictor eta = Xb; fitted rate = logistic(eta)
# Holding other predictors at their sample means:
other_terms <- setdiff(tbl_primary$term, c("(Intercept)", "public_share"))

# Build representative mean row
xbar <- dat |>
  summarise(
    lat_sc        = mean(lat_sc),
    logdens_sc    = mean(logdens_sc),
    pct_urban_notif = mean(pct_urban_notif),
    lognotif_sc   = mean(lognotif_sc)
  )

# Most common region for representative prediction
modal_region <- names(sort(table(dat$region), decreasing = TRUE))[1]

newdat_lo <- data.frame(
  public_share    = 0.5,
  lat_sc          = xbar$lat_sc,
  logdens_sc      = xbar$logdens_sc,
  pct_urban_notif = xbar$pct_urban_notif,
  region          = factor(modal_region, levels = levels(dat$region)),
  lognotif_sc     = xbar$lognotif_sc
)
newdat_hi <- newdat_lo
newdat_hi$public_share <- 0.9

pred_lo <- predict(fit_primary, newdata = newdat_lo, type = "response")
pred_hi <- predict(fit_primary, newdata = newdat_hi, type = "response")
effect_size <- pred_hi - pred_lo

cat("\n=== PUBLIC_SHARE EFFECT SIZE ===\n")
cat("Primary model: public_share coefficient (logit scale):", round(ps_coef, 4), "\n")
cat("Predicted surge_artefact_rate at public_share = 0.5:", round(pred_lo, 4), "\n")
cat("Predicted surge_artefact_rate at public_share = 0.9:", round(pred_hi, 4), "\n")
cat("Effect (change in artefact rate, 0.5 -> 0.9):", round(effect_size, 4), "\n")
cat("Direction: public_share coefficient is", ifelse(ps_coef > 0, "POSITIVE", "NEGATIVE"),
    "(p =", round(tbl_primary$p[tbl_primary$term == "public_share"], 4), ")\n")

ps_coef_sec <- tbl_secondary$estimate[tbl_secondary$term == "public_share"]
ps_p_sec    <- tbl_secondary$p[tbl_secondary$term == "public_share"]
cat("\nSecondary model (sign_flip_rate) public_share coef:", round(ps_coef_sec, 4),
    "(p =", round(ps_p_sec, 4), ")\n")

# ── 6. SPATIAL CHECK ──────────────────────────────────────────────────────────
cat("\n=== SPATIAL CHECK ===\n")
sf_use_s2(FALSE)
poly <- st_read(GEO, quiet = TRUE)
cat("Polygon rows:", nrow(poly), "\n")

# Align poly to dat
poly_matched <- poly[match(dat$dist_lgd, poly$dist_lgd_anl), ]
cat("Matched polygons:", sum(!is.na(poly_matched$dist_lgd_anl)), "of", nrow(dat), "\n")

# Queen contiguity neighbours
nb <- poly2nb(poly_matched, queen = TRUE)
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# Moran's I on primary model residuals
resid_primary <- residuals(fit_primary, type = "response")
moran_res <- moran.test(resid_primary, lw, zero.policy = TRUE)
cat("\nMoran's I on primary model residuals:\n")
print(moran_res)
moran_stat <- moran_res$statistic
moran_p    <- moran_res$p.value
cat("Moran I =", round(moran_stat, 4), ", p =", round(moran_p, 4), "\n")

# Spatial error model if Moran p < 0.05
spatial_model_fitted <- FALSE
if (moran_p < 0.05) {
  cat("\nMoran p < 0.05: fitting spatial error model (errorsarlm)...\n")
  # errorsarlm uses the untransformed response on logit scale via OLS surrogate;
  # use the logit of y_surge as the response for the spatial model
  dat$y_surge_logit <- log(dat$y_surge / (1 - dat$y_surge))
  f_sar <- y_surge_logit ~ public_share + lat_sc + logdens_sc +
           pct_urban_notif + region + lognotif_sc
  fit_sar <- errorsarlm(f_sar, data = dat, listw = lw, zero.policy = TRUE)
  cat("\n--- Spatial Error Model summary ---\n")
  print(summary(fit_sar))
  ps_sar   <- coef(fit_sar)["public_share"]
  ps_p_sar <- summary(fit_sar)$Coef["public_share", "Pr(>|z|)"]
  cat("\nSpatial error model: public_share coef =", round(ps_sar, 4),
      ", p =", round(ps_p_sar, 4), "\n")
  cat("Public_share association",
      ifelse(ps_p_sar < 0.05, "SURVIVES", "does NOT survive"),
      "spatial adjustment.\n")
  spatial_model_fitted <- TRUE
  # Add spatial model to output table
  cf_sar <- summary(fit_sar)$Coef
  tbl_sar <- tibble(
    term     = rownames(cf_sar),
    estimate = cf_sar[, "Estimate"],
    se       = cf_sar[, "Std. Error"],
    z        = cf_sar[, "z value"],
    p        = cf_sar[, "Pr(>|z|)"],
    ci_lo    = cf_sar[, "Estimate"] - 1.96 * cf_sar[, "Std. Error"],
    ci_hi    = cf_sar[, "Estimate"] + 1.96 * cf_sar[, "Std. Error"],
    model    = "surge_artefact_spatial_error"
  )
} else {
  cat("\nMoran p >= 0.05: no spatial autocorrelation in residuals; spatial model not needed.\n")
}

# ── 7. Save coefficient table ─────────────────────────────────────────────────
out_tbl <- bind_rows(tbl_primary, tbl_secondary)
if (spatial_model_fitted) out_tbl <- bind_rows(out_tbl, tbl_sar)
write_csv(out_tbl, file.path(OUT, "o3_heterogeneity.csv"))
cat("\nSaved:", file.path(OUT, "o3_heterogeneity.csv"), "\n")

# ── 8. Plot ───────────────────────────────────────────────────────────────────
# Panel A: Forest plot of primary model (mean coefficients, excluding intercept)
plot_dat <- tbl_primary |>
  filter(term != "(Intercept)") |>
  mutate(
    term_label = dplyr::recode(term,
      "public_share"      = "Public share",
      "lat_sc"            = "Latitude (scaled)",
      "logdens_sc"        = "Log density (scaled)",
      "pct_urban_notif"   = "% urban notif.",
      "lognotif_sc"       = "Log total notif. (scaled)",
      "regionEast"        = "Region: East",
      "regionNorth"       = "Region: North",
      "regionNortheast"   = "Region: Northeast",
      "regionOther"       = "Region: Other",
      "regionSouth"       = "Region: South",
      "regionWest"        = "Region: West"
    ),
    sig = ifelse(p < 0.05, "p < 0.05", "p >= 0.05")
  )

p_forest <- ggplot(plot_dat, aes(x = estimate, y = reorder(term_label, estimate),
                                  colour = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.3, linewidth = 0.6) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c("p < 0.05" = "#c0392b", "p >= 0.05" = "#7f8c8d"),
                      name = NULL) +
  labs(
    title    = "Primary model: surge_artefact_rate (beta regression)",
    subtitle = "Coefficients on logit scale with 95% CIs",
    x        = "Coefficient (logit scale)",
    y        = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# Panel B: Partial-effect plot of artefact rate vs public_share
ps_seq <- seq(min(dat$public_share), max(dat$public_share), length.out = 100)
newdat_pe <- data.frame(
  public_share    = ps_seq,
  lat_sc          = xbar$lat_sc,
  logdens_sc      = xbar$logdens_sc,
  pct_urban_notif = xbar$pct_urban_notif,
  region          = factor(modal_region, levels = levels(dat$region)),
  lognotif_sc     = xbar$lognotif_sc
)
# Bootstrap CIs for the partial effect (parametric)
pred_pe    <- predict(fit_primary, newdata = newdat_pe, type = "response")
# Simulate from vcov for approximate CIs
set.seed(SEED)
mu_vec <- coef(fit_primary)[names(coef(fit_primary)) %in%
            c("(Intercept)", "public_share", "lat_sc", "logdens_sc",
              "pct_urban_notif", paste0("region", levels(dat$region)[-1]),
              "lognotif_sc")]
# Build model matrix for partial-effect newdat
mm_pe <- model.matrix(~ public_share + lat_sc + logdens_sc +
                         pct_urban_notif + region + lognotif_sc,
                       data = newdat_pe)
# Use only mean part of betareg coef
coef_mean <- coef(fit_primary, model = "mean")
eta_pe    <- as.numeric(mm_pe %*% coef_mean)
pred_pe2  <- plogis(eta_pe)

vc_mean <- vcov(fit_primary)[names(coef_mean), names(coef_mean)]
sims    <- MASS::mvrnorm(2000, mu = coef_mean, Sigma = vc_mean)
pred_sims <- plogis(sims %*% t(mm_pe))
ci_lo_pe  <- apply(pred_sims, 2, quantile, 0.025)
ci_hi_pe  <- apply(pred_sims, 2, quantile, 0.975)

pe_df <- data.frame(
  public_share = ps_seq,
  pred         = pred_pe2,
  ci_lo        = ci_lo_pe,
  ci_hi        = ci_hi_pe
)

p_partial <- ggplot() +
  geom_ribbon(data = pe_df, aes(x = public_share, ymin = ci_lo, ymax = ci_hi),
              fill = "#2980b9", alpha = 0.2) +
  geom_line(data = pe_df, aes(x = public_share, y = pred),
            colour = "#2980b9", linewidth = 1) +
  geom_rug(data = dat, aes(x = public_share), sides = "b",
           alpha = 0.3, colour = "grey40") +
  labs(
    title    = "Partial effect: artefact rate vs public share",
    subtitle = "Other predictors held at sample mean/modal region",
    x        = "Public share (fraction of public-sector notifications)",
    y        = "Predicted surge_artefact_rate"
  ) +
  theme_bw(base_size = 11)

# Combine panels
p_combined <- p_forest + p_partial +
  plot_annotation(
    title   = "District-level predictors of the TB reporting artefact",
    caption = "Beta regression (logit link). Squeeze: y' = (y*(n-1)+0.5)/n. n = 744 districts."
  )

ggsave(file.path(OUT, "o3_heterogeneity.png"), p_combined,
       width = 13, height = 6, dpi = 150)
cat("Saved:", file.path(OUT, "o3_heterogeneity.png"), "\n")

# ── 9. VERIFICATION ───────────────────────────────────────────────────────────
cat("\n", strrep("=", 60), "\n")
cat("VERIFICATION\n")
cat(strrep("=", 60), "\n")

cat("n districts in regression:", nrow(dat), "\n")

cat("\nPRIMARY model (surge_artefact_rate) key coefficients:\n")
print(as.data.frame(tbl_primary[, c("term","estimate","ci_lo","ci_hi","p")]),
      digits = 4)

cat("\nPUBLIC_SHARE association:\n")
cat("  Direction:", ifelse(ps_coef > 0, "POSITIVE", "NEGATIVE"), "\n")
cat("  Logit-scale coefficient:", round(ps_coef, 4), "\n")
cat("  Effect (artefact rate change, ps 0.5->0.9):", round(effect_size, 4), "\n")
cat("  p-value:", round(tbl_primary$p[tbl_primary$term == "public_share"], 4), "\n")

cat("\nMoran's I on primary model residuals:\n")
cat("  I =", round(moran_stat, 4), ", p =", round(moran_p, 4), "\n")
if (spatial_model_fitted) {
  cat("  Spatial error model fitted.\n")
  cat("  public_share in spatial model: coef =", round(ps_sar, 4),
      ", p =", round(ps_p_sar, 4), "\n")
  cat("  Association", ifelse(ps_p_sar < 0.05, "SURVIVES", "does NOT survive"),
      "spatial adjustment.\n")
} else {
  cat("  No spatial model needed (Moran p >= 0.05).\n")
}

# VERDICT
cat("\n--- VERDICT ---\n")
ps_sig  <- tbl_primary$p[tbl_primary$term == "public_share"] < 0.05
lat_sig <- tbl_primary$p[tbl_primary$term == "lat_sc"]       < 0.05
den_sig <- tbl_primary$p[tbl_primary$term == "logdens_sc"]   < 0.05
notif_sig <- tbl_primary$p[tbl_primary$term == "lognotif_sc"] < 0.05

sig_terms <- tbl_primary$term[tbl_primary$p < 0.05 & tbl_primary$term != "(Intercept)"]
sig_labels <- dplyr::recode(sig_terms,
  "public_share" = "public_share", "lat_sc" = "latitude",
  "logdens_sc" = "log_density", "pct_urban_notif" = "pct_urban_notif",
  "lognotif_sc" = "log_total_notif"
)

cat("Significant predictors (p < 0.05):", paste(sig_labels, collapse = ", "), "\n")
cat("Public-share hypothesis",
    ifelse(ps_sig & ps_coef > 0, "SUPPORTED: higher public share -> more artefact.",
    ifelse(ps_sig & ps_coef < 0, "SUPPORTED (inverse): higher public share -> LESS artefact.",
    "NOT supported (p >= 0.05).")), "\n")

cat("\nScript completed successfully.\n")

## ===========================================================================
## from 35_heterogeneity_adjusted.R
## ===========================================================================

# 35_heterogeneity_adjusted.R
# M5 peer-review point: does the public_share -> artefact association survive
# adjustment for district socio-demography?
#
# Covariates added (Census 2011 PCA via SHRUG/Harvard Dataverse):
#   literacy_pct   : literate persons / total persons
#   scst_pct       : (SC + ST) persons / total persons
#
# Source: SHRUG v2 pc11_pca_clean_pc11dist.tab (File ID 10742786)
# DOI: 10.7910/DVN/DPESAK
#
# Outputs:
#   derived/modifiers_extended.rds            -- modifiers + socio-demographic covariates
#   out/o3_heterogeneity_adjusted.csv         -- coefficient table (unadjusted + adjusted)



# ── Packages ─────────────────────────────────────────────────────────────────
for (pkg in c("betareg", "sf", "spdep", "spatialreg")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# ── 1. Load core inputs ───────────────────────────────────────────────────────
payoff    <- read_csv(file.path(OUT, "o3_payoff.csv"), show_col_types = FALSE)
modifiers <- readRDS(file.path(DER, "modifiers.rds"))

cat("Payoff rows:", nrow(payoff), "\n")
cat("Modifiers rows:", nrow(modifiers), "\n")

# ── 2. Build socio-demographic covariates from Census 2011 PCA ─────────────────
cat("\n=== STEP 2: Loading Census 2011 PCA data ===\n")

pca_file <- file.path(DER, "pc11_pca_district.tab")
if (!file.exists(pca_file)) {
  stop("pc11_pca_district.tab not found in derived/. Run download step first.")
}

pca <- read.delim(pca_file, sep = "\t", header = TRUE,
                  colClasses = "character", check.names = FALSE)
cat("PCA rows loaded:", nrow(pca), "\n")

# Convert key columns to numeric
pca$pc11_state_id    <- as.integer(pca$pc11_state_id)
pca$pc11_district_id <- as.integer(pca$pc11_district_id)
pca$pc11_pca_tot_p   <- as.numeric(pca$pc11_pca_tot_p)
pca$pc11_pca_p_sc    <- as.numeric(pca$pc11_pca_p_sc)
pca$pc11_pca_p_st    <- as.numeric(pca$pc11_pca_p_st)
pca$pc11_pca_p_lit   <- as.numeric(pca$pc11_pca_p_lit)

# Compute proportions
pca_cov <- pca |>
  mutate(
    literacy_pct = pc11_pca_p_lit / pc11_pca_tot_p,
    scst_pct     = (pc11_pca_p_sc + pc11_pca_p_st) / pc11_pca_tot_p
  ) |>
  select(pc11_state_id, pc11_district_id, literacy_pct, scst_pct)

cat("PCA covariate summary:\n")
print(summary(pca_cov[, c("literacy_pct", "scst_pct")]))

# ── 3. Match covariates to modifiers via stcode11 + dtcode11 ─────────────────
cat("\n=== STEP 3: Matching covariates to modifiers ===\n")

modifiers_ext <- modifiers |>
  left_join(pca_cov,
            by = c("stcode11" = "pc11_state_id",
                   "dtcode11" = "pc11_district_id"))

n_matched <- sum(!is.na(modifiers_ext$literacy_pct))
cat("Matched:", n_matched, "of", nrow(modifiers_ext), "districts\n")
cat("Unmatched:", nrow(modifiers_ext) - n_matched,
    "(likely post-2011 district splits)\n")

# Save extended modifiers
saveRDS(modifiers_ext, file.path(DER, "modifiers_extended.rds"))
cat("Saved: derived/modifiers_extended.rds\n")

# ── 4. Join payoff + extended modifiers ───────────────────────────────────────
dat_all <- inner_join(payoff, modifiers_ext, by = "dist_lgd")
cat("\nAfter inner join:", nrow(dat_all), "districts\n")

# Key variables for base model
key_vars_base <- c("public_share", "latitude", "density_2024",
                   "pct_urban_notif", "region", "total_notif_2022_2025",
                   "surge_artefact_rate", "sign_flip_rate")

# Key variables for adjusted model (base + socio-demo)
key_vars_adj <- c(key_vars_base, "literacy_pct", "scst_pct")

dat_base <- dat_all[complete.cases(dat_all[, key_vars_base]), ]
dat_adj  <- dat_all[complete.cases(dat_all[, key_vars_adj]), ]

cat("n districts (base model - all covariates present):", nrow(dat_base), "\n")
cat("n districts (adjusted model - also has socio-demo):", nrow(dat_adj), "\n")
cat("Districts dropped due to missing socio-demo:", nrow(dat_base) - nrow(dat_adj), "\n")

# ── 5. Construct predictors ───────────────────────────────────────────────────
prepare_dat <- function(d) {
  N <- nrow(d)
  d |>
    mutate(
      region          = factor(region),
      lat_sc          = as.numeric(scale(latitude)),
      logdens_sc      = as.numeric(scale(log(density_2024 + 1))),
      lognotif_sc     = as.numeric(scale(log(total_notif_2022_2025 + 1))),
      literacy_sc     = as.numeric(scale(literacy_pct)),
      scst_sc         = as.numeric(scale(scst_pct)),
      y_surge         = (surge_artefact_rate * (N - 1) + 0.5) / N,
      y_flip          = (sign_flip_rate      * (N - 1) + 0.5) / N
    )
}

dat_base <- prepare_dat(dat_base)
dat_adj  <- prepare_dat(dat_adj)

# ── 6. Extract coefficient helper ─────────────────────────────────────────────
extract_betareg <- function(fit, model_name) {
  cf      <- coef(summary(fit))$mean
  ci      <- confint(fit)
  mean_terms <- rownames(cf)
  ci_mean <- ci[mean_terms, , drop = FALSE]
  tibble(
    term     = mean_terms,
    estimate = cf[, "Estimate"],
    se       = cf[, "Std. Error"],
    z        = cf[, "z value"],
    p        = cf[, "Pr(>|z|)"],
    ci_lo    = ci_mean[, 1],
    ci_hi    = ci_mean[, 2],
    model    = model_name
  )
}

# ── 7. BASE MODEL (replication of 30_heterogeneity.R primary model) ───────────
cat("\n=== BASE MODEL (unadjusted, n =", nrow(dat_base), ") ===\n")
f_base <- y_surge ~ public_share + lat_sc + logdens_sc +
          pct_urban_notif + region + lognotif_sc

fit_base <- betareg(f_base, data = dat_base, link = "logit")
tbl_base <- extract_betareg(fit_base, "surge_unadjusted")

ps_base_est <- tbl_base$estimate[tbl_base$term == "public_share"]
ps_base_lo  <- tbl_base$ci_lo[tbl_base$term == "public_share"]
ps_base_hi  <- tbl_base$ci_hi[tbl_base$term == "public_share"]
ps_base_p   <- tbl_base$p[tbl_base$term == "public_share"]

cat("public_share (unadjusted): est =", round(ps_base_est, 4),
    ", 95%CI [", round(ps_base_lo, 4), ",", round(ps_base_hi, 4), "]",
    ", p =", round(ps_base_p, 4), "\n")

# ── 8. ADJUSTED MODEL (adds literacy_sc + scst_sc) ────────────────────────────
cat("\n=== ADJUSTED MODEL (+ literacy_pct + scst_pct, n =", nrow(dat_adj), ") ===\n")
f_adj <- y_surge ~ public_share + lat_sc + logdens_sc +
         pct_urban_notif + region + lognotif_sc +
         literacy_sc + scst_sc

fit_adj <- betareg(f_adj, data = dat_adj, link = "logit")
tbl_adj <- extract_betareg(fit_adj, "surge_adjusted_sociodem")

ps_adj_est <- tbl_adj$estimate[tbl_adj$term == "public_share"]
ps_adj_lo  <- tbl_adj$ci_lo[tbl_adj$term == "public_share"]
ps_adj_hi  <- tbl_adj$ci_hi[tbl_adj$term == "public_share"]
ps_adj_p   <- tbl_adj$p[tbl_adj$term == "public_share"]

cat("public_share (adjusted): est =", round(ps_adj_est, 4),
    ", 95%CI [", round(ps_adj_lo, 4), ",", round(ps_adj_hi, 4), "]",
    ", p =", round(ps_adj_p, 4), "\n")

# ── 9. ADJUSTED MODEL also on the restricted sample (same n as adj) ───────────
# Re-run base model on the dat_adj sample to ensure fair comparison
cat("\n=== BASE MODEL on RESTRICTED SAMPLE (same districts as adjusted, n =",
    nrow(dat_adj), ") ===\n")
fit_base_restr <- betareg(f_base, data = dat_adj, link = "logit")
tbl_base_restr <- extract_betareg(fit_base_restr, "surge_unadjusted_restricted")

ps_base_r_est <- tbl_base_restr$estimate[tbl_base_restr$term == "public_share"]
ps_base_r_lo  <- tbl_base_restr$ci_lo[tbl_base_restr$term == "public_share"]
ps_base_r_hi  <- tbl_base_restr$ci_hi[tbl_base_restr$term == "public_share"]
ps_base_r_p   <- tbl_base_restr$p[tbl_base_restr$term == "public_share"]

cat("public_share (base, restricted sample): est =", round(ps_base_r_est, 4),
    ", 95%CI [", round(ps_base_r_lo, 4), ",", round(ps_base_r_hi, 4), "]",
    ", p =", round(ps_base_r_p, 4), "\n")

# ── 10. Save coefficient table ─────────────────────────────────────────────────
out_tbl <- bind_rows(tbl_base, tbl_base_restr, tbl_adj)
write_csv(out_tbl, file.path(OUT, "o3_heterogeneity_adjusted.csv"))
cat("\nSaved:", file.path(OUT, "o3_heterogeneity_adjusted.csv"), "\n")

# ── 11. VERIFICATION BLOCK ────────────────────────────────────────────────────
cat("\n", strrep("=", 70), "\n")
cat("VERIFICATION\n")
cat(strrep("=", 70), "\n")

cat("\nCOVARIATE SOURCES:\n")
cat("  FETCHED:   SHRUG v2 Census 2011 PCA (Harvard Dataverse DOI 10.7910/DVN/DPESAK)\n")
cat("             Variables: literacy_pct (literate / total), scst_pct ((SC+ST) / total)\n")
cat("  FAILED:    NFHS-5 district data - no downloadable district-level file found\n")
cat("             on Harvard Dataverse (only state-level tab); data.gov.in API key required;\n")
cat("             Zenodo/Kaggle search did not yield direct CSV access.\n")

cat("\nMATCH RATE:\n")
cat("  n_matched:", n_matched, "of", nrow(modifiers_ext), "districts (747)\n")
cat("  Unmatched:", nrow(modifiers_ext) - n_matched,
    "(post-2011 district splits lacking Census 2011 codes)\n")
cat("  Districts in adjusted model:", nrow(dat_adj),
    "(base model uses", nrow(dat_base), ")\n")

cat("\nPUBLIC_SHARE COEFFICIENT:\n")
cat("  BEFORE adjustment (full sample, n =", nrow(dat_base), "):\n")
cat("    est =", round(ps_base_est, 4),
    ", 95%CI [", round(ps_base_lo, 4), ",", round(ps_base_hi, 4), "]",
    ", p =", round(ps_base_p, 4), "\n")

cat("  BEFORE adjustment (restricted sample, n =", nrow(dat_adj), "):\n")
cat("    est =", round(ps_base_r_est, 4),
    ", 95%CI [", round(ps_base_r_lo, 4), ",", round(ps_base_r_hi, 4), "]",
    ", p =", round(ps_base_r_p, 4), "\n")

cat("  AFTER  adjustment (+ literacy + SC/ST, n =", nrow(dat_adj), "):\n")
cat("    est =", round(ps_adj_est, 4),
    ", 95%CI [", round(ps_adj_lo, 4), ",", round(ps_adj_hi, 4), "]",
    ", p =", round(ps_adj_p, 4), "\n")

# Verdict
same_sign  <- sign(ps_adj_est) == sign(ps_base_est)
adj_sig    <- ps_adj_p < 0.05
verdict    <- if (same_sign && adj_sig) "SURVIVES" else if (same_sign && !adj_sig) "ATTENUATED (same sign, p >= 0.05)" else "DOES-NOT-SURVIVE (sign change)"

cat("\n--- VERDICT ---\n")
cat(verdict, "\n")
cat("public_share association with surge_artefact_rate",
    ifelse(same_sign && adj_sig,
           "remains the same sign and statistically significant after adding",
           "is ATTENUATED or reversed after adding"),
    "district literacy and SC/ST share.\n")

cat("\nADJUSTED MODEL - all coefficients:\n")
print(as.data.frame(tbl_adj[, c("term","estimate","ci_lo","ci_hi","p")]), digits = 4)

cat("\nScript completed successfully.\n")
