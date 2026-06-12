# 06_figures_and_tables.R — Assemble manuscript figures and tables
#
# Renumbered, path-normalized reorganization of the working pipeline scripts: 40_manuscript_assets.R.
# Computation is unchanged from the scripts that produced the shipped results in out/.

source(here::here("code", "00_setup.R"))


## ===========================================================================
## from 40_manuscript_assets.R
## ===========================================================================

## 40_manuscript_assets.R - build the O1 seasonal figure + O2 coefficient table for the manuscript

suppressMessages({library(feasts); library(fabletools); library(ggplot2); library(glmmTMB)})
set.seed(SEED)

## ---- Part 1: O1 national + sector seasonal figure ----
dm  <- readRDS(file.path(DER, "panel_dm.rds")) |> dplyr::filter(complete)
agg <- dm |> index_by(month) |> summarise(public = sum(public), private = sum(private), total = sum(total))

seas_by_month <- function(var) {
  ag1  <- agg |> dplyr::transmute(y = .data[[var]])   # 1-measure tsibble, fixed name 'y'
  comp <- ag1 |> model(STL(log(y + 1) ~ trend() + season(window = 11), robust = TRUE)) |> components()
  comp |> tibble::as_tibble() |>
    mutate(m = lubridate::month(as.Date(month))) |>
    group_by(m) |> summarise(sf = mean(season_year), .groups = "drop") |>
    mutate(series = var, pct = 100 * (exp(sf) - 1))
}
sdf <- dplyr::bind_rows(seas_by_month("total"), seas_by_month("public"), seas_by_month("private"))
sdf$series <- factor(sdf$series, levels = c("total","public","private"), labels = c("Total","Public","Private"))

p <- ggplot(sdf, aes(m, pct, colour = series)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_line(linewidth = 0.9) + geom_point(size = 1.4) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_colour_manual(values = c(Total = "#222222", Public = "#08519c", Private = "#b2182b"), name = NULL) +
  labs(x = NULL, y = "Seasonal deviation from trend (%)",
       title = "Seasonal pattern of TB notification, by sector",
       subtitle = "STL multiplicative seasonal factor; May peak; public amplitude exceeds private") +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"), legend.position = "top")
ggsave(file.path(OUT, "o1_seasonal_plot.png"), p, width = 8, height = 4.8, dpi = 150)
cat("saved o1_seasonal_plot.png\n")

## ---- Part 2: O2 sector-model coefficients (rate ratios + 95% CI) ----
fp <- readRDS(file.path(DER, "fit_public.rds"))
fr <- readRDS(file.path(DER, "fit_private_ar1.rds"))
get_coefs <- function(m, sector) {
  cf <- summary(m)$coefficients$cond
  terms <- intersect(c("q_end","fy_end","campaign","t2m_l3","solar_l3","rh_l0"), rownames(cf))
  est <- cf[terms, "Estimate"]; se <- cf[terms, "Std. Error"]; pv <- cf[terms, 4]
  labs <- c(q_end = "Fiscal quarter-end", fy_end = "Fiscal year-end (March)", campaign = "100-day campaign",
            t2m_l3 = "Temperature (lag 3m)", solar_l3 = "Solar radiation (lag 3m)", rh_l0 = "Relative humidity")
  tibble::tibble(Sector = sector, Term = labs[terms],
                 Block = ifelse(terms %in% c("q_end","fy_end","campaign"), "Reporting calendar / campaign", "Climate"),
                 `Rate ratio` = round(exp(est), 3),
                 `95% CI` = sprintf("%.3f-%.3f", exp(est - 1.96*se), exp(est + 1.96*se)),
                 p = signif(pv, 2))
}
o2c <- dplyr::bind_rows(get_coefs(fp, "Public"), get_coefs(fr, "Private"))
readr::write_csv(o2c, file.path(OUT, "o2_coefficients.csv"))
print(as.data.frame(o2c))
cat("saved o2_coefficients.csv\n")
