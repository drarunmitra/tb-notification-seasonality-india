# ----------------------------------------------------------------------------
# PROVENANCE STEP. This script rebuilds inputs from TU-level Ni-kshay data that
# is NOT redistributable. The shipped files in data/ and files/ are the
# quick-start entry point; you only need this script to rebuild from source.
# To run it, set the environment variable TU_SOURCE_ROOT to your local source
# tree (the parent of side_quest_subdistrict/) before sourcing 00_setup.R.
# ----------------------------------------------------------------------------

# 01_build_district_panel.R — Build the district-month notification panel
#
# Renumbered, path-normalized reorganization of the working pipeline scripts: 00_crosswalk.R, 00b_parent_remap.R, 01_panel.R.
# Computation is unchanged from the scripts that produced the shipped results in out/.

source(here::here("code", "00_setup.R"))


## ===========================================================================
## from 00_crosswalk.R
## ===========================================================================


map <- read_csv(file.path(SQ,"tu_subdistrict_map_google.csv"), show_col_types=FALSE) |> select(tu_unique, sd_id)
sdm <- read_csv(file.path(SQ,"subdistrict_analysis.csv"), show_col_types=FALSE) |> select(sd_id, dist_lgd, stname, region)
ur  <- read_csv(file.path(SQ,"tu_urbanrural.csv"), show_col_types=FALSE) |> select(tu_unique, urban_rural)
excl<- read_csv(file.path(SQ,"excluded_tus.csv"), show_col_types=FALSE)$tu_unique
tu2d <- map |> left_join(sdm, by="sd_id") |> left_join(ur, by="tu_unique") |>
  filter(!tu_unique %in% excl, !is.na(dist_lgd))
frame <- tu2d |> distinct(dist_lgd, stname, region)
saveRDS(tu2d, file.path(DER,"tu2district.rds"))
saveRDS(frame, file.path(DER,"district_frame.rds"))
cat("TUs mapped:", nrow(tu2d), "| districts:", nrow(frame), "\n")
unmapped <- map |> left_join(sdm,by="sd_id") |> filter(is.na(dist_lgd)) |> nrow()
cat("rows with no dist_lgd (pre-exclusion):", unmapped, "\n")

## ===========================================================================
## from 00b_parent_remap.R
## ===========================================================================

# 00b_parent_remap.R
# Purpose: Fold 18 recently-created districts (no TUs) into their parent districts
# so that the population denominator aligns with where notifications are recorded.
# Also excludes Pakistan-administered Kashmir (dist_lgd == 0).
#
# Outputs:
#   derived/parent_remap.csv        -- 18-row child->parent lookup
#   derived/dist_lgd_anl_lookup.rds -- full code->analysis code tibble (all codes)
#   derived/district_pop_anl.rds    -- long pop table: dist_lgd_anl, yr, pop
#   derived/district_anl.gpkg       -- dissolved district polygons
#
# Population shape: LONG (dist_lgd_anl, yr, pop) for easy panel joins.


library(sf)

# ── 1. Read inputs ─────────────────────────────────────────────────────────────
sd_csv <- read_csv(
  file.path(SQ, "subdistrict_analysis.csv"),
  show_col_types = FALSE
)

dist_frame <- readRDS(file.path(DER, "district_frame.rds"))

# ── 2. Build normalisation helper ──────────────────────────────────────────────
norm <- function(x) {
  x |>
    tolower() |>
    trimws() |>
    gsub(pattern = "-", replacement = " ", fixed = TRUE) |>
    gsub(pattern = "\\s+", replacement = " ")
}

# District-level name<->code table (one row per unique dtname/stname/dist_lgd)
dist_codes <- sd_csv |>
  distinct(dtname, stname, dist_lgd) |>
  mutate(
    norm_dt = norm(dtname),
    norm_st = norm(stname)
  )

# ── 3. Hardcoded lineage table ─────────────────────────────────────────────────
# child_name, child_state, parent_name, note
# Names exactly as they appear in the specification; normalised for matching
lineage_raw <- tribble(
  ~child_name,             ~child_state,         ~parent_name,        ~note,
  "malerkotla",            "punjab",              "sangrur",           "",
  "noney",                 "manipur",             "tamenglong",        "",
  "pherzawl",              "manipur",             "churachandpur",     "",
  "kamjong",               "manipur",             "ukhrul",            "",
  "kra daadi",             "arunachal pradesh",   "kurung kumey",      "",
  "shi yomi",              "arunachal pradesh",   "west siang",        "",
  "pakke kessang",         "arunachal pradesh",   "east kameng",       "",
  "siang",                 "arunachal pradesh",   "west siang",        "split from West+East Siang; folded to West Siang",
  "lower siang",           "arunachal pradesh",   "west siang",        "split from West+East Siang; folded to West Siang",
  "kamle",                 "arunachal pradesh",   "lower subansiri",   "split from Lower+Upper Subansiri; folded to Lower Subansiri",
  "khawzawl",              "mizoram",             "champhai",          "",
  "saitual",               "mizoram",             "aizawl",            "split from Aizawl+Champhai; folded to Aizawl",
  "noklak",                "nagaland",            "tuensang",          "",
  "shamator",              "nagaland",            "tuensang",          "",
  "niuland",               "nagaland",            "dimapur",           "",
  "tseminyu",              "nagaland",            "kohima",            "",
  "east jaintia hills",    "meghalaya",           "west jaintia hills","parent 'Jaintia Hills' now split; folded to West Jaintia Hills",
  "south west khasi hills","meghalaya",           "west khasi hills",  ""
)

# ── 4. Resolve child and parent dist_lgd codes ────────────────────────────────
resolve_code <- function(norm_name, norm_state, dc) {
  hits <- dc |>
    filter(norm_dt == norm_name, norm_st == norm_state)
  if (nrow(hits) == 0) return(list(code = NA_integer_, dtname = NA_character_))
  if (nrow(hits) > 1) warning("Multiple matches for: ", norm_name, " / ", norm_state)
  list(code = hits$dist_lgd[1], dtname = hits$dtname[1])
}

lookup_rows <- lineage_raw |>
  rowwise() |>
  mutate(
    child_res  = list(resolve_code(child_name,  child_state, dist_codes)),
    parent_res = list(resolve_code(parent_name, child_state, dist_codes)),
    child_dist_lgd  = child_res$code,
    child_dtname    = child_res$dtname,
    parent_dist_lgd = parent_res$code,
    parent_dtname   = parent_res$dtname
  ) |>
  ungroup() |>
  select(
    child_dist_lgd, child_dtname, stname = child_state,
    parent_dist_lgd, parent_dtname, note
  )

# ── 5. Validate: check for unresolved codes ────────────────────────────────────
unresolved_child  <- lookup_rows |> filter(is.na(child_dist_lgd))
unresolved_parent <- lookup_rows |> filter(is.na(parent_dist_lgd))

if (nrow(unresolved_child) > 0) {
  cat("WARNING - could not resolve CHILD codes:\n")
  print(unresolved_child)
} else {
  cat("OK: all 18 child dist_lgd codes resolved\n")
}

if (nrow(unresolved_parent) > 0) {
  cat("WARNING - could not resolve PARENT codes:\n")
  print(unresolved_parent)
} else {
  cat("OK: all 18 parent dist_lgd codes resolved\n")
}

# ── 6. Validate: every parent must be in district_frame (747 frame) ────────────
parent_not_in_frame <- lookup_rows |>
  filter(!is.na(parent_dist_lgd)) |>
  filter(!parent_dist_lgd %in% dist_frame$dist_lgd)

if (nrow(parent_not_in_frame) > 0) {
  cat("WARNING - parent dist_lgd NOT in district_frame (747 frame):\n")
  print(parent_not_in_frame)
} else {
  cat("OK: all resolved parent codes exist in district_frame.rds\n")
}

# ── 7. Write derived/parent_remap.csv ──────────────────────────────────────────
write_csv(lookup_rows, file.path(FILES, "parent_remap.csv"))
cat("\nWritten: derived/parent_remap.csv (", nrow(lookup_rows), "rows)\n")

# ── 8. Build full dist_lgd_anl lookup (all codes present in subdistrict CSV) ──
all_codes <- dist_codes |>
  distinct(dist_lgd) |>
  mutate(dist_lgd_anl = case_when(
    dist_lgd == 0L ~
      NA_integer_,                         # PoK → excluded
    dist_lgd %in% lookup_rows$child_dist_lgd ~
      lookup_rows$parent_dist_lgd[match(dist_lgd, lookup_rows$child_dist_lgd)],
    TRUE ~
      dist_lgd
  ))

saveRDS(all_codes, file.path(FILES, "dist_lgd_anl_lookup.rds"))
cat("Written: derived/dist_lgd_anl_lookup.rds (", nrow(all_codes), "rows, covering all dist_lgd codes)\n")

# ── 9. Build analysis-vintage population ──────────────────────────────────────
pop_cols <- paste0("pop_", YEARS)

district_pop_anl <- sd_csv |>
  left_join(all_codes, by = "dist_lgd") |>
  filter(!is.na(dist_lgd_anl)) |>          # drop PoK
  group_by(dist_lgd_anl) |>
  summarise(across(all_of(pop_cols), sum, .names = "{.col}"), .groups = "drop") |>
  pivot_longer(
    cols      = all_of(pop_cols),
    names_to  = "yr",
    names_prefix = "pop_",
    values_to = "pop"
  ) |>
  mutate(yr = as.integer(yr))

saveRDS(district_pop_anl, file.path(DER, "district_pop_anl.rds"))
cat("Written: derived/district_pop_anl.rds  shape = long (dist_lgd_anl, yr, pop)\n")
cat("         rows:", nrow(district_pop_anl), " | years:", paste(sort(unique(district_pop_anl$yr)), collapse=","), "\n")

# ── 10. Re-dissolve geometry ──────────────────────────────────────────────────
# Use GEOS (sf_use_s2 = FALSE) for the dissolve: s2 rejects self-intersecting
# sub-district polygons that exist in the source data. We validate each row
# before union, then validate the dissolved result.
cat("\nReading subdistrict_analysis.gpkg ...\n")
sf_use_s2(FALSE)   # switch to GEOS; avoids s2 "Loop 0 not valid" on bad input polygons

sd_gpkg <- st_read(file.path(SQ, "subdistrict_analysis.gpkg"), quiet = TRUE)

# Validate individual sub-district geometries before dissolve
cat("Validating sub-district geometries ...\n")
sd_gpkg <- sd_gpkg |> st_make_valid()

district_anl_sf <- sd_gpkg |>
  left_join(all_codes, by = "dist_lgd") |>
  filter(!is.na(dist_lgd_anl)) |>
  group_by(dist_lgd_anl) |>
  summarise(geometry = st_union(geom), .groups = "drop") |>
  st_make_valid()

st_write(district_anl_sf, GEO,
         delete_dsn = TRUE, quiet = TRUE)
cat("Written: derived/district_anl.gpkg (", nrow(district_anl_sf), "features)\n")

# ── 11. VERIFICATION ──────────────────────────────────────────────────────────
cat("\n=== VERIFICATION ===\n")

# V1: child codes must NOT appear as keys in district_pop_anl
child_codes <- lookup_rows$child_dist_lgd[!is.na(lookup_rows$child_dist_lgd)]

pop_anl_wide <- district_pop_anl |>
  filter(yr == 2024L) |>
  select(dist_lgd_anl, pop)

child_in_pop <- child_codes[child_codes %in% pop_anl_wide$dist_lgd_anl]
child_in_gpkg <- child_codes[child_codes %in% district_anl_sf$dist_lgd_anl]

if (length(child_in_pop) == 0) {
  cat("OK: No child codes appear as keys in district_pop_anl\n")
} else {
  cat("FAIL: Child codes found in district_pop_anl:", paste(child_in_pop, collapse=", "), "\n")
}

if (length(child_in_gpkg) == 0) {
  cat("OK: No child codes appear as keys in district_anl.gpkg\n")
} else {
  cat("FAIL: Child codes found in district_anl.gpkg:", paste(child_in_gpkg, collapse=", "), "\n")
}

# V2: Population conservation
total_pop_anl <- district_pop_anl |>
  filter(yr == 2024L) |>
  summarise(pop = sum(pop, na.rm = TRUE)) |>
  pull(pop)

total_pop_raw <- sd_csv |>
  filter(dist_lgd != 0) |>
  summarise(pop = sum(pop_2024, na.rm = TRUE)) |>
  pull(pop)

cat("\nPOPULATION CONSERVATION (pop_2024):\n")
cat("  Raw sub-district total (excl. PoK):     ", format(round(total_pop_raw), big.mark=","), "\n")
cat("  Analysis district_pop_anl total:        ", format(round(total_pop_anl), big.mark=","), "\n")
diff_pop <- abs(total_pop_anl - total_pop_raw)
if (diff_pop < 1) {
  cat("  OK: Difference =", diff_pop, "(< 1, conserved)\n")
} else {
  cat("  FAIL: Difference =", diff_pop, "(> 1, population NOT conserved)\n")
}

# V3: Per-parent pop delta vs pre-remap pop
cat("\nPER-PARENT pop_2024 delta (affected parents only):\n")

# pre-remap: sum of sub-districts belonging to parent dist_lgd only (no children yet)
pre_remap_parent <- sd_csv |>
  filter(dist_lgd %in% unique(lookup_rows$parent_dist_lgd[!is.na(lookup_rows$parent_dist_lgd)])) |>
  group_by(dist_lgd) |>
  summarise(pop_pre = sum(pop_2024, na.rm = TRUE), .groups = "drop")

# pop that moved (children only)
child_pop_moved <- sd_csv |>
  filter(dist_lgd %in% child_codes) |>
  left_join(lookup_rows |> select(child_dist_lgd, parent_dist_lgd),
            by = c("dist_lgd" = "child_dist_lgd")) |>
  group_by(parent_dist_lgd) |>
  summarise(pop_added = sum(pop_2024, na.rm = TRUE), .groups = "drop")

# post-remap from analysis frame
post_remap_parent <- district_pop_anl |>
  filter(yr == 2024L,
         dist_lgd_anl %in% unique(lookup_rows$parent_dist_lgd[!is.na(lookup_rows$parent_dist_lgd)])) |>
  select(dist_lgd_anl, pop_post = pop)

delta_tbl <- pre_remap_parent |>
  left_join(child_pop_moved, by = c("dist_lgd" = "parent_dist_lgd")) |>
  left_join(post_remap_parent, by = c("dist_lgd" = "dist_lgd_anl")) |>
  left_join(distinct(dist_codes, dist_lgd, dtname), by = "dist_lgd") |>
  mutate(
    pop_added       = replace_na(pop_added, 0),
    expected_post   = pop_pre + pop_added,
    actual_post     = pop_post,
    matches         = abs(actual_post - expected_post) < 1
  ) |>
  select(dist_lgd, dtname, pop_pre, pop_added, expected_post, actual_post, matches)

print(delta_tbl, n = Inf)

if (all(delta_tbl$matches, na.rm = TRUE)) {
  cat("OK: All parent populations match expected (pre + child added)\n")
} else {
  cat("FAIL: Some parent populations do not match expected\n")
}

# V4: Geometry row count
n_distinct_anl <- n_distinct(all_codes$dist_lgd_anl[!is.na(all_codes$dist_lgd_anl)])
n_gpkg_rows    <- nrow(district_anl_sf)

cat("\nGEOMETRY ROW COUNT:\n")
cat("  Distinct dist_lgd_anl values:  ", n_distinct_anl, "\n")
cat("  Rows in district_anl.gpkg:     ", n_gpkg_rows, "\n")
if (n_gpkg_rows == n_distinct_anl) {
  cat("  OK: One polygon per analysis district\n")
} else {
  cat("  FAIL: Mismatch\n")
}

# V5: Expected district count
# 766 distinct dist_lgd in subdistrict_analysis
# minus 18 folded children
# minus 1 PoK (code 0)
# = 747 expected analysis districts
n_all_lgd  <- length(unique(sd_csv$dist_lgd))           # 766
n_expected <- n_all_lgd - 18L - 1L                      # 747

cat("\nDISTRICT COUNT CHECK:\n")
cat("  Distinct dist_lgd in subdistrict_analysis.csv:  ", n_all_lgd, "\n")
cat("  minus 18 folded child codes:                    -18\n")
cat("  minus 1 PoK code (dist_lgd == 0):               -1\n")
cat("  Expected analysis districts:                    ", n_expected, "\n")
cat("  Actual district_anl.gpkg row count:             ", n_gpkg_rows, "\n")
if (n_gpkg_rows == n_expected) {
  cat("  OK: Count matches expected (", n_expected, ")\n")
} else {
  cat("  NOTE: Count is", n_gpkg_rows, "vs expected", n_expected,
      "(review lookup_rows for duplicate parent mappings)\n")
}

cat("\n=== VERIFICATION COMPLETE ===\n")

## ===========================================================================
## from 01_panel.R
## ===========================================================================


suppressMessages({library(duckdb); library(DBI)})
tu2d <- readRDS(file.path(DER,"tu2district.rds"))
con <- dbConnect(duckdb::duckdb())
duckdb::duckdb_register(con, "tu2d", tu2d[,c("tu_unique","dist_lgd")])
p <- file.path(ROOT,"data/derived/tu_panel.parquet")
dm <- dbGetQuery(con, paste0(
  "SELECT t.dist_lgd, p.period, sum(p.public) public, sum(p.private) private, sum(p.total) total
   FROM read_parquet('",p,"') p JOIN tu2d t ON p.tu_unique=t.tu_unique
   GROUP BY t.dist_lgd, p.period ORDER BY t.dist_lgd, p.period"))
dbDisconnect(con, shutdown=TRUE)
dm <- dm |> mutate(month = yearmonth(as.Date(paste0(period,"-01")))) |>
  select(-period) |>
  as_tsibble(index=month, key=dist_lgd) |>
  fill_gaps(public=0L, private=0L, total=0L)
# reporting-lag completeness flag: trailing months whose national total < 85% of prior-12 median
natl <- dm |> index_by(month) |> summarise(total=sum(total)) |> as_tibble()
natl <- natl |> mutate(roll12 = sapply(seq_len(n()), function(i){j<-max(1,i-12):max(1,i-1); median(total[j])}), ratio=total/roll12)
inc <- logical(nrow(natl)); for(i in nrow(natl):1){ if(!is.na(natl$ratio[i]) && natl$ratio[i]<0.85) inc[i]<-TRUE else break }
incomplete_months <- natl$month[inc]
dm <- dm |> mutate(complete = !(month %in% incomplete_months))
saveRDS(dm, file.path(DER,"panel_dm.rds"))
cat("districts:", n_distinct(dm$dist_lgd), "| months:", length(unique(dm$month)),
    "| incomplete trailing months:", length(incomplete_months),
    paste(as.character(incomplete_months), collapse=","), "\n")
cat("rows:", nrow(dm), "| min total:", min(dm$total), "| any NA total:", any(is.na(dm$total)), "\n")
