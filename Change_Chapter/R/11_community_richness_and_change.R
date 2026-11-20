# ============================================================
# 11_community_richness_and_change.R
#
# Purpose
#   Community-level summaries across all modelled species:
#     Figure 1  Expected species richness per 1-km pixel in each atlas, and its
#               change.
#     Figure 2  Province-wide change for every species, ranked, by change class.
#     Figure 3  Number of species in each change class, province and per BCR.
#     Figure 4  Per 1-km pixel: species pool, and the number and share of it with
#               a well-supported increase / decrease (P >= PIXEL_P_PRIMARY).
#               One figure per pool threshold.
#     Figure 5  Sensitivity of Figure 4's shares: restricted to BCRs where 08
#               finds the species sufficiently observed, and with the
#               approximate P >= PIXEL_P_APPROX threshold.
#     Figure 6  Gains, losses and total change in species-level detection
#               probabilities, and a Bray-Curtis-style relative change.
#     Figure 7  Per 25-km hexagon: species with well-supported changes of at
#               least HEX_FOLD-fold (exact, from posterior draws).
#     Figure 8  Per 25-km hexagon: community offset (geometric-mean change
#               across species, with its interval) and species departures
#               from it (species doing well-supported better or worse than
#               the community by at least HEX_FOLD-fold).
#
# Self-contained: sources nothing from the workflow. All input locations are
# set in CONFIG; outputs go to their own subfolder of the output folder.
#
# Inputs (written by the main workflow)
#   predictions_<model>/<species>_1km.rds        07  per-pixel change summaries;
#                                                    25-km hex posterior draws;
#                                                    hyperparameters (fallback)
#   summaries_<model>/model_summaries.rds        07  hyperparameters (preferred)
#   rasters_<model>/<species>_a2.tif, _a3.tif    09  water-corrected, safe-date-
#                                                    trimmed posterior medians
#   change_estimates_<model>/regional_change_estimates_<model>.rds   08
#   birds/data_ready_for_analysis.rds            06  grid_OBBA2, bcr_sf, hex_grid_25km
#
# Expected richness (Figure 1)
#   PObs = P(>= 1 detection across POBS_N_SURVEYS 5-min point counts), computed
#   as on 09's page 14 (median expected count integrated over the omitted iid
#   variance under the fitted error family). Expected richness = sum of PObs
#   over species (Calabrese et al. 2014). From posterior medians: no uncertainty.
#
# Change classes (Figures 2 and 3), from 08's province / BCR estimates.
#   CLASS_PRECEDENCE = "direction": very likely > little change > likely >
#     uncertain (a precisely estimated small change is "very likely").
#   CLASS_PRECEDENCE = "magnitude": little change > very likely > likely >
#     uncertain (a change whose CI lies within the little-change band is
#     "little change" however certain its direction).
#   08's province-wide estimate for each species sums only the BCRs where that
#   species was sufficiently observed, so footprints differ among species.
#
# Pixel change (Figures 4 and 5)
#   07 saves, per pixel, the posterior median and 90% interval (5th / 95th
#   percentiles) of the change in expected count; the draws are not saved.
#   P(increase) is therefore approximated by a split normal (the half of the
#   interval facing zero sets the scale). This is exact at P = 0.95 (interval
#   bound at zero), so PIXEL_P_PRIMARY = 0.95 is equivalent to "the 90%
#   interval excludes zero". PIXEL_P_APPROX (0.80) relies on the normal shape
#   and is shown only as a sensitivity.
#   Several grid pixels can fall in one raster cell; P(increase) is averaged
#   over them before classifying.
#   Species pool: a species is present in a cell (counts 1) if its PObs exceeds
#   the threshold in either atlas AND it has a change summary there, so pool
#   and counts use the same cells. Shares = count / pool where the pool has
#   >= MIN_POOL_FOR_PROP species. "Changed" = increased + decreased.
#   Support-mask version: a species counts only in BCRs where 08 marks it as
#   sufficiently observed (sufficient_data), in both pool and counts.
#   For display, 1-km cells are aggregated by PLOT_AGG_FACT; shares are
#   recomputed from summed counts and pools (not averaged).
#   Counts reflect evidence of change, which depends on data: sparsely surveyed
#   regions show fewer well-supported changes for the same true change (see
#   script 16 for where comparable data in both atlases exist).
#
# Gains and losses (Figure 6), from PObs medians:
#   G = sum_s max(p3 - p2, 0), L = sum_s max(p2 - p3, 0), D = G + L,
#   richness change = G - L (Figure 1), relative change = D / sum_s (p2 + p3)
#   (Bray-Curtis form on detection probabilities; not compositional turnover).
#   Like Figure 1, these weight species by how much their PObs can move, so
#   saturated common species and very rare species contribute little.
#
# Hexagon magnitude assessment (Figure 7)
#   From each species' 25-km hex posterior draws (water-corrected; the
#   safe-date trim scales both atlases' draws almost equally, so ratios are
#   unaffected): P(A3 / A2 >= HEX_FOLD) and P(A3 / A2 <= 1 / HEX_FOLD), exact
#   over the draws. Supported if >= HEX_P. Hex pool: hex-mean PObs (over the
#   hex's modelled cells) > HEX_POOL_POBS in either atlas. Direction-only
#   counts (P(A3 > A2) >= HEX_P) are saved alongside for comparison.
#
# Community offset and species departures (Figure 8)
#   Offset: per posterior draw, the mean over a hex's pool species of
#   log(A3 / A2) (each species' log ratio capped at +/- OFFSET_LR_CAP, so near-
#   zero denominators cannot dominate); summarised as the posterior median and
#   90% interval of the geometric-mean ratio. Species are fit separately, so
#   pairing their draws by index samples the mean under independence; the
#   interval reflects species-level estimation uncertainty only. A bias shared
#   by all species (protocol, survey type, observers, timing) moves the offset
#   directly and is NOT in its interval: the offset needs external support
#   (e.g. BBS, a point-count-only refit) before it is read as real change.
#   Departures: a species departs upward if P(log ratio - offset >= log
#   HEX_FOLD) >= HEX_P, i.e. its (1 - HEX_P) quantile is >= offset +
#   log(HEX_FOLD) (downward symmetrically), with the offset fixed at its
#   posterior median (its uncertainty is small next to a single species').
#   Any multiplicative bias shared by all species cancels from departures, so
#   they show where species fared differently from the community. Like Figure
#   7's counts they depend on precision, so sparse hexes show fewer.
#   Hexes, BCRs: pooled offsets average species-hex log ratios per draw, which
#   keeps each species' spatial correlation across hexes.
#
# Caveats
#   - PObs conflates abundance and per-individual detectability.
#   - Only MODELLED species are counted.
#   - Species are classified one at a time; biases shared by all species
#     (protocol, observers, timing) are not represented.
#
# Outputs (OUT_DIR; <pool> = pobs010, pobs025, ...)
#   fig1_expected_richness_pixels.png, fig2_species_change_ranked.png,
#   fig3_change_classes.png, fig4_species_change_pixels_<pool>.png,
#   fig5_species_change_sensitivity_<pool>.png, fig6_gains_losses_pixels.png,
#   fig7_species_change_hex.png, fig8_community_offset_hex.png
#   richness_A2.tif, richness_A3.tif, richness_change_A3_minus_A2.tif
#   species_pool_change_<pool>.tif   (per variant: pool, increase, decrease,
#                                     share_increase, share_decrease, share_changed)
#   gains_losses.tif                 (gains, losses, total, relative)
#   hex_species_change_summary.csv (with majority BCR per hex),
#   hex_species_change_by_bcr.csv, species_hex_change.rds
#   species_community_table.csv, community_summaries.rds (incl. settings)
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# ============================================================
# CONFIG
# ============================================================

# Locations (all required).
PRED_DIR       <- "C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/Analysis/Ontario_BBA3_draft5/data_clean/model_output/predictions_PC_ARU_CL_SS"
RAST_DIR       <- "C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/Analysis/Ontario_BBA3_draft5/data_clean/model_output/rasters_PC_ARU_CL_SS"
CHANGE_TABLE   <- "C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/Analysis/Ontario_BBA3_draft5/data_clean/model_output/change_estimates_PC_ARU_CL_SS/regional_change_estimates_PC_ARU_CL_SS.rds"
DATA_READY_RDS <- "C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/Analysis/Ontario_BBA3_draft5/data_clean/birds/data_ready_for_analysis.rds"
OUTPUT_ROOT    <- "C:/Users/IlesD/Documents/Projects_Local/Ontario_Atlas/Change_Chapter/output"
OUT_DIR        <- file.path(OUTPUT_ROOT, "community_change")

# Probability of observation. Keep equal to 09's cfg$pobs_* so the maps agree
# with page 14 of the species PDFs.
POBS_N_SURVEYS         <- 25L
POBS_MARGINALIZE_TERMS <- c("kappa_diff", "site_iid")
POBS_GH_NODES          <- 41L

# Maps: pixels where neither atlas reaches this many expected species are not
# drawn (mostly open water), so the colour scale spans the land.
MIN_MAP_RICHNESS <- 10

# Change classes (Figures 2-3; see header).
P_VERY_LIKELY    <- 0.95
P_LIKELY         <- 0.80
STABLE_FOLD      <- 1.2       # little-change band: -17% to +20%
CLASS_PRECEDENCE <- "direction"   # "direction" or "magnitude"

# Pixel change (Figures 4-5). 07 summarises pixel posteriors with a 90%
# interval (5th / 95th percentiles).
PIXEL_CI_Z           <- stats::qnorm(0.95)
PIXEL_P_PRIMARY      <- 0.95      # exact from the saved quantiles
PIXEL_P_APPROX       <- 0.80      # split-normal approximation; sensitivity only
POOL_POBS_THRESHOLDS <- c(0.10, 0.25)   # species pool: PObs > threshold in either atlas
MIN_POOL_FOR_PROP    <- 5         # shares shown only where the pool has >= this many species

# Hexagon magnitude assessment (Figure 7).
HEX_FOLD          <- 1.2      # increase >= +20%; decrease <= 1/1.2 (-17%)
HEX_P             <- 0.95     # posterior support required (exact, from draws)
HEX_POOL_POBS     <- 0.25     # hex pool: hex-mean PObs > this in either atlas
HEX_MIN_CELL_PROP <- 0.25     # drop hexes with < this share of the fullest hex's modelled cells

# Community offset (Figure 8): cap on each species' |log(A3 / A2)| per draw, so a
# near-zero Atlas 2 (or Atlas 3) value cannot dominate a hex's mean.
OFFSET_LR_CAP <- log(10)

# BCR order for tables (north to south): Hudson Plains, Boreal Shield W / E,
# Temperate Mixed, Lower Great Lakes / St. Lawrence Plain.
BCR_ORDER <- c("74", "76", "77", "12", "13")

# Exclude species whose 08 province total fell back to the untrimmed (whole-
# province) aggregate because no BCR cleared min_sq_det.
DROP_UNTRIMMED_PROVINCE <- TRUE

CHANGE_LEVELS  <- c("Very likely decrease", "Likely decrease", "Little change",
                    "Uncertain", "Likely increase", "Very likely increase")
# Orange-purple diverging ramp (safe for red-green colour deficiency).
CHANGE_COLOURS <- setNames(c("#B35806", "#F1A340", "grey45",
                             "grey80", "#998EC3", "#542788"), CHANGE_LEVELS)

PLOT_AGG_FACT <- 2L      # display-only aggregation of 1-km rasters
FIG_DPI       <- 300

stopifnot(P_VERY_LIKELY > P_LIKELY, P_LIKELY > 0.5,
          CLASS_PRECEDENCE %in% c("direction", "magnitude"),
          PIXEL_P_PRIMARY > 0.5, PIXEL_P_APPROX > 0.5, HEX_FOLD > 1, HEX_P > 0.5,
          all(POOL_POBS_THRESHOLDS > 0 & POOL_POBS_THRESHOLDS < 1))

pool_tag <- function(t) sprintf("pobs%03d", round(100 * t))

# Figure 4-5 variants: P threshold and whether 08's support mask applies.
VARIANTS <- list(
  primary   = list(p = PIXEL_P_PRIMARY, support = FALSE,
                   label = sprintf("P \u2265 %.2f (all BCRs)", PIXEL_P_PRIMARY)),
  supported = list(p = PIXEL_P_PRIMARY, support = TRUE,
                   label = sprintf("P \u2265 %.2f, 08-supported BCRs only", PIXEL_P_PRIMARY)),
  approx    = list(p = PIXEL_P_APPROX, support = FALSE,
                   label = sprintf("P \u2265 %.2f (approximate)", PIXEL_P_APPROX))
)
METRICS <- c("pool", "increase", "decrease")
ACC_COLS <- as.vector(outer(names(VARIANTS), METRICS, paste, sep = "_"))

# ============================================================
# Helpers
# ============================================================

# Same slug rule as sp_filename() in inla_model_utils.R.
sp_filename <- function(x) gsub("^_|_$", "", gsub("[^a-z0-9]+", "_", tolower(x)))

# Median variance (1 / precision) of an iid term; 0 if the term was not fitted.
hyper_var <- function(hyper, term) {
  if (is.null(hyper)) return(0)
  i <- match(paste0("Precision for ", term), rownames(hyper))
  if (is.na(i)) return(0)
  v <- 1 / hyper[i, "0.5quant"]
  if (is.finite(v) && v > 0) v else 0
}

# Median nbinomial size; Inf for a Poisson fit.
nb_size <- function(hyper) {
  if (is.null(hyper)) return(Inf)
  i <- grep("size for the nbinomial", rownames(hyper), fixed = TRUE)
  if (length(i) == 0) return(Inf)
  s <- hyper[i[1], "0.5quant"]
  if (is.finite(s) && s > 0) s else Inf
}

# P(>= 1 detection in n_surveys) from median expected counts (vector form of
# marginal_pobs() in model_product_utils.R). NA stays NA; lambda <= 0 -> 0.
pobs_from_lambda <- function(lambda, sigma2, size, n_surveys, n_nodes) {
  out <- ifelse(is.finite(lambda), 0, NA_real_)
  pos <- which(is.finite(lambda) & lambda > 0)
  if (length(pos) == 0) return(out)
  l  <- lambda[pos]
  P0 <- if (is.finite(size)) function(m) (size / (size + m))^size else function(m) exp(-m)
  if (sigma2 > 0) {
    z  <- seq(from = -6, to = 6, length.out = n_nodes)
    w  <- stats::dnorm(z) / sum(stats::dnorm(z))
    p1 <- numeric(length(l))
    for (j in seq_along(z)) p1 <- p1 + w[j] * (1 - P0(l * exp(sqrt(sigma2) * z[j])))
  } else {
    p1 <- 1 - P0(l)
  }
  out[pos] <- 1 - (1 - p1)^n_surveys
  out
}

# P(change > 0) from a posterior median and interval (split normal: the half of
# the interval facing zero sets the scale; exact when a bound sits at zero).
p_increase_split <- function(q50, lower, upper, z = PIXEL_CI_Z) {
  s_lo <- pmax((q50 - lower) / z, 1e-12)
  s_hi <- pmax((upper - q50) / z, 1e-12)
  ifelse(q50 >= 0, stats::pnorm(q50 / s_lo), 1 - stats::pnorm(-q50 / s_hi))
}

# Unique raster cells and each pixel's group (index into the unique cells).
cell_groups <- function(cells) {
  u <- unique(cells[!is.na(cells)])
  list(u = u, g = match(cells, u))
}

# Single-layer raster -> data.frame(x, y, value), masked and aggregated (mean).
raster_df <- function(r, fact = PLOT_AGG_FACT) {
  r <- terra::mask(r, keep_rast)
  if (fact > 1) r <- terra::aggregate(r, fact = fact, fun = "mean", na.rm = TRUE)
  d <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(d)[3] <- "value"
  d
}

# Share = count / pool, aggregated by summing counts and pools over each display
# block; shown where the mean pool per cell is >= MIN_POOL_FOR_PROP.
share_df <- function(num, den, fact = PLOT_AGG_FACT) {
  num <- terra::mask(num, keep_rast)
  den <- terra::mask(den, keep_rast)
  n   <- !is.na(den)
  if (fact > 1) {
    n   <- terra::aggregate(n,   fact = fact, fun = "sum", na.rm = TRUE)
    num <- terra::aggregate(num, fact = fact, fun = "sum", na.rm = TRUE)
    den <- terra::aggregate(den, fact = fact, fun = "sum", na.rm = TRUE)
  }
  r <- terra::ifel(den / n >= MIN_POOL_FOR_PROP, num / den, NA)
  d <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(d)[3] <- "value"
  d
}

map_panel <- function(df, title, fill_scale) {
  ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    fill_scale +
    coord_equal(expand = FALSE) +
    labs(title = title) +
    theme_void() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

fill_diff <- function(dlim, name = "Atlas 3 \u2212\nAtlas 2", labels = waiver()) {
  scale_fill_gradient2(name = name, low = "#B35806", mid = "grey97", high = "#542788",
                       midpoint = 0, limits = c(-dlim, dlim), labels = labels, oob = scales::squish)
}

seq_fill <- function(high, name, limits = NULL, labels = waiver()) {
  scale_fill_gradient(name = name, low = "grey97", high = high, limits = limits, labels = labels,
                      oob = scales::squish)
}

pct_label <- function(x) {
  pct <- round(100 * (exp(x) - 1))
  ifelse(pct > 0, paste0("+", pct, "%"), paste0(pct, "%"))
}

classify_change <- function(p_increase, p_decrease, qlow, qhigh) {
  in_band <- log1p(qlow) >= -log(STABLE_FOLD) & log1p(qhigh) <= log(STABLE_FOLD)
  cls <- if (CLASS_PRECEDENCE == "magnitude") {
    dplyr::case_when(
      in_band                     ~ "Little change",
      p_increase >= P_VERY_LIKELY ~ "Very likely increase",
      p_decrease >= P_VERY_LIKELY ~ "Very likely decrease",
      p_increase >= P_LIKELY      ~ "Likely increase",
      p_decrease >= P_LIKELY      ~ "Likely decrease",
      TRUE                        ~ "Uncertain")
  } else {
    dplyr::case_when(
      p_increase >= P_VERY_LIKELY ~ "Very likely increase",
      p_decrease >= P_VERY_LIKELY ~ "Very likely decrease",
      in_band                     ~ "Little change",
      p_increase >= P_LIKELY      ~ "Likely increase",
      p_decrease >= P_LIKELY      ~ "Likely decrease",
      TRUE                        ~ "Uncertain")
  }
  factor(cls, levels = CHANGE_LEVELS)
}

# ============================================================
# Locations
# ============================================================

inputs  <- c(PRED_DIR = PRED_DIR, RAST_DIR = RAST_DIR,
             CHANGE_TABLE = CHANGE_TABLE, DATA_READY_RDS = DATA_READY_RDS)
missing <- inputs[!file.exists(inputs)]
if (length(missing) > 0) {
  stop("Not found (check CONFIG):\n  ",
       paste(names(missing), missing, sep = " = ", collapse = "\n  "), call. = FALSE)
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

model_out      <- dirname(PRED_DIR)
model_name     <- sub("^predictions_", "", basename(PRED_DIR))
summaries_path <- file.path(model_out, paste0("summaries_", model_name), "model_summaries.rds")

message("Predictions:  ", PRED_DIR, "\nRasters:      ", RAST_DIR,
        "\nChange table: ", CHANGE_TABLE, "\nOutput:       ", OUT_DIR)

# ============================================================
# Species inventory, hyperparameters, 08 change table
# ============================================================

species <- tibble(sp_file = sub("_a2\\.tif$", "",
                                list.files(RAST_DIR, pattern = "_a2\\.tif$"))) |>
  mutate(path_a2   = file.path(RAST_DIR, paste0(sp_file, "_a2.tif")),
         path_a3   = file.path(RAST_DIR, paste0(sp_file, "_a3.tif")),
         path_pred = file.path(PRED_DIR, paste0(sp_file, "_1km.rds")))

incomplete <- species |> filter(!file.exists(path_a3) | !file.exists(path_pred))
if (nrow(incomplete) > 0) {
  message("Skipping ", nrow(incomplete), " species without both rasters and a ",
          "prediction file: ", paste(incomplete$sp_file, collapse = ", "))
}
species <- species |> filter(file.exists(path_a3), file.exists(path_pred))
if (nrow(species) == 0) stop("No species with A2 + A3 rasters and predictions.", call. = FALSE)
message("Species with rasters + predictions: ", nrow(species))

model_summaries <- if (file.exists(summaries_path)) readRDS(summaries_path) else list()
ms_lookup <- setNames(model_summaries, sp_filename(names(model_summaries)))

species_hyper <- function(sp_file, preds) {
  ms <- ms_lookup[[sp_file]]
  if (!is.null(ms$summary_hyperpar)) {
    return(list(sp_name = ms$sp_name, hyper = ms$summary_hyperpar, vc = ms$variance_correction))
  }
  list(sp_name = preds$sp_name, hyper = preds$summary_hyperpar, vc = preds$variance_correction)
}

chg_all <- readRDS(CHANGE_TABLE)
needed  <- c("sp_english", "sp_file", "region_type", "Region_Name", "Region_Number",
             "prop_change_median", "prop_change_qlow", "prop_change_qhigh",
             "p_increase", "p_decrease", "sufficient_data", "province_trimmed", "ci_level")
missing_cols <- setdiff(needed, names(chg_all))
if (length(missing_cols) > 0) {
  stop("Change table lacks: ", paste(missing_cols, collapse = ", "), ". Re-run 08.", call. = FALSE)
}
# 08's support: BCRs where each species was sufficiently observed.
supported_bcrs <- chg_all |>
  filter(region_type == "BCR", sufficient_data %in% TRUE) |>
  group_by(sp_file) |>
  summarise(bcrs = list(as.character(Region_Number)), .groups = "drop")
supported_bcrs <- setNames(supported_bcrs$bcrs, supported_bcrs$sp_file)

# ============================================================
# Spatial lookups: pixel -> cell, cell -> BCR, cell -> hex
# ============================================================

template <- terra::rast(species$path_a2[1])[["mu_q50"]]
n_cell   <- terra::ncell(template)
cell_km2 <- terra::values(terra::cellSize(template, unit = "km"), mat = FALSE)
atlases  <- c(A2 = "path_a2", A3 = "path_a3")
tmpl_crs <- sf::st_crs(terra::crs(template))

ready    <- readRDS(DATA_READY_RDS)
grid2    <- ready$grid_OBBA2
bcr_sf   <- ready$bcr_sf
hex_grid <- ready$hex_grid_25km
rm(ready)
for (nm in c("grid2", "bcr_sf", "hex_grid")) {
  if (is.null(get(nm))) stop(nm, " missing from ", DATA_READY_RDS, call. = FALSE)
}

grid_geom <- sf::st_geometry(grid2)
if (!all(sf::st_geometry_type(grid_geom) == "POINT")) grid_geom <- suppressWarnings(sf::st_centroid(grid_geom))
grid_xy       <- sf::st_coordinates(sf::st_transform(grid_geom, crs = tmpl_crs))
grid_pixel_id <- grid2$pixel_id
grid_cell     <- terra::cellFromXY(template, grid_xy[, 1:2])
grid_groups   <- cell_groups(grid_cell)
n_dup <- length(grid_cell) - length(grid_groups$u) - sum(is.na(grid_cell))
message(sprintf("Prediction pixels: %d; raster cells holding >1 pixel: %d (their P(increase) is averaged).",
                length(grid_cell), n_dup))
rm(grid2, grid_geom, grid_xy)

bcr_v     <- terra::vect(sf::st_transform(sf::st_make_valid(bcr_sf), tmpl_crs))
bcr_v$row <- seq_len(nrow(bcr_v))
cell_bcr  <- as.character(bcr_sf$BCR)[terra::values(terra::rasterize(bcr_v, template, field = "row"),
                                                    mat = FALSE)]

hex_grid  <- sf::st_transform(hex_grid, tmpl_crs)
cell_hex  <- terra::values(terra::rasterize(terra::vect(hex_grid), template, field = "hex_id"), mat = FALSE)
hex_cells <- which(!is.na(cell_hex))
hex_ids   <- sort(unique(cell_hex[hex_cells]))
hex_row   <- match(cell_hex[hex_cells], hex_ids)
n_hex     <- length(hex_ids)

# Majority BCR per hex (by cell count) and the share of its cells in that BCR.
hex_bcr <- tibble(row = hex_row, bcr = cell_bcr[hex_cells]) |>
  filter(!is.na(bcr)) |>
  count(row, bcr) |>
  group_by(row) |>
  mutate(bcr_prop = n / sum(n)) |>
  slice_max(n, n = 1, with_ties = FALSE) |>
  ungroup()
hex_bcr <- tibble(hex_id   = hex_ids,
                  bcr      = hex_bcr$bcr[match(seq_len(n_hex), hex_bcr$row)],
                  bcr_prop = hex_bcr$bcr_prop[match(seq_len(n_hex), hex_bcr$row)])
bcr_names <- chg_all |>
  filter(region_type == "BCR") |>
  distinct(bcr = as.character(Region_Number), bcr_name = Region_Name)

# ============================================================
# Per-species loop
# ============================================================

soft     <- lapply(atlases, function(x) numeric(n_cell))
modelled <- lapply(atlases, function(x) logical(n_cell))
eligible <- logical(n_cell)                      # cells with any species' change summary
acc      <- lapply(setNames(POOL_POBS_THRESHOLDS, pool_tag(POOL_POBS_THRESHOLDS)),
                   function(t) matrix(0L, nrow = n_cell, ncol = length(ACC_COLS),
                                      dimnames = list(NULL, ACC_COLS)))
GL       <- matrix(0, nrow = n_cell, ncol = 3, dimnames = list(NULL, c("gains", "losses", "sum_pobs")))
hex_acc  <- matrix(0L, nrow = n_hex, ncol = 5,
                   dimnames = list(NULL, c("pool", "increase", "decrease", "increase_dir", "decrease_dir")))
# Figure 8: per hex and draw, the sum over pool species of capped log(A3 / A2),
# and the number of species summed (sized on the first species with draws).
hex_lr_sum <- NULL
hex_lr_n   <- integer(n_hex)
sp_dep     <- vector("list", nrow(species))   # per species: departure quantiles by hex
n_draw_mismatch <- 0L
sp_records <- vector("list", nrow(species))
sp_hex     <- vector("list", nrow(species))
n_no_change_summary <- 0L
n_no_hex_draws      <- 0L

for (i in seq_len(nrow(species))) {
  
  sp_file <- species$sp_file[i]
  rasts   <- lapply(atlases, function(col) terra::rast(species[[col]][i]))
  
  geom_ok <- all(vapply(rasts, function(r) {
    "mu_q50" %in% names(r) && terra::compareGeom(r, template, stopOnError = FALSE)
  }, logical(1)))
  if (!geom_ok) {
    warning(sp_file, ": missing mu_q50 layer or raster geometry differs; skipped.", call. = FALSE)
    next
  }
  
  preds <- readRDS(species$path_pred[i])
  hp    <- species_hyper(sp_file = sp_file, preds = preds)
  
  # Integrate only iid terms NOT already baked into q50 by 07's optional
  # variance correction (same rule as 09).
  applied <- if (is.null(hp$vc$present)) character(0) else names(which(hp$vc$present))
  sigma2  <- sum(vapply(setdiff(POBS_MARGINALIZE_TERMS, applied),
                        function(tm) hyper_var(hyper = hp$hyper, term = tm), numeric(1)))
  size    <- nb_size(hyper = hp$hyper)
  
  rec <- list(sp_file = sp_file, sp_english = if (is.null(hp$sp_name)) sp_file else hp$sp_name,
              sigma2_omit = sigma2, nb_size = size)
  
  # ---- PObs per cell, both atlases (Figure 1) ----
  p_at <- list()
  for (a in names(atlases)) {
    p  <- pobs_from_lambda(lambda    = terra::values(rasts[[a]][["mu_q50"]], mat = FALSE),
                           sigma2    = sigma2,
                           size      = size,
                           n_surveys = POBS_N_SURVEYS,
                           n_nodes   = POBS_GH_NODES)
    ok <- is.finite(p)
    soft[[a]][ok] <- soft[[a]][ok] + p[ok]
    modelled[[a]] <- modelled[[a]] | ok
    p_at[[a]]     <- p
    
    # PObs-weighted area: threshold-free analogue of occupied area.
    rec[[paste0("pobs_area_km2_", a)]] <- sum(cell_km2[ok] * p[ok])
    rec[[paste0("max_pobs_", a)]]      <- if (any(ok)) max(p[ok]) else NA_real_
  }
  p_max <- pmax(p_at$A2, p_at$A3, na.rm = TRUE)    # PObs in the higher atlas
  
  # ---- Gains and losses (Figure 6) ----
  both <- which(is.finite(p_at$A2) & is.finite(p_at$A3))
  dp   <- p_at$A3[both] - p_at$A2[both]
  GL[both, "gains"]    <- GL[both, "gains"]    + pmax(dp, 0)
  GL[both, "losses"]   <- GL[both, "losses"]   + pmax(-dp, 0)
  GL[both, "sum_pobs"] <- GL[both, "sum_pobs"] + p_at$A2[both] + p_at$A3[both]
  
  # ---- Pixel change classification (Figures 4-5) ----
  pchg <- preds$abs_change_Corrected_for_Water
  if (is.null(pchg) || !all(c("pixel_id", "abs_change_q50", "abs_change_lower", "abs_change_upper") %in% names(pchg))) {
    n_no_change_summary <- n_no_change_summary + 1L
  } else {
    ug <- if (identical(pchg$pixel_id, grid_pixel_id)) grid_groups else
      cell_groups(grid_cell[match(pchg$pixel_id, grid_pixel_id)])
    ok   <- !is.na(ug$g)
    p_px <- p_increase_split(q50   = pchg$abs_change_q50[ok],
                             lower = pchg$abs_change_lower[ok],
                             upper = pchg$abs_change_upper[ok])
    # One P(increase) per raster cell: mean over the pixels in it.
    keep_g <- is.finite(p_px)
    s      <- rowsum(p_px[keep_g], ug$g[ok][keep_g])
    uc     <- ug$u[as.integer(rownames(s))]
    p_inc  <- s[, 1] / tabulate(ug$g[ok][keep_g], nbins = length(ug$u))[as.integer(rownames(s))]
    
    eligible[uc] <- TRUE
    pm  <- p_max[uc]
    sup <- cell_bcr[uc] %in% supported_bcrs[[sp_file]]
    for (k in seq_along(POOL_POBS_THRESHOLDS)) {
      in_pool <- (pm > POOL_POBS_THRESHOLDS[k]) %in% TRUE
      for (v in names(VARIANTS)) {
        member <- if (VARIANTS[[v]]$support) in_pool & sup else in_pool
        thr    <- VARIANTS[[v]]$p - 1e-9
        acc[[k]][uc, paste0(v, "_pool")]     <- acc[[k]][uc, paste0(v, "_pool")]     + member
        acc[[k]][uc, paste0(v, "_increase")] <- acc[[k]][uc, paste0(v, "_increase")] + (member & p_inc >= thr)
        acc[[k]][uc, paste0(v, "_decrease")] <- acc[[k]][uc, paste0(v, "_decrease")] + (member & (1 - p_inc) >= thr)
      }
    }
  }
  
  # ---- Hexagon magnitude assessment (Figure 7) ----
  hd <- preds$hex_draws_Corrected_for_Water
  if (is.null(hd)) {
    n_no_hex_draws <- n_no_hex_draws + 1L
  } else {
    pv  <- p_max[hex_cells]
    okv <- is.finite(pv)
    hm  <- numeric(n_hex)
    s   <- rowsum(pv[okv], hex_row[okv])
    hm[as.integer(rownames(s))] <- s[, 1]
    hm  <- hm / pmax(tabulate(hex_row[okv], nbins = n_hex), 1)   # hex-mean PObs (higher atlas)
    
    j  <- match(hex_ids, hd$hex_id)
    m2 <- do.call(rbind, hd$mu_OBBA2)[j, , drop = FALSE]
    m3 <- do.call(rbind, hd$mu_OBBA3)[j, , drop = FALSE]
    r  <- m3 / m2                                                 # 0/0 -> NaN (dropped); x/0 -> Inf
    p_up  <- rowMeans(r >= HEX_FOLD, na.rm = TRUE)
    p_dn  <- rowMeans(r <= 1 / HEX_FOLD, na.rm = TRUE)
    p_dir <- rowMeans(r > 1, na.rm = TRUE)
    
    in_pool <- (hm > HEX_POOL_POBS) %in% TRUE
    hex_acc[, "pool"]         <- hex_acc[, "pool"]         + in_pool
    hex_acc[, "increase"]     <- hex_acc[, "increase"]     + (in_pool & (p_up >= HEX_P) %in% TRUE)
    hex_acc[, "decrease"]     <- hex_acc[, "decrease"]     + (in_pool & (p_dn >= HEX_P) %in% TRUE)
    hex_acc[, "increase_dir"] <- hex_acc[, "increase_dir"] + (in_pool & (p_dir >= HEX_P) %in% TRUE)
    hex_acc[, "decrease_dir"] <- hex_acc[, "decrease_dir"] + (in_pool & ((1 - p_dir) >= HEX_P) %in% TRUE)
    
    # ---- Community offset and departures (Figure 8) ----
    # Capped log ratio per draw; capping turns x/0 (+/-Inf) into the cap, while
    # 0/0 (NaN) stays non-finite and drops the species from that hex.
    lr     <- pmin(pmax(log(r), -OFFSET_LR_CAP), OFFSET_LR_CAP)
    use_lr <- in_pool & (rowSums(!is.finite(lr)) == 0)
    if (is.null(hex_lr_sum)) hex_lr_sum <- matrix(0, nrow = n_hex, ncol = ncol(lr))
    if (ncol(lr) != ncol(hex_lr_sum)) {
      n_draw_mismatch <- n_draw_mismatch + 1L      # cannot pair draws; left out of Figure 8
    } else if (any(use_lr)) {
      hex_lr_sum[use_lr, ] <- hex_lr_sum[use_lr, ] + lr[use_lr, , drop = FALSE]
      hex_lr_n[use_lr]     <- hex_lr_n[use_lr] + 1L
      # Quantiles used for departures once the offset is known (after the loop).
      lr_q <- t(apply(X = lr[use_lr, , drop = FALSE], MARGIN = 1, FUN = stats::quantile,
                      probs = c(1 - HEX_P, HEX_P), names = FALSE))
      sp_dep[[i]] <- tibble(sp_file = sp_file, hex_row = which(use_lr),
                            lr_q_low = lr_q[, 1], lr_q_high = lr_q[, 2])
    }
    
    if (any(in_pool)) {
      sp_hex[[i]] <- tibble(sp_file = sp_file, hex_id = hex_ids[in_pool], hex_mean_pobs = hm[in_pool],
                            ratio_q50 = apply(r[in_pool, , drop = FALSE], 1, stats::median, na.rm = TRUE),
                            p_increase_fold = p_up[in_pool], p_decrease_fold = p_dn[in_pool],
                            p_increase = p_dir[in_pool])
    }
  }
  rm(preds)
  
  sp_records[[i]] <- tibble::as_tibble(rec)
  if (i %% 25 == 0) message("  processed ", i, "/", nrow(species), " species")
}

if (n_no_change_summary > 0) message(n_no_change_summary, " species had no pixel change summary (absent from Figures 4-5).")
if (n_no_hex_draws > 0)      message(n_no_hex_draws, " species had no hex draws (absent from Figure 7).")
if (n_draw_mismatch > 0)     message(n_draw_mismatch, " species had a different number of hex draws (absent from Figure 8).")

sp_tab <- bind_rows(sp_records)
n_sp_richness <- nrow(sp_tab)
n_sp_pixel    <- n_sp_richness - n_no_change_summary
n_sp_hex      <- n_sp_richness - n_no_hex_draws
message("Richness built from ", n_sp_richness, " species.")

# ============================================================
# Rasters
# ============================================================

rich_vals <- lapply(names(atlases), function(a) ifelse(modelled[[a]], soft[[a]], NA_real_))
names(rich_vals) <- names(atlases)
richness <- lapply(rich_vals, function(v) {
  r <- terra::setValues(template, v)
  names(r) <- "expected_n_spp"
  r
})
richness_change <- richness$A3 - richness$A2

# Display mask: one mask for all maps. Saved rasters stay unmasked.
keep_cells <- pmax(rich_vals$A2, rich_vals$A3, na.rm = TRUE) >= MIN_MAP_RICHNESS
keep_cells <- keep_cells %in% TRUE
keep_rast  <- terra::setValues(template, ifelse(keep_cells, 1, NA))
any_modelled <- modelled$A2 | modelled$A3
message(sprintf("Maps show %.1f%% of modelled pixels (>= %d expected species in either atlas).",
                100 * sum(keep_cells) / sum(any_modelled), MIN_MAP_RICHNESS))

# Species pool / change rasters: counts and shares per variant (pool-eligible cells only).
pool_rasts <- lapply(acc, function(m) {
  layers <- list()
  for (v in names(VARIANTS)) {
    pool <- m[, paste0(v, "_pool")]; inc <- m[, paste0(v, "_increase")]; dec <- m[, paste0(v, "_decrease")]
    shr  <- function(n) ifelse(pool >= MIN_POOL_FOR_PROP, n / pmax(pool, 1), NA_real_)
    layers[[paste0(v, "_pool")]]           <- pool
    layers[[paste0(v, "_increase")]]       <- inc
    layers[[paste0(v, "_decrease")]]       <- dec
    layers[[paste0(v, "_share_increase")]] <- shr(inc)
    layers[[paste0(v, "_share_decrease")]] <- shr(dec)
    layers[[paste0(v, "_share_changed")]]  <- shr(inc + dec)
  }
  vals <- do.call(cbind, layers)
  vals[!eligible, ] <- NA
  r <- terra::rast(template, nlyrs = ncol(vals))
  terra::values(r) <- vals
  names(r) <- names(layers)
  r
})

GL_vals <- cbind(gains = GL[, "gains"], losses = GL[, "losses"],
                 total = GL[, "gains"] + GL[, "losses"],
                 relative = ifelse(GL[, "sum_pobs"] > 0, (GL[, "gains"] + GL[, "losses"]) / GL[, "sum_pobs"], NA))
GL_vals[!any_modelled, ] <- NA
gl_rast <- terra::rast(template, nlyrs = 4)
terra::values(gl_rast) <- GL_vals
names(gl_rast) <- colnames(GL_vals)
sum_pobs_rast <- terra::setValues(template, ifelse(any_modelled, GL[, "sum_pobs"], NA))
rm(soft, rich_vals, acc, GL, GL_vals)

richness_tags <- c(pobs_n_surveys     = as.character(POBS_N_SURVEYS),
                   pobs_marginalized  = paste(POBS_MARGINALIZE_TERMS, collapse = ","),
                   n_species          = as.character(n_sp_richness),
                   model_name         = model_name)
for (a in names(richness)) {
  terra::metags(richness[[a]]) <- c(richness_tags, atlas = a)
  terra::writeRaster(richness[[a]], filename = file.path(OUT_DIR, paste0("richness_", a, ".tif")),
                     overwrite = TRUE)
}
terra::metags(richness_change) <- c(richness_tags, atlas = "A3 minus A2")
terra::writeRaster(richness_change, filename = file.path(OUT_DIR, "richness_change_A3_minus_A2.tif"),
                   overwrite = TRUE)
for (k in seq_along(pool_rasts)) {
  terra::metags(pool_rasts[[k]]) <- c(richness_tags, pool_pobs_threshold = as.character(POOL_POBS_THRESHOLDS[k]),
                                      p_primary = as.character(PIXEL_P_PRIMARY),
                                      p_approx = as.character(PIXEL_P_APPROX),
                                      min_pool_for_prop = as.character(MIN_POOL_FOR_PROP))
  terra::writeRaster(pool_rasts[[k]], overwrite = TRUE,
                     filename = file.path(OUT_DIR, sprintf("species_pool_change_%s.tif", names(pool_rasts)[k])))
}
terra::metags(gl_rast) <- richness_tags
terra::writeRaster(gl_rast, filename = file.path(OUT_DIR, "gains_losses.tif"), overwrite = TRUE)

# ============================================================
# Figure 1: expected richness
# ============================================================

richness_caption <- sprintf(paste0("%d modelled species; sum over species of P(\u2265 1 detection in %d ",
                                   "five-minute point counts). Areas with < %d expected species in ",
                                   "both atlases (mostly open water) not shown."),
                            n_sp_richness, POBS_N_SURVEYS, MIN_MAP_RICHNESS)

pix <- lapply(list(A2 = richness$A2, A3 = richness$A3, chg = richness_change), raster_df)
fill_rich <- scale_fill_viridis_c(name = "Expected\nspecies",
                                  limits = range(c(pix$A2$value, pix$A3$value)))
dlim_pix  <- max(1, stats::quantile(abs(pix$chg$value), probs = 0.99))

fig1 <- (map_panel(df = pix$A2, title = "Atlas 2", fill_scale = fill_rich) |
           map_panel(df = pix$A3, title = "Atlas 3", fill_scale = fill_rich) |
           map_panel(df = pix$chg, title = "Change", fill_scale = fill_diff(dlim_pix))) +
  plot_annotation(title    = "Expected number of species detected",
                  subtitle = "1-km pixels, from posterior medians",
                  caption  = richness_caption,
                  theme    = theme(plot.title = element_text(size = 16, face = "bold")))
ggsave(filename = file.path(OUT_DIR, "fig1_expected_richness_pixels.png"),
       plot = fig1, width = 18, height = 7, dpi = FIG_DPI, bg = "white")

# ============================================================
# Figures 2-3: province and BCR change classes
# ============================================================

ci_level <- unique(chg_all$ci_level)[1]
ci_pct   <- round(100 * ci_level)

chg <- chg_all |>
  filter(region_type %in% c("province", "BCR")) |>
  mutate(change_class = classify_change(p_increase = p_increase,
                                        p_decrease = p_decrease,
                                        qlow       = prop_change_qlow,
                                        qhigh      = prop_change_qhigh))

prov <- chg |> filter(region_type == "province")
if (DROP_UNTRIMMED_PROVINCE) {
  untrimmed <- prov |> filter(!province_trimmed %in% TRUE)
  if (nrow(untrimmed) > 0) {
    message("Dropping ", nrow(untrimmed), " species with an untrimmed province total ",
            "(no BCR cleared min_sq_det in 08): ",
            paste(untrimmed$sp_english, collapse = ", "))
  }
  prov <- prov |> filter(province_trimmed %in% TRUE)
}
bcr  <- chg |> filter(region_type == "BCR", sufficient_data)
n_sp <- nrow(prov)

class_def <- sprintf(paste0("Very likely: P(direction) \u2265 %.2f;  likely: \u2265 %.2f;  ",
                            "little change: %d%% CI within %+.0f%% to %+.0f%%;  uncertain: otherwise."),
                     P_VERY_LIKELY, P_LIKELY, ci_pct,
                     100 * (1 / STABLE_FOLD - 1), 100 * (STABLE_FOLD - 1))
precedence_txt <- if (CLASS_PRECEDENCE == "magnitude") {
  "'Little change' takes precedence over direction."
} else {
  "Direction takes precedence: a precisely estimated small change can be 'very likely'."
}
footprint_txt <- paste0("Province-wide change for each species sums only the BCRs where it was sufficiently ",
                        "observed (script 08), so the province-wide footprint differs among species.")

# ---- Figure 2: ranked province-wide change, every species ----
cat_df <- prov |>
  mutate(lr_med  = log1p(prop_change_median),
         lr_low  = log1p(prop_change_qlow),
         lr_high = log1p(prop_change_qhigh)) |>
  arrange(lr_med) |>
  mutate(sp_english = factor(sp_english, levels = sp_english))

pct_breaks <- c(-80, -50, -33, 0, 50, 100, 200, 400)

fig2 <- ggplot(cat_df, aes(y = sp_english, x = lr_med, xmin = lr_low, xmax = lr_high,
                           colour = change_class)) +
  annotate(geom = "rect", xmin = -log(STABLE_FOLD), xmax = log(STABLE_FOLD),
           ymin = -Inf, ymax = Inf, fill = "grey60", alpha = 0.15) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_pointrange(size = 0.12, linewidth = 0.35) +
  scale_colour_manual(values = CHANGE_COLOURS, drop = FALSE, name = NULL) +
  scale_x_continuous(breaks = log1p(pct_breaks / 100), labels = pct_label) +
  coord_cartesian(xlim = log(c(0.15, 6))) +
  labs(x = sprintf("Province-wide population change, Atlas 2 to Atlas 3 (median, %d%% CI)", ci_pct),
       y = NULL, title = "Province-wide change by species",
       subtitle = paste0(gsub(";  little", ";\nlittle", class_def), "\n", precedence_txt,
                         "\nShaded band = 'little change' zone. Intervals beyond the axis are truncated."),
       caption = gsub("observed \\(", "observed\n(", footprint_txt)) +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_text(size = 5.5), legend.position = "top",
        panel.grid.major.y = element_blank(), plot.title = element_text(face = "bold")) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 0.5)))
ggsave(filename = file.path(OUT_DIR, "fig2_species_change_ranked.png"),
       plot = fig2, width = 9, height = max(6, 1.5 + 0.1 * n_sp),
       dpi = FIG_DPI, bg = "white", limitsize = FALSE)

# ---- Figure 3: class counts, province + sufficiently observed BCRs ----
region_levels <- rev(unique(c(prov$Region_Name, bcr$Region_Name)))
bar_df <- bind_rows(prov, bcr) |>
  count(Region_Name, change_class, .drop = FALSE) |>
  mutate(Region_Name = factor(Region_Name, levels = region_levels))
bar_tot <- bar_df |> group_by(Region_Name) |> summarise(n = sum(n), .groups = "drop")

fig3 <- ggplot(bar_df, aes(x = n, y = Region_Name, fill = change_class)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.7,
           colour = "grey30", linewidth = 0.2) +
  geom_text(data = bar_tot, aes(x = n, y = Region_Name, label = paste0("n = ", n)),
            inherit.aes = FALSE, hjust = -0.15, size = 3.5) +
  scale_fill_manual(values = CHANGE_COLOURS, drop = FALSE, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Number of species", y = NULL,
       title = "Species population change between atlases",
       subtitle = paste0(class_def, "\n", precedence_txt,
                         "\nBCRs show only species detected in \u2265 the 08 square threshold there."),
       caption = footprint_txt) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold")) +
  guides(fill = guide_legend(nrow = 1))
ggsave(filename = file.path(OUT_DIR, "fig3_change_classes.png"),
       plot = fig3, width = 12, height = 2 + 0.6 * length(region_levels),
       dpi = FIG_DPI, bg = "white")

message(sprintf("\nProvince-wide change classes (%d species, precedence = %s):", n_sp, CLASS_PRECEDENCE))
print(table(prov$change_class))

# ============================================================
# Figures 4-5: species pool and well-supported pixel change
# ============================================================

pixel_caption <- function(extra = "") {
  sprintf(paste0("Species pool: P(\u2265 1 detection in %d counts) above the threshold in either atlas, with a ",
                 "pixel change summary; each species counts as one. %d modelled species. Shares shown where the ",
                 "pool has \u2265 %d species%s.\nP \u2265 %.2f = the 90%% interval of the pixel change excludes ",
                 "zero (exact from 07's saved quantiles); P \u2265 %.2f is approximated from them. Counts show ",
                 "evidence of change, which depends on survey coverage (see script 16).\nSpecies are classified ",
                 "one at a time; biases shared by all species are not represented. Areas with < %d expected species ",
                 "in both atlases not shown."),
          POBS_N_SURVEYS, n_sp_pixel, MIN_POOL_FOR_PROP, extra, PIXEL_P_PRIMARY, PIXEL_P_APPROX, MIN_MAP_RICHNESS)
}

pool_summary <- list()
for (k in seq_along(POOL_POBS_THRESHOLDS)) {
  t   <- POOL_POBS_THRESHOLDS[k]
  pr  <- pool_rasts[[k]]
  tag <- names(pool_rasts)[k]
  L   <- function(v, m) pr[[paste0(v, "_", m)]]
  
  # ---- Figure 4: primary variant ----
  cnt <- lapply(c(pool = "pool", increase = "increase", decrease = "decrease"),
                function(m) raster_df(L("primary", m)))
  shr <- list(changed  = share_df(L("primary", "increase") + L("primary", "decrease"), L("primary", "pool")),
              increase = share_df(L("primary", "increase"), L("primary", "pool")),
              decrease = share_df(L("primary", "decrease"), L("primary", "pool")))
  n_lim <- c(0, max(cnt$increase$value, cnt$decrease$value))
  s_lim <- c(0, max(shr$increase$value, shr$decrease$value))
  
  fig4 <- (map_panel(df = cnt$pool, title = "Species pool",
                     fill_scale = scale_fill_viridis_c(name = "Species")) |
             map_panel(df = cnt$increase, title = "Increased", fill_scale = seq_fill("#542788", "Species", n_lim)) |
             map_panel(df = cnt$decrease, title = "Decreased", fill_scale = seq_fill("#B35806", "Species", n_lim))) /
    (map_panel(df = shr$changed, title = "Share of pool changed (either direction)",
               fill_scale = seq_fill("grey15", "Share of\npool", c(0, max(shr$changed$value)),
                                     scales::label_percent())) |
       map_panel(df = shr$increase, title = "Share of pool increased",
                 fill_scale = seq_fill("#542788", "Share of\npool", s_lim, scales::label_percent())) |
       map_panel(df = shr$decrease, title = "Share of pool decreased",
                 fill_scale = seq_fill("#B35806", "Share of\npool", s_lim, scales::label_percent()))) +
    plot_annotation(
      title    = sprintf("Species with a well-supported increase or decrease in each 1-km pixel (pool: PObs > %.2f)", t),
      subtitle = sprintf(paste0("A species changed if P(increase) or P(decrease) \u2265 %.2f in the pixel ",
                                "(direction only; no magnitude criterion; see Figure 7 for \u2265 %+.0f%% changes)."),
                         PIXEL_P_PRIMARY, 100 * (HEX_FOLD - 1)),
      caption  = pixel_caption(),
      theme    = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  ggsave(filename = file.path(OUT_DIR, sprintf("fig4_species_change_pixels_%s.png", tag)),
         plot = fig4, width = 18, height = 14, dpi = FIG_DPI, bg = "white")
  
  # ---- Figure 5: sensitivity of the shares ----
  sens <- lapply(names(VARIANTS), function(v) {
    list(changed  = share_df(L(v, "increase") + L(v, "decrease"), L(v, "pool")),
         increase = share_df(L(v, "increase"), L(v, "pool")),
         decrease = share_df(L(v, "decrease"), L(v, "pool")))
  })
  names(sens) <- names(VARIANTS)
  lim_of <- function(m) c(0, max(unlist(lapply(sens, function(x) x[[m]]$value))))
  cols   <- c(changed = "grey15", increase = "#542788", decrease = "#B35806")
  heads  <- c(changed = "Share changed", increase = "Share increased", decrease = "Share decreased")
  panels <- list()
  for (v in names(VARIANTS)) for (m in names(cols)) {
    panels[[paste(v, m)]] <- map_panel(df = sens[[v]][[m]],
                                       title = paste0(heads[[m]], "\n", VARIANTS[[v]]$label),
                                       fill_scale = seq_fill(cols[[m]], "Share of\npool", lim_of(m),
                                                             scales::label_percent()))
  }
  fig5 <- wrap_plots(panels, ncol = 3) +
    plot_annotation(
      title    = sprintf("Sensitivity of the pixel change shares (pool: PObs > %.2f)", t),
      subtitle = paste0("Rows: primary (Figure 4); restricted to BCRs where 08 finds the species sufficiently ",
                        "observed (pool and counts); approximate lower probability threshold."),
      caption  = pixel_caption(", per variant"),
      theme    = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  ggsave(filename = file.path(OUT_DIR, sprintf("fig5_species_change_sensitivity_%s.png", tag)),
         plot = fig5, width = 18, height = 21, dpi = FIG_DPI, bg = "white")
  
  v <- terra::values(terra::mask(pr, keep_rast), mat = TRUE)
  pool_summary[[k]] <- bind_rows(lapply(names(VARIANTS), function(vn) {
    tibble(pool_pobs_threshold = t, variant = vn,
           median_pool           = stats::median(v[, paste0(vn, "_pool")], na.rm = TRUE),
           median_increase       = stats::median(v[, paste0(vn, "_increase")], na.rm = TRUE),
           median_decrease       = stats::median(v[, paste0(vn, "_decrease")], na.rm = TRUE),
           median_share_changed  = stats::median(v[, paste0(vn, "_share_changed")], na.rm = TRUE),
           median_share_increase = stats::median(v[, paste0(vn, "_share_increase")], na.rm = TRUE),
           median_share_decrease = stats::median(v[, paste0(vn, "_share_decrease")], na.rm = TRUE))
  }))
}
pool_summary <- bind_rows(pool_summary)
message("\nFigures 4-5, medians across mapped pixels:")
print(as.data.frame(pool_summary), digits = 3, row.names = FALSE)

# ============================================================
# Figure 6: gains, losses, total and relative change
# ============================================================

gl <- lapply(c(gains = "gains", losses = "losses", total = "total"), function(l) raster_df(gl_rast[[l]]))
# Relative change for display: summed D over summed (p2 + p3) per display block.
gl$relative <- {
  num <- terra::mask(gl_rast[["total"]], keep_rast); den <- terra::mask(sum_pobs_rast, keep_rast)
  if (PLOT_AGG_FACT > 1) {
    num <- terra::aggregate(num, fact = PLOT_AGG_FACT, fun = "sum", na.rm = TRUE)
    den <- terra::aggregate(den, fact = PLOT_AGG_FACT, fun = "sum", na.rm = TRUE)
  }
  d <- terra::as.data.frame(terra::ifel(den > 0, num / den, NA), xy = TRUE, na.rm = TRUE)
  names(d)[3] <- "value"
  d
}
gl_lim <- c(0, max(stats::quantile(gl$gains$value, 0.995), stats::quantile(gl$losses$value, 0.995)))

fig6 <- (map_panel(df = gl$gains, title = "Gains: \u03a3 increases in PObs",
                   fill_scale = seq_fill("#542788", "Species", gl_lim)) |
           map_panel(df = gl$losses, title = "Losses: \u03a3 decreases in PObs",
                     fill_scale = seq_fill("#B35806", "Species", gl_lim))) /
  (map_panel(df = gl$total, title = "Total change: gains + losses",
             fill_scale = seq_fill("grey15", "Species", c(0, stats::quantile(gl$total$value, 0.995)))) |
     map_panel(df = gl$relative, title = "Relative change: total / \u03a3 (PObs A2 + PObs A3)",
               fill_scale = seq_fill("grey15", "Bray\u2013Curtis",
                                     c(0, stats::quantile(gl$relative$value, 0.995))))) +
  plot_annotation(
    title    = "Gains, losses and total change in species detection probabilities",
    subtitle = sprintf(paste0("From posterior-median PObs (\u2265 1 detection in %d counts). Gains \u2212 losses = ",
                              "the richness change in Figure 1; gains + losses also counts changes that cancel."),
                       POBS_N_SURVEYS),
    caption  = paste0("Relative change is a Bray\u2013Curtis form on detection probabilities, not compositional ",
                      "turnover. Like Figure 1, species count by how much their PObs can move: saturated common ",
                      "species and very rare species\ncontribute little. No posterior uncertainty (plug-in medians). ",
                      "Colour scales capped at the 99.5th percentile; full values in gains_losses.tif."),
    theme    = theme(plot.title = element_text(size = 16, face = "bold"))
  )
ggsave(filename = file.path(OUT_DIR, "fig6_gains_losses_pixels.png"),
       plot = fig6, width = 14, height = 14, dpi = FIG_DPI, bg = "white")

# ============================================================
# Figure 7: hexagon-level changes of at least HEX_FOLD-fold
# ============================================================

# ---- Community offset and species departures per hex (used by Figure 8) ----
if (is.null(hex_lr_sum)) stop("No species contributed hex draws; cannot build Figure 8.", call. = FALSE)

# Offset per draw: mean capped log ratio over the hex's pool species.
off_draws <- hex_lr_sum / ifelse(hex_lr_n > 0, hex_lr_n, NA)
off_q     <- t(apply(X = off_draws, MARGIN = 1, FUN = stats::quantile,
                     probs = c(0.05, 0.5, 0.95), names = FALSE, na.rm = TRUE))

# Departures: species whose change beats (or trails) the hex's offset by at
# least HEX_FOLD-fold with probability >= HEX_P, offset fixed at its median.
dep <- bind_rows(sp_dep)
dep$offset         <- off_q[dep$hex_row, 2]
dep$departs_up     <- dep$lr_q_low  >= dep$offset + log(HEX_FOLD)
dep$departs_down   <- dep$lr_q_high <= dep$offset - log(HEX_FOLD)
dep$hex_id         <- hex_ids[dep$hex_row]
hex_offset <- tibble(hex_id        = hex_ids,
                     offset_n      = hex_lr_n,
                     offset_lr_q05 = off_q[, 1],
                     offset_lr_q50 = off_q[, 2],
                     offset_lr_q95 = off_q[, 3],
                     depart_up     = tabulate(dep$hex_row[dep$departs_up],   nbins = n_hex),
                     depart_down   = tabulate(dep$hex_row[dep$departs_down], nbins = n_hex))

hex_n <- tabulate(hex_row[any_modelled[hex_cells]], nbins = n_hex)
hex_summary <- tibble(hex_id = hex_ids, n_cells = hex_n) |>
  bind_cols(as.data.frame(hex_acc)) |>
  left_join(hex_offset, by = "hex_id") |>
  mutate(offset_pct_q50     = 100 * (exp(offset_lr_q50) - 1),
         offset_pct_q05     = 100 * (exp(offset_lr_q05) - 1),
         offset_pct_q95     = 100 * (exp(offset_lr_q95) - 1),
         offset_ci_width    = offset_lr_q95 - offset_lr_q05,      # log scale
         departed           = depart_up + depart_down,
         share_depart_up    = ifelse(offset_n >= MIN_POOL_FOR_PROP, depart_up / offset_n, NA_real_),
         share_depart_down  = ifelse(offset_n >= MIN_POOL_FOR_PROP, depart_down / offset_n, NA_real_),
         share_departed     = ifelse(offset_n >= MIN_POOL_FOR_PROP, departed / offset_n, NA_real_),
         changed        = increase + decrease,
         share_increase = ifelse(pool >= MIN_POOL_FOR_PROP, increase / pool, NA_real_),
         share_decrease = ifelse(pool >= MIN_POOL_FOR_PROP, decrease / pool, NA_real_),
         share_changed  = ifelse(pool >= MIN_POOL_FOR_PROP, changed / pool, NA_real_),
         share_changed_dir = ifelse(pool >= MIN_POOL_FOR_PROP, (increase_dir + decrease_dir) / pool, NA_real_)) |>
  filter(n_cells >= HEX_MIN_CELL_PROP * max(n_cells)) |>
  left_join(hex_bcr, by = "hex_id") |>
  left_join(bcr_names, by = "bcr") |>
  relocate(bcr, bcr_name, bcr_prop, .after = n_cells)
utils::write.csv(hex_summary, file.path(OUT_DIR, "hex_species_change_summary.csv"), row.names = FALSE)
dep_tbl <- dep |>
  transmute(sp_file, hex_id,
            departure = ifelse(departs_up, "up", ifelse(departs_down, "down", "none")))
saveRDS(bind_rows(sp_hex) |>
          left_join(hex_bcr |> select(hex_id, bcr), by = "hex_id") |>
          left_join(dep_tbl, by = c("sp_file", "hex_id")),
        file.path(OUT_DIR, "species_hex_change.rds"))

# Per-BCR summary over mapped hexes (pool >= MIN_POOL_FOR_PROP): pooled shares
# (summed counts / summed pools) and hex medians.
hex_by_bcr <- hex_summary |>
  filter(pool >= MIN_POOL_FOR_PROP, !is.na(bcr)) |>
  group_by(bcr, bcr_name) |>
  summarise(n_hex                    = n(),
            median_pool              = stats::median(pool),
            median_increase          = stats::median(increase),
            median_decrease          = stats::median(decrease),
            median_share_changed     = stats::median(share_changed),
            median_share_changed_dir = stats::median(share_changed_dir),
            pooled_share_increase    = sum(increase) / sum(pool),
            pooled_share_decrease    = sum(decrease) / sum(pool),
            pooled_share_changed     = sum(changed) / sum(pool),
            pooled_share_changed_dir = sum(increase_dir + decrease_dir) / sum(pool),
            increase_to_decrease     = sum(increase) / sum(decrease),
            fold_of_dir              = sum(changed) / sum(increase_dir + decrease_dir),
            # Figure 8: pooled departure shares (summed counts / summed species).
            pooled_share_depart_up   = sum(depart_up) / sum(offset_n),
            pooled_share_depart_down = sum(depart_down) / sum(offset_n),
            pooled_share_departed    = sum(departed) / sum(offset_n),
            .groups = "drop")

# Figure 8: BCR community offset per draw, averaging all species-hex log ratios in
# the BCR's mapped hexes (keeps each species' correlation across hexes; the
# interval excludes biases shared by all species).
bcr_offset <- tibble(bcr = hex_by_bcr$bcr, offset_pct_q50 = NA_real_,
                     offset_pct_q05 = NA_real_, offset_pct_q95 = NA_real_)
for (b in seq_len(nrow(bcr_offset))) {
  ids   <- hex_summary$hex_id[hex_summary$bcr %in% bcr_offset$bcr[b] & hex_summary$pool >= MIN_POOL_FOR_PROP]
  rows  <- match(ids, hex_ids)
  draws <- colSums(hex_lr_sum[rows, , drop = FALSE]) / sum(hex_lr_n[rows])
  q     <- stats::quantile(draws, probs = c(0.05, 0.5, 0.95), names = FALSE)
  bcr_offset$offset_pct_q05[b] <- 100 * (exp(q[1]) - 1)
  bcr_offset$offset_pct_q50[b] <- 100 * (exp(q[2]) - 1)
  bcr_offset$offset_pct_q95[b] <- 100 * (exp(q[3]) - 1)
}
hex_by_bcr <- hex_by_bcr |>
  left_join(bcr_offset, by = "bcr") |>
  arrange(match(bcr, BCR_ORDER))
utils::write.csv(hex_by_bcr, file.path(OUT_DIR, "hex_species_change_by_bcr.csv"), row.names = FALSE)
message(sprintf("\nFigure 7 by BCR (hexes assigned to their majority BCR; fold_of_dir = share of direction-supported changes also >= %.2f-fold):",
                HEX_FOLD))
print(as.data.frame(hex_by_bcr), digits = 3, row.names = FALSE)

# Map only hexes with a usable pool (others are mostly open water).
hex_sf <- dplyr::inner_join(hex_grid, hex_summary |> filter(pool >= MIN_POOL_FOR_PROP), by = "hex_id")
hex_panel <- function(col, title, fill_scale) {
  ggplot(hex_sf) +
    geom_sf(aes(fill = .data[[col]]), colour = NA) +
    fill_scale +
    labs(title = title) +
    theme_void() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}
hn_lim <- c(0, max(hex_sf$increase, hex_sf$decrease))
hs_lim <- c(0, max(hex_sf$share_increase, hex_sf$share_decrease, na.rm = TRUE))

fig7 <- (hex_panel("pool", "Species pool", scale_fill_viridis_c(name = "Species")) |
           hex_panel("increase", sprintf("Increased \u2265 %+.0f%%", 100 * (HEX_FOLD - 1)),
                     seq_fill("#542788", "Species", hn_lim)) |
           hex_panel("decrease", sprintf("Decreased \u2265 %.0f%%", 100 * (1 - 1 / HEX_FOLD)),
                     seq_fill("#B35806", "Species", hn_lim))) /
  (hex_panel("share_changed", "Share of pool changed (either direction)",
             seq_fill("grey15", "Share of\npool", c(0, max(hex_sf$share_changed, na.rm = TRUE)),
                      scales::label_percent())) |
     hex_panel("share_increase", "Share of pool increased",
               seq_fill("#542788", "Share of\npool", hs_lim, scales::label_percent())) |
     hex_panel("share_decrease", "Share of pool decreased",
               seq_fill("#B35806", "Share of\npool", hs_lim, scales::label_percent()))) +
  plot_annotation(
    title    = sprintf("Species with well-supported changes of at least %+.0f%% / %.0f%% (25-km hexagons)",
                       100 * (HEX_FOLD - 1), -100 * (1 - 1 / HEX_FOLD)),
    subtitle = sprintf(paste0("A species changed if P(Atlas 3 / Atlas 2 \u2265 %.2f) or P(\u2264 %.3f) \u2265 %.2f, ",
                              "computed exactly from its posterior draws of hex mean abundance.\n",
                              "Species pool: hex-mean P(\u2265 1 detection in %d counts) > %.2f in either atlas; ",
                              "each counts as one."),
                       HEX_FOLD, 1 / HEX_FOLD, HEX_P, POBS_N_SURVEYS, HEX_POOL_POBS),
    caption  = sprintf(paste0("%d modelled species with hex draws. Hexes shown where the pool has \u2265 %d ",
                              "species. Hexes with < %.0f%% of the fullest hex's modelled cells dropped. Direction-only ",
                              "counts (P(A3 > A2) \u2265 %.2f)\nare in hex_species_change_summary.csv. Counts show ",
                              "evidence of change, which depends on survey coverage; biases shared by all species are ",
                              "not represented."),
                       n_sp_hex, MIN_POOL_FOR_PROP, 100 * HEX_MIN_CELL_PROP, HEX_P),
    theme    = theme(plot.title = element_text(size = 16, face = "bold"))
  )
ggsave(filename = file.path(OUT_DIR, "fig7_species_change_hex.png"),
       plot = fig7, width = 18, height = 14, dpi = FIG_DPI, bg = "white")

message(sprintf(paste0("\nFigure 7 (hexes, pool PObs > %.2f, P >= %.2f): median pool %.0f; median species ",
                       "changed >= %.2f-fold %.1f (share %.2f); direction only %.1f."),
                HEX_POOL_POBS, HEX_P, stats::median(hex_summary$pool), HEX_FOLD,
                stats::median(hex_summary$changed), stats::median(hex_summary$share_changed, na.rm = TRUE),
                stats::median(hex_summary$increase_dir + hex_summary$decrease_dir)))

# ============================================================
# Figure 8: community offset and species departures (hexagons)
# ============================================================

# Offset: diverging on the log scale, labelled as percent change.
dlim_off <- stats::quantile(abs(hex_sf$offset_lr_q50), probs = 0.99, na.rm = TRUE)
dep_lim  <- c(0, max(hex_sf$share_depart_up, hex_sf$share_depart_down, na.rm = TRUE))

fig8 <- (hex_panel("offset_lr_q50", "Community offset: geometric-mean change across species",
                   fill_diff(dlim = dlim_off, name = "Atlas 3 vs\nAtlas 2", labels = pct_label)) |
           hex_panel("offset_ci_width", "Offset precision: ratio of 90% interval bounds",
                     seq_fill("grey15", "Upper /\nlower",
                              c(0, stats::quantile(hex_sf$offset_ci_width, probs = 0.99, na.rm = TRUE)),
                              function(x) sprintf("×%.2f", exp(x))))) /
  (hex_panel("share_depart_up",
             sprintf("Share of species ≥ %+.0f%% above the community offset", 100 * (HEX_FOLD - 1)),
             seq_fill("#542788", "Share of\nspecies", dep_lim, scales::label_percent())) |
     hex_panel("share_depart_down",
               sprintf("Share of species ≥ %.0f%% below the community offset", 100 * (1 - 1 / HEX_FOLD)),
               seq_fill("#B35806", "Share of\nspecies", dep_lim, scales::label_percent()))) +
  plot_annotation(
    title    = "Community-wide change and species departures from it (25-km hexagons)",
    subtitle = sprintf(paste0("Top: per posterior draw, the mean over pool species of log(Atlas 3 / Atlas 2) ",
                              "(each capped at ×%.0f); posterior median and 90%% interval width. A bias shared ",
                              "by all species (survey type, protocol,\nobservers, timing) moves the offset and is ",
                              "not in its interval. Bottom: species whose change differs from the offset by at least ",
                              "%.2f-fold with P ≥ %.2f; shared biases cancel here."),
                       exp(OFFSET_LR_CAP), HEX_FOLD, HEX_P),
    caption  = sprintf(paste0("%d modelled species with hex draws. Species pool: hex-mean P(≥ 1 detection in %d ",
                              "counts) > %.2f in either atlas. Hexes shown where the pool has ≥ %d species. ",
                              "Intervals assume independence among species (separate fits).\nDepartures treat the ",
                              "offset as fixed at its median. Like Figure 7, departure counts depend on survey ",
                              "coverage: sparsely surveyed hexes show fewer. Colour scales capped at the 99th ",
                              "percentile."),
                       n_sp_hex, POBS_N_SURVEYS, HEX_POOL_POBS, MIN_POOL_FOR_PROP),
    theme    = theme(plot.title = element_text(size = 16, face = "bold"))
  )
ggsave(filename = file.path(OUT_DIR, "fig8_community_offset_hex.png"),
       plot = fig8, width = 14, height = 14, dpi = FIG_DPI, bg = "white")

message("\nFigure 8 by BCR: community offset (%, median and 90% interval) and pooled departure shares:")
print(as.data.frame(hex_by_bcr |>
                      select(bcr_name, offset_pct_q50, offset_pct_q05, offset_pct_q95,
                             pooled_share_depart_up, pooled_share_depart_down)),
      digits = 3, row.names = FALSE)

# ============================================================
# Species table and run record
# ============================================================

species_table <- sp_tab |>
  full_join(prov |>
              transmute(sp_file, sp_english_chg = sp_english,
                        pct_change_median = 100 * prop_change_median,
                        pct_change_qlow   = 100 * prop_change_qlow,
                        pct_change_qhigh  = 100 * prop_change_qhigh,
                        p_increase, p_decrease, change_class),
            by = "sp_file") |>
  mutate(sp_english           = coalesce(sp_english, sp_english_chg),
         pobs_area_change_pct = ifelse(pobs_area_km2_A2 > 0,
                                       100 * (pobs_area_km2_A3 / pobs_area_km2_A2 - 1), NA_real_)) |>
  select(-sp_english_chg) |>
  relocate(sp_english, sp_file) |>
  arrange(sp_english)

utils::write.csv(species_table, file.path(OUT_DIR, "species_community_table.csv"), row.names = FALSE)

saveRDS(
  list(
    settings = list(model_name = model_name, pred_dir = PRED_DIR, rast_dir = RAST_DIR,
                    change_table = CHANGE_TABLE, data_ready_rds = DATA_READY_RDS,
                    pobs_n_surveys = POBS_N_SURVEYS,
                    pobs_marginalize_terms = POBS_MARGINALIZE_TERMS,
                    pobs_gh_nodes = POBS_GH_NODES, min_map_richness = MIN_MAP_RICHNESS,
                    p_very_likely = P_VERY_LIKELY, p_likely = P_LIKELY,
                    stable_fold = STABLE_FOLD, class_precedence = CLASS_PRECEDENCE,
                    ci_level = ci_level, pixel_ci_z = PIXEL_CI_Z,
                    pixel_p_primary = PIXEL_P_PRIMARY, pixel_p_approx = PIXEL_P_APPROX,
                    pool_pobs_thresholds = POOL_POBS_THRESHOLDS, min_pool_for_prop = MIN_POOL_FOR_PROP,
                    hex_fold = HEX_FOLD, hex_p = HEX_P, hex_pool_pobs = HEX_POOL_POBS,
                    hex_min_cell_prop = HEX_MIN_CELL_PROP, offset_lr_cap = OFFSET_LR_CAP,
                    bcr_order = BCR_ORDER,
                    drop_untrimmed_province = DROP_UNTRIMMED_PROVINCE,
                    date_run = Sys.time()),
    species_table  = species_table,
    pool_summary   = pool_summary,
    hex_summary    = hex_summary,
    hex_by_bcr     = hex_by_bcr,
    change_classes = chg |> select(sp_english, sp_file, region_type, Region_Name, change_class)
  ),
  file.path(OUT_DIR, "community_summaries.rds")
)

message("\n11_community_richness_and_change.R complete -> ", OUT_DIR)