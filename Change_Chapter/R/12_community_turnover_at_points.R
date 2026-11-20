# ============================================================
# 12_community_turnover_at_points.R
#
# Purpose
#   For a region chosen by the user (a point and a radius), set the model's
#   predicted change between Atlas 2 and Atlas 3 beside the raw survey data
#   collected there, so that reviewers can judge how far the model estimates
#   are reflected by the data. Also diagnoses which species drive the
#   community-level change mapped by 11 (Figures 1 and 4).
#
# Run AFTER 11 (reads its community_summaries.rds and rasters).
#
# Regional model change (used for the arrows and "supported" labels)
#   From the posterior draws of 25-km hexagon mean abundance saved by 07
#   (hex_draws_Corrected_for_Water, as in 11's Figure 7). The region is the set
#   of hexagons whose centres lie within MODEL_RADIUS_KM of the point (the hex
#   containing the point when MODEL_RADIUS_KM = 0). Per draw, hexagon means are
#   averaged, each weighted by its number of modelled 1-km pixels, giving the
#   regional mean in each atlas; their ratio A3 / A2 per draw is the posterior
#   of regional change. Each draw is a joint sample across hexagons, so the
#   spatial correlation is kept and the interval is valid for the region.
#   Reported: posterior medians of the regional means (scaled by the omitted
#   iid variance, as below, so they are expected counts per survey), median
#   and REGION_CI interval of the change, P(increase), and P(change at least
#   REGION_FOLD-fold up / down). "Supported increase" = P(A3 / A2 >=
#   REGION_FOLD) >= REGION_P (direction only when REGION_FOLD = 1).
#
# Pixel-level model values per location in POINTS (longitude / latitude, WGS84)
# come from the 1-km pixel containing the point (MODEL_RADIUS_KM = 0) or the
# mean over pixels within MODEL_RADIUS_KM of it:
#   Expected count per survey, from 09's posterior-median rasters with the
#     omitted iid variance integrated back in (a mean across sites); plug-in
#     values without uncertainty, used where the hex draws are missing and in D1.
#   PObs in each atlas across POBS_N_SURVEYS counts and per single survey
#     (p1). dpobs = PObs(A3) - PObs(A2); summed over species it is the
#     richness change in 11's Figure 1.
#   P(increase), from 07's pixel posterior summary of absolute change (split
#     normal, as in 11), averaged over the pixels in each raster cell before
#     classifying (as in 11). "Supported" increase / decrease = P(direction) >=
#     11's pixel_p_primary. direction_pixel = the majority of cells in the
#     region individually supported (kept to compare with 11's Figure 4; it
#     depends on pixel precision, so it is not a test of regional change).
#     fig4_pool / fig4_increase / fig4_decrease = the
#     species' contribution to 11's Figure 4 primary variant for the species
#     pool POOL_POBS_THRESHOLD; summed over species they reproduce
#     species_pool_change_<pool>.tif at the point.
#   Province-wide change class from 08 (as in 11's Figures 2-3), for context.
# Raw data within RADIUS_KM, per species and atlas x survey type: surveys,
#   squares, detections, proportion of surveys with a detection, and mean count
#   per survey (pooled, and as the mean of per-square means). Naive: not
#   adjusted for time of day, date, duration, habitat or observer.
#
# Reviewer diagnostics (what each tests, and why)
#   D1 Like-for-like comparison in shared squares (fig: _shared_squares.png).
#      Raw point counts restricted to squares point-counted in BOTH atlases,
#      beside the model evaluated at the cells holding those same point counts
#      (squares weighted equally on both sides). Same places, same survey type
#      in both atlases, so the comparison is free of the ARU-vs-point-count
#      adjustment and of differences in where each atlas surveyed. If the
#      circle-wide model change differs from this, the circle-wide change is
#      coming from other data (ARUs, other squares) or from spatial smoothing.
#   D2 Agreement across ALL species (fig: _agreement.png; csv: _agreement*.csv).
#      Model log ratio vs raw log ratio for every species with enough raw
#      detections, for the like-for-like comparison and for the whole region
#      (with the model's interval from the hexagon draws), with correlation,
#      share agreeing in direction, medians, and the share of species whose
#      raw change lies inside the model interval. Also the least-squares slope
#      of raw on model change (D2b, csv: _slopes.csv): a slope near 1 means the
#      model's changes are about the size the data support; clearly above 1
#      among well-sampled species would mean the model compresses change
#      (over-shrinkage). Species with a zero in one atlas's raw data are drawn
#      hollow and reported separately, because their raw ratio is set by
#      RATIO_OFFSET rather than by the data. The drivers
#      figure shows only the TOP_N largest modelled changes; this shows
#      whether those are representative. The model is fit to these data, so
#      agreement is not independent validation: it shows whether a species'
#      modelled change is consistent with the point counts. Disagreement
#      points to change driven by ARU data, by the survey-type adjustment or
#      by smoothing.
#   D3 Survey timing (fig: _timing.png; csv: _timing.csv). Day of year and
#      time relative to sunrise by atlas x survey type, within the radius and
#      in the shared squares. Raw means are not adjusted for timing; if Atlas 3
#      surveys were later in the season or the morning, part of a raw change
#      is timing rather than abundance. The shaded band is the window used to
#      select ARU recordings comparable with point counts.
#   D4 Data behind each raw mean (drivers and shared-squares figures).
#      Detection counts printed beside each species, and raw points from fewer
#      than FEW_DETECTIONS detections drawn hollow: a fall to zero means little
#      when it rests on two detections.
#   D5 Species detected per point count (fig: _species_per_pc.png; csv:
#      _species_per_survey.csv, _species_per_pc_change.csv,
#      _species_per_pc_squares.csv). Observed number of modelled species
#      detected per point count by atlas, pooled and square-weighted, beside
#      the model's expected number at the same survey locations (sum over
#      species of P(>= 1 detection in one standardised point count), as 11's
#      PObs with n = 1), and the number of species each survey was analysed
#      for. Tests whether the community-level change in the raw data matches
#      the model, and whether it could be an artefact of surveys entering
#      fewer species' analysed data in one atlas (e.g. safe dates).
#
# Outputs (OUT_DIR; <name> = POINTS$name)
#   turnover_<name>_species.csv   one row per species, sorted by dpobs
#   turnover_<name>_effort.csv    surveys near the point by atlas x survey type
#   turnover_<name>_species_per_survey.csv
#                                 species detected per survey by atlas x survey
#                                 type (all squares, and squares point-counted in
#                                 both atlases): observed pooled and square-
#                                 weighted, species analysed per survey, and model
#                                 expected species per point count (D5)
#   turnover_<name>_species_per_pc_change.csv   D5 Atlas 2 -> 3 change, observed vs model
#   turnover_<name>_species_per_pc_squares.csv  D5 per-square point-count means
#   turnover_<name>_species_per_pc.png          D5 figure
#   turnover_<name>_drivers.png   maps (top): point-count and ARU locations in
#                                 each atlas, province-wide and zoomed to the
#                                 focal region (RADIUS_KM circle); below, top
#                                 species by |dpobs|: model (left) and raw mean
#                                 counts near the point (right): point counts,
#                                 plus ARU_WINDOW_DURATION-min ARU recordings
#                                 within the point counts' temporal window
#                                 (WINDOW_COLS, WINDOW_PROBS) as "x" (D4)
#   turnover_<name>_shared_squares.png   D1, same species rows as the drivers figure
#   turnover_<name>_agreement.png        D2
#   turnover_<name>_agreement.csv        D2 summary statistics per comparison
#   turnover_<name>_agreement_species.csv D2 per-species log ratios
#   turnover_<name>_slopes.csv           D2b slopes of raw on model change, by subset
#   turnover_<name>_timing.png, turnover_<name>_timing.csv   D3
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# ============================================================
# CONFIG
# ============================================================

# ---- Input and output locations ----
COMMUNITY_DIR  <- "C:/Users/IlesD/Documents/Projects_Local/Ontario_Atlas/Change_Chapter/output/community_change"   # 11's OUT_DIR
PRED_DIR       <- "C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/Analysis/Ontario_BBA3_draft5/data_clean/model_output/predictions_PC_ARU_CL_SS"
RAST_DIR       <- "C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/Analysis/Ontario_BBA3_draft5/data_clean/model_output/rasters_PC_ARU_CL_SS"
DATA_USED_DIR  <- "C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/Analysis/Ontario_BBA3_draft5/data_clean/model_output/data_used_PC_ARU_CL_SS"
DATA_READY_RDS <- "C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/Analysis/Ontario_BBA3_draft5/data_clean/birds/data_ready_for_analysis.rds"
OUT_DIR        <- file.path(COMMUNITY_DIR, "turnover")

# ---- Region of interest ----
# Locations to query (decimal degrees). Names are used in output file names.
POINTS <- data.frame(
  name = c("northern decline blob"),
  lon  = c(-83.5865),
  lat  = c(51.3161)
)

MODEL_RADIUS_KM <- 125              # model: 0 = the pixel containing the point; > 0 = mean over pixels within this distance
RADIUS_KM       <- MODEL_RADIUS_KM  # raw data: surveys summarised within this distance of each point

# ---- Model values ----
POOL_POBS_THRESHOLD <- 0.25     # 11's Figure 4 species pool to reproduce; one of 11's POOL_POBS_THRESHOLDS

# Regional change from the 25-km hexagon posterior draws (see header).
REGION_FOLD <- 1.2    # "supported" change must be at least this fold (+20% / -17%); 1 = direction only
REGION_P    <- 0.95   # posterior probability required for a supported change
REGION_CI   <- 0.90   # width of the credible interval reported for regional change
MODEL_METRIC        <- "count"  # drivers figure, left panel: "count" = model expected count per
# survey; "pobs" = P(>= 1 detection in pobs_n_surveys counts)

# ---- Raw data summaries ----
MIN_GROUP_SURVEYS <- 20         # plot a raw mean only if its atlas x survey type has this many surveys
MIN_GROUP_SQUARES <- 3          # ... and covers at least this many atlas squares
RAW_MEAN_METHOD   <- "square"   # "square" = mean of per-square mean counts (each square weighted
# equally); "survey" = mean pooled over surveys

# ARU recordings shown with the point counts ("x"): only ARU_WINDOW_DURATION-min
# recordings inside the point counts' temporal window near each point, i.e.
# within the central WINDOW_PROBS range of every WINDOW_COLS variable over that
# point's Atlas 2 + Atlas 3 point counts.
ARU_WINDOW_DURATION <- 5
WINDOW_COLS         <- c("DayOfYear", "Sunrise_Relative_Hours")   # all_surveys columns defining the window
WINDOW_LABELS       <- c(DayOfYear = "Day of year", Sunrise_Relative_Hours = "Hours relative to sunrise")
WINDOW_PROBS        <- c(0.025, 0.975)
ARU_WINDOW_TYPE     <- "ARU, point-count window"           # label used in the outputs
SHARED_PC_TYPE      <- "Point_Count, shared squares"       # label used in the outputs (D1)

# ---- Reviewer diagnostics (D1-D4; see header) ----
RATIO_OFFSET   <- 0.01   # D2: added to both atlases' mean counts before log(A3 / A2), so zeros give finite ratios
MIN_DETECTIONS <- 5      # D2: a species enters a comparison only with >= this many raw detections (both atlases summed)
FEW_DETECTIONS <- 5      # D4: raw point-count means from fewer detections than this are drawn hollow
N_LABEL_OUTLIERS <- 8    # D2: label this many species furthest from the 1:1 line in each comparison
SLOPE_DETECTION_TIERS <- c(20, 50)   # D2b: also fit slopes for species with at least these many raw detections

# ---- Figures ----
TOP_N            <- 30           # species shown in the drivers and shared-squares figures
COUNT_AXIS       <- "pseudo_log" # count axes: "pseudo_log", "sqrt" or "linear"
COUNT_AXIS_SIGMA <- 0.01         # pseudo_log: roughly linear below this count, logarithmic above
COUNT_BREAKS     <- c(0, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300)   # pseudo_log axis labels
ARROW_GAP_PROP   <- 0.012        # arrows stop this share of the x-axis range short of the Atlas 3 point

ATLAS_LABELS  <- c(OBBA2 = "Atlas 2", OBBA3 = "Atlas 3")
ATLAS_COLOURS <- c("Atlas 2" = "#E69F00", "Atlas 3" = "#0072B2")   # orange / blue
DIR_COLOURS   <- c("Supported increase" = "#542788", "Supported decrease" = "#B35806",
                   "Not supported" = "grey65")
FIG_DPI       <- 300

# Maps of survey locations (top of the drivers figure).
MAP_SURVEY_TYPES <- c("Point_Count", "ARU")           # point-based survey types mapped
SURVEY_COLOURS   <- c("Point count" = "grey25", "ARU" = "#CC79A7")
MAP_ZOOM_KM      <- MODEL_RADIUS_KM * 1.25            # half-width of the zoomed maps around the point
MAP_HEIGHT_IN    <- 4.5                               # height of the map row in the saved figure

# ============================================================
# Helpers
# ============================================================

# P(>= 1 detection in n_surveys) from median expected counts; identical to 11.
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

# P(change > 0) from a posterior median and interval; identical to 11.
p_increase_split <- function(q50, lower, upper, z) {
  s_lo <- pmax((q50 - lower) / z, 1e-12)
  s_hi <- pmax((upper - q50) / z, 1e-12)
  ifelse(q50 >= 0, stats::pnorm(q50 / s_lo), 1 - stats::pnorm(-q50 / s_hi))
}

# End of an Atlas 2 -> Atlas 3 arrow, stopped `gap` short of the Atlas 3 point so
# the arrowhead is not hidden under it; NA (no arrow) where the change is too
# small to draw.
arrow_end <- function(x2, x3, gap) {
  ifelse(abs(x3 - x2) > 2 * gap, x3 - sign(x3 - x2) * gap, NA_real_)
}

# Count summaries per atlas x survey type: pooled over surveys, and the mean of
# per-square mean counts (each atlas square weighted equally, so heavily
# surveyed squares do not dominate). Surveys without a square_id are left out of
# the square-level mean only.
summarise_counts <- function(d) {
  per_square <- d |>
    filter(!is.na(square_id)) |>
    group_by(Atlas, Survey_Type, square_id) |>
    summarise(square_mean = mean(count), .groups = "drop") |>
    group_by(Atlas, Survey_Type) |>
    summarise(n_squares = n(), mean_count_square = mean(square_mean), .groups = "drop")
  d |>
    group_by(Atlas, Survey_Type) |>
    summarise(n_surveys = n(), n_det = sum(count > 0), total_count = sum(count), .groups = "drop") |>
    left_join(per_square, by = c("Atlas", "Survey_Type"))
}

# Axis transformation for counts (COUNT_AXIS). On a log-type axis an arrow's
# length shows proportional change, comparable across common and rare species.
count_trans <- switch(COUNT_AXIS,
                      pseudo_log = scales::pseudo_log_trans(sigma = COUNT_AXIS_SIGMA, base = 10),
                      sqrt       = scales::sqrt_trans(),
                      linear     = scales::identity_trans())

# x scale with an optional transformation (ggplot2 >= 3.5.0 renamed trans to transform).
x_scale <- function(tr, limits, breaks = waiver()) {
  args <- list(limits = limits, breaks = breaks, labels = function(x) as.character(x))   # no 1e+01 style
  if (utils::packageVersion("ggplot2") >= "3.5.0") args$transform <- tr else args$trans <- tr
  do.call(scale_x_continuous, args)
}

# Arrow end on a transformed axis: the gap is taken in transformed units, so
# arrowheads stop the same visual distance short of the Atlas 3 point everywhere.
arrow_end_t <- function(x2, x3, limits, tr) {
  gap <- ARROW_GAP_PROP * diff(tr$transform(limits))
  tr$inverse(arrow_end(x2 = tr$transform(x2), x3 = tr$transform(x3), gap = gap))
}

# Position a share `prop` of the axis range past `x` on a transformed axis.
beyond <- function(x, prop, limits, tr) tr$inverse(tr$transform(x) + prop * diff(tr$transform(limits)))

# log(A3 / A2) of mean counts, with RATIO_OFFSET added so zeros give finite ratios (D2).
log_ratio <- function(a2, a3) log((a3 + RATIO_OFFSET) / (a2 + RATIO_OFFSET))

# Regional change for one species from its hexagon draws (see header).
#   hd      the species' hex_draws_Corrected_for_Water (hex_id; mu_OBBA2 and
#           mu_OBBA3: one vector of draws per hexagon)
#   region  list(hex_id, weight) for the region
#   scale   multiplier putting hexagon means on the expected-count scale
# Returns NA values when the species has no draws for the region.
region_change <- function(hd, region, scale) {
  out <- list(region_A2 = NA_real_, region_A3 = NA_real_, region_ratio_q50 = NA_real_,
              region_ratio_low = NA_real_, region_ratio_high = NA_real_, region_p_increase = NA_real_,
              region_p_up_fold = NA_real_, region_p_down_fold = NA_real_, region_n_hex = 0L)
  if (is.null(hd)) return(out)
  j    <- match(region$hex_id, hd$hex_id)
  keep <- !is.na(j)
  if (!any(keep)) return(out)
  m2 <- do.call(rbind, hd$mu_OBBA2[j[keep]])        # hexagons x draws
  m3 <- do.call(rbind, hd$mu_OBBA3[j[keep]])
  if (ncol(m2) != ncol(m3)) stop("Atlas 2 and Atlas 3 hexagon draws differ in number.", call. = FALSE)
  w   <- region$weight[keep] / sum(region$weight[keep])
  r2  <- colSums(w * m2, na.rm = TRUE)              # regional mean per draw (matrix x weight recycles by row)
  r3  <- colSums(w * m3, na.rm = TRUE)
  rat <- r3 / r2
  rat <- rat[is.finite(rat)]
  if (length(rat) == 0) return(out)
  q <- stats::quantile(rat, probs = c((1 - REGION_CI) / 2, 0.5, 1 - (1 - REGION_CI) / 2), names = FALSE)
  list(region_A2          = scale * stats::median(r2),
       region_A3          = scale * stats::median(r3),
       region_ratio_q50   = q[2],
       region_ratio_low   = q[1],
       region_ratio_high  = q[3],
       region_p_increase  = mean(rat > 1),
       region_p_up_fold   = mean(rat >= REGION_FOLD),
       region_p_down_fold = mean(rat <= 1 / REGION_FOLD),
       region_n_hex       = sum(keep))
}

# Least-squares fit of raw on model log ratio (D2b): slope with its 95%
# confidence interval, and the intercept (for drawing the line). NA with < 5 species.
slope_fit <- function(d) {
  out <- list(n_species = nrow(d), slope = NA_real_, slope_low = NA_real_,
              slope_high = NA_real_, intercept = NA_real_)
  if (nrow(d) < 5) return(out)
  fit <- stats::lm(formula = raw_lr ~ model_lr, data = d)
  ci  <- stats::confint(object = fit)["model_lr", ]
  out$slope      <- unname(stats::coef(fit)[["model_lr"]])
  out$slope_low  <- unname(ci[1])
  out$slope_high <- unname(ci[2])
  out$intercept  <- unname(stats::coef(fit)[["(Intercept)"]])
  out
}

# Axis labels for a log ratio, as percent change ("+50%", "-33%").
pct_lab <- function(x) {
  pct <- round(100 * (exp(x) - 1))
  ifelse(pct > 0, paste0("+", pct, "%"), paste0(pct, "%"))
}

# "Atlas 2 / Atlas 3" detection counts for a figure row (D4); "-" where no data.
det_label <- function(d2, d3) {
  f <- function(v) ifelse(is.na(v), "-", format(v, big.mark = ","))
  paste0(f(d2), " / ", f(d3))
}

# Point shapes: model and well-supported raw means solid; raw means from fewer
# than FEW_DETECTIONS detections hollow (D4); ARU recordings as "x".
FEW_LAB <- sprintf("Point count, < %d detections", FEW_DETECTIONS)
ARU_LAB <- sprintf("ARU %d-min, point-count window", ARU_WINDOW_DURATION)
SHAPES  <- setNames(c(16, 16, 1, 4), c("Model", "Point count", FEW_LAB, ARU_LAB))

# One panel of Atlas 2 -> Atlas 3 arrows, one row per species, on the count axis.
#   d               data frame with sp_english (factor, sets row order), a2, a3;
#                   with_detections = TRUE also needs det_a2, det_a3 and det_lab
#   with_detections raw-data panel: hollow points for thin means, and each
#                   row's detection counts printed on the right (D4)
#   x_values        further values the x-axis must reach (e.g. ARU means)
#   y_expand        vertical padding (in rows), e.g. room for a text box
arrow_panel <- function(d, title, xlab, with_detections = FALSE, x_values = NULL,
                        y_expand = expansion(add = 0.6)) {
  v_max <- max(c(d$a2, d$a3, x_values), na.rm = TRUE)
  if (!is.finite(v_max) || v_max <= 0) v_max <- 1
  room  <- if (with_detections) 0.22 else 0.05          # right-hand space for detection labels
  lim   <- c(0, beyond(x = v_max, prop = room, limits = c(0, v_max), tr = count_trans))
  lab_x <- beyond(x = v_max, prop = 0.03, limits = c(0, v_max), tr = count_trans)
  if (with_detections) {
    d$shape_a2 <- ifelse(d$det_a2 %in% NA | d$det_a2 >= FEW_DETECTIONS, "Point count", FEW_LAB)
    d$shape_a3 <- ifelse(d$det_a3 %in% NA | d$det_a3 >= FEW_DETECTIONS, "Point count", FEW_LAB)
  } else {
    d$shape_a2 <- "Model"
    d$shape_a3 <- "Model"
  }
  p <- ggplot(d, aes(y = sp_english)) +
    geom_segment(aes(x = a2, xend = arrow_end_t(x2 = a2, x3 = a3, limits = lim, tr = count_trans),
                     yend = sp_english),
                 colour = "grey40", arrow = grid::arrow(length = unit(0.15, "cm"), type = "closed"),
                 na.rm = TRUE) +
    geom_point(aes(x = a2, colour = "Atlas 2", shape = shape_a2), size = 2, stroke = 0.9, na.rm = TRUE) +
    geom_point(aes(x = a3, colour = "Atlas 3", shape = shape_a3), size = 2, stroke = 0.9, na.rm = TRUE) +
    scale_colour_manual(values = ATLAS_COLOURS, name = NULL) +
    scale_shape_manual(values = SHAPES, name = NULL, guide = if (with_detections) "legend" else "none") +
    scale_y_discrete(drop = FALSE, expand = y_expand) +
    x_scale(tr = count_trans, limits = lim, breaks = if (COUNT_AXIS == "pseudo_log") COUNT_BREAKS else waiver()) +
    labs(x = xlab, y = NULL, title = title) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top", legend.box = "vertical",       # colour and shape legends stacked
          panel.grid.major.y = element_line(colour = "grey92"))
  if (with_detections) {
    p <- p + geom_text(aes(x = lab_x, label = det_lab), hjust = 0, size = 2.6, colour = "grey30")
  }
  p
}

km_units <- function(crs) {
  u <- crs$units_gdal
  if (isTRUE(u %in% c("metre", "meter")))         return(1000)
  if (isTRUE(u %in% c("kilometre", "kilometer"))) return(1)
  stop("Unsupported CRS units: ", u, call. = FALSE)
}

# ============================================================
# Inputs
# ============================================================
# Check that every input exists and that the settings agree with 11's run,
# then read the survey table and the prediction grid.

inputs  <- c(COMMUNITY_SUMMARIES = file.path(COMMUNITY_DIR, "community_summaries.rds"),
             RICHNESS_CHANGE     = file.path(COMMUNITY_DIR, "richness_change_A3_minus_A2.tif"),
             POOL_CHANGE         = file.path(COMMUNITY_DIR, sprintf("species_pool_change_pobs%03d.tif",
                                                                    round(100 * POOL_POBS_THRESHOLD))),
             PRED_DIR = PRED_DIR, RAST_DIR = RAST_DIR, DATA_USED_DIR = DATA_USED_DIR,
             DATA_READY_RDS = DATA_READY_RDS)
missing <- inputs[!file.exists(inputs)]
if (length(missing) > 0) {
  stop("Not found (check CONFIG; COMMUNITY_DIR must hold the revised 11's outputs):\n  ",
       paste(names(missing), missing, sep = " = ", collapse = "\n  "), call. = FALSE)
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cs       <- readRDS(inputs[["COMMUNITY_SUMMARIES"]])
settings <- cs$settings
for (nm in c("pobs_n_surveys", "pobs_gh_nodes", "pixel_p_primary", "pool_pobs_thresholds", "pixel_ci_z")) {
  if (is.null(settings[[nm]])) stop("community_summaries.rds lacks settings$", nm, "; re-run the revised 11.", call. = FALSE)
}
if (REGION_FOLD < 1 || REGION_P <= 0.5 || REGION_CI <= 0 || REGION_CI >= 1) {
  stop("Need REGION_FOLD >= 1, REGION_P > 0.5 and 0 < REGION_CI < 1.", call. = FALSE)
}
if (!all(WINDOW_COLS %in% names(WINDOW_LABELS))) stop("Give every WINDOW_COLS entry a label in WINDOW_LABELS.", call. = FALSE)
if (!RAW_MEAN_METHOD %in% c("square", "survey")) stop("RAW_MEAN_METHOD must be \"square\" or \"survey\".", call. = FALSE)
if (!COUNT_AXIS %in% c("pseudo_log", "sqrt", "linear")) stop("COUNT_AXIS must be \"pseudo_log\", \"sqrt\" or \"linear\".", call. = FALSE)
if (!MODEL_METRIC %in% c("count", "pobs")) stop("MODEL_METRIC must be \"count\" or \"pobs\".", call. = FALSE)
if (!any(abs(settings$pool_pobs_thresholds - POOL_POBS_THRESHOLD) < 1e-9)) {
  stop("POOL_POBS_THRESHOLD = ", POOL_POBS_THRESHOLD, " is not one of 11's pool thresholds (",
       paste(settings$pool_pobs_thresholds, collapse = ", "), ").", call. = FALSE)
}
sp_info <- cs$species_table |>
  filter(!is.na(sigma2_omit)) |>
  transmute(sp_english, sp_file, sigma2_omit, nb_size,
            prov_change_class = change_class, prov_pct_change_median = pct_change_median)

template <- terra::rast(file.path(RAST_DIR, paste0(sp_info$sp_file[1], "_a2.tif")))[["mu_q50"]]
rast_crs <- sf::st_crs(terra::crs(template))
u        <- km_units(rast_crs)

dat         <- readRDS(DATA_READY_RDS)
missing_window <- setdiff(WINDOW_COLS, names(dat$all_surveys))
if (length(missing_window) > 0) {
  stop("all_surveys lacks WINDOW_COLS ", paste(missing_window, collapse = ", "), ". Available columns: ",
       paste(setdiff(names(dat$all_surveys), attr(dat$all_surveys, "sf_column")), collapse = ", "),
       call. = FALSE)
}
all_surveys <- dat$all_surveys |>
  select(any_of(c("survey_id", "Atlas", "Survey_Type", "Survey_Duration_Minutes", "square_id")),
         all_of(WINDOW_COLS))

# Prediction-grid pixels -> raster cells (as in 11).
grid2 <- dat$grid_OBBA2
if (is.null(grid2)) stop("No grid_OBBA2 in ", DATA_READY_RDS, call. = FALSE)
grid_geom <- sf::st_geometry(grid2)
if (!all(sf::st_geometry_type(grid_geom) == "POINT")) grid_geom <- suppressWarnings(sf::st_centroid(grid_geom))
grid_xy       <- sf::st_coordinates(sf::st_transform(grid_geom, crs = rast_crs))[, 1:2]
grid_pixel_id <- grid2$pixel_id
grid_cell     <- terra::cellFromXY(template, grid_xy)
# Province and BCR outlines for the maps, in the raster CRS.
if (is.null(dat$bcr_sf)) stop("No bcr_sf in ", DATA_READY_RDS, call. = FALSE)
bcr_map      <- sf::st_transform(sf::st_make_valid(dat$bcr_sf), crs = rast_crs)
prov_outline <- sf::st_union(bcr_map)

# 25-km hexagons (as in 11): each raster cell's hexagon, and each hexagon's
# weight = its number of modelled 1-km pixels (so partly modelled hexagons on
# the coast or the province edge count for less in the regional mean).
if (is.null(dat$hex_grid_25km)) stop("No hex_grid_25km in ", DATA_READY_RDS, call. = FALSE)
hex_grid    <- sf::st_transform(dat$hex_grid_25km, crs = rast_crs)
cell_hex    <- terra::values(terra::rasterize(terra::vect(hex_grid), template, field = "hex_id"), mat = FALSE)
hex_weight  <- table(cell_hex[grid_cell[!is.na(grid_cell)]])
hex_centres <- sf::st_coordinates(suppressWarnings(sf::st_centroid(sf::st_geometry(hex_grid))))[, 1:2]
rm(dat, grid2, grid_geom)

pts    <- sf::st_as_sf(POINTS, coords = c("lon", "lat"), crs = 4326)
pts_xy <- sf::st_coordinates(sf::st_transform(pts, crs = rast_crs))

# Model pixels per point: the pixel containing it, or all pixels within MODEL_RADIUS_KM.
model_px <- lapply(seq_len(nrow(POINTS)), function(p) {
  k <- if (MODEL_RADIUS_KM > 0) {
    which(sqrt((grid_xy[, 1] - pts_xy[p, 1])^2 + (grid_xy[, 2] - pts_xy[p, 2])^2) / u <= MODEL_RADIUS_KM)
  } else {
    which(grid_cell == terra::cellFromXY(template, pts_xy[p, , drop = FALSE]))
  }
  if (length(k) == 0) stop(POINTS$name[p], " is outside the prediction grid.", call. = FALSE)
  list(cell = grid_cell[k], pixel_id = grid_pixel_id[k])
})
rm(grid_xy)
model_lab <- if (MODEL_RADIUS_KM > 0) sprintf("mean of pixels within %g km", MODEL_RADIUS_KM) else "1-km pixel"

# Region hexagons per point: hexagons whose centres lie within MODEL_RADIUS_KM
# (or the hexagon containing the point), with their pixel-count weights.
region_hex <- lapply(seq_len(nrow(POINTS)), function(p) {
  ids <- if (MODEL_RADIUS_KM > 0) {
    hex_grid$hex_id[sqrt((hex_centres[, 1] - pts_xy[p, 1])^2 + (hex_centres[, 2] - pts_xy[p, 2])^2) / u <= MODEL_RADIUS_KM]
  } else {
    cell_hex[terra::cellFromXY(template, pts_xy[p, , drop = FALSE])]
  }
  w <- as.numeric(hex_weight[as.character(ids)])
  keep <- !is.na(w) & w > 0                       # hexagons with modelled pixels only
  if (!any(keep)) stop(POINTS$name[p], ": no modelled 25-km hexagons in the region.", call. = FALSE)
  message(sprintf("%s: region = %d hexagons (%s modelled 1-km pixels).",
                  POINTS$name[p], sum(keep), format(sum(w[keep]), big.mark = ",")))
  list(hex_id = ids[keep], weight = w[keep])
})
region_lab <- if (MODEL_RADIUS_KM > 0) sprintf("25-km hexagons centred within %g km", MODEL_RADIUS_KM) else
  "the 25-km hexagon containing the point"

# Surveys within RADIUS_KM of each point.
near <- sf::st_is_within_distance(sf::st_transform(pts, crs = sf::st_crs(all_surveys)),
                                  all_surveys, dist = units::set_units(RADIUS_KM, "km"))
# Survey locations (x, y in the raster CRS) are kept so model cells can be
# matched to the places where surveys were done (D1).
near_surveys <- lapply(near, function(k) {
  xy <- sf::st_coordinates(sf::st_transform(all_surveys[k, ], crs = rast_crs))
  sf::st_drop_geometry(all_surveys[k, ]) |>
    mutate(Survey_Type = ifelse(Survey_Type == "ARU",       # ARU durations have separate intercepts
                                paste0("ARU ", Survey_Duration_Minutes, "-min"), Survey_Type),
           x = xy[, 1],
           y = xy[, 2])
})
# ARU recordings within the point counts' temporal window, per point.
# Why: ARUs record across the whole season and day, point counts mostly on
# June mornings; restricting ARUs to the point counts' window makes the "x"
# means comparable in timing with the point-count means (D3 shows the window).
pc_window <- vector("list", length(near_surveys))
for (p in seq_along(near_surveys)) {
  x      <- near_surveys[[p]]
  is_pc  <- x$Survey_Type == "Point_Count"
  in_win <- x$Survey_Type == paste0("ARU ", ARU_WINDOW_DURATION, "-min")
  lims   <- matrix(NA_real_, nrow = 2, ncol = length(WINDOW_COLS),
                   dimnames = list(c("low", "high"), WINDOW_COLS))
  for (col in WINDOW_COLS) {
    lims[, col] <- stats::quantile(x[[col]][is_pc], probs = WINDOW_PROBS, na.rm = TRUE, names = FALSE)
    in_win      <- in_win & (x[[col]] >= lims["low", col] & x[[col]] <= lims["high", col]) %in% TRUE
  }
  near_surveys[[p]]$aru_pc_window <- in_win
  pc_window[[p]] <- lims
}

# ---- D1 set-up: squares point-counted in both atlases ----
# These are the only places where the raw data compare like with like: same
# squares, same survey type in both atlases. For each point we record the
# shared squares, flag their point counts, and list the model cells holding
# those point counts (both atlases pooled, so the model is evaluated at the
# same cells in Atlas 2 and Atlas 3).
shared_sq    <- vector("list", length(near_surveys))
shared_cells <- vector("list", length(near_surveys))
for (p in seq_along(near_surveys)) {
  x   <- near_surveys[[p]]
  pcs <- x[x$Survey_Type == "Point_Count" & !is.na(x$square_id), ]
  shared_sq[[p]] <- intersect(pcs$square_id[pcs$Atlas == "OBBA2"], pcs$square_id[pcs$Atlas == "OBBA3"])
  in_shared      <- pcs$square_id %in% shared_sq[[p]]
  shared_cells[[p]] <- tibble(square_id = pcs$square_id[in_shared],
                              cell      = terra::cellFromXY(template, cbind(pcs$x[in_shared], pcs$y[in_shared]))) |>
    filter(!is.na(cell)) |>
    distinct()
  near_surveys[[p]]$shared_pc <- x$Survey_Type == "Point_Count" & (x$square_id %in% shared_sq[[p]])
  message(sprintf("%s: %d squares point-counted in both atlases (%d model cells with point counts).",
                  POINTS$name[p], length(shared_sq[[p]]), nrow(shared_cells[[p]])))
}

# Survey locations for the maps: point-based survey types, in the raster CRS,
# one row per location (rounded to 100 m) x atlas x survey type.
map_surveys <- all_surveys |>
  filter(Survey_Type %in% MAP_SURVEY_TYPES) |>
  sf::st_transform(crs = rast_crs)
map_xy  <- sf::st_coordinates(map_surveys)
grid_m  <- 0.1 * u                                     # 100 m in CRS units
map_pts <- tibble(Atlas       = unname(ATLAS_LABELS[map_surveys$Atlas]),
                  Survey_Type = ifelse(map_surveys$Survey_Type == "ARU", "ARU", "Point count"),
                  x           = round(map_xy[, 1] / grid_m) * grid_m,
                  y           = round(map_xy[, 2] / grid_m) * grid_m) |>
  distinct()
rm(map_surveys, map_xy)
rm(all_surveys)

# ============================================================
# Per-species values at every point
# ============================================================
# For each species and point: model values over the region (circle or pixel)
# and at the shared-square point counts (D1), and raw count summaries by atlas x
# survey type. One pass over species, because each species' prediction and
# data files are large and are read once here.

message("Extracting pixel predictions and nearby detections for ", nrow(sp_info), " species ...")

pix_rows  <- vector("list", nrow(sp_info))
data_rows <- vector("list", nrow(sp_info))

# Per-survey tallies near each point, filled species by species:
#   sv_det        number of modelled species detected on the survey
#   sv_used       analysed for at least one species (else it is not a real zero)
#   sv_n_analysed number of species whose analysed data include the survey (D5)
#   sv_model_all  model: expected number of species detected by one standardised
#                 point count at the survey's cell, summed over ALL modelled species (D5)
#   sv_model_an   the same, summed only over the species the survey was analysed
#                 for, i.e. the same species set the observed count can include (D5)
sv_det        <- lapply(near_surveys, function(x) integer(nrow(x)))
sv_used       <- lapply(near_surveys, function(x) logical(nrow(x)))
sv_n_analysed <- lapply(near_surveys, function(x) integer(nrow(x)))
sv_model_all  <- lapply(near_surveys, function(x) numeric(nrow(x)))
sv_model_an   <- lapply(near_surveys, function(x) numeric(nrow(x)))

# Raster cell of each nearby survey, and which atlas raster applies to it.
sv_cell     <- lapply(near_surveys, function(x) terra::cellFromXY(template, cbind(x$x, x$y)))
ATLAS_RAST  <- c(OBBA2 = "A2", OBBA3 = "A3")

for (k in seq_len(nrow(sp_info))) {
  
  s     <- sp_info[k, ]
  rasts <- lapply(c(A2 = "_a2.tif", A3 = "_a3.tif"),
                  function(suf) terra::rast(file.path(RAST_DIR, paste0(s$sp_file, suf)))[["mu_q50"]])
  pobs  <- function(l, n) pobs_from_lambda(lambda = l, sigma2 = s$sigma2_omit, size = s$nb_size,
                                           n_surveys = n, n_nodes = settings$pobs_gh_nodes)
  pred_path <- file.path(PRED_DIR, paste0(s$sp_file, "_1km.rds"))
  preds     <- if (file.exists(pred_path)) readRDS(pred_path) else NULL
  pchg      <- preds$abs_change_Corrected_for_Water   # pixel change summary (Figure 4 classification)
  hd        <- preds$hex_draws_Corrected_for_Water    # hexagon draws (regional change)
  rm(preds)
  
  pix_rows[[k]] <- bind_rows(lapply(seq_len(nrow(POINTS)), function(p) {
    mp    <- model_px[[p]]
    cells <- unique(mp$cell)                  # raster cells; several pixels can share one
    g     <- match(mp$cell, cells)            # each pixel's cell
    lam   <- lapply(rasts, function(r) terra::extract(r, cells)[, 1])
    P     <- lapply(lam, pobs, n = settings$pobs_n_surveys)
    P1    <- lapply(lam, pobs, n = 1L)
    
    # P(increase) per pixel, then averaged over the pixels in each cell (as in 11).
    p_cell <- rep(NA_real_, length(cells))
    if (!is.null(pchg)) {
      p_px <- rep(NA_real_, length(mp$cell))
      j    <- match(mp$pixel_id, pchg$pixel_id)
      ok   <- !is.na(j)
      p_px[ok] <- p_increase_split(q50 = pchg$abs_change_q50[j[ok]], lower = pchg$abs_change_lower[j[ok]],
                                   upper = pchg$abs_change_upper[j[ok]], z = settings$pixel_ci_z)
      # Sum and count of finite pixel values per cell in one pass (rowsum), then divide.
      keep    <- is.finite(p_px)
      sums    <- rowsum(p_px[keep], g[keep])
      cell_ix <- as.integer(rownames(sums))
      p_cell[cell_ix] <- sums[, 1] / tabulate(g[keep], nbins = length(cells))[cell_ix]
    }
    
    # Classify as 11's Figure 4 primary variant (same small tolerance as 11).
    has_chg <- is.finite(p_cell)
    thr_p   <- settings$pixel_p_primary - 1e-9
    sup_inc <- has_chg & p_cell >= thr_p
    sup_dec <- has_chg & (1 - p_cell) >= thr_p
    in_pool <- has_chg & (pmax(P$A2, P$A3, na.rm = TRUE) > POOL_POBS_THRESHOLD) %in% TRUE
    
    # Expected count per survey: posterior-median count with the omitted iid
    # variance integrated back in (lognormal mean), the same terms PObs integrates.
    cnt <- lapply(lam, function(l) l * exp(s$sigma2_omit / 2))
    
    # D1: model expected count at the cells holding the shared-square point
    # counts; cells averaged within each square, then squares averaged, which
    # mirrors the square-weighted raw means.
    shared_A2 <- NA_real_
    shared_A3 <- NA_real_
    sc <- shared_cells[[p]]
    if (nrow(sc) > 0) {
      sc_cnt  <- lapply(rasts, function(r) terra::extract(r, sc$cell)[, 1] * exp(s$sigma2_omit / 2))
      sq_mean <- function(v) mean(tapply(v, sc$square_id, mean, na.rm = TRUE), na.rm = TRUE)
      shared_A2 <- sq_mean(sc_cnt$A2)
      shared_A3 <- sq_mean(sc_cnt$A3)
    }
    
    # Regional change from the hexagon draws, on the expected-count scale.
    rc <- region_change(hd = hd, region = region_hex[[p]], scale = exp(s$sigma2_omit / 2))
    
    tibble(point         = POINTS$name[p],
           sp_file       = s$sp_file,
           !!!rc,
           count_A2      = mean(cnt$A2, na.rm = TRUE),
           count_A3      = mean(cnt$A3, na.rm = TRUE),
           model_shared_A2 = shared_A2,
           model_shared_A3 = shared_A3,
           pobs_A2       = mean(P$A2, na.rm = TRUE),
           pobs_A3       = mean(P$A3, na.rm = TRUE),
           dpobs         = mean(P$A3 - P$A2, na.rm = TRUE),
           p1_A2         = mean(P1$A2, na.rm = TRUE),
           p1_A3         = mean(P1$A3, na.rm = TRUE),
           abund_ratio   = mean(lam$A3, na.rm = TRUE) / mean(lam$A2, na.rm = TRUE),
           p_increase    = mean(p_cell, na.rm = TRUE),
           share_supported_increase = mean(sup_inc),
           share_supported_decrease = mean(sup_dec),
           fig4_pool     = mean(in_pool),
           fig4_increase = mean(in_pool & sup_inc),
           fig4_decrease = mean(in_pool & sup_dec))
  }))
  
  # D5: model probability that ONE standardised 5-min point count detects this
  # species (p1, as 11's PObs with n = 1), at each nearby survey's cell, using
  # the raster for that survey's atlas. Cells outside the species' predictions
  # contribute 0. Computed before the data-used check, so species without an
  # analysed-data file still count towards the all-species model total.
  sv_p1 <- vector("list", nrow(POINTS))
  for (p in seq_len(nrow(POINTS))) {
    p1 <- numeric(length(sv_cell[[p]]))
    for (a in names(ATLAS_RAST)) {
      rows <- which(near_surveys[[p]]$Atlas == a & !is.na(sv_cell[[p]]))
      if (length(rows) == 0) next
      lam_sv   <- terra::extract(rasts[[ATLAS_RAST[[a]]]], sv_cell[[p]][rows])[, 1]
      p1[rows] <- pobs(lam_sv, n = 1L)
    }
    p1[!is.finite(p1)] <- 0
    sv_p1[[p]]        <- p1
    sv_model_all[[p]] <- sv_model_all[[p]] + p1
  }
  
  du_path <- file.path(DATA_USED_DIR, paste0(s$sp_file, "_1km.rds"))
  if (!file.exists(du_path)) next
  counts <- readRDS(du_path)$survey_counts[, c("survey_id", "count")]
  
  data_rows[[k]] <- bind_rows(lapply(seq_len(nrow(POINTS)), function(p) {
    near_counts <- near_surveys[[p]] |> inner_join(counts, by = "survey_id")
    bind_rows(
      summarise_counts(near_counts),
      # ARU recordings inside the point counts' temporal window (an extra "survey type").
      summarise_counts(near_counts |> filter(aru_pc_window) |> mutate(Survey_Type = ARU_WINDOW_TYPE)),
      # D1: point counts in squares point-counted in both atlases (another extra "survey type").
      summarise_counts(near_counts |> filter(shared_pc) |> mutate(Survey_Type = SHARED_PC_TYPE))
    ) |>
      mutate(point = POINTS$name[p], sp_file = s$sp_file)
  }))
  
  for (p in seq_len(nrow(POINTS))) {
    m <- match(near_surveys[[p]]$survey_id, counts$survey_id)
    sv_used[[p]] <- sv_used[[p]] | !is.na(m)
    sv_det[[p]]  <- sv_det[[p]] + (!is.na(m) & counts$count[m] > 0)   # FALSE & NA is FALSE
    sv_n_analysed[[p]] <- sv_n_analysed[[p]] + !is.na(m)
    sv_model_an[[p]]   <- sv_model_an[[p]] + ifelse(is.na(m), 0, sv_p1[[p]])
  }
  
  if (k %% 50 == 0) message("  ", k, "/", nrow(sp_info))
}

# Surveys outside the prediction grid have no model value.
for (p in seq_len(nrow(POINTS))) {
  sv_model_all[[p]][is.na(sv_cell[[p]])] <- NA_real_
  sv_model_an[[p]][is.na(sv_cell[[p]])]  <- NA_real_
}

pix_tab  <- bind_rows(pix_rows) |>
  # direction: regional change from the hexagon draws (used in all figures).
  # direction_pixel: majority of cells individually supported (as 11's Figure 4).
  mutate(direction = factor(case_when(region_p_up_fold   >= REGION_P ~ "Supported increase",
                                      region_p_down_fold >= REGION_P ~ "Supported decrease",
                                      TRUE                           ~ "Not supported"),
                            levels = names(DIR_COLOURS)),
         direction_pixel = factor(case_when(share_supported_increase >= 0.5 ~ "Supported increase",
                                            share_supported_decrease >= 0.5 ~ "Supported decrease",
                                            TRUE                            ~ "Not supported"),
                                  levels = names(DIR_COLOURS)))
data_tab <- bind_rows(data_rows) |>
  mutate(rate       = n_det / n_surveys,
         mean_count = total_count / n_surveys,
         group = paste(unname(ATLAS_LABELS[Atlas]), Survey_Type))

# Raw means for one survey type at one point, one column per atlas:
# <prefix>_OBBA2 / _OBBA3 (mean count, by RAW_MEAN_METHOD) and <prefix>_det_OBBA2 /
# _OBBA3 (detections). With apply_minimums, a mean is NA unless its atlas x survey
# type meets MIN_GROUP_SURVEYS and MIN_GROUP_SQUARES; detection counts are kept
# either way, so figures can still show how little data there was (D4).
raw_col <- if (RAW_MEAN_METHOD == "square") "mean_count_square" else "mean_count"
raw_by_atlas <- function(nm, type, prefix, apply_minimums = TRUE) {
  d   <- data_tab[data_tab$point == nm & data_tab$Survey_Type == type, ]
  sp  <- unique(data_tab$sp_file)
  out <- tibble(sp_file = sp)
  for (a in names(ATLAS_LABELS)) {
    da    <- d[d$Atlas == a, ]
    i     <- match(sp, da$sp_file)
    means <- da[[raw_col]][i]
    if (apply_minimums) {
      enough <- (da$n_surveys[i] >= MIN_GROUP_SURVEYS & da$n_squares[i] >= MIN_GROUP_SQUARES) %in% TRUE
      means[!enough] <- NA_real_
    }
    out[[paste0(prefix, "_", a)]]     <- means
    out[[paste0(prefix, "_det_", a)]] <- da$n_det[i]
  }
  out
}

# 11's community rasters at the points (same cells as the species values).
pool_rast <- terra::rast(inputs[["POOL_CHANGE"]])[[c("primary_pool", "primary_increase", "primary_decrease")]]
fig_vals <- lapply(seq_len(nrow(POINTS)), function(p) {
  cells <- unique(model_px[[p]]$cell)
  list(richness_change = mean(terra::extract(terra::rast(inputs[["RICHNESS_CHANGE"]]), cells)[, 1], na.rm = TRUE),
       fig4 = colMeans(terra::extract(pool_rast, cells), na.rm = TRUE))
})

# ============================================================
# Outputs per point
# ============================================================
# For each point: tables, a console summary, the drivers figure (maps plus
# model vs raw change for the TOP_N species), then the reviewer diagnostics
# D1-D3 (D4 is built into the drivers and shared-squares figures).

for (p in seq_len(nrow(POINTS))) {
  
  nm <- POINTS$name[p]
  
  data_wide <- data_tab |>
    filter(point == nm) |>
    select(sp_file, group, n_surveys, n_squares, n_det, rate, mean_count, mean_count_square) |>
    pivot_wider(names_from = group, values_from = c(rate, mean_count, mean_count_square, n_det,
                                                    n_surveys, n_squares), names_sep = ": ")
  
  sp_pt <- sp_info |>
    select(sp_english, sp_file, prov_change_class, prov_pct_change_median) |>
    left_join(pix_tab |> filter(point == nm) |> select(-point), by = "sp_file") |>
    relocate(prov_change_class, prov_pct_change_median, .after = last_col()) |>
    left_join(data_wide, by = "sp_file") |>
    arrange(dpobs)
  
  utils::write.csv(sp_pt, file.path(OUT_DIR, paste0("turnover_", nm, "_species.csv")), row.names = FALSE)
  
  effort <- near_surveys[[p]] |>
    group_by(Atlas = unname(ATLAS_LABELS[Atlas]), Survey_Type) |>
    summarise(n_surveys      = n(),
              n_squares      = n_distinct(square_id),
              median_minutes = stats::median(Survey_Duration_Minutes),
              .groups = "drop")
  utils::write.csv(effort, file.path(OUT_DIR, paste0("turnover_", nm, "_effort.csv")), row.names = FALSE)
  
  # ==========================================================
  # D5: species detected per survey, observed vs model
  # ==========================================================
  # Tests: does the number of species detected per point count change between
  # atlases as the model predicts, and is any observed change an artefact of
  # how surveys enter each species' analysed data?
  # Why: the observed number of species per point count is the most direct raw
  # summary of community change, but three things can distort it:
  #   1. Unequal effort among squares: a pooled mean over surveys is dominated
  #      by heavily surveyed squares, so we also report the mean of per-square
  #      means (each square weighted equally), as for the species-level means.
  #   2. Which species a survey counts for: a detection is counted only if the
  #      survey is in that species' analysed data (data-used file), which can
  #      exclude surveys, e.g. outside the species' safe dates. If Atlas 3
  #      surveys were analysed for fewer species, they detect fewer species for
  #      reasons unrelated to birds. n_analysed (species each survey was
  #      analysed for) shows whether that happened.
  #   3. What the model expects: the model's expected number of species per
  #      standardised point count at each survey's location (sum over species
  #      of p1). model_an sums only the species the survey was analysed for,
  #      matching the observed count; model_all sums all modelled species.
  # Reading it: observed change close to model_an change = the data and the
  # model agree at these surveys. model_an change differing from model_all
  # change = the analysed-species filter alone shifts the comparison. Observed
  # counts are not adjusted for time, date, duration or observer; the model is a
  # standardised 5-min point count, so model columns are given for point counts
  # only.
  sv <- near_surveys[[p]] |>
    mutate(n_species  = sv_det[[p]],
           n_analysed = sv_n_analysed[[p]],
           model_an   = ifelse(Survey_Type == "Point_Count", sv_model_an[[p]],  NA_real_),
           model_all  = ifelse(Survey_Type == "Point_Count", sv_model_all[[p]], NA_real_)) |>
    filter(sv_used[[p]])
  sv_sub <- bind_rows(
    sv |> mutate(subset = "All squares"),
    sv |> filter(shared_pc) |> mutate(subset = "Squares point-counted in both atlases")
  ) |>
    mutate(Atlas = unname(ATLAS_LABELS[Atlas]))
  
  # Per-square means (surveys without a square are left out of these only).
  per_square <- sv_sub |>
    filter(!is.na(square_id)) |>
    group_by(subset, Atlas, Survey_Type, square_id) |>
    summarise(n_surveys  = n(),
              n_species  = mean(n_species),
              n_analysed = mean(n_analysed),
              model_an   = mean(model_an),
              model_all  = mean(model_all),
              .groups = "drop")
  sq_weighted <- per_square |>
    group_by(subset, Atlas, Survey_Type) |>
    summarise(mean_species_sq  = mean(n_species),
              model_an_sq      = mean(model_an),
              model_all_sq     = mean(model_all),
              .groups = "drop")
  
  per_survey <- sv_sub |>
    group_by(subset, Atlas, Survey_Type) |>
    summarise(n_surveys        = n(),
              n_squares        = n_distinct(square_id),
              mean_species     = mean(n_species),                       # pooled over surveys
              se               = stats::sd(n_species) / sqrt(n()),
              n_analysed_mean  = mean(n_analysed),
              n_analysed_min   = min(n_analysed),
              n_analysed_max   = max(n_analysed),
              model_an         = mean(model_an),
              model_all        = mean(model_all),
              .groups = "drop") |>
    left_join(sq_weighted, by = c("subset", "Atlas", "Survey_Type"))
  utils::write.csv(per_survey, file.path(OUT_DIR, paste0("turnover_", nm, "_species_per_survey.csv")),
                   row.names = FALSE)
  utils::write.csv(per_square |> filter(Survey_Type == "Point_Count"),
                   file.path(OUT_DIR, paste0("turnover_", nm, "_species_per_pc_squares.csv")), row.names = FALSE)
  
  # Atlas 2 -> Atlas 3 change in species per point count, observed vs model.
  pc_change <- per_survey |>
    filter(Survey_Type == "Point_Count") |>
    select(subset, Atlas, n_surveys, n_squares, mean_species, mean_species_sq,
           n_analysed_mean, model_an_sq, model_all_sq) |>
    tidyr::pivot_wider(names_from = Atlas, values_from = -c(subset, Atlas), names_sep = "_")
  # An atlas with no point counts in a subset gives missing columns; add them as NA.
  for (col in c("n_surveys", "n_squares", "mean_species", "mean_species_sq", "n_analysed_mean",
                "model_an_sq", "model_all_sq")) {
    for (a in c("Atlas 2", "Atlas 3")) {
      if (!paste0(col, "_", a) %in% names(pc_change)) pc_change[[paste0(col, "_", a)]] <- NA_real_
    }
  }
  pc_change <- pc_change |>
    mutate(observed_pct_pooled = 100 * (`mean_species_Atlas 3`    / `mean_species_Atlas 2`    - 1),
           observed_pct_sq     = 100 * (`mean_species_sq_Atlas 3` / `mean_species_sq_Atlas 2` - 1),
           model_an_pct_sq     = 100 * (`model_an_sq_Atlas 3`     / `model_an_sq_Atlas 2`     - 1),
           model_all_pct_sq    = 100 * (`model_all_sq_Atlas 3`    / `model_all_sq_Atlas 2`    - 1),
           n_analysed_pct      = 100 * (`n_analysed_mean_Atlas 3` / `n_analysed_mean_Atlas 2` - 1))
  utils::write.csv(pc_change, file.path(OUT_DIR, paste0("turnover_", nm, "_species_per_pc_change.csv")),
                   row.names = FALSE)
  
  # ---- Console summary ----
  # Checks that the species values reproduce 11's rasters at the point
  # (richness change and Figure 4 counts), then lists the largest changes.
  fv <- fig_vals[[p]]
  message(sprintf("\n==== %s (%s) ====", nm, model_lab))
  message(sprintf("Expected richness: Atlas 2 %.1f, Atlas 3 %.1f; change %+.1f (11's Figure 1 raster: %+.1f)",
                  sum(sp_pt$pobs_A2, na.rm = TRUE), sum(sp_pt$pobs_A3, na.rm = TRUE),
                  sum(sp_pt$dpobs, na.rm = TRUE), fv$richness_change))
  message(sprintf(paste0("Figure 4 (pool PObs > %.2f, P >= %.2f): pool %.1f, increased %.1f, decreased %.1f, ",
                         "net %+.1f (11's raster: %.1f, %.1f, %.1f, %+.1f)"),
                  POOL_POBS_THRESHOLD, settings$pixel_p_primary,
                  sum(sp_pt$fig4_pool, na.rm = TRUE),
                  sum(sp_pt$fig4_increase, na.rm = TRUE), sum(sp_pt$fig4_decrease, na.rm = TRUE),
                  sum(sp_pt$fig4_increase - sp_pt$fig4_decrease, na.rm = TRUE),
                  fv$fig4[["primary_pool"]], fv$fig4[["primary_increase"]], fv$fig4[["primary_decrease"]],
                  fv$fig4[["primary_increase"]] - fv$fig4[["primary_decrease"]]))
  message("Surveys within ", RADIUS_KM, " km:")
  print(as.data.frame(effort), row.names = FALSE)
  message("Species detected per survey (modelled species only; observed naive, model standardised):")
  print(as.data.frame(per_survey |> select(subset, Atlas, Survey_Type, n_surveys, n_squares, mean_species,
                                           mean_species_sq, n_analysed_mean, model_an_sq, model_all_sq)),
        digits = 3, row.names = FALSE)
  message("D5 species per point count, Atlas 2 -> Atlas 3 (sq = mean of per-square means):")
  for (i in seq_len(nrow(pc_change))) {
    x <- pc_change[i, ]
    message(sprintf("  %s", x$subset))
    message(sprintf("    observed  pooled %.2f -> %.2f (%+.0f%%); per square %.2f -> %.2f (%+.0f%%)",
                    x$`mean_species_Atlas 2`, x$`mean_species_Atlas 3`, x$observed_pct_pooled,
                    x$`mean_species_sq_Atlas 2`, x$`mean_species_sq_Atlas 3`, x$observed_pct_sq))
    message(sprintf("    model, analysed species  %.2f -> %.2f (%+.0f%%);  model, all species  %.2f -> %.2f (%+.0f%%)",
                    x$`model_an_sq_Atlas 2`, x$`model_an_sq_Atlas 3`, x$model_an_pct_sq,
                    x$`model_all_sq_Atlas 2`, x$`model_all_sq_Atlas 3`, x$model_all_pct_sq))
    message(sprintf("    species each point count was analysed for  %.1f -> %.1f (%+.0f%%)",
                    x$`n_analysed_mean_Atlas 2`, x$`n_analysed_mean_Atlas 3`, x$n_analysed_pct))
  }
  show_cols <- c("sp_english", "dpobs", "region_ratio_q50", "region_ratio_low", "region_ratio_high",
                 "direction", "direction_pixel", "prov_change_class")
  n_no_draws <- sum(is.na(sp_pt$region_ratio_q50))
  if (n_no_draws > 0) message(n_no_draws, " species lack hexagon draws for the region: ",
                              "their arrows use pixel plug-in values and they count as 'Not supported'.")
  message(sprintf("Regional change (%s): %d species supported increases, %d supported decreases (P >= %.2f, >= %.2f-fold).",
                  region_lab, sum(sp_pt$direction == "Supported increase"),
                  sum(sp_pt$direction == "Supported decrease"), REGION_P, REGION_FOLD))
  message("Largest declines in PObs:")
  print(as.data.frame(sp_pt |> filter(dpobs < 0) |> select(all_of(show_cols)) |> head(15)),
        digits = 2, row.names = FALSE)
  message("Largest gains in PObs:")
  print(as.data.frame(sp_pt |> filter(dpobs > 0) |> arrange(desc(dpobs)) |>
                        select(all_of(show_cols)) |> head(15)),
        digits = 2, row.names = FALSE)
  
  # ---- Maps of survey locations ----
  # Shows: where each atlas surveyed, by survey type, province-wide and around
  # the point. Why: in much of the north Atlas 2 is point counts and Atlas 3 is
  # mostly ARUs in different places, so reviewers need to see how much the two
  # atlases' data overlap before reading the change estimates.
  focal_pt  <- sf::st_sfc(sf::st_point(pts_xy[p, 1:2]), crs = rast_crs)
  focal_buf <- sf::st_buffer(focal_pt, dist = RADIUS_KM * u)
  zoom_x    <- pts_xy[p, 1] + c(-1, 1) * MAP_ZOOM_KM * u
  zoom_y    <- pts_xy[p, 2] + c(-1, 1) * MAP_ZOOM_KM * u
  
  # Survey-type legend drawn inside the last map's bottom-right corner, so it
  # takes no space from the panels ("inside" needs ggplot2 >= 3.5.0).
  legend_inside <- if (utils::packageVersion("ggplot2") >= "3.5.0") {
    theme(legend.position = "inside", legend.position.inside = c(0.98, 0.02))
  } else {
    theme(legend.position = c(0.98, 0.02))
  }
  
  survey_map <- function(atlas, zoom, show_legend) {
    d <- map_pts[map_pts$Atlas == atlas, ]
    if (zoom) d <- d[d$x >= zoom_x[1] & d$x <= zoom_x[2] & d$y >= zoom_y[1] & d$y <= zoom_y[2], ]
    ggplot() +
      geom_sf(data = bcr_map, fill = "grey97", colour = "grey75", linewidth = 0.2) +
      geom_sf(data = prov_outline, fill = NA, colour = "grey40", linewidth = 0.3) +
      geom_point(data = d, aes(x = x, y = y, colour = Survey_Type),
                 size = if (zoom) 0.8 else 0.15, alpha = 0.7, stroke = 0, na.rm = TRUE) +
      geom_sf(data = focal_buf, fill = NA, colour = "white", linewidth = 1.6) +   # halo
      geom_sf(data = focal_buf, fill = NA, colour = "black", linewidth = 0.7) +
      geom_sf(data = focal_pt, shape = 3, size = 2, colour = "black") +
      scale_colour_manual(values = SURVEY_COLOURS, limits = names(SURVEY_COLOURS), name = NULL) +
      guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
      coord_sf(crs = rast_crs, xlim = if (zoom) zoom_x else NULL, ylim = if (zoom) zoom_y else NULL,
               expand = FALSE, datum = NA) +
      labs(title = if (zoom) sprintf("%s, %d-km circle", atlas, RADIUS_KM) else paste(atlas, "(province)"),
           x = NULL, y = NULL) +
      theme_void(base_size = 11) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5),
            panel.border = if (zoom) element_rect(colour = "grey40", fill = NA) else element_blank(),
            legend.position = "none",
            legend.justification = c(1, 0),
            legend.background = element_rect(fill = "white", colour = "grey60")) +
      (if (show_legend) legend_inside else theme())
  }
  map_row <- survey_map(atlas = "Atlas 2", zoom = FALSE, show_legend = FALSE) |
    survey_map(atlas = "Atlas 3", zoom = FALSE, show_legend = FALSE) |
    survey_map(atlas = "Atlas 2", zoom = TRUE,  show_legend = FALSE) |
    survey_map(atlas = "Atlas 3", zoom = TRUE,  show_legend = TRUE)
  
  # ---- Drivers figure ----
  # Shows: the TOP_N species with the largest modelled change in PObs at the
  # point (which drive 11's community maps), model change on the left and the
  # raw data near the point on the right, row by row.
  top <- sp_pt |>
    filter(is.finite(dpobs)) |>
    slice_max(abs(dpobs), n = TOP_N, with_ties = FALSE) |>
    arrange(dpobs) |>
    mutate(sp_english = factor(sp_english, levels = rev(sp_english)))   # largest decline on top
  
  # Extra room below the bottom species for the right panel's text box; the same
  # expansion on both panels keeps their rows aligned.
  Y_EXPAND <- expansion(add = c(2.6, 0.6))
  
  # Left-panel metric (MODEL_METRIC). "count": posterior medians of the regional
  # mean expected count from the hexagon draws (pixel plug-in means where a
  # species lacks draws). "pobs": pixel PObs (plug-in; no regional draws).
  if (MODEL_METRIC == "count") {
    top$model_A2 <- ifelse(is.na(top$region_A2), top$count_A2, top$region_A2)
    top$model_A3 <- ifelse(is.na(top$region_A3), top$count_A3, top$region_A3)
    model_xlab   <- sprintf("Expected count per survey (standardised), regional posterior median")
  } else {
    top$model_A2 <- top$pobs_A2
    top$model_A3 <- top$pobs_A3
    model_xlab   <- sprintf("P(\u2265 1 detection in %d counts), %s", settings$pobs_n_surveys, model_lab)
  }
  
  # Regional change and its interval, printed at the right of each row, with
  # a triangle where the change is supported: e.g. "\u25bc -42% [-60, -20]".
  pct <- function(r) sprintf("%+.0f", 100 * (r - 1))
  top$change_lab <- ifelse(is.na(top$region_ratio_q50), "no draws",
                           paste0(ifelse(top$direction == "Supported increase", "\u25b2 ",
                                         ifelse(top$direction == "Supported decrease", "\u25bc ", "")),
                                  pct(top$region_ratio_q50), "% [", pct(top$region_ratio_low), ", ",
                                  pct(top$region_ratio_high), "]"))
  
  m_max <- if (MODEL_METRIC == "count") max(c(top$model_A2, top$model_A3), na.rm = TRUE) else 1
  if (!is.finite(m_max) || m_max <= 0) m_max <- 1
  m_trans  <- if (MODEL_METRIC == "count") count_trans else scales::identity_trans()
  m_breaks <- if (MODEL_METRIC == "count" && COUNT_AXIS == "pseudo_log") COUNT_BREAKS else waiver()
  m_lim    <- c(0, beyond(x = m_max, prop = 0.42, limits = c(0, m_max), tr = m_trans))   # room for change labels
  lab_x    <- beyond(x = m_max, prop = 0.03, limits = c(0, m_max), tr = m_trans)
  
  p_model <- ggplot(top, aes(y = sp_english)) +
    geom_segment(aes(x = model_A2, xend = arrow_end_t(x2 = model_A2, x3 = model_A3, limits = m_lim, tr = m_trans),
                     yend = sp_english, linetype = direction),
                 colour = "grey40", arrow = grid::arrow(length = unit(0.15, "cm"), type = "closed"),
                 na.rm = TRUE) +
    geom_point(aes(x = model_A2, colour = "Atlas 2"), size = 2, na.rm = TRUE) +
    geom_point(aes(x = model_A3, colour = "Atlas 3"), size = 2, na.rm = TRUE) +
    geom_text(aes(x = lab_x, label = change_lab), hjust = 0, size = 2.6, colour = "grey20") +
    scale_colour_manual(values = ATLAS_COLOURS, name = NULL) +
    scale_linetype_manual(values = c("Supported increase" = "solid", "Supported decrease" = "solid",
                                     "Not supported" = "dotted"), guide = "none") +
    scale_y_discrete(expand = Y_EXPAND) +
    x_scale(tr = m_trans, limits = m_lim, breaks = m_breaks) +
    labs(x = model_xlab, y = NULL,
         title = sprintf("Model: %s", region_lab)) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top", panel.grid.major.y = element_blank())
  
  # Right panel: raw data within RADIUS_KM for the same species rows, drawn like
  # the left. Point counts as arrows; windowed ARU recordings as "x" (Atlas 3
  # only in practice). D4: detections printed per row, thin means hollow.
  pc <- top |>
    select(sp_file, sp_english) |>                                          # every top species keeps its row
    left_join(raw_by_atlas(nm = nm, type = "Point_Count",   prefix = "PC"),  by = "sp_file") |>
    left_join(raw_by_atlas(nm = nm, type = ARU_WINDOW_TYPE, prefix = "ARU"), by = "sp_file")
  d_raw <- tibble(sp_english = pc$sp_english,
                  a2 = pc$PC_OBBA2, a3 = pc$PC_OBBA3,
                  det_a2 = pc$PC_det_OBBA2, det_a3 = pc$PC_det_OBBA3,
                  aru2 = pc$ARU_OBBA2, aru3 = pc$ARU_OBBA3,
                  det_lab = paste0(det_label(d2 = pc$PC_det_OBBA2, d3 = pc$PC_det_OBBA3),
                                   ifelse(is.na(pc$ARU_det_OBBA3), "",
                                          paste0("  x ", format(pc$ARU_det_OBBA3, big.mark = ",")))))
  
  # Point counts and squares behind this panel: point counts within RADIUS_KM
  # analysed for at least one species (from the per-survey summary above).
  pc_n <- per_survey |> filter(subset == "All squares", Survey_Type == "Point_Count")
  pc_line <- function(a) {
    i <- which(pc_n$Atlas == a)
    n_surveys <- if (length(i) == 0) 0 else pc_n$n_surveys[i]
    n_squares <- if (length(i) == 0) 0 else pc_n$n_squares[i]
    sprintf("%s: %s point counts in %s squares", a,
            format(n_surveys, big.mark = ","), format(n_squares, big.mark = ","))
  }
  aru_n <- sv |>
    filter(aru_pc_window) |>
    group_by(Atlas = unname(ATLAS_LABELS[Atlas])) |>
    summarise(n_surveys = n(), n_squares = n_distinct(square_id), .groups = "drop")
  aru_lines <- sprintf("%s: %s ARU %d-min recordings (x) in %s squares", aru_n$Atlas,
                       format(aru_n$n_surveys, big.mark = ","), ARU_WINDOW_DURATION,
                       format(aru_n$n_squares, big.mark = ","))
  pc_label <- paste(c(pc_line("Atlas 2"), pc_line("Atlas 3"), aru_lines), collapse = "\n")
  
  p_data <- arrow_panel(d = d_raw, with_detections = TRUE, x_values = c(d_raw$aru2, d_raw$aru3),
                        y_expand = Y_EXPAND,
                        title = sprintf("Point counts and ARU recordings within %d km (raw)", RADIUS_KM),
                        xlab = if (RAW_MEAN_METHOD == "square") "Mean count per survey (raw; mean of square means)" else
                          "Mean count per survey (raw; pooled over surveys)") +
    geom_point(aes(x = aru2, colour = "Atlas 2", shape = ARU_LAB), size = 2.5, stroke = 1, na.rm = TRUE) +
    geom_point(aes(x = aru3, colour = "Atlas 3", shape = ARU_LAB), size = 2.5, stroke = 1, na.rm = TRUE) +
    annotate(geom = "label", x = Inf, y = -Inf, label = pc_label, hjust = 1, vjust = 0,
             size = 3, fill = "white", label.r = unit(0, "pt"))
  
  drivers_h <- 2 + 0.25 * nrow(top)
  fig <- (map_row / (p_model | p_data) + plot_layout(heights = c(MAP_HEIGHT_IN, drivers_h))) +
    plot_annotation(
      title    = sprintf("Species driving predicted change at %s (%.2f, %.2f)", nm, POINTS$lon[p], POINTS$lat[p]),
      subtitle = sprintf(paste0("Top %d species by |change in PObs|, sorted from largest decline. Arrows run from ",
                                "Atlas 2 to Atlas 3.\nLeft labels: regional change from the 25-km hexagon draws, ",
                                "posterior median [%.0f%% interval]. \u25b2 / \u25bc and solid arrows = supported ",
                                "change, P(change \u2265 %+.0f%% / \u2264 %.0f%%) \u2265 %.2f; dotted = not supported."),
                         nrow(top), 100 * REGION_CI, 100 * (REGION_FOLD - 1), -100 * (1 - 1 / REGION_FOLD), REGION_P),
      caption  = paste0("Right: raw mean count per survey among surveys analysed for each species, ",
                        if (RAW_MEAN_METHOD == "square") "averaging per-square means (each atlas square weighted equally)"
                        else "pooled over surveys",
                        " (unadjusted for time, date, duration,\nhabitat or observer), shown where an atlas has \u2265 ",
                        MIN_GROUP_SURVEYS, " such surveys in \u2265 ", MIN_GROUP_SQUARES, " squares. ",
                        "The box gives point counts within the ",
                        "radius analysed for at least one species, and the squares they fall in.\n",
                        sprintf(paste0("Numbers at right: point-count detections, Atlas 2 / Atlas 3 (x = detections in ",
                                       "windowed ARU recordings); hollow points rest on < %d detections.\n"),
                                FEW_DETECTIONS),
                        sprintf(paste0("x = mean count per %d-min ARU recording within the point counts' window ",
                                       "(central %.0f%% of the Atlas 2 + 3 point counts near the point): %s."),
                                ARU_WINDOW_DURATION, 100 * diff(WINDOW_PROBS),
                                paste(sprintf("%s %.1f to %.1f", WINDOW_COLS, pc_window[[p]]["low", ],
                                              pc_window[[p]]["high", ]), collapse = "; "))),
      theme    = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  ggsave(filename = file.path(OUT_DIR, paste0("turnover_", nm, "_drivers.png")), plot = fig,
         width = 14, height = MAP_HEIGHT_IN + drivers_h, dpi = FIG_DPI, bg = "white")
  
  # ==========================================================
  # D1: like-for-like comparison in shared squares
  # ==========================================================
  # Tests: does the model's change at the places surveyed by point counts in
  # BOTH atlases match the raw point counts there?
  # Why: these squares compare the same places with the same survey type, so
  # the comparison does not depend on the model's ARU-vs-point-count adjustment
  # or on differences in where each atlas surveyed. Model and raw values are
  # both averaged over the same squares, each square weighted equally. A
  # species whose circle-wide model change differs from its shared-square
  # change owes that change to other data (ARUs, other squares) or smoothing.
  # Note: the model is standardised to reference survey conditions and the raw
  # means are not, so compare directions and relative sizes, not levels (see D3
  # for timing).
  sh_n <- per_survey |> filter(subset == "Squares point-counted in both atlases")
  if (length(shared_sq[[p]]) == 0 || nrow(sh_n) == 0) {
    message("D1: no squares point-counted in both atlases near ", nm, "; shared-squares figure skipped.")
  } else {
    sh <- top |>
      select(sp_file, sp_english, model_shared_A2, model_shared_A3) |>
      left_join(raw_by_atlas(nm = nm, type = SHARED_PC_TYPE, prefix = "SH", apply_minimums = FALSE),
                by = "sp_file")
    d_sh_model <- tibble(sp_english = sh$sp_english, a2 = sh$model_shared_A2, a3 = sh$model_shared_A3)
    d_sh_raw   <- tibble(sp_english = sh$sp_english, a2 = sh$SH_OBBA2, a3 = sh$SH_OBBA3,
                         det_a2 = sh$SH_det_OBBA2, det_a3 = sh$SH_det_OBBA3,
                         det_lab = det_label(d2 = sh$SH_det_OBBA2, d3 = sh$SH_det_OBBA3))
    sh_line <- function(a) {
      i <- which(sh_n$Atlas == a)
      sprintf("%s: %s point counts", a, if (length(i) == 0) "0" else format(sh_n$n_surveys[i], big.mark = ","))
    }
    p_sh_model <- arrow_panel(d = d_sh_model, title = "Model at the shared squares' point counts",
                              xlab = "Expected count per survey (standardised)")
    p_sh_raw   <- arrow_panel(d = d_sh_raw, with_detections = TRUE,
                              title = "Raw point counts in the shared squares",
                              xlab = "Mean count per point count (raw; mean of square means)") +
      annotate(geom = "label", x = Inf, y = -Inf, hjust = 1, vjust = 0, size = 3, fill = "white",
               label.r = unit(0, "pt"),
               label = paste(c(sprintf("%d shared squares", length(shared_sq[[p]])),
                               sh_line("Atlas 2"), sh_line("Atlas 3")), collapse = "\n"))
    
    fig_sh <- (p_sh_model | p_sh_raw) +
      plot_annotation(
        title    = sprintf("Like-for-like check at %s: squares point-counted in both atlases", nm),
        subtitle = paste0("Same species rows as the drivers figure. Left: model expected count at the cells ",
                          "holding these point counts. Right: raw point counts in the same squares.\n",
                          "Both sides weight each square equally. Arrows run from Atlas 2 to Atlas 3."),
        caption  = sprintf(paste0("Numbers at right: detections, Atlas 2 / Atlas 3; hollow points rest on < %d ",
                                  "detections. Raw means are unadjusted for time, date, observer or weather;\n",
                                  "the model is standardised to reference conditions, so compare directions and ",
                                  "relative sizes rather than levels."), FEW_DETECTIONS),
        theme    = theme(plot.title = element_text(size = 14, face = "bold"))
      )
    ggsave(filename = file.path(OUT_DIR, paste0("turnover_", nm, "_shared_squares.png")), plot = fig_sh,
           width = 14, height = 2.5 + 0.25 * nrow(top), dpi = FIG_DPI, bg = "white")
  }
  
  # ==========================================================
  # D2: agreement between model and raw change across ALL species
  # ==========================================================
  # Tests: across every species with enough raw detections (not only the
  # TOP_N shown in the figures), how closely does the model's log ratio
  # (Atlas 3 / Atlas 2) track the raw log ratio?
  # Why: the drivers figure is selected for large modelled change, which can
  # make agreement look better or worse than it is. Two comparisons:
  #   like-for-like : model at the shared squares' point counts vs raw point
  #                   counts in those squares (as D1)
  #   whole region  : regional model change from the hexagon draws (as the
  #                   drivers figure), with its REGION_CI interval, vs raw
  #                   point counts in all squares within RADIUS_KM. Also
  #                   reported: the share of species whose raw change lies
  #                   inside the model's interval.
  # The model is fit to these data, so agreement is not independent
  # validation; it shows whether a species' modelled change is consistent with
  # the point counts. Disagreement points to change driven by ARU data, the
  # survey-type adjustment or spatial smoothing.
  comp_lfl   <- "Like-for-like: shared squares"
  comp_whole <- "Whole region: hexagon draws vs all point-count squares"
  sh_all <- raw_by_atlas(nm = nm, type = SHARED_PC_TYPE, prefix = "SH", apply_minimums = FALSE)
  pc_all <- raw_by_atlas(nm = nm, type = "Point_Count",  prefix = "PC")
  agree <- bind_rows(
    # Like-for-like: plug-in model values at the shared cells (no draws at that
    # scale), with the same RATIO_OFFSET as the raw data.
    sp_pt |>
      transmute(sp_english, sp_file, direction,
                model_lr      = log_ratio(a2 = model_shared_A2, a3 = model_shared_A3),
                model_lr_low  = NA_real_,
                model_lr_high = NA_real_) |>
      left_join(sh_all |> transmute(sp_file, raw_A2 = SH_OBBA2, raw_A3 = SH_OBBA3,
                                    det_A2 = SH_det_OBBA2, det_A3 = SH_det_OBBA3), by = "sp_file") |>
      mutate(comparison = comp_lfl),
    # Whole region: exact posterior of the regional ratio (no offset needed).
    sp_pt |>
      transmute(sp_english, sp_file, direction,
                model_lr      = log(region_ratio_q50),
                model_lr_low  = log(region_ratio_low),
                model_lr_high = log(region_ratio_high)) |>
      left_join(pc_all |> transmute(sp_file, raw_A2 = PC_OBBA2, raw_A3 = PC_OBBA3,
                                    det_A2 = PC_det_OBBA2, det_A3 = PC_det_OBBA3), by = "sp_file") |>
      mutate(comparison = comp_whole)
  ) |>
    mutate(det_total    = coalesce(det_A2, 0L) + coalesce(det_A3, 0L),
           raw_lr       = log_ratio(a2 = raw_A2, a3 = raw_A3),
           lr_diff      = raw_lr - model_lr,
           raw_in_model_ci = raw_lr >= model_lr_low & raw_lr <= model_lr_high,
           # A zero mean in either atlas: the raw ratio is then set by RATIO_OFFSET.
           raw_zero     = (raw_A2 == 0 | raw_A3 == 0) %in% TRUE) |>
    filter(is.finite(model_lr), is.finite(raw_lr), det_total >= MIN_DETECTIONS)
  utils::write.csv(agree |> arrange(comparison, desc(abs(lr_diff))),
                   file.path(OUT_DIR, paste0("turnover_", nm, "_agreement_species.csv")), row.names = FALSE)
  
  # Summary statistics per comparison (correlations need at least 3 species).
  agree_stats <- agree |>
    group_by(comparison) |>
    summarise(n_species            = n(),
              pearson_r            = if (n() >= 3) stats::cor(model_lr, raw_lr) else NA_real_,
              spearman_rho         = if (n() >= 3) stats::cor(model_lr, raw_lr, method = "spearman") else NA_real_,
              share_same_direction = mean(sign(model_lr) == sign(raw_lr)),
              median_model_pct     = 100 * (exp(stats::median(model_lr)) - 1),
              median_raw_pct       = 100 * (exp(stats::median(raw_lr)) - 1),
              median_abs_lr_diff   = stats::median(abs(lr_diff)),
              share_raw_in_model_ci = if (all(is.na(raw_in_model_ci))) NA_real_ else mean(raw_in_model_ci, na.rm = TRUE),
              .groups = "drop")
  utils::write.csv(agree_stats, file.path(OUT_DIR, paste0("turnover_", nm, "_agreement.csv")), row.names = FALSE)
  message("D2 agreement between model and raw log ratios (species with >= ", MIN_DETECTIONS, " detections):")
  print(as.data.frame(agree_stats), digits = 3, row.names = FALSE)
  
  # ---------------------------------------------------------------
  # D2b: is the model's change the right size? Slope of raw on model
  # ---------------------------------------------------------------
  # Why: correlation says whether model and raw changes rank species alike, not
  # whether the model's changes are the right SIZE. A Bayesian model shrinks
  # its estimates, and a reviewer may read a steep cloud of points as the
  # model understating change. This block tests that directly.
  #
  # The regression: raw log ratio ~ model log ratio, ordinary least squares,
  # one point per species.
  #   - Raw ratios carry sampling error (few detections, few squares). Error in
  #     the response spreads points vertically, which makes the cloud LOOK
  #     steeper than 1:1, but it does not bias the least-squares slope. So the
  #     slope, not the visual impression, is the test.
  #   - A well-calibrated model is expected to give a slope near 1 even though
  #     it shrinks: calibrated shrinkage is exactly the amount at which the
  #     data, regressed on the estimates, follow a 1:1 line.
  # Interpretation:
  #   slope near 1 (interval includes 1)    model change is about the size the
  #                                         data support
  #   slope clearly > 1, persisting among   model compresses change
  #     well-sampled species without zeros  (over-shrinkage: look at the change-
  #                                         field prior and spatial range)
  #   slope clearly < 1                     model exaggerates change relative to
  #                                         the point counts
  #   > 1 only in "All species"             driven by zero-baseline species,
  #                                         whose raw ratio depends on RATIO_OFFSET
  #   ~1 like-for-like but > 1 whole region a data-source or spatial mismatch
  #                                         (ARUs, different squares), not shrinkage
  # Limits: the model is fit to these data (not independent), and with few
  # species the intervals are wide; "no evidence of compression" is not
  # "no compression". Subsets: all species; without raw zeros; and without
  # raw zeros at each SLOPE_DETECTION_TIERS level of raw detections.
  slope_rows <- list()
  for (comparison in c(comp_lfl, comp_whole)) {
    d_comp  <- agree[agree$comparison == comparison, ]
    subsets <- list("All species" = rep(TRUE, nrow(d_comp)),
                    "No raw zeros" = !d_comp$raw_zero)
    for (tier in SLOPE_DETECTION_TIERS) {
      subsets[[sprintf("No raw zeros, >= %d detections", tier)]] <- !d_comp$raw_zero & d_comp$det_total >= tier
    }
    for (subset_name in names(subsets)) {
      fit_vals <- slope_fit(d = d_comp[subsets[[subset_name]], ])
      slope_rows[[length(slope_rows) + 1]] <- tibble(comparison = comparison, subset = subset_name, !!!fit_vals)
    }
  }
  slopes <- bind_rows(slope_rows)
  utils::write.csv(slopes, file.path(OUT_DIR, paste0("turnover_", nm, "_slopes.csv")), row.names = FALSE)
  
  message("D2b slope of raw on model change (1 = model change the right size; > 1 = model compresses change):")
  for (comparison in unique(slopes$comparison)) {
    message("  ", comparison)
    d_sl <- slopes[slopes$comparison == comparison, ]
    for (i in seq_len(nrow(d_sl))) {
      if (is.na(d_sl$slope[i])) {
        message(sprintf("    %-34s n = %3d  (too few species)", d_sl$subset[i], d_sl$n_species[i]))
      } else {
        message(sprintf("    %-34s n = %3d  slope = %.2f [%.2f, %.2f]", d_sl$subset[i], d_sl$n_species[i],
                        d_sl$slope[i], d_sl$slope_low[i], d_sl$slope_high[i]))
      }
    }
  }
  
  if (nrow(agree) > 0) {
    # Facet labels carry the statistics, so the figure stands alone.
    # Slopes shown in the headers: all species, and without raw zeros.
    slope_txt <- function(comparison, subset) {
      x <- slopes[slopes$comparison == comparison & slopes$subset == subset, ]
      if (nrow(x) == 0 || is.na(x$slope)) return("slope n/a")
      sprintf("slope %.2f [%.2f, %.2f]", x$slope, x$slope_low, x$slope_high)
    }
    r_no_zero <- function(comparison) {
      x <- agree[agree$comparison == comparison & !agree$raw_zero, ]
      if (nrow(x) < 3) return(NA_real_)
      stats::cor(x$model_lr, x$raw_lr)
    }
    facet_lab <- agree_stats |>
      rowwise() |>
      mutate(lab = paste0(sprintf("%s\nn = %d species; same direction %.0f%%; median model %+.0f%%, raw %+.0f%%",
                                  comparison, n_species, 100 * share_same_direction,
                                  median_model_pct, median_raw_pct),
                          ifelse(is.na(share_raw_in_model_ci), "\n",
                                 sprintf("\nRaw change inside the model's %.0f%% interval: %.0f%% of species",
                                         100 * REGION_CI, 100 * share_raw_in_model_ci)),
                          sprintf("\nAll species: r = %.2f, %s", pearson_r, slope_txt(comparison, "All species")),
                          sprintf("\nWithout raw zeros (n = %d): r = %.2f, %s",
                                  sum(agree$comparison == comparison & !agree$raw_zero),
                                  r_no_zero(comparison), slope_txt(comparison, "No raw zeros")))) |>
      ungroup()
    facet_levels <- facet_lab$lab[match(c(comp_lfl, comp_whole), facet_lab$comparison)]
    agree$facet  <- factor(facet_lab$lab[match(agree$comparison, facet_lab$comparison)],
                           levels = facet_levels[!is.na(facet_levels)])
    # Label the species furthest from the 1:1 line in each comparison.
    outliers <- agree |>
      group_by(comparison) |>
      slice_max(abs(lr_diff), n = N_LABEL_OUTLIERS, with_ties = FALSE) |>
      ungroup()
    lr_lim <- 1.05 * max(abs(c(agree$model_lr, agree$raw_lr)), log(2))   # 5% margin keeps edge points
    # Interval bars are clipped to the axis range rather than widening it.
    agree_bars <- agree |>
      filter(!is.na(model_lr_low)) |>
      mutate(model_lr_low = pmax(model_lr_low, -lr_lim), model_lr_high = pmin(model_lr_high, lr_lim))
    lr_brk <- log(c(0.1, 0.25, 0.5, 1, 2, 4, 10))
    lr_brk <- lr_brk[abs(lr_brk) <= lr_lim]
    
    # Fitted least-squares line per panel (species without raw zeros), to set
    # against the dashed 1:1 line.
    fit_lines <- slopes[slopes$subset == "No raw zeros" & !is.na(slopes$slope), ]
    fit_lines$facet <- factor(facet_lab$lab[match(fit_lines$comparison, facet_lab$comparison)],
                              levels = levels(agree$facet))
    
    # Point fill: hollow where the raw ratio rests on a zero (set by RATIO_OFFSET).
    agree$raw_type <- ifelse(agree$raw_zero, "Zero in one atlas (position set by offset)", "Detected in both atlases")
    
    fig_ag <- ggplot(agree, aes(x = model_lr, y = raw_lr)) +
      geom_hline(yintercept = 0, colour = "grey70") +
      geom_vline(xintercept = 0, colour = "grey70") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey30") +
      geom_abline(data = fit_lines, aes(slope = slope, intercept = intercept), colour = "black", linewidth = 0.5) +
      geom_errorbarh(data = agree_bars, aes(xmin = model_lr_low, xmax = model_lr_high, colour = direction),
                     height = 0, alpha = 0.5, linewidth = 0.4) +
      geom_point(aes(size = det_total, colour = direction, shape = raw_type), alpha = 0.8, stroke = 0.9) +
      geom_text(data = outliers, aes(label = sp_english), size = 2.6, vjust = -0.9, check_overlap = TRUE) +
      scale_colour_manual(values = DIR_COLOURS,
                          name = sprintf("Regional model change\n(hexagon draws; P \u2265 %.2f,\n\u2265 %+.0f%% / \u2264 %.0f%%)",
                                         REGION_P, 100 * (REGION_FOLD - 1), -100 * (1 - 1 / REGION_FOLD))) +
      scale_size_area(name = "Raw detections\n(both atlases)", max_size = 5) +
      scale_shape_manual(values = c("Detected in both atlases" = 16,
                                    "Zero in one atlas (position set by offset)" = 1),
                         name = "Raw data") +
      scale_x_continuous(limits = c(-lr_lim, lr_lim), breaks = lr_brk, labels = pct_lab) +
      scale_y_continuous(limits = c(-lr_lim, lr_lim), breaks = lr_brk, labels = pct_lab) +
      coord_equal() +
      facet_wrap(~ facet) +
      labs(x = "Model change, Atlas 2 to Atlas 3", y = "Raw point-count change, Atlas 2 to Atlas 3",
           title = sprintf("Model vs raw change across all species at %s", nm),
           subtitle = sprintf(paste0("Each point is a species with \u2265 %d raw detections. Dashed line = perfect ",
                                     "agreement; points above it rose more (or fell less) in the raw data than in ",
                                     "the model.\nSolid line: least-squares fit of raw on model change, species ",
                                     "without raw zeros. Labelled: the %d species furthest from the 1:1 line.\n",
                                     "Horizontal bars (right): the model's %.0f%% interval for regional change."),
                              MIN_DETECTIONS, N_LABEL_OUTLIERS, 100 * REGION_CI),
           caption = sprintf(paste0("Log ratios of mean counts per survey; %.2f is added to both atlases' raw means ",
                                    "(and to the like-for-like model values) so zeros are defined. Colour shows ",
                                    "regional support in both panels.\nSlope near 1 = model change about the size the ",
                                    "data support; > 1 = model compresses change (95%% confidence intervals). ",
                                    "The model is fit to these data,\nso agreement shows consistency, not independent ",
                                    "validation."), RATIO_OFFSET)) +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold"), strip.text = element_text(size = 9))
    ggsave(filename = file.path(OUT_DIR, paste0("turnover_", nm, "_agreement.png")), plot = fig_ag,
           width = 14, height = 8, dpi = FIG_DPI, bg = "white")
  }
  
  # ==========================================================
  # D3: survey timing by atlas and survey type
  # ==========================================================
  # Tests: were Atlas 3 surveys done at a different time of season or of day
  # than Atlas 2 surveys, within the radius and in the shared squares?
  # Why: raw means are not adjusted for timing, and singing rates change
  # through the season and the morning, so part of a raw change can be timing
  # rather than abundance. The shaded band is the window used to pick ARU
  # recordings comparable with the point counts.
  timing <- bind_rows(
    sv |> mutate(subset = sprintf("All analysed surveys within %d km", RADIUS_KM)),
    sv |> filter(shared_pc) |> mutate(subset = "Shared squares (point counts)")
  ) |>
    mutate(group = case_when(Survey_Type == "Point_Count" ~ paste(unname(ATLAS_LABELS[Atlas]), "point counts"),
                             grepl("^ARU", Survey_Type)   ~ paste(unname(ATLAS_LABELS[Atlas]), "ARU recordings"),
                             TRUE                         ~ NA_character_),
           atlas_lab = unname(ATLAS_LABELS[Atlas])) |>
    filter(!is.na(group)) |>
    select(subset, group, atlas_lab, all_of(WINDOW_COLS)) |>
    pivot_longer(cols = all_of(WINDOW_COLS), names_to = "variable", values_to = "value") |>
    filter(is.finite(value)) |>
    mutate(variable = unname(WINDOW_LABELS[variable]),
           group    = factor(group, levels = rev(c("Atlas 2 point counts", "Atlas 2 ARU recordings",
                                                   "Atlas 3 point counts", "Atlas 3 ARU recordings"))))
  
  timing_tab <- timing |>
    group_by(subset, group, variable) |>
    summarise(n_surveys = n(),
              q05 = stats::quantile(value, 0.05), q25 = stats::quantile(value, 0.25),
              median = stats::median(value),
              q75 = stats::quantile(value, 0.75), q95 = stats::quantile(value, 0.95),
              .groups = "drop")
  utils::write.csv(timing_tab, file.path(OUT_DIR, paste0("turnover_", nm, "_timing.csv")), row.names = FALSE)
  message("D3 survey timing (medians):")
  print(as.data.frame(timing_tab |> select(subset, group, variable, n_surveys, median)), digits = 3, row.names = FALSE)
  
  if (nrow(timing) > 0) {
    win_df <- tibble(variable = unname(WINDOW_LABELS[WINDOW_COLS]),
                     low      = pc_window[[p]]["low", ],
                     high     = pc_window[[p]]["high", ]) |>
      filter(is.finite(low), is.finite(high))
    fig_tm <- ggplot(timing, aes(x = value, y = group, colour = atlas_lab)) +
      geom_rect(data = win_df, aes(xmin = low, xmax = high, ymin = -Inf, ymax = Inf),
                inherit.aes = FALSE, fill = "grey92") +
      geom_boxplot(outlier.size = 0.6, width = 0.6) +
      facet_grid(subset ~ variable, scales = "free", space = "free_y") +
      scale_colour_manual(values = ATLAS_COLOURS, name = NULL) +
      labs(x = NULL, y = NULL,
           title = sprintf("Survey timing at %s", nm),
           subtitle = paste0("Analysed surveys only. Shaded: the point counts' window used to select comparable ARU ",
                             "recordings (central ", round(100 * diff(WINDOW_PROBS)), "% of Atlas 2 + 3 point counts).")) +
      theme_bw(base_size = 11) +
      theme(legend.position = "none", plot.title = element_text(face = "bold"))
    ggsave(filename = file.path(OUT_DIR, paste0("turnover_", nm, "_timing.png")), plot = fig_tm,
           width = 12, height = 5.5, dpi = FIG_DPI, bg = "white")
  }
  
  # ---- D5 figure: species per point count, observed vs model ----
  # Top: species detected per point count by atlas. Grey points and lines are
  # individual squares (observed; lines join a square's two atlases where it
  # has point counts in both), so reviewers can see whether one square drives
  # the change. Black: observed, square-weighted. Coloured: model, square-
  # weighted, for the analysed species and for all modelled species.
  # Bottom: number of species each point count was analysed for; a drop from
  # Atlas 2 to Atlas 3 would lower observed counts for non-biological reasons.
  pc_sq  <- per_square |> filter(Survey_Type == "Point_Count")
  pc_sum <- per_survey |> filter(Survey_Type == "Point_Count")
  if (nrow(pc_sum) > 0) {
    SOURCE_COLOURS <- c("Observed (square-weighted)"         = "black",
                        "Model, species analysed per survey" = "#CC79A7",
                        "Model, all modelled species"        = "#009E73")
    SOURCE_SHAPES  <- c("Observed (square-weighted)"         = 18,
                        "Model, species analysed per survey" = 17,
                        "Model, all modelled species"        = 15)
    pc_long <- bind_rows(
      pc_sum |> transmute(subset, Atlas, source = names(SOURCE_COLOURS)[1], value = mean_species_sq),
      pc_sum |> transmute(subset, Atlas, source = names(SOURCE_COLOURS)[2], value = model_an_sq),
      pc_sum |> transmute(subset, Atlas, source = names(SOURCE_COLOURS)[3], value = model_all_sq)
    ) |>
      mutate(source = factor(source, levels = names(SOURCE_COLOURS)))
    dodge <- position_dodge(width = 0.35)
    
    p_rich <- ggplot() +
      geom_line(data = pc_sq, aes(x = Atlas, y = n_species, group = square_id), colour = "grey80") +
      geom_point(data = pc_sq, aes(x = Atlas, y = n_species, size = n_surveys), colour = "grey60",
                 shape = 1) +
      geom_line(data = pc_long, aes(x = Atlas, y = value, colour = source, group = source),
                position = dodge, linewidth = 0.8) +
      geom_point(data = pc_long, aes(x = Atlas, y = value, colour = source, shape = source),
                 position = dodge, size = 3.2) +
      facet_wrap(~ subset) +
      scale_colour_manual(values = SOURCE_COLOURS, name = NULL) +
      scale_shape_manual(values = SOURCE_SHAPES, name = NULL) +
      scale_size_area(name = "Point counts\nin square", max_size = 4) +
      expand_limits(y = 0) +
      labs(x = NULL, y = "Species detected per point count",
           title = "Species detected per point count: observed vs model") +
      theme_bw(base_size = 11) +
      theme(legend.position = "right", plot.title = element_text(face = "bold"))
    
    p_an <- ggplot(pc_sum, aes(x = Atlas, y = n_analysed_mean)) +
      geom_linerange(aes(ymin = n_analysed_min, ymax = n_analysed_max), colour = "grey60", linewidth = 1) +
      geom_point(size = 3) +
      geom_line(aes(group = 1)) +
      facet_wrap(~ subset) +
      expand_limits(y = 0) +
      labs(x = NULL, y = "Species analysed\nper point count",
           title = "Number of species each point count was analysed for (mean; bar = range)") +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 11))
    
    fig_rich <- (p_rich / p_an + plot_layout(heights = c(2.2, 1))) +
      plot_annotation(
        title    = sprintf("Community check at %s: species per point count", nm),
        subtitle = paste0("Modelled species only. Grey: individual squares (observed). Observed counts are naive ",
                          "(unadjusted for time, date, duration or observer);\nthe model is a standardised 5-min ",
                          "point count at each survey's location, summed over species (P of \u2265 1 detection)."),
        caption  = paste0("'Species analysed per survey': a survey counts towards a species only if it is in that ",
                          "species' analysed data. Fewer species analysed per\nsurvey in Atlas 3 would lower its ",
                          "observed counts for reasons unrelated to birds; the pink model line uses the same species ",
                          "set as the observed count."),
        theme    = theme(plot.title = element_text(size = 14, face = "bold"))
      )
    ggsave(filename = file.path(OUT_DIR, paste0("turnover_", nm, "_species_per_pc.png")), plot = fig_rich,
           width = 12, height = 8, dpi = FIG_DPI, bg = "white")
  }
}

message("\n12_community_turnover_at_points.R complete -> ", OUT_DIR)