# ============================================================
# 09_generate_model_products.R
#
# Purpose
#   Convert fitted-model prediction RDS files into map-ready rasters and
#   multi-page species assessment PDFs.
#
# Running order and dependencies
#   Runs AFTER 08_export_regional_change_estimates.R. The province-wide and
#   per-BCR posterior change numbers shown on page 3 are READ from the table 08
#   saves (regional_change_estimates_<model_name>.rds) rather than recomputed
#   here, so the figures and the exported table cannot disagree. The per-hex
#   change map (summarize_hex_draw_change) is computed here -- it is cheap and
#   species-specific.
#
#   Every species-invariant spatial operation is done ONCE before the species
#   loop (see "Precompute species-independent objects"): pixel -> raster-cell
#   index, region masks, hex -> cell lookups, survey -> region / hex membership
#   and clipped region boundaries. Inside the loop rasters are filled, clipped
#   and summarised by vector indexing.
#   Set cfg$report_timings = TRUE to log time per step for each species.
#
#   If available, 07b's paired_summaries.rds is used for the page-3 repeated-
#   survey comparison (change estimates, locations map, and summary table). It is
#   optional: species with full-model predictions still receive the full PDF when
#   no paired result exists.
#
# Main outputs
#   - model_output/rasters_<model_name>/<species>_{a2,a3,chg}.tif
#   - model_output/figures_<model_name>/<species>.pdf   (14 pages; honeycomb-only
#     species get 6 pages, one per region)
#
# Notes
#   - Both the medians and the CI widths of the mapped surfaces are taken from
#     the water-corrected summaries.
#   - CRS choices are inherited from the saved analysis-ready data; objects are
#     aligned to the prediction-grid CRS only.
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(terra)
  library(patchwork)
  library(magick)
  library(viridis)
})

source(here::here("R", "00_config_paths.R"))
source(file.path(paths$functions, "model_product_utils.R"))
source(file.path(paths$functions, "spatial_utils.R"))      # make_hex_grid()
source(file.path(paths$functions, "inla_model_utils.R"))   # get_var()


# ============================================================
# CONFIG
# ============================================================

model_name <- "PC_ARU_CL_SS"

# Zero predictions (pixel summaries and hex draws) in biological regions where
# the species has no safe dates, using each prediction file's own sp_safe_dates.
# Keep the same setting in 08 and 09 so page-3 numbers match the maps.
TRIM_TO_SAFE_DATES <- TRUE

cfg <- list(
  # Rasterization
  plot_res = 1001,            # raster resolution in map units (m in EPSG:3978)
  
  # Assessment hex grids / scaling
  n_hexagons    = 2000,
  rast_max_q    = 0.99,
  chg_max_q     = 0.99,
  count_max_q   = 0.99,
  max_surveys_q = 0.80,
  
  # Colour-scale transform for the relative-abundance (raster) panel. One of
  # "identity", "sqrt", "log", etc. "sqrt" compresses a long upper tail and
  # reveals low-abundance structure.
  relabund_transform = "sqrt",
  
  # Change summaries. ci_level MUST match the value 08 used to build the change
  # table (checked below), so the page-3 intervals match the exported table.
  # ci_level is the TWO-SIDED credible interval used for the page-3 interval
  # displays (map title "Overall change = ... [low to high]" and the comparison-
  # plot error bars).
  ci_level   = 0.90,
  # ONE-SIDED confidence for the hex change map's "changed by at least this much"
  # colour bands (classify_min_supported_change). A hex is coloured e.g. ">=1.1x"
  # only if P(change >= 1.1x) >= map_change_confidence, i.e. the (1 - value)
  # posterior quantile of the log-ratio is above the band threshold. This is a
  # DIRECTIONAL claim and is deliberately separate from the two-sided ci_level
  # above: reusing ci_level's 5th/95th percentiles as one-sided bounds would make
  # each band a 95%-confident claim while the legend reads "90% confident".
  map_change_confidence = 0.90,
  min_sq_det = 10,            # hide BCR-level change unless detected in >= this many squares
  
  # Relative-abundance display threshold. First find the density cutoff below
  # which the low-density tail contains at most (1 - relabund_coverage) = 1% of
  # the predicted population. Then CAP that cutoff at relabund_absent_limit =
  # 1/1000 expected birds per standardized 5-min count. Thus the final threshold
  # is min(T_99%, 1/1000): the 99% rule can make the map more inclusive, but can
  # never mask values above 1/1000.
  relabund_coverage = 0.99,
  relabund_absent_limit = 1/1000,
  
  
  # Probability-of-observation surface ----------------------------------------
  pobs_marginalize_terms = c("kappa_diff","site_iid"),  # omitted iid terms integrated back in
  pobs_gh_nodes          = 41L,              # Gaussian integration nodes
  pobs_floor             = 0,                # white out P(Obs) below this (0 = show all)
  pobs_n_surveys = 25L,
  
  # PDF page rendering. pdf_engine = "png" renders 300-dpi PNGs then binds them
  # with magick (keeps the PNGs). "cairo" writes the PDF directly (raster layers
  # embedded as images), which is much faster.
  page_width  = 20,
  page_height = 10,
  page_res    = 300,
  pdf_engine  = "png",
  
  # Log seconds per processing step for each species.
  report_timings = TRUE
)

# Survey columns used by the products. all_surveys is slimmed to these (plus
# geometry) once, so per-species joins and filters don't copy every covariate.
SURVEY_COLS <- c("survey_id", "Atlas", "BCR", "Survey_Type",
                 "Survey_Duration_Minutes", "square_id")

# Bird Conservation Regions used for the sub-region panels.
# Codes come from bcr_regions$BCR (character); each BCR gets its own pair of
# assessment pages. make_region_geom() stops if a code listed here is absent
# from bcr_regions, catching a stale data_ready_for_analysis.rds.
# List order sets the order of the region pages in the output PDF.
region_defs <- list(
  province = list(bcr = NULL,          suffix = ""),   # whole study area
  bcr12    = list(bcr = c("12"),       suffix = "Temperate Mixed"),
  bcr13    = list(bcr = c("13"),       suffix = "Lower Great Lakes"),
  bcr76    = list(bcr = c("76"),       suffix = "Boreal Shield West"),
  bcr77    = list(bcr = c("77"),       suffix = "Boreal Shield East"),
  hb       = list(bcr = c("74"),       suffix = "Hudson Plains")
)

# Display order (south -> north) for the regional change-comparison plot and the
# repeated-survey summary table on page 3. Defined once so the two panels cannot
# disagree, and checked against bcr_regions$BCR_Label below. These are BCR_Label
# values, not the shorter page-title suffixes used in region_defs above.
region_label_order <- c(
  "Lower Great Lakes / St. Lawrence Plain",
  "Temperate Mixed",
  "Boreal Shield West",
  "Boreal Shield East",
  "Hudson Plains"
)

# ------------------------------------------------------------
# Paths and directories
# ------------------------------------------------------------

in_data  <- file.path(paths$data_clean, "birds", "data_ready_for_analysis.rds")
pred_dir <- file.path(paths$model_output, paste0("predictions_", model_name))
rast_dir <- file.path(paths$model_output, paste0("rasters_", model_name))
fig_dir  <- file.path(paths$model_output, paste0("figures_", model_name))
png_dir  <- file.path(fig_dir, "PNGs")
data_used_dir <- file.path(paths$model_output, paste0("data_used_", model_name))
out_dir  <- paths$model_output

paired_analysis_path <- file.path(out_dir, "paired_summaries", "paired_summaries.rds")

# Regional change table produced by script 08 (read, not recomputed).
change_table_path <- file.path(
  out_dir, paste0("change_estimates_", model_name),
  paste0("regional_change_estimates_", model_name, ".rds")
)

for (dir_i in c(rast_dir, fig_dir, png_dir)) {
  dir.create(dir_i, recursive = TRUE, showWarnings = FALSE)
}

stopifnot(file.exists(in_data))

if (!file.exists(change_table_path)) {
  stop("Regional change table not found at:\n  ", change_table_path,
       "\nRun 08_export_regional_change_estimates.R first (same model_name).",
       call. = FALSE)
}

# ============================================================
# Load analysis-ready data and the change table
# ============================================================

dat                 <- readRDS(in_data)
regional_change_all <- readRDS(change_table_path)

# The exported table carries the ci_level 08 used. If it differs from cfg here,
# the page-3 intervals would silently disagree with the numbers plotted, so stop.
if ("ci_level" %in% names(regional_change_all)) {
  tbl_ci <- unique(regional_change_all$ci_level)
  if (length(tbl_ci) == 1 && !isTRUE(all.equal(tbl_ci, cfg$ci_level))) {
    stop("cfg$ci_level (", cfg$ci_level, ") does not match the change table's ci_level (",
         tbl_ci, "). Re-run 08 or align cfg$ci_level.", call. = FALSE)
  }
}

# Geometry only. Prediction rows in 07 align with the FULL grids, so the grids
# are not na.omit()'d (doing so would misalign values if any row were dropped).
grid2          <- dat$grid_OBBA2["pixel_id"]
grid3          <- dat$grid_OBBA3["pixel_id"]
study_boundary <- dat$study_boundary
hex_grid       <- dat$hex_grid_25km

# Align all map layers to the Atlas 2 prediction-grid CRS.
crs_use        <- sf::st_crs(grid2)
study_boundary <- sf::st_transform(study_boundary, crs_use)
grid2          <- sf::st_transform(grid2, crs_use)
grid3          <- sf::st_transform(grid3, crs_use)

# ------------------------------------------------------------
# BCRs (panel boundaries; also used to mask/annotate change)
# ------------------------------------------------------------

bcr_regions <- dat$bcr_sf %>%
  sf::st_transform(crs_use) %>%
  sf::st_make_valid() %>%
  dplyr::arrange(BCR)

stopifnot(all(c("BCR", "BCR_Label") %in% names(bcr_regions)))

bcr_regions$BCR       <- as.character(bcr_regions$BCR)
bcr_regions$BCR_Label <- as.character(bcr_regions$BCR_Label)

# Guard against a data_ready_for_analysis.rds whose BCR set differs from region_label_order.
if (!setequal(region_label_order, bcr_regions$BCR_Label)) {
  stop(
    "region_label_order does not match the BCR labels in dat$bcr_sf.\n",
    "  region_label_order: ", paste(region_label_order, collapse = ", "), "\n",
    "  bcr_sf BCR_Label:   ", paste(sort(unique(bcr_regions$BCR_Label)), collapse = ", "), "\n",
    "Re-run 06_filter_and_finalize_surveys.R so bcr_sf matches region_label_order.",
    call. = FALSE
  )
}

# Clipped copy reused by the repeated-survey locations map (built once).
bcr_regions_clipped <- suppressWarnings(
  bcr_regions %>% sf::st_intersection(study_boundary)
)

# Shared survey table. Per-species data_used files store only survey_id + count;
# geometry and attributes live here once and are rejoined by survey_id.
all_surveys <- dat$all_surveys %>%
  dplyr::select(dplyr::all_of(SURVEY_COLS)) %>%
  sf::st_transform(crs_use)

stopifnot("BCR" %in% names(all_surveys))
all_surveys$BCR <- as.character(all_surveys$BCR)

if (any(is.na(all_surveys$BCR))) {
  message(sum(is.na(all_surveys$BCR)),
          " surveys have no BCR assigned; they will be dropped from regional ",
          "change masking.")
}

# Pixel -> Biol_Region lookups for safe-date trimming (straight off 06's grids).
pixel_region_OBBA2 <- grid_region_lookup(dat$grid_OBBA2)
pixel_region_OBBA3 <- grid_region_lookup(dat$grid_OBBA3)

# ============================================================
# Precompute species-independent objects (built ONCE)
# ============================================================

# Relative-abundance colour palette (deuteranomaly-friendly green ramp).
colpal_relabund <- grDevices::colorRampPalette(c(
  "#FBF7E2", "#FCF8D0", "#EEF7C2",
  "#CEF2B0",
  "#94E5A0", "#51C987", "#18A065", "#008C59",
  "#007F53", "#006344"
))(100)

# Probability-of-observation palette.
colpal_pobs <- c("white", viridis::mako(6, direction = -1))

message("Precomputing species-independent spatial indices ...")

# Pixel -> raster-cell index, built once so each species' rasters are filled by
# vector indexing (rasterize_from_index()).
rast_index2 <- build_raster_index(grid2, res = cfg$plot_res)
rast_index3 <- if (identical(sf::st_geometry(grid2), sf::st_geometry(grid3))) {
  rast_index2
} else {
  build_raster_index(grid3, res = cfg$plot_res)
}
rm(grid2, grid3)

# Region geometries, assessment hex grids, clipped BCR boundaries, and the
# raster / survey indices that let every per-region panel skip vector geometry.
region_specs <- lapply(names(region_defs), function(nm) {
  d    <- region_defs[[nm]]
  geom <- make_region_geom(d$bcr, bcr_regions, study_boundary)
  hex  <- make_hex_grid(study_area = geom, n_hexagons = cfg$n_hexagons)
  list(
    geom   = geom,
    suffix = d$suffix,
    hex    = hex,
    bounds = suppressWarnings(sf::st_intersection(bcr_regions, geom)),
    index  = c(
      list(in_col = paste0("in_", nm), hex_col = paste0("hex_", nm)),
      build_region_raster_index(template = rast_index2$template,
                                region_geom = geom, hex_grid = hex)
    )
  )
})
names(region_specs) <- names(region_defs)

# Survey -> region membership and region hex id, once for all species.
all_surveys <- add_region_survey_index(all_surveys, region_specs)

# ============================================================
# Per-species products
#
# Iterate EVERY species that has a saved honeycomb record (data_used). Species
# with fitted-model predictions get the full multi-page assessment PDF. Paired-
# survey results are optional: when absent, page 3 shows the full-model change
# summaries only. Species without model predictions get honeycomb-only PDFs.
# ============================================================

dat_used_files <- list.files(data_used_dir, pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(dat_used_files) > 0)
message("Species honeycomb records found: ", length(dat_used_files))

# Paired-survey change estimates (raw shared points live here; the page-3 panel
# needs them, so they are read directly rather than from 08's export).
if (file.exists(paired_analysis_path)) {
  paired_summaries <- readRDS(paired_analysis_path)
} else {
  message("Paired summaries not found at: ", paired_analysis_path,
          "; full model products will still be generated, with paired panels omitted.")
  paired_summaries <- list()
}

for (i in seq_along(dat_used_files)) {
  
  dat_used <- readRDS(dat_used_files[i])
  
  # sp_file is taken from the filename so it always matches the names script 07
  # used for the prediction and raster files.
  sp_file    <- basename(dat_used_files[i]) |> stringr::str_remove("_1km\\.rds$")
  sp_english <- dat_used$sp_english
  
  message("\nGenerating products for: ", sp_english)
  
  pdf_path  <- file.path(fig_dir, paste0(sp_file, ".pdf"))
  pred_path <- file.path(pred_dir, paste0(sp_file, "_1km.rds"))
  
  # ----------------------------------------------------------
  # Decide what, if anything, needs to be generated
  #   1. PDF already exists            -> skip.
  #   2. Modelable but no predictions  -> skip (don't make a honeycomb-only PDF).
  #   3. No predictions / not modelable -> honeycomb-only PDF.
  #   4. Predictions exist              -> full package (paired summary optional).
  # ----------------------------------------------------------
  
   if (file.exists(pdf_path)) {
     message("Skipping ", sp_english, "; PDF already exists at: ", pdf_path)
     next
   }
  
  if (!"is_modelable" %in% names(dat_used)) {
    dat_used$is_modelable <- file.exists(pred_path)
  }
  
  is_modelable <- isTRUE(dat_used$is_modelable)
  has_model    <- file.exists(pred_path)
  paired       <- paired_summaries[[sp_english]]
  has_paired   <- !is.null(paired)
  
  if (is_modelable && !has_model) {
    message("Skipping ", sp_english,
            "; species is modelable, but model predictions were not found.")
    next
  }
  
  if (has_model && !has_paired) {
    message("  Paired-survey summaries not found; generating full model products ",
            "and omitting paired-analysis elements from page 3.")
  }
  
  # Reconstruct this species' survey sf (geometry + attributes + count) by
  # rejoining the lightweight count table to the shared all_surveys table.
  sp_surveys <- load_species_surveys(dat_used, all_surveys)
  
  # ----------------------------------------------------------
  # Honeycomb-only branch for species without model predictions
  # ----------------------------------------------------------
  
  if (!has_model) {
    message("Honeycomb-only maps for ", sp_english, " (no model predictions available).")
    
    honey_prep <- prepare_honeycomb_surveys(sp_surveys)
    A2_obs     <- honey_prep$data %>% dplyr::filter(Atlas == "OBBA2")
    A3_obs     <- honey_prep$data %>% dplyr::filter(Atlas == "OBBA3")
    A2_special <- honey_prep$special %>% dplyr::filter(Atlas == "OBBA2")
    A3_special <- honey_prep$special %>% dplyr::filter(Atlas == "OBBA3")
    
    honey_size_label <- honeycomb_size_label(honey_prep$per_minute)
    page_subtitle    <- if (honey_prep$per_minute) {
      "Observed survey effort and mean counts per survey-minute (no model fitted)"
    } else {
      "Observed survey effort and mean counts (no model fitted)"
    }
    
    honey_pages <- lapply(region_specs, function(spec) {
      h <- build_period_honeycombs(
        region_geom      = spec$geom,
        region_hex       = spec$hex,
        region_bounds    = spec$bounds,
        A2_dat           = A2_obs,
        A3_dat           = A3_obs,
        A2_special       = A2_special,
        A3_special       = A3_special,
        sp_english       = sp_english,
        title_suffix     = spec$suffix,
        max_surveys_q    = cfg$max_surveys_q,
        count_max_q      = cfg$count_max_q,
        count_size_label = honey_size_label,
        region_index     = spec$index
      )
      
      page_title <- if (nzchar(spec$suffix)) paste0(sp_english, " - ", spec$suffix) else sp_english
      
      (h$A3 + margin_theme) + (h$A2 + margin_theme) +
        patchwork::plot_layout(ncol = 2) +
        patchwork::plot_annotation(
          title    = page_title,
          subtitle = page_subtitle,
          theme = ggplot2::theme(
            plot.title    = ggplot2::element_text(size = 18, face = "bold"),
            plot.subtitle = ggplot2::element_text(size = 12, colour = "grey30")
          )
        )
    })
    
    save_pdf_pages(
      plots     = honey_pages,
      pdf_path  = pdf_path,
      png_paths = file.path(png_dir, paste0(sp_file, "_honey_page",
                                            seq_along(honey_pages), ".png")),
      width     = cfg$page_width,
      height    = cfg$page_height,
      res       = cfg$page_res,
      engine    = cfg$pdf_engine
    )
    next
  }
  
  # ----------------------------------------------------------
  # Full model products
  # ----------------------------------------------------------
  
  tick  <- make_step_timer(cfg$report_timings)
  preds <- readRDS(pred_path)
  tick("read predictions")
  
  # Trim out predictions for areas with no safe dates
  # NOTE: this ensures a species cannot occur in a region where it has no safe dates
  if (TRIM_TO_SAFE_DATES) {
    preds <- trim_predictions_to_safe_dates(
      preds                 = preds,
      lookup_obba2          = pixel_region_OBBA2,
      lookup_obba3          = pixel_region_OBBA3,
      hex_no_safe_threshold = 0.5,   # zero a hex when >= this fraction of its pixels are no-safe
      trim_unclassified     = FALSE
    )
  }
  
  # Change numbers for this species come from 08's saved table.
  rc_sp <- regional_change_all %>% dplyr::filter(sp_english == !!sp_english)
  if (nrow(rc_sp) == 0) {
    message("Skipping ", sp_english,
            "; no rows in the regional change table. Re-run ",
            "08_export_regional_change_estimates.R for model ", model_name, ".")
    next
  }
  
  sp_code <- dat$species_detection_summaries$species_id[
    which(dat$species_detection_summaries$sp_english == sp_english)
  ][1]
  
  if (length(sp_code) == 0) {
    warning("No species_id found for ", sp_english, "; using NA in raster metadata.")
    sp_code <- NA_character_
  }
  

  rast_path_a2  <- file.path(rast_dir, paste0(sp_file, "_a2.tif"))
  rast_path_a3  <- file.path(rast_dir, paste0(sp_file, "_a3.tif"))
  rast_path_chg <- file.path(rast_dir, paste0(sp_file, "_chg.tif"))
  
  # ----------------------------------------------------------
  # Raster products: Atlas 2, Atlas 3, absolute change
  # ----------------------------------------------------------
  
  # Filled from the precomputed pixel -> cell index (no per-species geometry).
  # Both the median and the interval width come from the water-corrected summaries.
  # CI_90_width = 95th - 5th posterior percentile (the interval 07 summarises).
  p2   <- preds$OBBA2_Corrected_for_Water
  p3   <- preds$OBBA3_Corrected_for_Water
  pchg <- preds$abs_change_Corrected_for_Water
  meta <- c(species_name = sp_english, species_id = sp_code,
            species_filename = sp_file, model_name = model_name)
  
  r2 <- rasterize_from_index(
    index    = rast_index2,
    values   = list(mu_q50 = p2$OBBA2_q50, CI_90_width = p2$OBBA2_upper - p2$OBBA2_lower),
    metadata = c(meta, units = "Expected count (Atlas 2)")
  )
  r3 <- rasterize_from_index(
    index    = rast_index3,
    values   = list(mu_q50 = p3$OBBA3_q50, CI_90_width = p3$OBBA3_upper - p3$OBBA3_lower),
    metadata = c(meta, units = "Expected count (Atlas 3)")
  )
  rchg <- rasterize_from_index(
    index    = rast_index3,
    values   = list(chg_q50 = pchg$abs_change_q50,
                    CI_90_width = pchg$abs_change_upper - pchg$abs_change_lower),
    metadata = c(meta, units = "Change in expected counts (Atlas 2 to Atlas 3)")
  )
  tick("rasterize")
  
  # ----------------------------------------------------------
  # Relative-abundance rasters (shared absence floor + scale)
  # ----------------------------------------------------------
  
  rasters_relabund_prepared <- prepare_relative_abundance_rasters(
    Atlas2            = r2[["mu_q50"]],
    Atlas3            = r3[["mu_q50"]],
    absent_method     = "cumulative",
    coverage          = cfg$relabund_coverage,
    coverage_basis    = "each",     # retain >= coverage in BOTH periods (never erase a period's signal)
    max_absent_limit  = cfg$relabund_absent_limit,  # cap: final threshold = min(T_99%, 1/1000)
    rast_max_quantile = 0.99
  )
  
  r2_clamped <- rasters_relabund_prepared$rasters$Atlas2
  r3_clamped <- rasters_relabund_prepared$rasters$Atlas3
  
  rchg_prepared <- prepare_change_raster(
    chg_q50   = rchg[["chg_q50"]],
    ci_width  = rchg[["CI_90_width"]],
    footprint = terra::cover(r2_clamped, r3_clamped),
    max_q     = cfg$chg_max_q
  )
  
  # Save rasters. Below-floor pixels INSIDE the modelled area are written as 0
  # ("effectively absent"); NA is reserved for pixels that were never modelled
  # (outside the study area, large waterbodies). cover() fills the NAs created by
  # thresholding with 0 wherever the unthresholded surface has a value. The
  # NA-thresholded r2_clamped / r3_clamped are left untouched for the PDF maps.
  r2_to_save   <- r2
  r3_to_save   <- r3
  rchg_to_save <- rchg
  
  r2_to_save$mu_q50    <- terra::cover(x = r2_clamped, y = r2[["mu_q50"]] * 0)
  r3_to_save$mu_q50    <- terra::cover(x = r3_clamped, y = r3[["mu_q50"]] * 0)
  rchg_to_save$chg_q50 <- terra::cover(x = rchg_prepared$chg_q50, y = rchg[["chg_q50"]] * 0)
  
  # NA footprint must be unchanged by thresholding (only modelled-pixel NAs -> 0)
  stopifnot(
    identical(is.na(terra::values(r2_to_save$mu_q50, mat = FALSE)),
              is.na(terra::values(r2[["mu_q50"]], mat = FALSE))),
    identical(is.na(terra::values(r3_to_save$mu_q50, mat = FALSE)),
              is.na(terra::values(r3[["mu_q50"]], mat = FALSE))),
    identical(is.na(terra::values(rchg_to_save$chg_q50, mat = FALSE)),
              is.na(terra::values(rchg[["chg_q50"]], mat = FALSE)))
  )
  
  absent_tags <- c(
    display_absent_limit = format(rasters_relabund_prepared$rast_absent_limit),
    display_absent_value = "0"
  )
  
  terra::metags(r2_to_save) <- absent_tags
  terra::metags(r3_to_save) <- absent_tags
  terra::metags(rchg_to_save) <- c(
    absent_tags,
    display_max_abs_change = format(rchg_prepared$zlim[2]),
    display_max_ci_width   = format(rchg_prepared$width_max),
    display_max_q          = format(cfg$chg_max_q)
  )
  
  terra::writeRaster(r2_to_save,   filename = rast_path_a2,  overwrite = TRUE)
  terra::writeRaster(r3_to_save,   filename = rast_path_a3,  overwrite = TRUE)
  terra::writeRaster(rchg_to_save, filename = rast_path_chg, overwrite = TRUE)
  tick("prepare + write rasters")
  
  # ----------------------------------------------------------
  # Hex-level change map (Page 3, left)
  # ----------------------------------------------------------
  
  # Province annotation comes from 08's table (no polygon x hex intersection here).
  prov_change <- rc_sp %>% dplyr::filter(region_type == "province")
  stopifnot(nrow(prov_change) == 1)
  
  # Per-hex change surface IS still computed here (cheap; species-specific).
  hex_change_sf <- summarize_hex_draw_change(
    hex_grid            = hex_grid,
    hex_draws           = preds$hex_draws_Corrected_for_Water,
    ci_level            = cfg$ci_level,
    classify_confidence = cfg$map_change_confidence
  ) |>
    classify_min_supported_change()
  
  Hex_Change_Map <- make_hex_abs_change_map(
    species_name      = sp_english,
    hex_change_sf     = hex_change_sf,
    region            = study_boundary,
    region_boundaries = bcr_regions_clipped,
    clip_boundaries   = FALSE,
    prov_change       = prov_change,
    map_confidence    = cfg$map_change_confidence,
    zlim              = rasters_relabund_prepared$zlim,
    water             = NULL,
    water_fill        = "gray97"
  )
  tick("hex change map")
  
  # ----------------------------------------------------------
  # Regional change estimates: full model vs paired surveys (Page 3)
  # ----------------------------------------------------------
  
  # ---- Full spatial model, per BCR (read from 08's table) ----
  regional_estimates_FullModel <- rc_sp %>%
    dplyr::filter(region_type == "BCR") %>%
    dplyr::transmute(
      Region_Name, Region_Number,
      pct_change_median, pct_change_qlow, pct_change_qhigh,
      direction
    ) %>%
    dplyr::relocate(Region_Name, Region_Number)
  
  # ---- Paired analysis (surveys within shared 100 m footprint) ----
  # Optional. If 07b has no result for this species, all full-model products are
  # still generated and only the paired-analysis elements on page 3 are omitted.
  if (has_paired) {
    # NOTE: 07b must summarize by BCR; require_bcr_key() fails loudly on
    # Biol_Region-keyed files.
    paired_change_summary <- paired$shared_change_summary %>%
      require_bcr_key("shared_change_summary", sp_english) %>%
      dplyr::rename(
        pct_change_median = pct_change_q50,
        pct_change_qlow   = pct_change_q05,
        pct_change_qhigh  = pct_change_q95
      )
    
    paired_data <- paired$shared_data %>%
      require_bcr_key("shared_data", sp_english)
  } else {
    paired_change_summary <- NULL
    paired_data           <- NULL
  }
  
  # Hide change estimates where a species was detected in too few squares. The
  # detection counts are read from 08's table, keyed by BCR; 08 counts
  # point-count and ARU detections only.
  data_availability <- rc_sp %>%
    dplyr::filter(region_type == "BCR") %>%
    dplyr::transmute(
      BCR = Region_Number,
      n_sq_det_OBBA2, n_sq_det_OBBA3,
      sufficient_data = n_sq_det_OBBA2 >= cfg$min_sq_det |
        n_sq_det_OBBA3 >= cfg$min_sq_det
    )
  
  change_estimate_cols <- c("pct_change_median", "pct_change_qlow", "pct_change_qhigh")
  
  regional_estimates_FullModel_masked <- regional_estimates_FullModel %>%
    dplyr::left_join(data_availability, by = c("Region_Number" = "BCR")) %>%
    dplyr::mutate(sufficient_data = dplyr::coalesce(sufficient_data, FALSE)) %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(change_estimate_cols),
        ~ dplyr::if_else(sufficient_data, .x, NA_real_)
      ),
      direction = dplyr::if_else(sufficient_data, direction, NA_character_)
    )
  
  if (has_paired) {
    paired_change_summary_masked <- paired_change_summary %>%
      dplyr::left_join(data_availability, by = c("BCR" = "BCR")) %>%
      dplyr::mutate(sufficient_data = dplyr::coalesce(sufficient_data, FALSE)) %>%
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(change_estimate_cols),
          ~ dplyr::if_else(sufficient_data, .x, NA_real_)
        )
      )
  } else {
    paired_change_summary_masked <- NULL
  }
  
  change_comparison_plot <- make_change_comparison_plot(
    regional_estimates_FullModel = regional_estimates_FullModel_masked,
    paired_change_summary        = paired_change_summary_masked,
    region_order = region_label_order,
    title    = "Regional population change estimates",
    subtitle = if (has_paired) {
      "Comparison of full atlas model vs repeated surveys only"
    } else {
      "Full atlas model"
    }
  )
  
  if (has_paired) {
    paired_locations_map <-
      ggplot2::ggplot() +
      ggplot2::geom_sf(data = study_boundary, colour = "black") +
      ggplot2::geom_sf(data = paired_data, size = 0.1, col = "#1FAA59") +
      ggplot2::geom_sf(data = bcr_regions_clipped, fill = "transparent",
                       colour = "black", linewidth = 0.3) +
      ggplot2::ggtitle("Repeated survey locations") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size = 13),
        plot.margin = ggplot2::margin(5, 5, 5, 5)
      )
    
    paired_data_summary <- paired_data %>%
      as.data.frame() %>%
      dplyr::group_by(BCR, Atlas) %>%
      dplyr::summarize(
        mean_count = mean(count),
        PObs       = mean(count > 0),
        n_surveys  = dplyr::n(),
        n_squares  = length(unique(square_id)),
        .groups    = "drop"
      ) %>%
      dplyr::left_join(
        as.data.frame(bcr_regions)[, c("BCR", "BCR_Label")],
        by = c("BCR" = "BCR")
      ) %>%
      dplyr::mutate(BCR_Label = factor(BCR_Label, levels = region_label_order))
    
    paired_summary_table <- make_paired_summary_table(paired_data_summary)
  }
  
  # ----------------------------------------------------------
  # Survey data prepared once, then assessed per region
  # ----------------------------------------------------------
  
  honey_prep     <- prepare_honeycomb_surveys(sp_surveys)
  A2_dat_to_plot <- honey_prep$data %>% dplyr::filter(Atlas == "OBBA2")
  A3_dat_to_plot <- honey_prep$data %>% dplyr::filter(Atlas == "OBBA3")
  A2_special     <- honey_prep$special %>% dplyr::filter(Atlas == "OBBA2")
  A3_special     <- honey_prep$special %>% dplyr::filter(Atlas == "OBBA3")
  honey_size_label <- honeycomb_size_label(honey_prep$per_minute)
  
  # One unified call per region defined in region_defs.
  assessments <- lapply(region_specs, function(spec) {
    build_period_assessments(
      region_geom       = spec$geom,
      region_hex        = spec$hex,
      r2                = r2_clamped,
      r3                = r3_clamped,
      A2_dat            = A2_dat_to_plot,
      A3_dat            = A3_dat_to_plot,
      A2_special        = A2_special,
      A3_special        = A3_special,
      region_boundaries = spec$bounds,   # pre-clipped to the region
      region_index      = spec$index,
      sp_english        = sp_english,
      title_suffix      = spec$suffix,
      colpal            = colpal_relabund,
      rast_max_q        = cfg$rast_max_q,
      count_max_q       = cfg$count_max_q,
      max_surveys_q     = cfg$max_surveys_q,
      transform         = cfg$relabund_transform,
      count_size_label  = honey_size_label
    )
  })
  names(assessments) <- names(region_specs)
  tick("regional assessment panels")
  
  # ----------------------------------------------------------
  # Probability-of-observation maps (Page 14)
  # ----------------------------------------------------------
  # Per-survey detection probability for one 5-min count at a random square in
  # the pixel's neighbourhood. Built directly on the single-layer median surface
  # (NOT via prepare_relative_abundance_rasters, whose density floor and 0.99
  # clamp are for abundance) so the full 0-1 colour scale is available.
  
  # Integrate only the omitted iid terms NOT already baked into q50 (none when
  # vc_terms was empty in 07), so nothing is double-counted.
  vc_present    <- preds$variance_correction$present
  applied_terms <- if (is.null(vc_present)) character(0) else names(which(vc_present))
  marg_terms    <- setdiff(cfg$pobs_marginalize_terms, applied_terms)
  
  sigma2_omit <- sum(vapply(marg_terms,
                            function(tm) get_var(hyper = preds$summary_hyperpar, name = tm),
                            numeric(1)))
  nb_size     <- nb_size_from_hyper(preds$summary_hyperpar)
  
  pobs2 <- marginal_pobs(
    r2[["mu_q50"]],
    sigma2_omit = sigma2_omit,
    size = nb_size,
    n_nodes = cfg$pobs_gh_nodes,
    n_surveys = cfg$pobs_n_surveys
  )
  
  pobs3 <- marginal_pobs(
    r3[["mu_q50"]],
    sigma2_omit = sigma2_omit,
    size = nb_size,
    n_nodes = cfg$pobs_gh_nodes,
    n_surveys = cfg$pobs_n_surveys
  )
  if (cfg$pobs_floor > 0) {
    pobs2[pobs2 < cfg$pobs_floor] <- NA
    pobs3[pobs3 < cfg$pobs_floor] <- NA
  }
  
  # Clip to the province with the precomputed mask (raster-raster).
  pobs2 <- clip_to_region_mask(pobs2, region_specs$province$index$mask)
  pobs3 <- clip_to_region_mask(pobs3, region_specs$province$index$mask)
  
  pobs_Map_Atlas2 <- make_map(
    species_name = sp_english, 
    subtitle = "Probability of ≥1 Detection Across 25 Independent 5-min Point Counts - Atlas 2",
    legend_title = "P(≥1 detection)\nacross 25\npoint counts",
    rast = pobs2, region = study_boundary, clipped = TRUE,
    water = NULL, colpal = colpal_pobs, water_fill = "white",
    transform = "identity", zlim = c(0, 1), zbreaks = seq(0, 1, length.out = 5)
  )
  
  pobs_Map_Atlas3 <- make_map(
    species_name = sp_english, 
    subtitle = "Probability of ≥1 Detection Across 25 Independent 5-min Point Counts - Atlas 3",
    legend_title = "P(≥1 detection)\nacross 25\npoint counts",
    rast = pobs3, region = study_boundary, clipped = TRUE,
    water = NULL, colpal = colpal_pobs, water_fill = "white",
    transform = "identity", zlim = c(0, 1), zbreaks = seq(0, 1, length.out = 5)
  )
  
  tick("probability-of-observation maps")
  
  # ----------------------------------------------------------
  # Compose pages (BCR 76 and 77 get separate page pairs)
  # ----------------------------------------------------------
  
  # Page 3 layout. When paired results exist, retain the comparison plot, repeated-
  # survey locations, and repeated-survey summary table. When they do not, keep
  # the full-model regional estimates and omit only the paired-analysis elements.
  if (has_paired) {
    bottom_right_panel <-
      paired_locations_map |
      paired_summary_table +
      patchwork::plot_layout(widths = c(1, 1))
    
    right_panel <-
      change_comparison_plot /
      bottom_right_panel +
      patchwork::plot_layout(heights = c(0.5, 1))
  } else {
    right_panel <- change_comparison_plot
  }
  
  population_change_page <-
    Hex_Change_Map |
    right_panel +
    patchwork::plot_layout(widths = c(1.5, 1))
  
  pobs_page <-
    pobs_Map_Atlas3 + pobs_Map_Atlas2 +
    patchwork::plot_layout(ncol = 2, widths = c(1, 1))
  
  # Page plots in final PDF order.
  page_plots <- list(
    assessments$province$A3$plot_combined,   # page 1
    assessments$province$A2$plot_combined,   # page 2
    population_change_page,                  # page 3
    assessments$bcr12$A3$plot_combined,      # page 4   (BCR 12, Temperate Mixed)
    assessments$bcr12$A2$plot_combined,      # page 5
    assessments$bcr13$A3$plot_combined,      # page 6   (BCR 13, Lower Great Lakes)
    assessments$bcr13$A2$plot_combined,      # page 7
    assessments$bcr76$A3$plot_combined,      # page 8   (BCR 76, Boreal Shield West)
    assessments$bcr76$A2$plot_combined,      # page 9
    assessments$bcr77$A3$plot_combined,      # page 10  (BCR 77, Boreal Shield East)
    assessments$bcr77$A2$plot_combined,      # page 11
    assessments$hb$A3$plot_combined,         # page 12  (BCR 74, Hudson Plains)
    assessments$hb$A2$plot_combined,         # page 13
    pobs_page                                # page 14
  )
  
  save_pdf_pages(
    plots     = page_plots,
    pdf_path  = pdf_path,
    png_paths = file.path(png_dir, paste0(sp_file, "_page", seq_along(page_plots), ".png")),
    width     = cfg$page_width,
    height    = cfg$page_height,
    res       = cfg$page_res,
    engine    = cfg$pdf_engine
  )
  tick("render + write PDF")
}

message("09_generate_model_products.R complete")
