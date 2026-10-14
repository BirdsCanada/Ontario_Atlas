# ============================================================
# 05_extract_covariates.R
#
# Purpose
#   Extract covariate values for prediction-grid cells (per atlas) and for
#   survey locations (1-km square footprints), assign each survey to an atlas
#   square, and crop the prediction grids to the exact study boundary.
#
# Output
#   data_clean/birds/analysis_data_covariates.rds
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(terra)
  library(exactextractr)
  library(here)
})

source(here::here("R", "00_config_paths.R"))
source(file.path(paths$functions, "covariate_processing_utils.R"))

birds_dir   <- file.path(paths$data_clean, "birds")
spatial_dir <- file.path(paths$data_clean, "spatial")
dir.create(birds_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Config and inputs
# ------------------------------------------------------------

crs_rast <- sf::st_crs(3978)   # extraction CRS (matches the processed rasters)

out_file <- file.path(birds_dir, "analysis_data_covariates.rds")

in_grid         <- file.path(spatial_dir, "prediction_grid_1000_m.rds")
in_surveys      <- file.path(paths$data_clean, "surveys",  "surveys_raw.rds")
in_count_matrix <- file.path(paths$data_clean, "surveys",  "count_matrix_raw.rds")
in_species      <- file.path(paths$data_clean, "metadata", "species_list.rds")
in_study_area   <- file.path(spatial_dir, "study_area.rds")

# Raw inputs from ONTARIO_ATLAS_DATA; stops here if any are missing.
inputs <- require_inputs(c("atlas_squares", "lcc_2010", "lcc_2020"))

# ------------------------------------------------------------
# Load inputs
# ------------------------------------------------------------

all_surveys  <- readRDS(in_surveys)
count_matrix <- readRDS(in_count_matrix)
all_species  <- readRDS(in_species)

# Re-assert the survey <-> count row alignment BY KEY, so the count_matrix
# saved below shares all_surveys' row order. 02 enforces this at the source;
# re-asserting here stops a stale raw file reintroducing a misalignment
# silently. (06 also re-keys by survey_id before subsetting counts.)
stopifnot(
  setequal(rownames(count_matrix), all_surveys$survey_id),
  anyDuplicated(all_surveys$survey_id) == 0
)
count_matrix <- count_matrix[all_surveys$survey_id, , drop = FALSE]
stopifnot(identical(rownames(count_matrix), all_surveys$survey_id))

grid_obj       <- readRDS(in_grid)
grid_cells     <- st_transform(grid_obj$grid_cells, crs_rast)
grid_centroids <- st_transform(grid_obj$grid_centroids, crs_rast)

study_area   <- readRDS(in_study_area)
boundary_buf <- study_area$boundary_buffer_5km %>%
  st_transform(crs_rast) %>%
  st_make_valid()

# Stable ordering field for surveys.
if (!("obs_idx" %in% names(all_surveys))) {
  all_surveys$obs_idx <- seq_len(nrow(all_surveys))
}
all_surveys   <- all_surveys %>% relocate(obs_idx)
all_surveys_m <- st_transform(all_surveys, crs_rast)

# 1-km square footprints centred on each survey (for area-based extraction).
survey_pixels <- make_square_buffer(all_surveys_m, half_width = 500)
survey_pixels <- st_set_geometry(all_surveys_m, survey_pixels)

# ------------------------------------------------------------
# Load covariate rasters
# ------------------------------------------------------------

# Rasters written by 04; stops if any are missing.
read_04_rast <- function(file) {
  read_required_rast(file.path(spatial_dir, file),
                     hint = "Run 04_process_vectors_as_rasters.R first.")
}
roads_2005   <- read_04_rast("roads_2005_buf_125m.tif")
roads_2025   <- read_04_rast("roads_2025_buf_125m.tif")
rivers_large <- read_04_rast("rivers_large_buf_500m.tif")
rivers_small <- read_04_rast("rivers_small_buf_250m.tif")
lakes_large  <- read_04_rast("lakes_large_buf_500m.tif")
lakes_small  <- read_04_rast("lakes_small_buf_250m.tif")
coastline    <- read_04_rast("coastline_buf_3km.tif")
great_lakes  <- read_04_rast("great_lakes_buf_3km.tif")
open_water   <- read_04_rast("waterbodies.tif")

# Land cover.
# NOTE: this should eventually be supplemented with additional rasters
# (e.g. SCANFI, CANLABS, Far North Land Cover).
LCC2020 <- terra::rast(inputs$lcc_2020)
LCC2010 <- terra::rast(inputs$lcc_2010)

# ------------------------------------------------------------
# Build grid + survey covariate tables (OBBA2 = 2010 LCC, OBBA3 = 2020 LCC)
# ------------------------------------------------------------

surveys2 <- survey_pixels %>% filter(Atlas == "OBBA2")
surveys3 <- survey_pixels %>% filter(Atlas == "OBBA3")

grid_OBBA2 <- grid_centroids
grid_OBBA3 <- grid_centroids

# Land-cover class fractions -> LCC_<class> columns.
LCC2010_frac <- extract_frac(LCC2010, grid_cells, boundary_buf, prefix = "LCC", drop0 = TRUE, clean_names = FALSE)
LCC2020_frac <- extract_frac(LCC2020, grid_cells, boundary_buf, prefix = "LCC", drop0 = TRUE, clean_names = FALSE)
grid_OBBA2 <- bind_cols(grid_OBBA2, LCC2010_frac)
grid_OBBA3 <- bind_cols(grid_OBBA3, LCC2020_frac)

surveys2 <- bind_cols(surveys2, extract_frac(LCC2010, surveys2, boundary_buf, prefix = "LCC", drop0 = TRUE, clean_names = FALSE))
surveys3 <- bind_cols(surveys3, extract_frac(LCC2020, surveys3, boundary_buf, prefix = "LCC", drop0 = TRUE, clean_names = FALSE))

# Roads (atlas-specific: 2005 vs 2025 networks).
grid_OBBA2$road <- extract_mean(roads_2005, grid_cells, boundary_buf)
grid_OBBA3$road <- extract_mean(roads_2025, grid_cells, boundary_buf)
surveys2$road   <- extract_mean(roads_2005, surveys2, boundary_buf)
surveys3$road   <- extract_mean(roads_2025, surveys3, boundary_buf)

# Hydro covariates are time-invariant, so OBBA3 grid values are copied from OBBA2
# (surveys are still extracted separately per atlas because their footprints differ).
grid_OBBA2$rivers_large <- extract_mean(rivers_large, grid_cells, boundary_buf)
grid_OBBA3$rivers_large <- grid_OBBA2$rivers_large
surveys2$rivers_large   <- extract_mean(rivers_large, surveys2, boundary_buf)
surveys3$rivers_large   <- extract_mean(rivers_large, surveys3, boundary_buf)

grid_OBBA2$rivers_small <- extract_mean(rivers_small, grid_cells, boundary_buf)
grid_OBBA3$rivers_small <- grid_OBBA2$rivers_small
surveys2$rivers_small   <- extract_mean(rivers_small, surveys2, boundary_buf)
surveys3$rivers_small   <- extract_mean(rivers_small, surveys3, boundary_buf)

grid_OBBA2$lakes_large  <- extract_mean(lakes_large, grid_cells, boundary_buf)
grid_OBBA3$lakes_large  <- grid_OBBA2$lakes_large
surveys2$lakes_large    <- extract_mean(lakes_large, surveys2, boundary_buf)
surveys3$lakes_large    <- extract_mean(lakes_large, surveys3, boundary_buf)

grid_OBBA2$lakes_small  <- extract_mean(lakes_small, grid_cells, boundary_buf)
grid_OBBA3$lakes_small  <- grid_OBBA2$lakes_small
surveys2$lakes_small    <- extract_mean(lakes_small, surveys2, boundary_buf)
surveys3$lakes_small    <- extract_mean(lakes_small, surveys3, boundary_buf)

grid_OBBA2$great_lakes  <- extract_mean(great_lakes, grid_cells, boundary_buf)
grid_OBBA3$great_lakes  <- grid_OBBA2$great_lakes
surveys2$great_lakes    <- extract_mean(great_lakes, surveys2, boundary_buf)
surveys3$great_lakes    <- extract_mean(great_lakes, surveys3, boundary_buf)

grid_OBBA2$coastline    <- extract_mean(coastline, grid_cells, boundary_buf)
grid_OBBA3$coastline    <- grid_OBBA2$coastline
surveys2$coastline      <- extract_mean(coastline, surveys2, boundary_buf)
surveys3$coastline      <- extract_mean(coastline, surveys3, boundary_buf)

grid_OBBA2$open_water   <- extract_mean(open_water, grid_cells, boundary_buf)
grid_OBBA3$open_water   <- grid_OBBA2$open_water
surveys2$open_water     <- extract_mean(open_water, surveys2, boundary_buf)
surveys3$open_water     <- extract_mean(open_water, surveys3, boundary_buf)

# Keep geometry last.
grid_OBBA2 <- grid_OBBA2 %>% relocate(geometry, .after = last_col())
grid_OBBA3 <- grid_OBBA3 %>% relocate(geometry, .after = last_col())

# Recombine surveys into a single point sf, ordered by obs_idx.
all_surveys_with_covs <- bind_rows(surveys2, surveys3) %>%
  arrange(obs_idx) %>%
  st_centroid()

# ------------------------------------------------------------
# Assign each survey to an atlas square
# ------------------------------------------------------------

atlas_squares <- st_read(inputs$atlas_squares, quiet = TRUE) %>%
  st_make_valid() %>%
  st_transform(st_crs(all_surveys_with_covs)) %>%
  select(square_id)

hit       <- st_within(all_surveys_with_covs, atlas_squares)
square_id <- rep(NA_character_, nrow(all_surveys_with_covs))
has_hit   <- lengths(hit) > 0L
if (any(has_hit)) {
  idx <- vapply(hit[has_hit], function(x) x[1], integer(1))
  square_id[has_hit] <- atlas_squares$square_id[idx]
}
all_surveys_with_covs$square_id <- square_id

# ------------------------------------------------------------
# Crop prediction grids to the exact study boundary
# ------------------------------------------------------------

bndry      <- study_area$boundary %>% st_transform(st_crs(grid_OBBA2))
grid_OBBA2 <- st_filter(grid_OBBA2, bndry, .predicate = st_intersects)
grid_OBBA3 <- st_filter(grid_OBBA3, bndry, .predicate = st_intersects)

# ------------------------------------------------------------
# Save (surveys returned to km-unit CRS for downstream scripts)
# ------------------------------------------------------------

all_surveys_with_covs <- st_transform(all_surveys_with_covs, st_crs(all_surveys)) %>%
  relocate(geometry, .after = last_col())

saveRDS(
  list(
    all_surveys_with_covs = all_surveys_with_covs,
    count_matrix          = count_matrix,
    grid_OBBA2            = grid_OBBA2,   # EPSG:3978
    grid_OBBA3            = grid_OBBA3,   # EPSG:3978
    grid_centroids        = grid_centroids,
    grid_cells            = grid_cells,
    all_species           = all_species,
    boundary              = study_area$boundary,
    boundary_buffer_5km   = study_area$boundary_buffer_5km,
    date_created          = Sys.time()
  ),
  file = out_file
)

message("05_extract_covariates.R complete.")
