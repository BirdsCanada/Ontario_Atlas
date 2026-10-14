# ============================================================
# 00_config_paths.R
#
# Purpose
#   Centralize the paths used by the Distribution_Models workflow. Every
#   numbered script sources this file first and then uses `paths` and
#   require_inputs().
#
# Project root
#   here::here() resolves to this workflow folder when Distribution_Models.Rproj
#   is open. Outputs are written to data_clean/ inside it (git-ignored).
#
# Raw data
#   Raw inputs live OUTSIDE the repository, in the folder named by the
#   ONTARIO_ATLAS_DATA environment variable (set it in .Renviron; see README).
#   Every raw input is listed once in RAW_INPUTS below, relative to that folder.
#   Scripts get full paths from require_inputs(), which stops with a list of
#   anything missing, so a wrong data location fails immediately and loudly.
# ============================================================

suppressPackageStartupMessages({
  library(here)
})

# ------------------------------------------------------------
# Raw inputs, relative to ONTARIO_ATLAS_DATA
# ------------------------------------------------------------

DATA_ENV_VAR <- "ONTARIO_ATLAS_DATA"

RAW_INPUTS <- list(
  bcr             = file.path("Spatial", "National", "BCR_2026", "bcr2026_statprov.shp"),
  atlas_squares   = file.path("Spatial", "National", "AtlasSquares", "NationalSquares_FINAL.shp"),
  lcc_2010        = file.path("Spatial", "National", "LandCoverCanada2010", "landcover-2010-classification.tif"),
  lcc_2020        = file.path("Spatial", "National", "LandCoverCanada2020", "landcover-2020-classification.tif"),
  roads_2005      = file.path("Spatial", "Provincial", "Ontario", "RoadNetwork", "2005", "grnf035r05a_e.shp"),
  roads_2025      = file.path("Spatial", "Provincial", "Ontario", "RoadNetwork", "2025", "lrnf000r25a_e.shp"),
  ohn_waterbody   = file.path("Spatial", "Provincial", "Ontario", "Ontario_Hydro_Network_(OHN)_-_Waterbody",
                              "Ontario_Hydro_Network_(OHN)_-_Waterbody.shp"),
  ohn_watercourse = file.path("Spatial", "Provincial", "Ontario", "Ontario_Hydro_Network_(OHN)_-_Watercourse",
                              "Ontario_Hydro_Network_(OHN)_-_Watercourse.shp"),
  biol_regions    = file.path("Spatial", "Provincial", "Ontario", "Ontario_Atlas_biol_regions",
                              "Atlas_biol_regions.shp"),
  safe_dates      = file.path("Bird_Data_Raw", "OBBA3",
                              "OnAtlasSafeDates_update_for_NatureCounts_2026-07-20.xlsx")
)

# Shapefiles also need these sidecar files to be readable.
SHAPEFILE_SIDECARS <- c("shx", "dbf")

# ------------------------------------------------------------
# Project paths
# ------------------------------------------------------------

paths <- list(
  functions    = here::here("R", "functions"),
  data_clean   = here::here("data_clean"),
  model_output = here::here("data_clean", "model_output"),
  data         = Sys.getenv(DATA_ENV_VAR, unset = "")
)

if (!dir.exists(paths$functions)) {
  stop("R/functions not found under ", here::here(), ".\n",
       "Open Distribution_Models.Rproj (not the repository-level project) so ",
       "here::here() resolves to the workflow folder.", call. = FALSE)
}

# ------------------------------------------------------------
# require_inputs(): full paths to raw inputs, or a loud failure
# ------------------------------------------------------------
# names  Keys of RAW_INPUTS to check and return; NULL (default) checks all.
# Returns a named list of absolute paths. Stops, listing every missing file,
# if ONTARIO_ATLAS_DATA is unset or does not exist, or if any requested input
# (or a shapefile's sidecar) is absent.
require_inputs <- function(names = NULL) {

  if (is.null(names)) names <- names(RAW_INPUTS)
  unknown <- setdiff(names, names(RAW_INPUTS))
  if (length(unknown) > 0) {
    stop("Unknown raw input(s): ", paste(unknown, collapse = ", "),
         ". Add them to RAW_INPUTS in 00_config_paths.R.", call. = FALSE)
  }

  if (!nzchar(paths$data)) {
    stop(DATA_ENV_VAR, " is not set. Add a line such as\n  ", DATA_ENV_VAR,
         "=C:/path/to/Data\nto your .Renviron and restart R (see README).", call. = FALSE)
  }
  data_dir <- normalizePath(paths$data, winslash = "/", mustWork = FALSE)
  if (!dir.exists(data_dir)) {
    stop(DATA_ENV_VAR, " points to a folder that does not exist:\n  ", data_dir,
         call. = FALSE)
  }

  full <- lapply(RAW_INPUTS[names], function(p) file.path(data_dir, p))

  needed <- unlist(lapply(full, function(p) {
    if (grepl("\\.shp$", p)) c(p, paste0(sub("shp$", "", p), SHAPEFILE_SIDECARS)) else p
  }), use.names = FALSE)
  missing <- needed[!file.exists(needed)]

  if (length(missing) > 0) {
    stop(length(missing), " required raw input file(s) not found under ",
         DATA_ENV_VAR, " = ", data_dir, ":\n  ",
         paste(substring(missing, nchar(data_dir) + 2), collapse = "\n  "),
         "\nSee README for the expected data structure.", call. = FALSE)
  }

  full
}
