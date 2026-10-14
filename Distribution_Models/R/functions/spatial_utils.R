# ============================================================
# spatial_utils.R
#
# Shared spatial helpers for the Distribution_Models workflow.
#
# Contents
#   get_aea_km_crs   Default kilometre-based projected CRS for spatial modelling.
#   make_hex_grid    Hexagon grid over a study area (by cell width or by count).
#
# Sourced by 01, 06 and 09.
# ============================================================

# Return the default kilometre-based projected CRS used for spatial modelling.
# A Lambert conformal conic in km units keeps SPDE mesh edge lengths, buffers,
# and simplification tolerances all expressed directly in kilometres.
get_aea_km_crs <- function() {
  # Single-line PROJ string to avoid newline/whitespace surprises.
  proj <- paste(
    "+proj=lcc",
    "+lat_1=49",
    "+lat_2=77",
    "+lat_0=49",
    "+lon_0=-95",
    "+datum=NAD83",
    "+units=km",
    "+no_defs"
  )
  sf::st_crs(proj)
}

# Hexagon grid over a study area, keeping only hexagons that intersect it.
#
# Size the hexagons with EXACTLY ONE of:
#   width_km    Cell size passed straight to st_make_grid(), in CRS units (km
#               for the study CRS). hex_id is assigned BEFORE non-intersecting
#               hexagons are dropped, so ids are not contiguous. 06 uses this
#               for the 25-km aggregation grid; that hex_id is carried into every
#               prediction file (07) and region summary (08), so keep it stable.
#   n_hexagons  Target number of hexagons intersecting the study area. The area
#               is validated and unioned first (projected to EPSG:3978 if it is
#               geographic, then returned in the input CRS), and hex_id is
#               renumbered 1..n after filtering. 09 uses this for the
#               per-region assessment grids.
#
# Returns an sf with columns hex_id (integer) and geometry.
make_hex_grid <- function(study_area, width_km = NULL, n_hexagons = NULL) {
  
  stopifnot(inherits(study_area, "sf"))
  if (is.null(width_km) == is.null(n_hexagons)) {
    stop("make_hex_grid(): supply exactly one of width_km or n_hexagons.")
  }
  
  # ---- Fixed cell width ----
  if (!is.null(width_km)) {
    hex <- sf::st_make_grid(study_area, cellsize = width_km, square = FALSE)
    
    hex_sf <- sf::st_sf(
      hex_id   = seq_along(hex),
      geometry = hex
    )
    
    return(hex_sf[sf::st_intersects(hex_sf, study_area, sparse = FALSE), , drop = FALSE])
  }
  
  # ---- Target hexagon count ----
  old_s2 <- sf::sf_use_s2()
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  sf::sf_use_s2(FALSE)
  
  original_crs <- sf::st_crs(study_area)
  
  if (sf::st_is_longlat(study_area)) {
    study_area_work <- sf::st_transform(study_area, 3978)
  } else {
    study_area_work <- study_area
  }
  
  study_area_work <- study_area_work |>
    sf::st_make_valid() |>
    sf::st_collection_extract("POLYGON") |>
    sf::st_buffer(0)
  
  study_area_union <- sf::st_union(study_area_work)
  
  target_area  <- as.numeric(sf::st_area(study_area_union)) / n_hexagons
  hex_cellsize <- sqrt(target_area / (sqrt(3) / 2))
  
  hex_grid <- sf::st_make_grid(
    study_area_union,
    cellsize = hex_cellsize,
    square   = FALSE,
    what     = "polygons"
  )
  
  hex_sf <- sf::st_sf(
    hex_id   = seq_along(hex_grid),
    geometry = hex_grid,
    crs      = sf::st_crs(study_area_work)
  )
  
  keep   <- lengths(sf::st_intersects(hex_sf, study_area_union)) > 0
  hex_sf <- hex_sf[keep, ]
  hex_sf$hex_id <- seq_len(nrow(hex_sf))
  
  if (!is.na(original_crs)) {
    hex_sf <- sf::st_transform(hex_sf, original_crs)
  }
  
  hex_sf
}
