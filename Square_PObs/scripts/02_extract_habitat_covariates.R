# Extract habitat covariates for three-atlas PObs analysis #
# Author: Sarah Mueller
# Last updated: Nov 19, 2026

# Run script 1 first #

# Raw data available at:
# https://doi.org/10.6084/m9.figshare.35063872

# Set up --------------

rm(list = ls())

library(tidyverse)
library(sf)
library(terra)
library(exactextractr)
library(compositions)

# Settings ----------

# path to raw data folder
data_path <- Sys.getenv("ONTARIO_ATLAS_DATA", unset = "")

# path to output folder
output_path <- Sys.getenv("OUTPUT_DIR", unset = "")

# Read in data -----------

# Atlas square shapefile with BCRs 
atlas_squares <- read_sf(paste0(data_path, "atlas_squares_complete_with_BCR/atlas_squares_with_bcr.shp"))

# Landcover rasters
# derived from GLC_FCS30D, see https://zenodo.org/records/15063683
landcover1985 <- rast(paste0(data_path, "landcover/ON_Landcover_1985_clipped.tif"))
landcover2003 <- rast(paste0(data_path, "landcover/ON_Landcover_2003_clipped.tif"))
landcover2022 <- rast(paste0(data_path, "landcover/ON_Landcover_2022_clipped.tif"))

# read in the value-class keys
class_keys <- read_csv(paste0(data_path, "landcover/Landcover_class_values_key.csv"))

# Extract habitat covariates -------------

# match CRS of square polygons to rasters
atlas_squares <- st_transform(atlas_squares, crs = st_crs(landcover1985))

# Extract landcover raster values one square at a time
# (all at once takes too much memory)

# list to store results, same length as squares dataset
result_list <- vector("list", nrow(atlas_squares))

# loop over squares, extracting pixel fractions in each square
for(i in seq_len(nrow(atlas_squares))){
  
  # get the ith square
  square_i <- atlas_squares[i,]
  
  # get the spatial extent of the ith square
  square_ext <- ext(vect(square_i))
  
  # crop the landcover rasters to just the ith square
  lc_i_1985 <- crop(landcover1985, square_ext)
  lc_i_2003 <- crop(landcover2003, square_ext)
  lc_i_2022 <- crop(landcover2022, square_ext)
  
  # extract coverage fractions from each raster for just the ith square 
  # dataframes have one row per raster cell that overlaps square i, with:
  # value = landcover class of that pixel
  # coverage fraction = proportion of the pixel's area inside the square
  df_i_1985 <- exact_extract(lc_i_1985,
                             square_i,
                             fun = NULL, # returns per-pixel data
                             progress = F)[[1]] %>%
    # add a year column
    mutate(lc_year = 1985)
  df_i_2003 <- exact_extract(lc_i_2003,
                             square_i,
                             fun = NULL, # returns per-pixel data
                             progress = F)[[1]] %>%
    # add a year column
    mutate(lc_year = 2003)
  df_i_2022 <- exact_extract(lc_i_2022,
                             square_i,
                             fun = NULL, # returns per-pixel data
                             progress = F)[[1]] %>%
    # add a year column
    mutate(lc_year = 2022)
  
  # bind output dataframes together
  df_i <- df_i_1985 %>%
    bind_rows(df_i_2003, df_i_2022) %>%
    # drop rows with value = NA, 0, or 250
    drop_na(value) %>%
    filter(!value %in% c(0, 250))
  
  # calculate proportion of each landcover in the square and add to result list
  result_list[[i]] <- df_i %>%
    group_by(value, lc_year) %>%
    summarize(area = sum(coverage_fraction),
              .groups = "drop") %>%
    group_by(lc_year) %>%
    mutate(year_total_area = sum(area)) %>%
    ungroup() %>%
    # the proportion of each lc type in each year = area of the type / total area each year
    mutate(prop = area / year_total_area,
           square_id = square_i$utm_square) %>%
    select(square_id,
           lc_year,
           value,
           prop,
           year_total_area) 
  
  rm(df_i, df_i_1985, df_i_2003, df_i_2022, lc_i_1985, lc_i_2003, lc_i_2022)
  if(i %% 10 == 0) gc() # clean up garbage every 10th loop
  
}

# bind into one dataframe
prop_df <- bind_rows(result_list) %>%
  arrange(square_id, lc_year, value)

# join to proportions dataframe
class_props <- prop_df %>%
  left_join(class_keys, by = "value") %>%
  # combine proportions for simplified classes (sum proportions of component classes)
  group_by(square_id, lc_year, simplified_class, year_total_area) %>%
  summarize(simplified_prop = sum(prop)) %>%
  arrange(square_id, lc_year, simplified_class)

# pivot wider so that each class proportion is a column, one row per square per Atlas period
prop_wide <- class_props %>%
  # add a column with classes reformatted as variable names
  mutate(variable_name = str_to_lower(str_replace_all(simplified_class,
                                                      " ",
                                                      "_"))) %>%
  pivot_wider(id_cols = c(square_id, lc_year, year_total_area),
              names_from = variable_name,
              values_from = simplified_prop,
              values_fill = 0,
              names_prefix = "prop_") %>%
  # add an Atlas period variable
  mutate(atlas_period = case_when(
    lc_year == 1985 ~ "Atlas 1",
    lc_year == 2003 ~ "Atlas 2",
    lc_year == 2022 ~ "Atlas 3",
    T ~ NA))

# save the habitat dataframe
write_csv(prop_wide, file = paste0(output_path, "habitat.csv"))
