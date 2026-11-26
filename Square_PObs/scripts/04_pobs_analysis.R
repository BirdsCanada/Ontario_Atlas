# Estimate BCR-level probability of observation across Atlas 1, 2, and 3 #
# Author: Sarah Mueller
# Last updated: Nov 19, 2026

##### run first three scripts (download and prep data) first

# Set up --------------

rm(list = ls())

library(tidyverse)
library(sf)
# note: requires ggplot2 >= 4.0.0 to use `space` in facet_wrap() 

# for regression
library(glmmTMB)
library(emmeans)
library(DHARMa)

# Settings -----------

# path to raw data folder
data_path <- Sys.getenv("ONTARIO_ATLAS_DATA", unset = "")

# path to output folder
output_path <- Sys.getenv("OUTPUT_DIR", unset = "")

# output directories
final_output_dir <- paste0(output_path, "species_predictions")
model_dir <- paste0(output_path, "species_models")
# create directories if they do not exist
if(!dir.exists(final_output_dir)) dir.create(final_output_dir, recursive = TRUE)
if(!dir.exists(model_dir)) dir.create(model_dir, recursive = TRUE)

# BCR numbers in the order they should appear in factors and plots
bcr_number_order <- c("13", "12", "77", "76", "74")

# color palette for Atlas periods
atlas_colors <- RColorBrewer::brewer.pal(9, "Blues")[c(5,7,9)]
names(atlas_colors) <- c("Atlas 1", "Atlas 2", "Atlas 3")

# Read in data -----------

# Atlas species list
species_list <- read_csv(paste0(data_path, "atlas3_specieslist.csv"))

# Load the Atlas data prepped in the previous scripts
all_data_model_sub <- readRDS(paste0("output/final_PObs_dataset/final_PObs_dataset_", data_date, ".rds"))

# Load Atlas square shapefile
atlas_squares <- read_sf("output/atlas_squares_complete_with_BCR/atlas_squares_with_bcr.shp")

# Get number of squares in each BCR to use as area weights ---------

n_squares <- atlas_squares %>%
  st_drop_geometry() %>%
  mutate(bcr.f = factor(bcr_number, levels = bcr_number_order)) %>%
  count(bcr.f, name = "n_squares")

# Source functions -----------

# Functions:
# 1) get_species_data returns a subset of data for one species for the analysis
# 2) fit_species_model_auto fits model with appropriate terms for a species
# 3) predict_20h genenerates model predictions at 20 h effort
#       returns list(prediction table, ggplot, model)
# 4) save_species_plot saves model and species prediction plot 
#       returns prediction table

source("scripts/00_pobs_analysis_functions.R")

# Loop through species to fit models and save prediction plots and table --------------

species_tables <- vector("list", nrow(species_list))
failed_species <- character()

NoSleepR::with_nosleep({
  
  for(i in seq_len(nrow(species_list))){
    
    species_i <- species_list$ScientificName[i]
    
    # only attempt species that are in the dataset used for modelling
    if(!species_i %in% all_data_model_sub$scientific_name) next
    
    # one failure shouldn't kill the whole run
    species_tables[[i]] <- tryCatch(
      save_species_plot(species_i),
      error = function(e){
        message("FAILED: ", species_i, " - ", conditionMessage(e))
        failed_species <<- c(failed_species, species_i)
        NULL
      }
    )
    
    if(i %% 25 == 0) gc()
  }
  
})

all_species_table <- bind_rows(species_tables)

# save the species prediction table
write_csv(all_species_table, file = paste0(output_path, "all_species_predictions_table.csv"))

if(length(failed_species) > 0){
  message(length(failed_species), " species failed; see output/failed_species.txt")
  writeLines(failed_species, paste0(output_path, "failed_species.txt"))
}
