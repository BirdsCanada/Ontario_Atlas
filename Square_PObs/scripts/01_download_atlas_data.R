# Download and clean Atlas-3 data #
# Author: Sarah Mueller
# Last updated: Nov 19, 2026

# Get Atlas data from NatureCounts and format for easy use
# Update species' scientific and common names

# NatureCounts account is required, and must request access to the following Atlas datasets:
# ONATLAS3BE_DO
# ONATLAS3PC
# ONATLAS3SS
# OBBA2BE_RAW
# OBBA2PC
# OBBA1BE_RAW

# Set up -----------------

# load naturecounts and other packages
library(naturecounts)
library(NoSleepR) 
library(tidyverse)
library(janitor)

# update NatureCounts metadata
nc_metadata()

# Settings ------------

# path to raw data folder
data_path <- Sys.getenv("ONTARIO_ATLAS_DATA", unset = "")

# path to output folder
output_path <- Sys.getenv("OUTPUT_DIR", unset = "")

# directory for downloaded data
download_dir <- paste0(output_path, "download")

# create the directory if it doesn't exist already
if (!dir.exists(download_dir)) dir.create(download_dir, recursive = TRUE)

# Info for downloads ---------

my_username <- Sys.getenv("NATURECOUNTS_USER", unset = "")

# Short description of reason for download
my_info <- Sys.getenv("NATURECOUNTS_INFO", unset = "")

# fill in with numeric request IDs
request_3BE <- Sys.getenv("NATURECOUNTS_request_3BE", unset = "")
request_3PC <- Sys.getenv("NATURECOUNTS_request_3PC", unset = "")
request_3SS <- Sys.getenv("NATURECOUNTS_request_3SS", unset = "")

# Load taxonomy update table ----------

# Taxonomy update lookup table
taxonomy_lookup <- read_csv(paste0(data_path, "species_taxonomy_lookup_table.csv")) %>%
  distinct(original_scientific_name, .keep_all = T) %>%
  select(-original_common_name, -atlas_period)

# Helper functions to clean up dataframes --------

# Cleans up downloaded dataframes
# 1) removes empty columns
# 2) cleans up variable names
# 3) filters to only Ontario and Akimiski NU records
# 4) trim leading and trailing whitespace from all character variables
# 5) add a simplified survey event ID
# 6) update species names
clean_dataframe <- function(df) {
  
  df %>%
    remove_empty("cols") %>%
    clean_names() %>%
    filter(statprov_code == "ON" | subnational2_code == "CA.NU.QK") %>%
    # trim leading and trailing whitespace from all character variables
    mutate(across(where(is.character), str_trim)) %>%
    # add a simple sampling event ID with just the unique identifier number
    # e.g. take just the middle number from "ONATLAS-153157-0"
    mutate(event_id = str_split_i(sampling_event_identifier, "-", 2)) %>%
    rename(original_scientific_name = scientific_name,
           original_common_name = common_name) %>%
    left_join(taxonomy_lookup, by = "original_scientific_name") %>%
    # recode to Y or N (sometimes "Yes")
    mutate(all_species_reported = case_when(
      all_species_reported == "No" ~ "N",
      all_species_reported == "Yes" ~ "Y",
      T ~ NA
    ))
  
}

# Fix data types in Atlas 1 and 2 dataframes
fix_types12 <- function(df) {
  
  df %>%
    mutate(across(record_id, as.character),
           across(c(obs_count_at_least, decimal_longitude, decimal_latitude, 
                    year_collected, effort_measurement1),
                  as.numeric))
  
}

# Download Atlas-3 datasets -------------

# enter NatureCounts password when prompted for first download

##### Download and prep ONATLAS3BE_DO dataset ---------------

with_nosleep({
  
  atlas3_raw <- nc_data_dl(collections = "ONATLAS3BE_DO",
                       info = my_info,
                       request_id = request_3BE,
                       username = my_username,
                       timeout = 10000)
  
})

# clean up
atlas3 <- atlas3_raw %>%
  clean_dataframe()

##### Download and prep Atlas-3 point count dataset ---------------

point_counts_raw <- nc_data_dl(collections = "ONATLAS3PC",
                           info = my_info,
                           request_id = request_3PC,
                           username = my_username,
                           timeout = 10000)

# remove columns that contain only NAs and keep only important columns
point_counts <- point_counts_raw %>%
  clean_dataframe()

##### Download and prep Atlas-3 special survey dataset ---------------

special_surveys_raw <- nc_data_dl(collections = "ONATLAS3SS",
                           info = my_info,
                           request_id = request_3SS,
                           username = my_username,
                           timeout = 10000)

# remove columns that contain only NAs and keep only important columns
special_surveys <- special_surveys_raw %>%
  clean_dataframe()

##### Combine Atlas-3 datasets --------------

# Add data from other Atlas-3 datasets not already present in ONATLAS3BE_DO

# subset of point counts not in ONATLAS3BE_DO
point_counts_sub <- point_counts %>%
  filter(!event_id %in% atlas3$event_id)

# subset of special surveys not in ONATLAS3BE_DO 
special_surveys_sub <- special_surveys %>%
  filter(!event_id %in% atlas3$event_id)

# Add to ONATLAS3BE_DO
atlas3_complete <- atlas3 %>%
  bind_rows(point_counts_sub, special_surveys_sub) %>%
  mutate(survey_type = case_when(
    event_id %in% point_counts$event_id ~ "point count",
    event_id %in% special_surveys$event_id ~ "special survey",
    T ~ "checklist"
  )) %>%
  # fix data types
  mutate(across(record_id, as.character),
         across(c(decimal_longitude, decimal_latitude, 
                  year_collected, month_collected, day_collected, time_collected,
                  julian_day, time_observations_started, time_observations_ended,
                  duration_in_hours, number_of_observers,
                  observation_count, date_uncertainty_in_days),
                as.numeric))

# Download Atlas-2 datasets -------------

# download breeding evidence
atlas2_raw <- nc_data_dl(collections = "OBBA2BE_RAW",
                         info = my_info,
                         username = my_username,
                         timeout = 10000)

# download point counts
atlas2_pc_raw <- nc_data_dl(collections = "OBBA2PC",
                            info = my_info,
                            username = my_username,
                            timeout = 10000)

# prep Atlas-2 datasets
atlas2 <- atlas2_raw %>%
  clean_dataframe() %>%
  fix_types12()

atlas2_pc <- atlas2_pc_raw %>%
  clean_dataframe() %>%
  mutate(across(record_id, as.character),
         across(c(decimal_longitude, decimal_latitude, observation_count, 
                  year_collected, month_collected, day_collected,
                  time_collected, duration_in_hours, effort_measurement1),
                as.numeric))

# Download Atlas-1 dataset -------------

atlas1_raw <- nc_data_dl(collections = "OBBA1BE_RAW",
                     info = my_info,
                     username = my_username,
                     timeout = 10000)

atlas1 <- atlas1_raw %>%
  clean_dataframe() %>%
  fix_types12()

# save the Atlas data ----------------

saveRDS(atlas1, file = file.path(download_dir, "atlas1.rds"))
saveRDS(atlas2, file = file.path(download_dir, "atlas2.rds"))
saveRDS(atlas2_pc, file = file.path(download_dir, "atlas2_pc.rds"))
saveRDS(atlas3_complete, file = file.path(download_dir, "atlas3_complete.rds"))


