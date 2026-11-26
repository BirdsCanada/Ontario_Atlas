# Prep data for three-Atlas PObs analysis #
# Author: Sarah Mueller
# Last updated: Nov 20, 2026

# IF NOT USING THE PROVIDED ATLAS DATASETS: run scripts 1 and 2 first
# Otherwise, start here

# Output dataframes have one row per atlasser x square x year x species, with:
#   binary detected variable,
#   effort, habitat PCs, BCR
# 1) all_model_data has the complete dataset
# 2) all_model_data_sub has data from BCRs where a species does not breed removed

# Set up --------------

rm(list = ls())

library(tidyverse)
library(sf)

# Settings ---------

# path to raw data folder
data_path <- Sys.getenv("ONTARIO_ATLAS_DATA", unset = "")

# path to output folder
output_path <- Sys.getenv("OUTPUT_DIR", unset = "")

# directory for downloaded NatureCounts data
download_dir <- paste0(output_path, "download")

# path for output prepared dataframe
prepped_data_dir <- paste0(output_path, "final_PObs_dataset")
if(!dir.exists(prepped_data_dir)) dir.create(prepped_data_dir, recursive = TRUE)

# Breeding codes (for all three Atlases)
# "X" is the only non-breeding code (or NA), so any not listed are treated as non-breeding
breeding_codes <- c("SH", "SM", "H", "S", "M", "P", "T", "D", "V", "A", "B", "N",
                    "NB", "DD", "NU", "FY", "AE", "FS", "CF", "NE", "NY")

# Drop species with fewer detections than this
min_detections <- 10

# Number of habitat PCs to keep
n_pcs <- 4

# Warn if the kept PCs explain less variance than this threshold
min_var <- 0.80

# Share of highest effort ratios dropped for conversion factor
ratio_trim <- 0.1

# Helper functions -----------

# min/max functions that return NA (not Inf and a warning) when every value in a group is NA
safe_min <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)

# Read in data ---------

# Atlas square shapefile with BCRs 
atlas_squares <- read_sf(paste0(data_path, "atlas_squares_complete_with_BCR/atlas_squares_with_bcr.shp"))

# Atlas species list 
species_list <- read_csv(paste0(data_path, "atlas3_specieslist.csv"))

# Special survey species list
ss_species <- read_csv(paste0(data_path, "special_survey_species.csv"))

# read in safe dates spreadsheet and process
safe_dates <- read_csv(paste0(data_path, "OnAtlasSafeDates_update_2026-11-18.csv")) %>%
  # recode Redpoll (Hoary) and Redpoll (Common)
  mutate(english_name = if_else(str_detect(english_name, "Redpoll"), "Redpoll", english_name)) %>%
  # join scientific names
  left_join(species_list, by = c("english_name" = "CommonName"))
# levels: 3 = warning (non-breeding), 2 = caution (shoulder), 1 = breeding (safe)

# CWS waterfowl survey info
cws_surveys <- read_csv(paste0(data_path, "waterfowl_aerial_survey_for_analysis_good.csv"),
                        # match column types in Atlas-3 dataset
                        col_types = cols(survey_month = col_double(),
                                         survey_day = col_double(),
                                         MonthCollected = col_double(),
                                         DayCollected = col_double(),
                                         SpeciesCode = col_character())) %>%
  janitor::clean_names() %>%
  # add sampling_event_identifier to uniquely identify surveys = square + date + time started
  mutate(sampling_event_identifier = paste(utm_square, observation_date, time_collected, sep = ", "),
         # and protocol code
         protocol_code = "CWS waterfowl surveys",
         # add collector number by survey type
         collector_number = str_split_i(record_id, "_", 1)) %>%
  left_join(species_list, by = c("common_name" = "CommonName")) %>%
  mutate(scientific_name = coalesce(scientific_name, ScientificName)) %>%
  drop_na(scientific_name) # drops WRSA, REPH, and BBPL

# Habitat covariates
prop_wide <- read_csv(paste0(output_path, "habitat.csv"))

# Load the Atlas 1-3 data downloaded in script 1
# 3 dataframes: atlas1, atlas2, atlas3_complete
atlas1 <- readRDS(file.path(download_dir, "atlas1.rds"))
atlas2 <- readRDS(file.path(download_dir, "atlas2.rds"))
atlas3_complete <- readRDS(file.path(download_dir, "atlas3_complete.rds"))

# Define sets of target species for each survey type ----------

# General and point counts - all species
all_species <- species_list %>%
  rename(scientific_name = ScientificName) %>%
  pull(scientific_name)
# Owl species
all_owls <- ss_species %>%
  filter(species_type == "owl") %>%
  pull(scientific_name)
# Nightjar species
all_nightjars <- ss_species %>%
  filter(species_type == "nightjar") %>%
  pull(scientific_name)
# Marsh bird species
all_marshbirds <- ss_species %>%
  filter(species_type == "marshbird") %>%
  pull(scientific_name)

# Prep Atlas 1 and 2 data --------------------

# join Atlas 1 and Atlas 2 dataframes
# convert Effort to numeric type
# keep only records with breeding evidence
atlas12 <- atlas1 %>%
  bind_rows(atlas2) %>%
  mutate(total_hours = effort_measurement1) %>%
  filter(breeding_bird_atlas_code %in% breeding_codes)

# Summarize Atlas-3 data to level of Atlas-1 and 2 --------------

# one record atlasser x square x year x species
# drop incidental records
# maintain "no species observed" checklists, which have species_id 0 or NA

# prep Atlas-3 data for summarizing
atlas3_tmp1 <- atlas3_complete %>%
  # drop incidental surveys
  filter(all_species_reported != "N") %>%
  # join CWS waterfowl surveys
  bind_rows(cws_surveys) %>%
  # variable indicating no species observed on a survey
  mutate(no_species = species_id == 0 | is.na(species_id)) %>%
  # whether record has breeding evidence (any level)
  mutate(breeding = breeding_bird_atlas_code %in% breeding_codes) %>%
  # did the atlasser find breeding evidence (any level) for a species in a 
  # square in a year?
  group_by(scientific_name, utm_square, survey_year, collector_number) %>%
  mutate(breeding_in_square = any(breeding)) %>%
  ungroup() %>%
  # 1) clean durations and start/end times
  rename(time_start = time_observations_started,
         time_end = time_observations_ended,
         duration = duration_in_hours) %>%
  mutate(
    # durations >= 24 h are impossible, make NA
    duration = if_else(duration >= 24, NA_real_, duration),
    # durations > 20 h with end times > 30 are probably inaccurate
    duration = ifelse(duration > 20 & time_end > 30, NA_real_, duration),
    # 2) fix weird start times
    time_start_new = case_when(
      # values 500-2400 are probably clock times, convert to hours
      time_start >= 500 & time_start <= 2400 ~ time_start/100, 
      # 0s and remaining values > 24 are not usable 
      time_start > 24 | time_start == 0 ~ NA_real_, 
      # otherwise, keep original time
      T ~ time_start),
    # 3) fix weird end times
    time_end_new = case_when(
      # clock time -> hours
      time_end >= 500 & time_end <= 2400 ~ time_end/100, 
      # end times can legitimately be > 24 (span midnight) but not >= 26
      time_end >= 26 | time_end == 0 ~ NA_real_, 
      # end earlier than start: use start + duration
      time_end < time_start_new ~ time_start_new + duration,
      # otherwise, keep original time
      T ~ time_end)
  ) %>%
  # drop rows that have no start time, no end time, AND no usable dration:
  # there is no way to work out effort for these
  filter(!(is.na(time_end_new) & is.na(time_start_new) & (duration == 0 | is.na(duration)))) %>%
  # drop collector number 0
  filter(collector_number != 0)

# survey-level data (one row per checklist) to calculate daily survey time
atlas3_survey_level <- atlas3_tmp1 %>%
  # one row per survey
  distinct(sampling_event_identifier, .keep_all = T) %>%
  # if duration is 0 or NA, replace it with the end time - start time 
  # (as long as that is not also NA)
  mutate(missing_duration = is.na(duration) | duration == 0,
         new_duration = time_end_new - time_start_new,
         has_new_duration = !(is.na(new_duration)),
         duration_filled = case_when(
           # use new duration if original is missing and it has a new one
           missing_duration & has_new_duration ~ new_duration,
           # if it's not missing duration, use the original
           !missing_duration ~ duration,
           # if duration is 0 and wasn't filled in with the new one above, keep it
           duration == 0 ~ duration,
           # otherwise NA
           T ~ NA)) 
# duration_filled = original duration, unless it is missing or 0 and 
# end - start time is available, in which case end - start is used
summary(atlas3_survey_level$duration_filled)

##### Summarize times at daily level -----------

# calculate daily survey times
# one row per atlasser per square per day
daily_times <- atlas3_survey_level %>%
  # calculate daily duration and survey start and end times for each square-atlasser-day
  group_by(collector_number, utm_square, observation_date) %>%
  mutate(
    # daily start time is the earliest start time in the group
    daily_start_time = safe_min(time_start_new),
    # last daily start time is the latest start time in the group
    last_daily_start_time = safe_max(time_start_new),
    # max end time is the latest end time in the group
    max_time_end_new = safe_max(time_end_new),
    # many records are missing end times, so if missing, use last start time
    daily_end_time = coalesce(max_time_end_new, last_daily_start_time),
    # daily end minus start time
    daily_start_to_end = daily_end_time - daily_start_time,
    # sum of all durations in that group
    daily_duration = sum(duration_filled, na.rm = T),
    # paste together a string of all the remarks in that group
    all_remarks = paste(na.omit(remarks), collapse = "; ")
  ) %>%
  ungroup() %>%
  # keep one row per atlasser per square per day
  # NOTE: remarks comes from the first record of the day only
  distinct(collector_number, utm_square, observation_date, .keep_all = T) %>%
  select(collector_number, utm_square, observation_date, survey_year,
         daily_start_time, last_daily_start_time, daily_end_time,
         daily_start_to_end, daily_duration, all_remarks) %>%
  # get a final measure of daily square effort
  mutate(daily_effort = case_when(
    # summed duration > 22 h is probably wrong; use end - start if available
    daily_duration > 22 & !(is.na(daily_start_to_end) | daily_start_to_end == 0) ~ daily_start_to_end,
    # when summed duration is > 24 and daily_start_to_end is 0, 
    # use daily_start_to_end (only four records)
    daily_duration > 24 & daily_start_to_end == 0 ~ 0,
    # otherwise, use summed duration
    T ~ daily_duration))
# daily_effort is usually the summed durations for an atlasser-square-day,
# unless the daily duration is wrong; then it's daily end - start time
summary(daily_times$daily_effort)

##### Summarize times at yearly level ------------

# yearly effort = sum of daily effort for an atlasser-square-year
yearly_times <- daily_times %>%
  group_by(collector_number, utm_square, survey_year) %>%
  summarize(yearly_effort = sum(daily_effort, na.rm = T)) %>%
  ungroup() 
summary(yearly_times$yearly_effort)

##### Prepare atlasser-square-year-species level Atlas-3 dataset --------

# Summarize Atlas-3 dataset to atlasser-square-year-species level
atlas3_summ <- atlas3_tmp1 %>%
  # Add a survey-level species detection variable, which is 0 if it's a no-species survey or the 
  # observation count is 0, but otherwise 1
  mutate(detected_survey = if_else(no_species | observation_count == "0", 0, 1,
                                   # if observation_count is missing, make it 1
                                   missing = 1)) %>%
  # group by atlasser-square-year
  group_by(collector_number, utm_square, survey_year) %>%
  # Simplify protocol codes and collect all survey within an atlasser-square-year record
  # all surveys for an atlasser-square-year are special if none are "NORAC" or "Point Count 6 Interval"
  mutate(simple_protocol_code = case_when(
           grepl("OWL|LEOW", protocol_code) ~ "Owls",
           grepl("NORAC|Point Count|CWS", protocol_code) ~ "General",
           T ~ protocol_code),
         all_protocols = paste(sort(unique(simple_protocol_code)), collapse = "; "),
         all_special = !grepl("General", all_protocols)) %>%
  ungroup() %>%
  # Species was detected at the atlasser-square-year level if there were 
  # any survey-level detections within that group
  group_by(collector_number, utm_square, survey_year, scientific_name) %>%
  mutate(detected = as.numeric(any(detected_survey == 1))) %>%
  ungroup() %>%
  distinct(collector_number, utm_square, survey_year, scientific_name, detected,
           no_species, all_protocols, all_special, breeding_in_square) %>%
  # drop "no species observed rows" if there are other rows for the atlasser-square-year group
  group_by(collector_number, utm_square, survey_year) %>%
  filter(!(no_species & any(no_species == F))) %>%
  ungroup() %>%
  # join effort hours by atlasser-square-year
  left_join(yearly_times, 
            by = c("collector_number", "utm_square", "survey_year")) %>%
  # rename columns to match Atlas 1-2
  rename(total_hours = yearly_effort) %>%
  ungroup() 

# Join Atlas 1, 2, and 3 data -------------------

# Join all datasets and prep for analysis
all_data <- atlas12 %>%
  as_tibble() %>%
  # add detection variable with value 1 for all rows for Atlas 1-2
  # and fill in "General" for protocol
  mutate(detected = 1,
         all_protocols = "General",
         breeding_in_square = breeding_bird_atlas_code %in% breeding_codes) %>%
  # add Atlas 3 data
  bind_rows(atlas3_summ) %>%
  # keep only relevant columns
  select(names(atlas3_summ)) %>%
  # add a unique record identifier: atlasser ID + square ID + survey year
  mutate(record_id = paste(collector_number, utm_square, survey_year))
summary(all_data$total_hours) # NAs and 0s in total hours are still present
table(all_data$detected, useNA = "ifany") # small number of 0s 

##### Fix some incorrect square IDs in Atlas dataset -----------

# Fix wrong square IDs in Atlas data
all_data <- all_data %>%
  mutate(utm_square = case_when(
    # in Michigan as recorded, but lat/long is around River Aux Sables PP
    utm_square == "16THS73" ~ "17TMM12",
    # probably a typo, lat/long is near Lac Seul
    utm_square == "45UZS98" ~ "15UWS26",
    # otherwise don't change
    T ~ utm_square)) 

##### Drop duplicates -------------------------

# drop duplicates
all_data_unique <- all_data %>%
  # replace remaining NAs in scientific name with "no species observed"
  mutate(scientific_name = replace_na(scientific_name, "no species observed")) %>%
  # drop duplicates by atlasser-square-species-year and effort hours
  distinct(utm_square, collector_number, survey_year, scientific_name, 
           total_hours, .keep_all = T) %>%
  # if every record in an atlasser-square-year has NA hours, set hours to 0
  group_by(record_id) %>%
  mutate(new_total_hours = if (all(is.na(total_hours))) 0 else total_hours) %>%
  # does the group have another record with hours?
  mutate(group_has_hours = any(!is.na(new_total_hours)),
         # get the max hours in that group
         group_max_hours = safe_max(new_total_hours)) %>%
  ungroup() %>%
  # if a species is missing hours but its group has records with hours, use the group max hours
  mutate(new_total_hours = if_else(is.na(new_total_hours) & group_has_hours,
                                   group_max_hours,
                                   new_total_hours)) %>%
  # for Atlas-1 records with 2 records per atlasser-square-year-species,
  # keep only the first record (with smaller effort) 
  arrange(collector_number, utm_square, survey_year, scientific_name, total_hours) %>%
  group_by(collector_number, utm_square,survey_year, scientific_name) %>%
  slice(1) %>%
  ungroup() %>%
  select(-group_has_hours, -group_max_hours) %>%
  # rename new total hours column
  rename(original_total_hours = total_hours,
         total_hours = new_total_hours)
summary(all_data_unique$total_hours) # no NAs, still has 0s

# check that all atlasser-square-year-species have 1 row
stopifnot(nrow(all_data_unique %>%
  count(collector_number, utm_square, survey_year, scientific_name) %>%
  filter(n > 1)) == 0)

##### Fill in non-detections --------------

# Fill in non-detections for the appropriate species, and only for the square-years an
# atlasser actually surveyed

# atlasser x square x year x species-specific hours, joined back after completing the data
hours_lookup <- all_data_unique %>%
  select(record_id, scientific_name, total_hours) 

# Atlasser-square-year records with general atlassing are complete for all species
# Records with only special surveys are complete only for their target species
# Fill in non-detections of appropriate species for each combo of survey types 
# all species

# function to complete records with appropriate species based on target groups
complete_records <- function(df = df,
                            groups = c("Owl", "Nightjars", "Marsh Birds", "General")){
  
  # search term and target species list
  terms <- character(0)
  target_species <- character(0)
  if("Owl" %in% groups) {
    terms <- c(terms, "OWLS", "LEOW")
    target_species <- c(target_species, all_owls)
    }
  if("Nightjars" %in% groups) {
    terms <- c(terms, "Nightjars")
    target_species <- c(target_species, all_nightjars)
    }
  if("Marsh Birds" %in% groups) {
    terms <- c(terms, "Marsh")
    target_species <- c(target_species, all_marshbirds)
  }
  if("General" %in% groups) {
    terms <- c(terms, "General")
    target_species <- c(target_species, all_species)
  }
  
  search_term <- paste(terms, collapse = "|")
  target_species <- unique(target_species)
  
  # complete records for the target species
  df %>%
    filter(grepl(search_term, all_protocols)) %>%
    select(-original_total_hours, -total_hours) %>%
    ungroup() %>%
    complete(nesting(utm_square, survey_year, collector_number, record_id, all_protocols),
             scientific_name = target_species,
             fill = list(detected = 0)) %>%
    arrange(record_id, scientific_name)

}

# Split dataset into list of separate dataframes based on protocol combinations
protocol_subsets <- split(all_data_unique, all_data_unique$all_protocols)

# Loop across each dataframe and complete them, based on target species
#protocols_complete <- tibble()
protocols_complete <- vector("list", 0)
for(i in seq_len(length(protocol_subsets))){
  
  groups_i <- str_split_1(protocol_subsets[[i]]$all_protocols[1], "; ")
  subset_complete <- complete_records(df = protocol_subsets[[i]], 
                                      groups = groups_i)
  protocols_complete[[i]] <- subset_complete
  # protocols_complete <- protocols_complete %>%
  #   bind_rows(subset_complete)
}

# Bind completed datasets together
# and drop detections without breeding evidence
all_data_complete <- bind_rows(protocols_complete) %>%
  arrange(record_id, scientific_name) %>%
  # join back total hours by atlasser-square-year group and scientific name
  left_join(hours_lookup, 
            by = c("record_id", "scientific_name")) %>%
  # non-detections (and NA hours) get the maximum hours in the atlasser-square-year
  group_by(record_id) %>%
  mutate(total_hours = coalesce(total_hours, max(total_hours, na.rm = T))) %>%
  ungroup() %>%
  # drop records where species was detected but wasn't recorded breeding
  # records added as non-detections in completion step have NA for breeding_in_square and are retained 
  filter(!(detected & !breeding_in_square))
summary(all_data_complete$total_hours)

# there should be one row per atlasser-square-year-species 
stopifnot(nrow(all_data_complete) == nrow(distinct(all_data_complete,
                                         utm_square, collector_number, survey_year, scientific_name)))

# # table number of detections and non-detections for each species
# with(all_data_complete, table(scientific_name, detected))

##### Add an Atlas period variable ----------------

all_data_tmp1 <- all_data_complete %>%
  # keep only rows within the Atlas data collection periods
  filter(survey_year %in% c(1981:1985, 2001:2005, 2021:2025)) %>%
  # add an "atlas period" variable
  mutate(atlas_period = case_when(survey_year %in% c(1981:1985) ~ "Atlas 1",
                                  survey_year %in% c(2001:2005) ~ "Atlas 2",
                                  survey_year %in% c(2021:2025) ~ "Atlas 3",
                                  T ~ NA_character_)) 

##### Join habitat covariates and BCR by square and atlas period ----------

atlas_squares_df <- st_drop_geometry(atlas_squares)
# a square must appear once, otherwise the join duplicates rows
stopifnot(!anyDuplicated(atlas_squares_df$utm_square))

# join habitat covariates and BCR
all_data_tmp2 <- all_data_tmp1 %>%
  left_join(prop_wide, 
            by = c("utm_square" = "square_id",
                   "atlas_period" = "atlas_period")) %>%
  left_join(atlas_squares_df,
            by = c("utm_square"))

# check new BCR variable
table(all_data_tmp2$bcr_name, useNA = "ifany")

# check habitat variables
summary(all_data_tmp2$prop_water)
table(is.na(all_data_tmp2$prop_water))

##### Prep dataframe for modelling and PCA --------------------

# For Atlas-3, total hours is the yearly (uncorrected) total duration

# Filtering and processing steps
# 1) Drop NAs in modelling variables
# 2) Drop records with total hours = 0
# 3) Drop species with few detections across all BCRs and years
all_data_forPCA <- all_data_tmp2 %>%
  # drop NAs in variables used in modeling
  drop_na(detected, scientific_name, utm_square, collector_number,
          total_hours, survey_year, bcr_name, atlas_period, starts_with("prop_")) %>%
  # drop zero-effort records
  filter(!total_hours == 0) %>%
  # drop species with few detections across all BCRs and years
  group_by(scientific_name) %>%
  mutate(n_total_detections = sum(detected)) %>%
  ungroup() %>%
  filter(n_total_detections >= min_detections)
  
# Run PCA on habitat variables ---------------

# PCA on raw proportions is problematic (bounded 0-1, skewed, parts sum to 1 so
# the covariance matrix is singular). Fix: log-ratio (ILR) transform first.
# ILR can't handle 0s, so add 0.5/n to every proportion (n = pixels in that
# square-period), e.g. 0 of 900 pixels becomes half a pixel out of 900.
# a zero out of a small n gets a larger nudge than a zero out of a large n

# get a dataframe with only one row per square per atlas period
square_period_data <- all_data_forPCA %>%
  distinct(utm_square, atlas_period, .keep_all = T)

habitat_mat <- square_period_data %>%
  mutate(across(starts_with("prop_"), ~ .x + 0.5/year_total_area)) %>%
  dplyr::select(starts_with("prop_")) %>%
  as.data.frame()
# habitat_mat has just the proportion columns with 0s replaced

# ILR transformation
habitat_acomp <- compositions::acomp(habitat_mat) # tells R this is compositional data
habitat_ilr <- compositions::ilr(habitat_acomp)

# PCA
pca_habitat <- prcomp(habitat_ilr, center = T, scale. = F)

# variance explained per axis
summary(pca_habitat)
screeplot(pca_habitat, type = "lines")
pca_habitat$rotation # loadings 

# decide how many PCs to keep
# cumulative variance threshold (e.g., 80-90%)
pca_var      <- pca_habitat$sdev^2
pca_prop_var <- pca_var / sum(pca_var)
pca_cum_var  <- cumsum(pca_prop_var)
pca_prop_var
pca_cum_var  # first 4 PCs explained ~81% when last checked

if (pca_cum_var[n_pcs] < min_var) {
  warning("First ", n_pcs, " PCs explain only ",
          round(100 * pca_cum_var[n_pcs]), "% of habitat variance.")
}

# extract PC scores as new predictors
pc_scores <- pca_habitat$x[, seq_len(n_pcs)] %>%
  as.data.frame() %>%
  setNames(paste0("pc", seq_len(n_pcs))) %>%
  bind_cols(select(square_period_data, utm_square, atlas_period))

# Join the PC scores to the detection data -----------------

# join PC scores to detection data by UTM square and atlas period
all_data_model <- all_data_forPCA %>%
  left_join(pc_scores, by = c("utm_square", "atlas_period"))

# Calculate conversion factor for Atlas-3 effort times -----------

# pattern to match ARU recordings in Remarks field
aru_pat <- paste0(
  "(",
  paste(c(
    "\\baru\\b",
    "aru_",
    "song\\s*meter", "songmeter",
    "wildlife\\s*acoustics",
    "audiomoth",
    "autonomous\\s*(record(ing)?|recorder)",
    "\\bwav\\b", "\\bmp3\\b"
  ), collapse = "|"),
  ")"
)

# ratio of time in square (end - start time) to summed effort, per atlasser-square-day 
ratios <- daily_times %>%
  # drop rows where one or the other time measure is 0
  filter(daily_effort != 0, !daily_start_to_end == 0) %>%
  # replace daily efforts = 0 with 1/60
  mutate(daily_effort = ifelse(daily_effort == 0, 
                               1/60,
                               daily_effort),
         is_aru = grepl(aru_pat, all_remarks, ignore.case = TRUE),
         ratio = daily_start_to_end / daily_effort) %>%
  # drop ARUs and ratios < 1 (end - start is less than daily effort)
  filter(!is_aru, ratio >= 1)

# drop the highest 10% of ratios (these seem inaccurate)
ratio_cutoff <- quantile(ratios$ratio, 1 - ratio_trim, na.rm = TRUE)
final_ratios <- filter(ratios, ratio <= ratio_cutoff)

summary(final_ratios$ratio) 
ggplot(final_ratios, aes(ratio)) +
  geom_histogram()
conversion_factor <- mean(final_ratios$ratio)

# Filter out regions where species don't breed ---------

# Based on safe dates: if a species has no safe/shoulder dates in a region, 
# drop all of that species' data from the region

# check if any missing
safe_dates %>%
  distinct(english_name, ScientificName) %>%
  filter(is.na(ScientificName))
# all fine, not species to be analyzed

# A species breeds in a region if its minimum level there is < 3
species_breeding_by_region <- safe_dates %>%
  group_by(english_name, ScientificName, BCR) %>%
  summarize(min_region_level = min(level), .groups = "drop") %>%
  mutate(breeds_in_region = min_region_level < 3) %>%
  rename(common_name = english_name,
         scientific_name = ScientificName) %>%
  # recode BCRs to match Atlas data
  # BCR 8E = 77, Boreal Softwood Shield East
  # BCR 8W = 76, Boreal Softwood Shield West
  # BCR 7(H) = 74, Hudson Plains
  # BCR 12 and 13 same
  mutate(bcr_number = case_when(
    BCR == "8e" ~ "77",
    BCR == "8w" ~ "76",
    BCR == "7" ~ "74",
    T ~ BCR
  ))
  
# add biological regions to Atlas data, then join breeding info
all_data_model <- all_data_model %>%
  # join breeding regions by BCR and species
  left_join(species_breeding_by_region, by = c("bcr_number", "scientific_name"))

# View species-BCR combinations with no breeding info (to be dropped)
to_drop <- all_data_model %>%
  filter(!breeds_in_region | is.na(breeds_in_region)) %>%
  distinct(scientific_name, bcr_number) 
print(to_drop, n = 100)

# Drop rows in regions where species doesn't breed
all_data_model_breeding <- all_data_model %>%
  filter(breeds_in_region)

# Round Atlas-2 and 3 effort and log-transform effort ------------

# 1) Correct Atlas-3 hours
# 2) If total hours < 1, round up to one, otherwise round to nearest whole number
# 3) Log-transform rounded effort
all_data_model_sub <- all_data_model_breeding %>%
  # correct Atlas-3 hours
  mutate(total_hours_corrected = if_else(atlas_period == "Atlas 3",
                                                total_hours * conversion_factor,
                                                total_hours),
         # then round
         total_hours_rounded = if_else(total_hours_corrected < 1,
                                       1,
                                       round(total_hours_corrected, 0)),
         # then take log
         log_effort_rounded = log(total_hours_rounded))

# Save the final prepped dataframe and conversion factor ----------------

saveRDS(all_data_model_sub, 
        file = file.path(prepped_data_dir, "final_PObs_dataset.rds"))
