library(tidyverse)
library(sf)
eccc_core <-
  "c:/users/hoped/OneDrive - EC-EC/Avian_Core_20251124.xlsx" |>
  readxl::read_excel()

atlas_groups <- read_csv("Data/species groupings.csv") |> 
  pivot_longer(cols = -c(1:3), names_to = "Group", values_to = "ingroup") |> 
  janitor::clean_names() |> 
  filter(!is.na(ingroup)) |> 
  mutate(species_name_model = str_to_lower(str_replace_all(english_common_name,
                                                           "\\'|-|\\s", "_")))
paired_sum <- read_rds(
  "C:/Users/HopeD/OneDrive - EC-EC/Projects/Atlas-3/Ontario_BBA3_draft5/data_clean/model_output/paired_summaries/paired_summaries.rds"
) |>
  transpose() |>
  pluck("shared_change_summary") |>
  list_rbind() |>
  mutate(model_name = "Paired Sites - Point Counts", sp_english = sp_name)

paired_sum_aru <-
  read_rds(
    "C:/Users/HopeD/OneDrive - EC-EC/Projects/Atlas-3/Ontario_BBA3_draft5/data_clean/model_output/paired_summaries/paired_summaries_with_ARU.rds") |> 
  # "C:/Users/HopeD/Documents/TMP_A3_SHINY_DEP/RAW_DATA/ATLAS_RAW_DATA/Ontario_BBA3_draft5/data_clean/model_output/paired_summaries/paired_summaries_with_ARU.rds"
  # ) |>
  transpose() |>
  pluck("shared_change_summary") |>
  list_rbind() |>
  mutate(model_name = "Paired Sites - PC_ARU", sp_english = sp_name)


paired <- bind_rows(paired_sum, paired_sum_aru) |>
  mutate(
    pct_change_median = pct_change_q50,
    pct_change_qlow = pct_change_q10,
    pct_change_qhigh = pct_change_q90,
    ci_level = 0.9) 

sg <- distinct(
  # soc_groups[soc_groups$include...4 == 'Y', ],
  atlas_groups,
  english_common_name,
  group, species_name_model
) |>
  mutate(
    sp_english = case_when(
      english_common_name == "Herring Gull" ~ "American Herring Gull",
      english_common_name == "Warbling Vireo" ~ "Eastern Warbling Vireo",
      english_common_name ==
        "Golden-winged x Blue-winged Warbler" ~ "Golden-winged/Blue-winged Warbler",
      english_common_name == "Whimbrel" ~ "Hudsonian Whimbrel",
      english_common_name == "House Wren" ~ "Northern House Wren",
      english_common_name == "Yellow Warbler" ~ "Northern Yellow Warbler",
      str_detect(english_common_name, "Rock Pigeon") ~ "Rock Pigeon (Feral Pigeon)",
      str_detect(english_common_name, "Redpoll") ~ "Redpoll",
      TRUE ~ english_common_name
    )
  )

d_change_models <-
  list.files(
    "C:/Users/HopeD/OneDrive - EC-EC/Projects/Atlas-3/Ontario_BBA3_draft5/data_clean/model_output/",
    "regional_change_estimates\\w+.csv",
    recursive = T,
    full.names = T
  ) |>
  str_subset("_PC_ARU_CL_SS", negate = F) |>
  map(read_csv) |>
  list_rbind() |> 
  # filter(model_name != "PC_ARU_CL") |> 
  mutate(pct_change_se=(pct_change_qhigh-pct_change_mean)/qnorm(ci_level)) |> 
  filter(sufficient_data)

sp_include <- distinct(d_change_models, sp_english, Region_Name,sufficient_data)

# mod_ss <- d_change_models[d_change_models$model_name=="PC_ARU_CL_SS",]
# d_change_models <- 
#   bind_rows(
#     mod_ss,
#     filter(d_change_models, model_name=="PC_ARU" & sp_english %notin% mod_ss$sp_english )
#   )|> mutate(pct_change_se=(pct_change_qhigh-pct_change_mean)/qnorm(ci_level))



d_change <-
  d_change_models |>
  bind_rows(
    paired |>
      left_join(
        distinct(d_change_models, Region_Number, Region_Name),
        by = join_by(BCR == Region_Number)
      ) |>
      left_join(sp_include, by =join_by(sp_english, Region_Name)) |> 
      filter(!is.na(sufficient_data) )|> 
      select(any_of(names(d_change_models)))
  ) |>
  mutate(
    fc_region = factor(
      Region_Name,
      levels = c(
        "Province-wide",
        "Lower Great Lakes / St. Lawrence Plain",
        "Temperate Mixed",
        "Boreal Shield West",
        "Boreal Shield East",
        "Hudson Plains"
      )),
      model_name = str_replace(model_name, "PC_ARU_CL_SS", "Spatial model"),
    
    change_group = case_when(
      (pct_change_qlow<0&pct_change_qhigh<0)~"Decline",
      ( pct_change_qlow<0&pct_change_qhigh>0)~"Increase",
      TRUE~"Uncertain or no change"
    )
  )

d_change_gr <- d_change |>
  left_join(
    sg,
    by = join_by(sp_english)
  ) |>
  replace_na(list(group = "None")) |> 
  mutate(english_common_name=ifelse(is.na(english_common_name), sp_english, english_common_name    ))
