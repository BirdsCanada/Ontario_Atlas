# Functions for PObs analysis #
# Author: Sarah Mueller
# Last updated: Nov 12, 2026

# Functions:
# 1) get_species_data returns a subset of data for one species for the analysis
# 2) fit_species_model_auto fits model with appropriate terms for a species
# 3) predict_20h genenerates model predictions at 20 h effort
#       returns list(prediction table, ggplot, model)
# 4) save_species_plot saves model and species prediction plot 
#       returns prediction table


# Get the subset of data for a species --------------

# Filters applied for data subset, in order:
#   1. drop outliers in total effort hours ( > mean + 4 SD)
#   2. drop BCRs with 0 detections
#   3. drop Atlas periods with 0 detections
#   4. drop BCR-Atlas period combinations with fewer than 2 detections
# Returns NULL (with a warning) if the species is not in the dataset
get_species_data <- function(species_name, 
                             df = all_data_model_sub) {
  
  # Return NULL with warning if species is not in dataset
  if(!species_name %in% df$scientific_name){
    warning("Species name as entered not in Atlas dataset: ", species_name)
    return(NULL)
  } 
  
  # get data for only the focal species
  species_df <- df %>%
    filter(scientific_name == species_name)
  
  # if the species is warbling vireo, drop BCR 74
  if(species_name == "Vireo gilvus"){
    species_df <- species_df %>%
      filter(!bcr_number == "74")
  }
  
  # apply data filters
  species_df %>%
    # drop square-atlasser-year records with outlying total effort hours
    filter(total_hours_rounded <= (mean(total_hours_rounded) + 4 * sd(total_hours_rounded))) %>%
    # drop BCRs with 0 detections
    group_by(bcr_number) %>%
    filter(sum(detected) > 0) %>%
    ungroup() %>%
    # drop Atlas periods with 0 detections
    group_by(atlas_period) %>%
    filter(sum(detected) > 0) %>%
    ungroup() %>%
    # drop BCR-Atlas combinations with fewer than 2 detections
    group_by(bcr_number, atlas_period) %>%
    filter(sum(detected) > 1) %>%
    ungroup() %>%
    # factor versions of variables used in models
    mutate(atlasser_id.f = factor(collector_number),
           square_id.f = factor(utm_square),
           bcr.f = droplevels(factor(bcr_number, levels = bcr_number_order)),
           atlas_period.f = factor(atlas_period)) 
}

# Functions to fit and check model -------------

# Logistic regression

# Response: binary, species detected or not in a square-year
# Fixed effects:
#   BCR
#   Atlas period
#   BCR x atlas period
#   Log effort
#   Habitat covariates (PC 1-4)
# Random effects:
#   Atlasser ID (intercept)
#   Square (intercept)

# fit full model with BCR x Atlas period interaction
fit_species_model <- function(species_data){
  
  glmmTMB(detected ~ bcr.f*atlas_period.f + log_effort_rounded + 
            pc1 + pc2 + pc3 + pc4 +
            (1|atlasser_id.f) + (1|square_id.f),
          family = binomial,
          data = species_data)
  
}

# fit model for species present in only one BCR (no BCR effects)
fit_species_model_noBCR <- function(species_data){
  
  glmmTMB(detected ~ atlas_period.f + log_effort_rounded + 
            pc1 + pc2 + pc3 + pc4 +
            (1|atlasser_id.f) + (1|square_id.f),
          family = binomial,
          data = species_data)
  
}

# fit model for species present in only one Atlas period (no Atlas period effects)
fit_species_model_noPeriod <- function(species_data){
  
  glmmTMB(detected ~ bcr.f + log_effort_rounded + 
            pc1 + pc2 + pc3 + pc4 +
            (1|atlasser_id.f) + (1|square_id.f),
          family = binomial,
          data = species_data)
  
}

# fit model for species present in only one Atlas period and one BCR
# (no Atlas period or BCR effects)
fit_species_model_none <- function(species_data){
  
  glmmTMB(detected ~ log_effort_rounded + 
            pc1 + pc2 + pc3 + pc4 +
            (1|atlasser_id.f) + (1|square_id.f),
          family = binomial,
          data = species_data)
  
}

# pick and fit the appropriate model based on the BCRs/Atlas periods present
fit_species_model_auto <- function(species_data){

  # whether or not to incude BCR and period
  include_bcr <- n_distinct(species_data$bcr_number) > 1
  include_period <- n_distinct(species_data$atlas_period) > 1
  
  if (include_bcr & include_period) {
    # fit model with all terms
    fit_species_model(species_data)
  } else if (!include_bcr & !include_period) {
    # fit model with no BCR and no period
    fit_species_model_none(species_data)
  } else if (!include_bcr) {
    # fit model with no BCR terms
    fit_species_model_noBCR(species_data)
  } else {
    # fit model with no Atlas period terms
    fit_species_model_noPeriod(species_data)
  }
  
}

# Predict detection probability after 20 h of effort and plot --------

# predictions generated at 20 h of effort,
# mean habitat covariate values (across the species dataset),
# and for a typical atlasser and square

# returns list(prediction table, ggplot, model)
predict_20h <- function(species_name){
  
  # get species data
  species_data <- get_species_data(species_name)
  
  # if there is no data for the species, warn and return NULL
  if(is.null(species_data) || nrow(species_data) == 0){
    warning("predict_20h(): no data for this species, returning NULL")
    return(NULL)
  }
    
  # make sure any unused levels of BCR are dropped
  species_data <- species_data %>%
    mutate(bcr.f = droplevels(bcr.f))
  
  # fit appropriate model for the species
  species_mod <- fit_species_model_auto(species_data)
  
  # get terms actually in model (atlas period and/or bcr)
  variable_names <- names(fixef(species_mod)$cond)
  has_bcr <- any(str_detect(variable_names, "bcr.f"))
  has_period <- any(str_detect(variable_names, "atlas_period.f"))
  
  # summarize sample size
  species_n <- species_data %>%
    group_by(atlas_period.f, bcr.f) %>%
    summarize(n_total_obs = n(), 
              n_detections = sum(detected),
              .groups = "drop")
  
  # define set of covariate values for predictions
  covariate_values <- list(log_effort_rounded = log(20))
  
  # significant pairwise Atlas period contrasts (stays NULL if no period term)
  pair_contrasts <- NULL
  
  # function to tidy contrasts by:
  #   keeping significant ones
  #   splitting "Atlas 1 - Atlas 2" into two columns
  tidy_contrasts <- function(x){
    x %>%
      as_tibble() %>%
      filter(p.value < 0.05) %>%
      mutate(significance = "*") %>%
      separate(contrast, into = c("period1", "period2"), sep = " [-]") %>%
      mutate(period1 = str_trim(period1),
             period2 = str_trim(period2))
  }
  
  # function to add plotting group (BCR-level vs. range-wide) and order BCRs
  format_bcr <- function(x){
    x %>%
      mutate(grouping = factor(if_else(bcr.f == "Range-wide", "province", "bcr")),
             bcr.f = droplevels(factor(as.character(bcr.f),
                                       levels = c(bcr_number_order, "Range-wide"))))
  }
  
  # Species with a single BCR and Atlas period (one prediction):
  if(!has_period & !has_bcr){
    
    # just one prediction (a single BCR and a single Atlas period)
    species_predictions <- emmeans(species_mod, 
                                 ~ 1, 
                                 at = covariate_values,
                                 type = "response") %>%
      as_tibble() %>%
      rename(predicted = prob,
             conf.low = asymp.LCL,
             conf.high = asymp.UCL) %>%
      # add sample size, BCR, and atlas period
      mutate(n_total_obs = species_n$n_total_obs,
             n_detections = species_n$n_detections,
             bcr.f = species_n$bcr.f,
             atlas_period.f = species_n$atlas_period.f)
    
    # Species present in multiple periods and only one BCR:
  } else if(has_period & !has_bcr){
    
    # get predictions by Atlas period
    species_predictions <- emmeans(species_mod, 
                                 ~ atlas_period.f, 
                                 at = covariate_values,
                                 type = "response") %>%
      as_tibble() %>%
      rename(predicted = prob,
             conf.low = asymp.LCL,
             conf.high = asymp.UCL) %>%
      # join sample size by atlas period
      left_join(species_n, by = c("atlas_period.f"))
    
    # get pairwise comparisons of Atlas periods in the single BCR
    pair_contrasts <- emmeans(species_mod, pairwise ~ atlas_period.f,
                              at = covariate_values)$contrasts %>%
      tidy_contrasts() %>%
      mutate(bcr.f = as.character(species_predictions$bcr.f[1]))
  
    # Species present in multiple BCRs and only one Atlas period: 
    # (only Lark Bunting, and uncertainty very high)
  } else if(has_bcr & !has_period){
    
    # get predictions by BCR
    species_predictions <- emmeans(species_mod, 
                                 ~ bcr.f, 
                                 at = covariate_values,
                                 type = "response") %>%
      as_tibble() %>%
      rename(predicted = prob,
             conf.low = asymp.LCL,
             conf.high = asymp.UCL) %>%
      # join sample size by BCR
      left_join(species_n, by = c("bcr.f"))
   
    # Species present in multiple Atlas periods and BCRs: 
  } else {
    
    # get predictions by BCR x Atlas period
    species_predictions <- emmeans(species_mod, 
                                ~ bcr.f*atlas_period.f, 
                                 at = covariate_values,
                                 type = "response") %>%
      as_tibble() %>%
      rename(predicted = prob,
             conf.low = asymp.LCL,
             conf.high = asymp.UCL) %>%
      # join sample size by atlas period and BCR
      left_join(species_n, by = c("atlas_period.f", "bcr.f"))
    
    # area weights: share of squares in each BCR (only BCRs where species is present)
    # (named by BCR so they can be matched to the emmeans grid)
    bcr_weights <- n_squares %>%
      filter(bcr.f %in% levels(species_data$bcr.f)) %>%
      mutate(bcr_weight = n_squares / sum(n_squares)) %>%
      with(setNames(bcr_weight, as.character(bcr.f)))
    # output is a vector of the BCR weights, named by BCR number
    
    # get sample size province-wide in each Atlas period
    species_n_overall <- species_n %>%
      group_by(atlas_period.f) %>%
      summarize(n_total_obs = sum(n_total_obs),
                n_detections = sum(n_detections))
    
    # Ontario-wide predictions: area-weighted means of BCR-specific probabilities
    # get BCR-specific probabilities first, then average across those to get the area-weighted mean prob
    
    # first get logit-scale predictions
    emm_full <- emmeans(species_mod,
                        ~ bcr.f * atlas_period.f,
                        at = covariate_values)
    
    # convert logit-scale grid to probability scale
    emm_resp <- regrid(emm_full)
    
    # Empty BCR-period cells (inestimable) get weight 0, i.e. they contribute
    # a detection probability of ~0 (species not detected there much/at all)
    # row order of the grid
    grid <- as.data.frame(emm_full)[, c("bcr.f", "atlas_period.f", "emmean")] 
    # get a vector of the weights that is the same length and Atlas period-BCR order as the grid
    w <- bcr_weights[as.character(grid$bcr.f)]
    # zero out the inestimable cell(s)
    w[is.na(grid$emmean)] <- 0
    
    # one coefficient vector per Atlas period
    periods <- levels(grid$atlas_period.f)
    coefs <- lapply(periods, function(p) ifelse(grid$atlas_period.f == p, w, 0))
    names(coefs) <- periods
    # coefs is a named list storing one vector of weights per Atlas period
    
    # get Ontario-wide averages
    emm_complete <- contrast(emm_resp, method = coefs)
    
    overall_predictions <- summary(emm_complete, infer = T) %>%
      as_tibble() %>%
      rename(atlas_period.f = contrast,
             predicted = estimate,
             conf.low = asymp.LCL,
             conf.high = asymp.UCL) %>%
        # join sample size by atlas period
        left_join(species_n_overall, by = c("atlas_period.f")) %>%
        # add "Range-wide" as BCR
        mutate(bcr.f = "Range-wide")

    species_predictions <- bind_rows(species_predictions, overall_predictions)
    
    # BCR-level pairwise comparisons between Atlas periods 
    pair_bcr_contrasts <- emmeans(species_mod, 
                                  pairwise ~ atlas_period.f | bcr.f,
                                  at = covariate_values)$contrasts %>%
      as_tibble() %>%
      mutate(bcr.f = as.character(bcr.f))
    
    # Range-wide pairwise comparisons between Atlas periods (area-weighted)
    pair_range_contrasts <- summary(pairs(emm_complete, adjust = "tukey")) %>%
      as_tibble() %>%
      mutate(bcr.f = "Range-wide")
    
    # bind together and format
    pair_contrasts <- bind_rows(pair_bcr_contrasts, pair_range_contrasts) %>%
      tidy_contrasts()
    
  }
  
  # add grouping variable for plotting (bcr-level vs. range-wide)
  species_predictions <- format_bcr(species_predictions)
  if(!is.null(pair_contrasts)) pair_contrasts <- format_bcr(pair_contrasts)
  
  # dodge constant
  dodge <- 0.7
  
  ##### Bracket geometry -------
  
  # offset of each Atlas period within a dodged group
  lev <- levels(factor(species_predictions$atlas_period.f))
  n_lev <- length(lev)
  offsets <- setNames((seq_len(n_lev) - (n_lev + 1) / 2) * (dodge / n_lev), lev)
    
  # position of each BCR within its own facet panel
  bcr_lookup <- species_predictions %>%
    distinct(grouping, bcr.f) %>%
    arrange(grouping, bcr.f) %>%
    group_by(grouping) %>%
    mutate(x_num = row_number()) %>%
    ungroup()
    
  # top of the tallest error bar in any BCR
  y_range <- max(species_predictions$conf.high, na.rm = T)
  
  
  char_w <- 0.015 # width of one character as a fraction of the panel (tune)
  step <- 0.04 # spacing between stacked brackets as fraction of panel
  tick_frac <- 0.015 # bracket tick length as fraction of panel
  
  # largest number of stacked brackets in any BCR
  n_max <- if(is.null(pair_contrasts) || nrow(pair_contrasts) == 0) {
    0
  } else if(nrow(pair_contrasts) == 0){
    0
  } else {
    max(count(pair_contrasts, grouping, bcr.f)$n, na.rm = T)
  }
  max_chr <- max(nchar(paste0("n = ", species_predictions$n_total_obs)), na.rm = T)
    
  # fraction of the axis reserved for label + gap + stacked brackets + margin
  reserve <- 0.02 + char_w * max_chr + 0.03 + step * n_max + 0.03
  # final axis length in data units
  L <- y_range / (1 - reserve)
    
  # height at which each BCR's "n" labels end
  ymax_by_bcr <- species_predictions %>%
    mutate(label_end = conf.high +
             L * (0.02 + char_w * nchar(paste0("n = ", n_total_obs)))) %>%
    group_by(bcr.f, grouping) %>%
    summarize(ymax = max(label_end, na.rm = T), 
              .groups = "drop")
    
  # Only add brackets if atlas period is in the model:
  if(has_period){
    
    # add variables to give parameters of brackets
    brackets <- pair_contrasts %>%
      left_join(bcr_lookup, by = c("bcr.f", "grouping")) %>%
      left_join(ymax_by_bcr, by = c("bcr.f", "grouping")) %>%
      group_by(bcr.f, grouping) %>%
      mutate(k = row_number()) %>%
      ungroup() %>%
      mutate(x1 = x_num + offsets[period1],
             x2 = x_num + offsets[period2],
             y = ymax + L * (0.03 + step * (k - 1)),
             tick = tick_frac * L)
    
  } 
  
  # Plot ------
  
  # y axis is fixed to L when brackets are drawn, otherwise just leave room for labels
  y_scale <- if(has_period){
    scale_y_continuous(limits = c(0, L), expand = expansion(mult = c(0, 0.02)))
  } else {
    scale_y_continuous(expand = expansion(mult = c(0, 0.17)))
  }
  
  species_plot <- ggplot(species_predictions, aes(bcr.f, predicted, 
                                                   color = atlas_period.f,
                                                   fill = atlas_period.f)) +
    geom_col(width = 0.7,
             position = position_dodge(width = dodge),
             alpha = 0.7,
             na.rm = T) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                  position = position_dodge(width = dodge),
                  width = 0.15) +
    geom_text(aes(y = conf.high,
                  label = paste0("n = ", n_total_obs)),
              position = position_dodge(width = dodge),
              hjust = -0.1,
              color = "black",
              size = 3,
              show.legend = F,
              na.rm = T) +
    cowplot::theme_cowplot() +
    scale_color_manual(values = atlas_colors,
                       aesthetics = c("color", "fill")) +
    y_scale +
    labs(y = "Square-level probability of observation",
         x = "BCR",
         color = "",
         fill = "") +
    guides(color = guide_legend(reverse = T),
           fill = guide_legend(reverse = T)) +
    coord_flip() +
    facet_wrap(~grouping,
               ncol = 1,
               scales = "free_y",
               space = "free_y") +
    theme(strip.text = element_blank(),
          strip.background = element_blank())
  
  
  ##### Add brackets to plot --------
  
  # add brackets only if atlas period is in the model
  if(has_period) {
    
    # gap between bracket line and the asterisks
    asterisk_nudge <- 0.008*L
    
    species_plot <- species_plot +
      # add bracket line
      geom_segment(data = brackets, aes(x = x1, xend = x2, y = y, yend = y),
                   inherit.aes = FALSE) +
      # add ticks at each end of the bracket
      geom_segment(data = brackets, aes(x = x1, xend = x1, y = y, yend = y - tick),
                   inherit.aes = FALSE) +
      geom_segment(data = brackets, aes(x = x2, xend = x2, y = y, yend = y - tick),
                   inherit.aes = FALSE) +
      # add significance label
      geom_text(data = brackets,
                aes(x = (x1 + x2) / 2, 
                    y = y + asterisk_nudge, 
                    label = significance),
                inherit.aes = FALSE,
                angle = 0,
                hjust = 0, 
                vjust = 0.75,
                size = 5)
  }
  
  return(list(species_predictions, species_plot, species_mod))
  
}

# Save the model, prediction plot, and table for one species -----------

# saves model and prediction plot
# returns prediction table (or NULL if species can't be modelled)
save_species_plot <- function(species_name){
  
  # get species data
  species_data <- get_species_data(species_name)
  
  # skip species without data
  if(is.null(species_data) || nrow(species_data) == 0){
    message("Skipping ", species_name, ": no usable data after filtering")
    return(NULL)
  }
  
  # fit model, generate predictions, and plot
  species_pred_list <- predict_20h(species_name)
  
  # skip and return NULL for species with no predictions
  if(is.null(species_pred_list)){
    message("Skipping ", species_name, ": no predictions")
    return(NULL)
  }
  
  # get common name
  common_name <- species_data$common_name[1]
  
  # common name for file name (spaces and hyphens -> underscore, apostrophes removed)
  common_name_file <- common_name %>%
    str_replace_all("[ -]", "_") %>%
    str_remove_all("'")
  
  # flag models whose Hessian is not positive definite (unreliable SEs/CIs)
  species_model <- species_pred_list[[3]]
  if(!isTRUE(species_model$sdr$pdHess)){
    warning(species_name, ": non-positive-definite Hessian, check this model")
  }
  
  # save the model object
  save(species_model, file = file.path(model_dir, paste0(common_name_file, ".RData")))
  
  # save plot
  ggsave(species_pred_list[[2]], 
         file = file.path(final_output_dir, paste0(common_name_file, "_PObs.png")),
         bg = "white",
         height = 5.5,
         width = 7,
         dpi = 400)
  
  # make a table of PObs by BCR and Atlas period
  # change_from_prev is only filled in when the previous row is the immediately
  # preceding Atlas (e.g. Atlas 3 - Atlas 2); it is NA if that period is missing
  species_table <- species_pred_list[[1]] %>%
    mutate(common_name = common_name) %>%
    select(common_name, bcr.f, atlas_period.f, predicted, conf.low, conf.high, 
           n_total_obs, n_detections) %>%
    # calculate change from previous Atlas
    mutate(period_num = as.integer(str_extract(as.character(atlas_period.f), "\\d+"))) %>%
    arrange(bcr.f, period_num) %>%
    group_by(bcr.f) %>%
    mutate(change_from_prev = ifelse(period_num - lag(period_num) == 1,
                                     predicted - lag(predicted),
                                     NA_real_)) %>%
    ungroup() %>%
    select(-period_num) %>%
    mutate(bcr.f = as.character(bcr.f),
           atlas_period.f = as.character(atlas_period.f))
  
  return(species_table)
}
