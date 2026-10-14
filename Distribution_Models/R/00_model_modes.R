# ============================================================
# 00_model_modes.R
#
# Model-mode definitions and the modelability rule, shared by
#   06_filter_and_finalize_surveys.R  (builds the species x mode modelability table)
#   07_fit_models_and_predict.R       (fits the mode selected by `survey_set`)
# Edit modes and thresholds HERE only, so the table in 06 and the gate in 07
# cannot drift apart. After editing, rerun 06 before 07.
#
# Each mode carries:
#   model_name        -> output folder names
#   fit_function      -> NAME of the fit function (resolved with match.fun() in 07,
#                        so this file does not need the fit functions loaded)
#   survey_types      -> survey types kept per species
#   survey_durations  -> allowed PC/ARU durations (checklists/special exempt)
#   include_site / include_square -> the two structural iid switches
#   special_*         -> SS modes only (see below)
# ============================================================

# ------------------------------------------------------------
# Modelability rule
# ------------------------------------------------------------
# A species is modelable in a mode iff it is detected in >= MIN_SQUARES_PER_ATLAS
# distinct atlas squares within at least ONE atlas, counted on exactly the rows
# that mode fits (its survey types, after the special-survey thresholds below and
# after counts in no-safe-date regions are zeroed in 06). No detection-count
# threshold. See assess_modelability(). In SS modes this is judged on the
# COMBINED survey set, whether or not any special type qualifies.
MIN_SQUARES_PER_ATLAS <- 20

# ------------------------------------------------------------
# Special (targeted single-protocol) surveys -- SS modes only
# ------------------------------------------------------------
# Each type gets its OWN intercept and NO effort term. A type is analysed only
# where the species is detected in >= its square threshold (squares pooled
# across atlases); types below threshold have their rows dropped. If NO special
# type qualifies, the SS fit carries no special terms and is identical to the
# corresponding non-SS mode. See select_mode_rows() / fit_PC_ARU_CL().
#
# Must be a subset of SPECIAL_SURVEY_TYPES retained by 06 (checked there).
MODEL_SPECIAL_SURVEY_TYPES <- c(
  "Eastern Screech-Owl Survey", "Long-eared Owl Survey", "Marshbird Survey",
  "Nightjar Survey Protocol",   "Central Ontario Owls",  "Northern Hawk Owl Survey"
)

# One threshold (atlas squares with a detection) for every special type, or a
# named vector to override per type, e.g.
#   c(default = 10, "Marshbird Survey" = 20, "Nightjar Survey Protocol" = 15)
MIN_SPECIAL_SQUARES <- 10

# FALSE: SS modes fit every modelable species; where no special type qualifies,
#        the fit is identical to the corresponding non-SS mode (duplicate compute).
# TRUE:  skip those species in SS modes (fit them in the non-SS mode instead).
# Read by 06 (species_modelability$will_fit) and 07 (the skip itself).
SS_SKIP_IF_NO_SPECIAL <- FALSE

# special_share_tod / special_share_doy control, PER TYPE, whether a special
# survey shares the diurnal TOD / seasonal DOY smooths. FALSE keeps a type out
# (timing/season absorbed by its intercept) -- the safe default for the
# deep-night owl protocols. Named TRUE entries opt individual types in.
SPECIAL_SHARE_TOD <- c("Northern Hawk Owl Survey" = TRUE,
                       "Marshbird Survey"         = TRUE)
SPECIAL_SHARE_DOY <- FALSE

# ------------------------------------------------------------
# Modes
# ------------------------------------------------------------
mode_settings <- list(

  # ---- Point counts + ARUs (1, 3, 5 min) ----
  PC_ARU = list(
    model_name       = "PC_ARU",
    fit_function     = "fit_PC_ARU",
    survey_types     = c("Point_Count", "ARU"),
    survey_durations = c(1, 3, 5),
    include_site     = TRUE,
    include_square   = TRUE
  ),

  # ---- Point counts + ARUs + checklists (BBA stationary + linear transects) ----
  # Checklists are OBBA3-only. One shared observation model and error family
  # across all surveys; linear transects keep their own effort terms.
  PC_ARU_CL = list(
    model_name       = "PC_ARU_CL",
    fit_function     = "fit_PC_ARU_CL",
    survey_types     = c("Point_Count", "ARU",
                         "Breeding Bird Atlas", "Linear transect"),
    survey_durations = c(1, 3, 5),
    include_site     = TRUE,
    include_square   = TRUE
  ),

  # ---- PC + ARU + special single-protocol surveys ----
  PC_ARU_SS = list(
    model_name           = "PC_ARU_SS",
    fit_function         = "fit_PC_ARU_CL",
    survey_types         = c("Point_Count", "ARU",
                             MODEL_SPECIAL_SURVEY_TYPES),
    survey_durations     = c(1, 3, 5),
    include_site         = TRUE,
    include_square       = TRUE,
    special_survey_types = MODEL_SPECIAL_SURVEY_TYPES,
    min_special_squares  = MIN_SPECIAL_SQUARES,
    special_share_tod    = SPECIAL_SHARE_TOD,
    special_share_doy    = SPECIAL_SHARE_DOY
  ),

  # ---- PC + ARU + checklists + special single-protocol surveys ----
  PC_ARU_CL_SS = list(
    model_name           = "PC_ARU_CL_SS",
    fit_function         = "fit_PC_ARU_CL",
    survey_types         = c("Point_Count", "ARU",
                             "Breeding Bird Atlas", "Linear transect",
                             MODEL_SPECIAL_SURVEY_TYPES),
    survey_durations     = c(1, 3, 5),
    include_site         = TRUE,
    include_square       = TRUE,
    special_survey_types = MODEL_SPECIAL_SURVEY_TYPES,
    min_special_squares  = MIN_SPECIAL_SQUARES,
    special_share_tod    = SPECIAL_SHARE_TOD,
    special_share_doy    = SPECIAL_SHARE_DOY
  )
)

# Mode names must match their model_name (06's table is keyed on the name, and
# 07's output folders on model_name).
stopifnot(identical(names(mode_settings),
                    vapply(mode_settings, `[[`, character(1), "model_name",
                           USE.NAMES = FALSE)))
