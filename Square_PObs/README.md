# Atlas_PObs

R workflow that fits logistic regression models for Atlas square-level 
probability of observation for the Ontario Breeding Bird Atlas. Produces 
plots of BCR-level and range-wide mean predicted probability of observation 
after 20 hours of effort in a square.

Models are fit with glmmTMB to square-level data (with each record representing
a single atlasser x square x year), using all survey types except incidental 
surveys. 

The Atlas dataset is not included in the repository, but is downloaded in the 
first script.

---

## Folder layout

```
Square_PObs/
├── Square_PObs.Rproj   open this; it sets the project root
├── README.md
├── scripts/
│   ├── 00_pobs_analysis_functions.R      functions for analysis sourced by 04
│   └── 01_ ... 04_*.R                    the pipeline, run in order
└── output/                 Atlas data, fitted model objects, prediction plots (git-ignored)
```

---

## Setup

### 1. R packages

CRAN:

```r
install.packages(c(
  "NoSleepR", "tidyverse", "janitor", "sf", "terra", 
  "exactextractr", "compositions", "glmmTMB", "emmeans", "DHARMa"
))
```

Not on CRAN:

```r
# NatureCounts (Birds Canada)
remotes::install_github("BirdsCanada/naturecounts")
```

### 2. Environment variables (`.Renviron`)

Some settings are specific to each user and are read from environment
variables, so no personal paths or usernames are stored in the code.

| Variable | Used by | Value |
|---|---|---|
| `ONTARIO_ATLAS_DATA` | 01-04 | Absolute path to your raw-data folder (structure below) |
| `NATURECOUNTS_USER`  | 01 | Your NatureCounts username |
| `NATURECOUNTS_INFO`  | 01 | A short description of your reason for downloading |
| `NATURECOUNTS_request_*`  | 01 | Your six-digit NatureCounts request ID numbers for the respective collections |
| `OUTPUT_DIR` | 01-04 | Absolute path to your output folder |

Create a file called `.Renviron` in this folder (next to the `.Rproj`):

```
ONTARIO_ATLAS_DATA=C:/Users/you/project/data/
NATURECOUNTS_USER=your_naturecounts_username
NATURECOUNTS_INFO=Short description of reason for download
NATURECOUNTS_request_3BE=******
NATURECOUNTS_request_3PC=******
NATURECOUNTS_request_3SS=******
OUTPUT_DIR=C:/Users/you/project/output/
```

Use forward slashes, no quotes, and end the file with a blank line. Restart R
(Session → Restart R) after editing it. Check with
`Sys.getenv("ONTARIO_ATLAS_DATA")`.

A project `.Renviron` replaces your user-level `~/.Renviron` for this project,
so if you keep other variables there (e.g. API keys), either copy them in or set
these two variables in `~/.Renviron` instead (`usethis::edit_r_environ()`).
`.Renviron` is git-ignored: never commit it.

The NatureCounts account must have been granted access to the atlas
collections listed below; `naturecounts` prompts for the password on the first
download.

---

### Raw inputs (`ONTARIO_ATLAS_DATA`)

Raw inputs are not distributed with this repository; download from:
https://doi.org/10.6084/m9.figshare.35063872. 

Each user keeps raw inputs in a folder of their choice and points 
`ONTARIO_ATLAS_DATA` at it. Inside that folder the files must be arranged 
exactly as follows (shapefiles need at least their `.shp`, `.shx` and `.dbf`; 
keep `.prj` too):

```
<ONTARIO_ATLAS_DATA>/   
├── atlas_squares_complete_with_BCR/  
│   └── atlas_squares_with_bcr.shp                        02, 03
├── atlas3_specieslist.csv                                03
├── landcover/                                              
│   ├── Landcover_class_values_key.csv                    02
│   ├── ON_Landcover_1985_clipped.tif                     02
│   ├── ON_Landcover_2003_clipped.tif                     02
│   └── ON_Landcover_2022_clipped.tif                     02
├── OnAtlasSafeDates_update_2026-11-18.csv                03
├── special_survey_species.csv                            03
├── species_taxonomy_lookup_table.csv                     01
└── waterfowl_aerial_survey_for_analysis_good.csv         03
```

The number on the right is the script that reads each file.

| Input | Description |
|---|---|
| atlas_squares_with_bcr | spatial Atlas squares with their assigned Bird Conservation Regions
| atlas3_specieslist | common and scientific names of all species with accounts in Atlas-3
| landcover | land cover from from GLC_FCS30D, 30 m (1985 for OBBA1, 2003 for OBBA2, 2022 for OBBA3)
| ONAtlasSafeDates_... | Species breeding safe dates by biological region
| special_survey_species | target species of Atlas special surveys
| species_taxonomy_...| lookup table to update species taxonomy between Atlases
| waterfowl_...| CWS waterfowl survey data

### Data downloaded at run time

- **NatureCounts** (01): `ONATLAS3BE_DO`, `ONATLAS3PC`, `ONATLAS3SS`,
  `OBBA2BE_RAW`, `OBBA2PC`, `OBBA1BE_RAW`.
  Downloads can take hours.
  
### Outputs (`output/`)

Everything the workflow writes goes to `output/` folder, created on first run:

```
output/
├── all_species_predictions_table.csv    model prediction table for all species (04)
├── download/                            contains datasets from OBBA1-OBBA3, downloaded from NatureCounts (01)
├── failed_species.txt                   species for which analysis failed (04)
├── final_PObs_dataset/      
│   └── final_PObs_dataset.rds  dataset for analysis (03)
├── habitat.csv                          square-level habitat covariates (02)
├── species_models/                      fitted model objects for each species (04)
└── species_predictions/                 pngs of predicted PObs (04)
```

The number on the right is the script that produces each file.

`output/` is git-ignored and **must never be committed**: it contains precise 
survey locations and NatureCounts data that may be subject to data-access 
agreements. To share results, share specific products (e.g. PNGs or prediction
tables) through other channels.

---

## Workflow

Run the numbered scripts in order, each in a fresh R session (each starts with
`rm(list = ls())`). Steps 01–03 build the data and only need rerunning when
inputs change; 04 fits models and generates outputs.

| Script | Does | Main outputs |
|---|---|---|
| `00_pobs_analysis_functions.R` | Functions used by 04 for analysis | None |
| `01_download_atlas_data.R` | Downloads and harmonizes OBBA1, OBBA2, and OBBA3 data from NatureCounts | `output/download/` |
| `02_extract_habitat_covariates.R` |  Extracts square-level habitat covariates | `output/habitat.csv` |
| `03_pobs_data_prep.R` | Generates analysis-ready dataset | `output/final_PObs_dataset/` |
| `04_pobs_analysis.R` | Fits logistic regression models and generates predictions | `output/species_models/`, `output/species_predictions`, `output/all_species_predictions_table.csv` |
