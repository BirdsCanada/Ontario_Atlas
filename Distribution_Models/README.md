# Distribution_Models

R workflow that fits joint OBBA2 / OBBA3 species distribution models for the
Ontario Breeding Bird Atlas, and turns them into relative-abundance maps,
change maps, regional change estimates and per-species assessment PDFs.

Models are fit with INLA / inlabru to point counts, ARU recordings, atlas
checklists and (optionally) targeted special surveys, using covariates derived
from land cover, roads and hydrology.

This folder is self-contained. **The repository holds code only.** Raw data and
everything the workflow produces live outside version control (see
[Data](#data)).

---

## Folder layout

```
Distribution_Models/
├── Distribution_Models.Rproj   open this; it sets the project root
├── README.md
├── R/
│   ├── 00_config_paths.R       paths + raw-input registry (sourced by every script)
│   ├── 00_model_modes.R        model modes and modelability rule (06, 07)
│   ├── 01_ ... 09_*.R          the pipeline, run in order
│   └── functions/              shared helpers sourced by the scripts
└── data_clean/                 ALL outputs (created on first run; git-ignored)
```

`here::here()` resolves paths relative to the folder holding
`Distribution_Models.Rproj`. Always work from that project, not from a project
at the repository root, or the scripts will stop with
`R/functions not found under ...`.

---

## Setup

### 1. R packages

CRAN:

```r
install.packages(c(
  "here", "dplyr", "tidyr", "purrr", "readr", "readxl", "stringr", "tibble",
  "lubridate", "sf", "terra", "exactextractr", "units", "smoothr", "dbscan",
  "fmesher", "inlabru", "matrixStats", "suncalc", "lutz",
  "rnaturalearth", "ggplot2", "ggspatial", "ggtext", "ggforce", "patchwork",
  "scales", "viridis", "colorspace", "RColorBrewer", "ragg", "magick", "remotes"
))
```

Not on CRAN:

```r
# INLA
install.packages("INLA",
  repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"),
  dep = TRUE)

# NatureCounts (Birds Canada) and high-resolution Natural Earth data
remotes::install_github("BirdsCanada/naturecounts")
remotes::install_github("ropensci/rnaturalearthhires")
```

### 2. Environment variables (`.Renviron`)

Two settings are specific to each user and are read from environment
variables, so no personal paths or usernames are stored in the code.

| Variable | Used by | Value |
|---|---|---|
| `ONTARIO_ATLAS_DATA` | 01, 04, 05, 06 | Absolute path to your raw-data folder (structure below) |
| `NATURECOUNTS_USER`  | 02 | Your NatureCounts username (needed only when a download is not already cached) |

Create a file called `.Renviron` in this folder (next to the `.Rproj`):

```
ONTARIO_ATLAS_DATA=C:/Users/you/Landbird-Distribution-Modeling-ECCC/Data
NATURECOUNTS_USER=your_naturecounts_username
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

## Data

### Raw inputs (`ONTARIO_ATLAS_DATA`)

Raw inputs are not distributed with this repository. Each user keeps them in a
folder of their choice and points `ONTARIO_ATLAS_DATA` at it. Inside that folder
the files must be arranged exactly as follows (shapefiles need at least their
`.shp`, `.shx` and `.dbf`; keep `.prj` too):

```
<ONTARIO_ATLAS_DATA>/
├── Spatial/
│   ├── National/
│   │   ├── BCR_2026/bcr2026_statprov.shp                                  01
│   │   ├── AtlasSquares/NationalSquares_FINAL.shp                         05
│   │   ├── LandCoverCanada2010/landcover-2010-classification.tif          05
│   │   └── LandCoverCanada2020/landcover-2020-classification.tif          05
│   └── Provincial/
│       └── Ontario/
│           ├── RoadNetwork/
│           │   ├── 2005/grnf035r05a_e.shp                                 04
│           │   └── 2025/lrnf000r25a_e.shp                                 04
│           ├── Ontario_Hydro_Network_(OHN)_-_Waterbody/
│           │   └── Ontario_Hydro_Network_(OHN)_-_Waterbody.shp            04
│           ├── Ontario_Hydro_Network_(OHN)_-_Watercourse/
│           │   └── Ontario_Hydro_Network_(OHN)_-_Watercourse.shp          04
│           └── Ontario_Atlas_biol_regions/Atlas_biol_regions.shp          06
└── Bird_Data_Raw/
    └── OBBA3/OnAtlasSafeDates_update_for_NatureCounts_2026-07-20.xlsx     06
```

The number on the right is the script that reads each file.

| Input | Description |
|---|---|
| BCR_2026 | Bird Conservation Regions intersected with provinces/territories; also supplies the Ontario boundary |
| AtlasSquares | National 10-km atlas squares (`square_id`) |
| LandCoverCanada2010 / 2020 | Land Cover of Canada, 30 m classification (2010 for OBBA2, 2020 for OBBA3) |
| RoadNetwork 2005 / 2025 | Statistics Canada road network files |
| OHN Waterbody / Watercourse | Ontario Hydro Network |
| Atlas_biol_regions | OBBA biological regions (used with safe dates) |
| OnAtlasSafeDates_... | Species breeding safe dates by biological region |

All of these are listed once, in `RAW_INPUTS` in `R/00_config_paths.R`. Scripts
obtain full paths with `require_inputs()`, which stops immediately, listing
every missing file, if `ONTARIO_ATLAS_DATA` is unset, points to a folder that
does not exist, or is missing any required file. Script 01 checks **all** raw
inputs before doing anything, so a misplaced data folder is caught at the start
of a run. To add or rename a raw input, edit `RAW_INPUTS` only.

### Data downloaded at run time

- **NatureCounts** (02): `ONATLAS3PC`, `OBBA2PC` (point counts),
  `ONATLAS3BE_DO` (checklists), `ONATLAS3SS`, `ONOWLS`, `NIGHTJAR`, `MMPBIRDS`
  (special surveys). Each download is cached in
  `data_clean/surveys/nc_downloads/` and reused on later runs; set
  `REFRESH_DOWNLOADS <- TRUE` in 02, or delete a cached file, to pull fresh data.
  Downloads can take hours.
- **Natural Earth** (04): coastline (via `rnaturalearthhires`) and Great Lakes
  polygons (`rnaturalearth::ne_download()`, needs internet access).

### Outputs (`data_clean/`)

Everything the workflow writes goes to `data_clean/` inside this folder,
created on first run:

```
data_clean/
├── spatial/       study area, prediction grid, covariate rasters (01, 03, 04, 06)
├── surveys/       raw harmonized surveys, count matrix, NatureCounts cache (02)
├── metadata/      species list (02)
├── birds/         analysis-ready data and per-species records (05, 06)
└── model_output/  predictions, summaries, change estimates, rasters, PDFs (07-09)
```

`data_clean/` is git-ignored and **must never be committed**: it contains
precise survey locations and NatureCounts data that may be subject to
data-access agreements. To share results, share
specific products (e.g. the PDFs or change-estimate tables) through other
channels.

---

## Workflow

Run the numbered scripts in order, each in a fresh R session (each starts with
`rm(list = ls())`). Steps 01–06 build the data and only need rerunning when
inputs change; 07–09 are run per model.

| Script | Does | Main outputs |
|---|---|---|
| `01_define_study_area.R` | Ontario boundary (+ Akimiski Island), 5-km buffer, simplified BCR polygons | `spatial/study_area.rds` |
| `02_load_format_all_data.R` | Downloads (or loads cached) NatureCounts data; harmonizes OBBA2/OBBA3 point counts, ARUs, checklists and special surveys; builds survey table and survey × species count matrix | `surveys/surveys_raw.rds`, `surveys/count_matrix_raw.rds`, `metadata/species_list.rds` |
| `03_define_prediction_grid.R` | 1-km prediction grid (EPSG:3978) | `spatial/prediction_grid_1000_m.rds` |
| `04_process_vectors_as_rasters.R` | Buffers roads, rivers, lakes, coastline and Great Lakes and rasterizes them at 30 m | `spatial/*.tif` |
| `05_extract_covariates.R` | Extracts covariates to grid cells and 1-km survey footprints; assigns surveys to atlas squares | `birds/analysis_data_covariates.rds` |
| `06_filter_and_finalize_surveys.R` | Filters surveys, applies safe dates, assesses modelability per species and mode, writes per-species records and covariate stacks | `birds/data_ready_for_analysis.rds`, `birds/species_data/`, `birds/species_modelability*.csv` |
| `07_fit_models_and_predict.R` | Fits the joint OBBA2/OBBA3 model for each modelable species and predicts on the 1-km grid | `model_output/predictions_<model>/`, `summaries_<model>/`, `data_used_<model>/` |
| `07b_fit_models_to_shared_footprint.R` | *Optional.* Paired analysis: change estimated from point counts in the footprint shared by both atlases; independent of 07 | `model_output/paired_summaries/` |
| `08_export_regional_change_estimates.R` | Province-wide and per-BCR posterior change estimates (and 07b's paired estimates, if present) | `model_output/change_estimates_<model>/` |
| `09_generate_model_products.R` | Map rasters and multi-page species assessment PDFs; reads 08's table | `model_output/rasters_<model>/`, `figures_<model>/` |

Shared helpers are in `R/functions/`; each script sources only the helper files
it uses.

### Choosing a model

Model modes are defined in `R/00_model_modes.R`:

| Mode | Survey types |
|---|---|
| `PC_ARU` | point counts + ARUs |
| `PC_ARU_CL` | + atlas checklists |
| `PC_ARU_SS` | + special surveys (own intercepts) |
| `PC_ARU_CL_SS` | + checklists and special surveys |

- Set `survey_set` at the top of **07** to the mode to fit.
- Set `model_name` at the top of **08** and **09** to the same value.
- `species_test` in 07 restricts a run to a few species for testing
  (`NULL` = all modelable species).

A species is fit in a mode if it is detected in at least
`MIN_SQUARES_PER_ATLAS` atlas squares in at least one atlas, counted on the
surveys that mode uses. Modes and thresholds are edited only in
`00_model_modes.R`; rerun 06 after changing them so its modelability table
matches what 07 fits.

### Reproducibility

Random steps use fixed seeds (06: `set.seed(123)` for the 25-km hex grid; 07: `prediction_seed`
for posterior prediction draws). NatureCounts data change as records are
added; the cached downloads in `data_clean/surveys/nc_downloads/` define the
data a given run used.
