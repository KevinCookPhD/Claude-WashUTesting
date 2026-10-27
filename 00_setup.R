# 00_setup.R
# R03 preliminary data: project setup (paths, packages, folders, helpers)
# Run once at the top of every analysis script: source("00_setup.R")

# ---- Paths -------------------------------------------------------------
# Forward slashes work on Windows and avoid backslash-escaping problems.
# Edit `base_dir` only if the project folder moves.
base_dir <- file.path(
  "C:/Users/Kevin/OneDrive - Childrens National Hospital",
  "ComputerSyncing/03_Funding/02_Drafting",
  "2025_NIH_R03_PTMED_NICUOpioidExposures",
  "A1 - November 2026/Preliminary Data"
)

dirs <- list(
  raw      = file.path(base_dir, "data_raw"),        # original files, never edited
  clean    = file.path(base_dir, "data_clean"),      # analysis-ready datasets
  code     = file.path(base_dir, "code"),
  tables   = file.path(base_dir, "output", "tables"),
  figures  = file.path(base_dir, "output", "figures"),
  logs     = file.path(base_dir, "output", "logs")
)

if (!dir.exists(base_dir)) {
  stop("Base directory not found. Check `base_dir` in 00_setup.R:\n", base_dir)
}
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# ---- Packages ----------------------------------------------------------
pkgs <- c(
  "tidyverse",   # data wrangling + ggplot2
  "janitor",     # clean_names(), tabyl()
  "gtsummary",   # Table 1 and model tables
  "gt",          # table export
  "naniar",      # missingness summaries/plots
  "pwr",         # basic power calculations
  "broom",       # tidy model output
  "readxl", "haven"  # Excel / SAS / SPSS / Stata import
)
missing_pkgs <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing_pkgs)) install.packages(missing_pkgs)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- Options -----------------------------------------------------------
set.seed(20261101)                       # reproducible simulations/bootstraps
options(scipen = 999, dplyr.summarise.inform = FALSE)
theme_set(theme_minimal(base_size = 12))

# ---- Helpers -----------------------------------------------------------
out_path <- function(type = c("tables", "figures", "logs", "clean", "raw"), file) {
  type <- match.arg(type)
  file.path(dirs[[type]], file)
}

message("Setup complete. Project root: ", base_dir)
