# 01_load_actigraphy.R
# Step 1: load the actigraphy CSV, check its structure, and print
# non-identifying diagnostics. Run from the repo root: source('R/01_load_actigraphy.R')
#
# Privacy rules for everything printed here:
#   - no IDs, no dates, no tables of sex/race/ethnicity
#   - any count between 1 and (minCell - 1) is shown as "<minCell"

source('config.R')

# ---- helpers ----------------------------------------------------------------
suppress <- function(n) ifelse(n > 0 & n < minCell, paste0('<', minCell), as.character(n))

show_counts <- function(x, label) {
  tab <- table(x, useNA = 'ifany')
  out <- data.frame(value = names(tab), n = suppress(as.integer(tab)))
  cat('\n', label, '\n', sep = '')
  print(out, row.names = FALSE)
}

show_quantiles <- function(x, label, digits = 2) {
  x <- x[!is.na(x)]
  cat('\n', label, ' (n = ', length(x), ')\n', sep = '')
  if (length(x) == 0) { cat('  no non-missing values\n'); return(invisible()) }
  print(round(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), digits))
}

to_num <- function(x) suppressWarnings(as.numeric(x))

# ---- load -------------------------------------------------------------------
if (!file.exists(actigraphyFile)) {
  csvs <- list.files(dataDir, pattern = '\\.csv$', ignore.case = TRUE)
  stop('Actigraphy file not found: ', actigraphyFile,
       '\nCSV files in dataDir: ', paste(csvs, collapse = ', '),
       '\nEdit actigraphyFile in config.R.', call. = FALSE)
}
# check.names = FALSE keeps the original headers (e.g. "datetime_vrrsb+srs_stop")
act <- read.csv(actigraphyFile, stringsAsFactors = FALSE, check.names = FALSE)

# ---- 1. expected columns ----------------------------------------------------
expectedCols <- c(
  'List', 'study_id', 'record_id', 'sex', 'ethnicity', 'race', 'kidokoro', 'timepoint',
  'wear_time_ave', 'total_mvt_time_ave', 'l_time_ave', 'l_magnitude_ave', 'r_time_ave',
  'r_magnitude_ave', 'l_magnitude_sd_ave', 'r_magnitude_sd_ave', 'bilateral_magnitude_ave',
  'l_jerk_avg_ave', 'r_jerk_avg_ave', 'l_sd_freq_ave', 'r_sd_freq_ave',
  'simultaneous_time_ave', 'agemonths', 'ageaccel', 'datetime_hx_start', 'date',
  'datetime_vrrsb+srs_stop', 'acceldatetime', 'acceldate', 'agemanualaccel'
)
cat('=== 1. Columns ===\n')
cat('Missing from file:', if (length(setdiff(expectedCols, names(act)))) paste(setdiff(expectedCols, names(act)), collapse = ', ') else 'none', '\n')
cat('Unexpected in file:', if (length(setdiff(names(act), expectedCols))) paste(setdiff(names(act), expectedCols), collapse = ', ') else 'none', '\n')
stopifnot('study_id' %in% names(act))

# ---- 2. rows and children ---------------------------------------------------
cat('\n=== 2. Rows and children ===\n')
cat('Collections (rows):', nrow(act), '\n')
cat('Children (unique study_id):', length(unique(act$study_id)), '\n')
rowsPerChild <- as.integer(table(act$study_id))
show_counts(rowsPerChild, 'Children by number of collections (rows per child)')

if ('record_id' %in% names(act)) {
  cat('\nDuplicated record_id values:', sum(duplicated(act$record_id) & !is.na(act$record_id)), '\n')
}

# ---- 3. timepoint -----------------------------------------------------------
if ('timepoint' %in% names(act)) {
  cat('\n=== 3. Timepoint ===')
  show_counts(act$timepoint, 'Collections by timepoint')
  # does a child ever appear twice at the same timepoint?
  dupTP <- sum(duplicated(act[, c('study_id', 'timepoint')]))
  cat('\nRows repeating a study_id + timepoint combination:', dupTP, '\n')
}

# ---- 4. age at collection ---------------------------------------------------
cat('\n=== 4. Age at collection ===\n')
ageDb <- to_num(act$ageaccel)
ageManual <- to_num(act$agemanualaccel)
show_quantiles(ageDb, 'ageaccel (database)')
show_quantiles(ageManual, 'agemanualaccel (manual)')

# agemanualaccel was described as a difference in days; if the values look like
# days (median well above any plausible age in months), convert for comparison.
manualMedian <- median(ageManual, na.rm = TRUE)
if (!is.na(manualMedian) && manualMedian > 150) {
  cat('\nNOTE: agemanualaccel median is ', round(manualMedian), ' -> looks like DAYS. ',
      'Converting to months (/30.4375) for the comparison below.\n', sep = '')
  ageManualMonths <- ageManual / 30.4375
} else {
  ageManualMonths <- ageManual
}
show_counts(floor(ageManualMonths), 'Collections by whole month of age (agemanualaccel)')

ageDiff <- ageDb - ageManualMonths
show_quantiles(ageDiff, 'ageaccel minus agemanualaccel (months)')
bothPresent <- sum(!is.na(ageDiff))
cat('\nCollections with both ages present:', bothPresent, '\n')
cat('Differing by more than 1 month:', suppress(sum(abs(ageDiff) > 1, na.rm = TRUE)), '\n')

# ---- 5. wear time -----------------------------------------------------------
if ('wear_time_ave' %in% names(act)) {
  cat('\n=== 5. Wear time ===')
  show_quantiles(to_num(act$wear_time_ave), 'wear_time_ave (units as stored)')
}

# ---- 6. missingness ---------------------------------------------------------
cat('\n=== 6. Missingness per column (% of rows blank or NA) ===\n')
isMissing <- function(x) is.na(x) | trimws(as.character(x)) == ''
miss <- data.frame(
  column = names(act),
  pct_missing = round(100 * sapply(act, function(x) mean(isMissing(x))), 1),
  row.names = NULL
)
print(miss, row.names = FALSE)

cat('\nDone. Nothing was written to disk.\n')
