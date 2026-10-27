# 01_medication_exposure.R
# Characterize medication exposures in preemiesMedicationSummary.csv
# Input:  data_raw/preemiesMedicationSummary.csv (one row per infant, presumably)
# Output: output/tables/medication_exposure_summary.csv

source("00_setup.R")

meds_raw <- read_csv(out_path("raw", "preemiesMedicationSummary.csv"),
                     show_col_types = FALSE)

# ---- 1. Structure checks ----------------------------------------------
cat("Rows:", nrow(meds_raw), " Columns:", ncol(meds_raw), "\n")
cat("Unique IDs:", n_distinct(meds_raw$ID),
    "| duplicated IDs:", sum(duplicated(meds_raw$ID)), "\n")
glimpse(meds_raw)

non_numeric <- meds_raw |> select(-ID) |> select(where(~ !is.numeric(.x))) |> names()
if (length(non_numeric)) warning("Non-numeric drug columns: ", paste(non_numeric, collapse = ", "))

# ---- 2. Resolve the fentanyl / fentaNYL duplicate columns -------------
# Renamed by hand (not clean_names) so the two columns stay distinguishable.
meds <- meds_raw |>
  rename(fentanyl_alt = fentaNYL, fentanyl_alt_byWgt = fentaNYL_byWgt)

# Do the two columns describe the same infants?
cat("\nFentanyl column overlap (dose > 0):\n")
meds |>
  mutate(fentanyl     = coalesce(fentanyl, 0) > 0,
         fentanyl_alt = coalesce(fentanyl_alt, 0) > 0) |>
  count(fentanyl, fentanyl_alt) |>
  print()

# DECISION: set TRUE only after confirming the two columns are the same drug
# under two spellings. Row-wise sum; stays NA only if both are NA.
combine_fentanyl <- FALSE
if (combine_fentanyl) {
  both_na <- function(a, b) is.na(a) & is.na(b)
  meds <- meds |>
    mutate(
      fentanyl        = if_else(both_na(fentanyl, fentanyl_alt), NA_real_,
                                coalesce(fentanyl, 0) + coalesce(fentanyl_alt, 0)),
      fentanyl_byWgt  = if_else(both_na(fentanyl_byWgt, fentanyl_alt_byWgt), NA_real_,
                                coalesce(fentanyl_byWgt, 0) + coalesce(fentanyl_alt_byWgt, 0))
    ) |>
    select(-fentanyl_alt, -fentanyl_alt_byWgt)
}

# ---- 3. Missing vs zero ------------------------------------------------
# Open question: does NA mean "not given" or "not recorded"?
cat("\nNA and zero counts per column:\n")
meds |>
  summarise(across(-ID, list(na = ~ sum(is.na(.x)), zero = ~ sum(.x == 0, na.rm = TRUE)),
                   .names = "{.col}__{.fn}")) |>
  pivot_longer(everything(), names_to = c("column", "stat"), names_sep = "__") |>
  pivot_wider(names_from = stat, values_from = value) |>
  print(n = Inf)

# ---- 4. Per-drug exposure summary -------------------------------------
# Exposed = total dose > 0. NA is treated as NOT exposed (revisit per item 3).
drugs <- meds |> select(-ID, -ends_with("_byWgt")) |> names()
n_infants <- n_distinct(meds$ID)
q <- function(x) sprintf("%.2f [%.2f, %.2f]", median(x), quantile(x, .25), quantile(x, .75))

exposure_summary <- map_dfr(drugs, function(d) {
  dose <- meds[[d]]
  wgt  <- meds[[paste0(d, "_byWgt")]]
  exposed <- !is.na(dose) & dose > 0
  tibble(
    drug            = d,
    n_missing       = sum(is.na(dose)),
    n_exposed       = sum(exposed),
    pct_exposed     = round(100 * sum(exposed) / n_infants, 1),
    total_dose_med_iqr  = if (any(exposed)) q(dose[exposed]) else NA_character_,
    dose_per_kg_med_iqr = if (any(exposed & !is.na(wgt))) q(wgt[exposed & !is.na(wgt)]) else NA_character_
  )
}) |> arrange(desc(n_exposed))

print(exposure_summary, n = Inf, width = Inf)
write_csv(exposure_summary, out_path("tables", "medication_exposure_summary.csv"))

# ---- 5. Opioid exposure flags (the aim-relevant exposure) -------------
# Edit class membership as needed; fentanyl_alt included until merged.
opioids <- intersect(c("fentanyl", "fentanyl_alt", "morphine", "hydromorphone"), drugs)

meds <- meds |>
  mutate(any_opioid = if_any(all_of(opioids), ~ !is.na(.x) & .x > 0),
         n_opioids  = rowSums(across(all_of(opioids), ~ !is.na(.x) & .x > 0)))

cat("\nInfants with any opioid exposure:", sum(meds$any_opioid), "of", n_infants,
    sprintf("(%.1f%%)\n", 100 * mean(meds$any_opioid)))
cat("Number of distinct opioids per infant:\n")
print(count(meds, n_opioids))
