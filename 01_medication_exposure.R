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

if (nrow(meds_raw) != 194) warning("Expected 194 participants, found ", nrow(meds_raw))

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

# fentaNYL is the tall-man spelling of fentanyl. Only 1 of 194 infants has both,
# so the columns are treated as one drug. Row-wise sum; NA only if both are NA.
combine_fentanyl <- TRUE
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
n_infants <- nrow(meds)   # ID was removed from the de-identified file, so count rows

# Units of the raw dose columns; _byWgt is the same unit per kg.
# CONFIRM each entry; these are guesses based on usual neonatal dosing.
dose_units <- c(caffeine = "mg", fentanyl = "mcg", fentanyl_alt = "mcg", midazolam = "mg",
                morphine = "mg", acetaminophen = "mg", dexmedetomidine = "mcg",
                hydromorphone = "mg", lorazepam = "mg")
q <-function(x) sprintf("%.2f [%.2f, %.2f]", median(x), quantile(x, .25), quantile(x, .75))

exposure_summary <- map_dfr(drugs, function(d) {
  dose <- meds[[d]]
  wgt  <- meds[[paste0(d, "_byWgt")]]
  exposed <- !is.na(dose) & dose > 0
  tibble(
    drug            = d,
    unit            = unname(dose_units[d]),
    n_missing      = sum(is.na(dose)),
    n_exposed       = sum(exposed),
    pct_exposed     = round(100 * sum(exposed) / n_infants, 1),
    total_dose_med_iqr  = if (any(exposed)) q(dose[exposed]) else NA_character_,
    dose_per_kg_med_iqr = if (any(exposed & !is.na(wgt))) q(wgt[exposed & !is.na(wgt)]) else NA_character_
  )
}) |> arrange(desc(n_exposed))

print(exposure_summary, n = Inf, width = Inf)
write_csv(exposure_summary, out_path("tables", "medication_exposure_summary.csv"))

# ---- 4b. Implied weight check -----------------------------------------
# dose / dose-per-kg = implied weight (kg). If _byWgt were total dose over one
# fixed weight, this would match across drugs within an infant. Spread within an
# infant suggests per-administration weights were used.
implied_wt <- map_dfr(drugs, function(d) {
  w <- meds[[d]] / meds[[paste0(d, "_byWgt")]]
  w <- w[is.finite(w) & meds[[d]] > 0]
  tibble(drug = d, n = length(w),
         min = min(w), median = median(w), max = max(w))
}) |> filter(n > 0)
cat("\nImplied weight (kg) per drug among exposed infants:\n")
print(implied_wt, n = Inf)

# Infant-drug pairs with implausible implied weight (term-equivalent infants are
# under ~4 kg): likely entry errors in total dose or per-kg dose. Check at source.
implausible_wt <- map_dfr(drugs, function(d) {
  tibble(ID = meds$ID, drug = d, dose = meds[[d]], perkg = meds[[paste0(d, "_byWgt")]]) |>
    filter(dose > 0, perkg > 0, dose / perkg > 4) |>
    mutate(implied_wt_kg = dose / perkg)
})
cat("\nInfant-drug pairs with implied weight > 4 kg:\n")
print(implausible_wt, n = Inf)

# Within-infant spread across drugs; large values support the per-administration idea
within_infant <- meds |>
  mutate(row = row_number()) |>
  select(row, all_of(drugs), ends_with("_byWgt")) |>
  pivot_longer(-row, names_to = "col", values_to = "v") |>
  mutate(drug = sub("_byWgt$", "", col), kind = if_else(grepl("_byWgt$", col), "perkg", "dose")) |>
  select(-col) |>
  pivot_wider(names_from = kind, values_from = v) |>
  filter(dose > 0, perkg > 0) |>
  mutate(wt = dose / perkg) |>
  group_by(row) |>
  filter(n() > 1) |>
  summarise(wt_range = max(wt) - min(wt))
cat("\nWithin-infant range of implied weight (kg), infants with 2+ drugs:\n")
print(summary(within_infant$wt_range))

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

# ---- 6. Counts and percentages table ----------------------------------
# Exposure = any administration (dose > 0). 95% CIs are exact (Clopper-Pearson).
flag <- function(x) !is.na(x) & x > 0
nb <- "  "   # non-breaking spaces to indent sub-rows in Word

meds <- meds |>
  mutate(multi_opioid = n_opioids > 1,
         any_benzo    = flag(midazolam) | flag(lorazepam))

flags <- list(
  "Any opioid"           = meds$any_opioid,
  "Fentanyl"             = flag(meds$fentanyl),
  "Morphine"             = flag(meds$morphine),
  "Hydromorphone"        = flag(meds$hydromorphone),
  "More than one opioid" = meds$multi_opioid,
  "Any benzodiazepine"   = meds$any_benzo,
  "Midazolam"            = flag(meds$midazolam),
  "Lorazepam"            = flag(meds$lorazepam),
  "Dexmedetomidine"      = flag(meds$dexmedetomidine),
  "Caffeine"             = flag(meds$caffeine),
  "Acetaminophen"        = flag(meds$acetaminophen)
)
indented <- c("Fentanyl", "Morphine", "Hydromorphone", "Midazolam", "Lorazepam")

exposure_counts <- imap_dfr(flags, function(x, nm) {
  ci <- binom.test(sum(x), length(x))$conf.int
  tibble(medication = if (nm %in% indented) paste0(nb, nm) else nm,
         n = sum(x), N = length(x),
         pct = 100 * mean(x), ci_low = 100 * ci[1], ci_high = 100 * ci[2])
})
print(exposure_counts, n = Inf)
write_csv(exposure_counts, out_path("tables", "exposure_counts.csv"))

exposure_counts |>
  transmute(Medication = medication,
            `n (%)`    = sprintf("%d (%.1f%%)", n, pct),
            `95% CI`   = sprintf("%.1f–%.1f", ci_low, ci_high)) |>
  gt() |>
  tab_header(title = sprintf("Medication exposure, birth to term-equivalent MRI (N = %d)", n_infants)) |>
  gtsave(out_path("tables", "exposure_counts.docx"))

# ---- 7. Feasibility of non-opioid exposures ---------------------------
# 7a. Co-exposure: can a drug's effect be separated from the opioids'?
cat("\nOpioid x benzodiazepine co-exposure:\n")
print(count(meds, any_opioid, any_benzo))
cat("\nOpioid x acetaminophen co-exposure:\n")
print(count(meds, any_opioid, acetaminophen_exposed = flag(acetaminophen)))
cat("\nBenzodiazepine or dexmedetomidine (any sedative):\n")
print(count(meds, any_sedative = any_benzo | flag(dexmedetomidine)))

# Spearman correlations of cumulative per-kg dose across all infants
cat("\nSpearman correlations, cumulative per-kg dose:\n")
meds |>
  select(caffeine = caffeine_byWgt, fentanyl = fentanyl_byWgt, morphine = morphine_byWgt,
         midazolam = midazolam_byWgt, acetaminophen = acetaminophen_byWgt) |>
  cor(method = "spearman") |>
  round(2) |>
  print()

# 7b. Minimum detectable effect (Cohen's d, exposed vs unexposed) at 80% power.
# Assumes the full sample (310) keeps the current exposure prevalence.
# Second alpha = 0.05 / 5 drugs, a Bonferroni correction for testing five drugs.
mde_tbl <- expand_grid(medication = names(flags), N_total = c(n_infants, 310)) |>
  mutate(p = map_dbl(medication, ~ mean(flags[[.x]])),
         n_exposed = round(p * N_total),
         n_unexposed = N_total - n_exposed) |>
  filter(n_exposed >= 2, n_unexposed >= 2) |>
  rowwise() |>
  mutate(mde_d_alpha05 = pwr.t2n.test(n1 = n_exposed, n2 = n_unexposed,
                                      power = 0.8, sig.level = 0.05)$d,
         mde_d_alpha01 = pwr.t2n.test(n1 = n_exposed, n2 = n_unexposed,
                                      power = 0.8, sig.level = 0.01)$d) |>
  ungroup() |>
  mutate(across(starts_with("mde"), ~ round(.x, 2)), p = round(100 * p, 1))
print(mde_tbl, n = Inf)
write_csv(mde_tbl, out_path("tables", "minimum_detectable_effects.csv"))

# Dose-response across all infants (continuous cumulative dose): detectable correlation
cat("\nMinimum detectable correlation at 80% power:\n")
for (n in c(n_infants, 310)) {
  cat(" N =", n, " r =", round(pwr.r.test(n = n, power = 0.8, sig.level = 0.05)$r, 2), "\n")
}
