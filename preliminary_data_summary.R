# ============================================================================
# Preliminary data summary: serial MRI scans across studies
# Combines the omnibus sheet and the preemie scan QA sheet, removes the
# redundancy between them, summarizes each study cohort, and draws a swimmer plot.
#
# Run the whole script top to bottom. Items marked CHECK need your eyes the
# first time you run it.
# ============================================================================

library(tidyverse)
library(lubridate)

# ---- 0. Settings -----------------------------------------------------------
dir <- "C:/Users/Kevin/OneDrive - Childrens National Hospital/ComputerSyncing/03_Funding/02_Drafting/2026_ARPA-H_ORBIT/Letter of Intent/Preliminary Data"

omni_file <- "omnibussheet_01252023.csv"          # CHECK extension
qa_file   <- "preemieScanQAUpdate_07022025.csv"   # CHECK extension
delim     <- NULL        # NULL = auto-detect (tab, comma or semicolon); or set "\t" / ","

TERM_PMA  <- 37          # weeks: scans before this = preterm period, at/after = term age

# Fixed color per cohort (first three slots of a colorblind-checked palette)
cohort_cols <- c("Preterm" = "#2a78d6", "Onesie" = "#eb6834", "PROTECT" = "#1baf7a")

# ---- 1. Read both sheets -----------------------------------------------------
# Detect the delimiter from the header line (tab, comma, semicolon) unless set above
detect_delim <- function(path) {
  header <- readLines(path, n = 1, warn = FALSE)
  counts <- c("\t" = str_count(header, "\t"), "," = str_count(header, ","), ";" = str_count(header, ";"))
  names(which.max(counts))
}
# Lowercase and trim names and strip any byte-order mark, so "Project" / "project" / " Project" all match
clean_names_simple <- function(df) {
  names(df) <- names(df) %>% str_remove("^﻿") %>% str_remove("^\\W+") %>% str_trim() %>% str_to_lower()
  df
}

if (is.null(delim)) {
  omni <- read_delim(file.path(dir, omni_file), delim = detect_delim(file.path(dir, omni_file)), show_col_types = FALSE)
  qa   <- read_delim(file.path(dir, qa_file),   delim = detect_delim(file.path(dir, qa_file)),   show_col_types = FALSE)
} else {
  omni <- read_delim(file.path(dir, omni_file), delim = delim, show_col_types = FALSE)
  qa   <- read_delim(file.path(dir, qa_file),   delim = delim, show_col_types = FALSE)
}
omni <- clean_names_simple(omni)
qa   <- clean_names_simple(qa)

cat("Omnibus columns:", paste(names(omni), collapse = ", "), "\n")
cat("QA columns:     ", paste(names(qa),   collapse = ", "), "\n")

# ---- 2. Harmonize names and key formats ---------------------------------------
# Trim/lowercase the key columns so "P01 " and "p01" match across sheets.
omni2 <- omni %>%
  rename(project = project, shortID = id, scanNo = scan, scanPMA = scanpma,
         scanner = scanner, sex = sex, birthPMA = birthga,
         kidokoro = kidokoro, vol = nvol, FullID = fullid) %>%
  mutate(across(c(project, shortID), ~ str_to_lower(str_trim(as.character(.x)))),
         scanNo = as.integer(scanNo))

qa2 <- qa %>%
  rename(shortID = shortid, scanNo = scanno, scanPMA = scanpma, vol = finalvol,
         scanDate = scandate, birthPMA = birthpma) %>%
  mutate(across(c(project, shortID), ~ str_to_lower(str_trim(as.character(.x)))),
         scanNo = as.integer(scanNo))

# ---- 3. Diagnose the redundancy between sheets --------------------------------
keys <- c("project", "shortID", "scanNo")

cat("\n--- Redundancy check ---\n")
cat("Duplicate scan keys within omnibus:", nrow(filter(count(omni2, across(all_of(keys))), n > 1)), "\n")
cat("Duplicate scan keys within QA:     ", nrow(filter(count(qa2,   across(all_of(keys))), n > 1)), "\n")
cat("Scans in both sheets:", nrow(inner_join(omni2, qa2, by = keys)), "\n")
cat("Omnibus only:        ", nrow(anti_join(omni2, qa2, by = keys)), "\n")
cat("QA only:             ", nrow(anti_join(qa2, omni2, by = keys)), "\n")
cat("Project labels, omnibus:", paste(sort(unique(omni2$project)), collapse = ", "), "\n")
cat("Project labels, QA:     ", paste(sort(unique(qa2$project)),   collapse = ", "), "\n")
# CHECK: if "in both" is much lower than expected, the sheets label projects or
# scan numbers differently; match on FullID or scanPMA instead.

cat("\n--- Why don't the scans match? ---\n")
cat("Infants present in BOTH sheets:", length(intersect(paste(omni2$project, omni2$shortID),
                                                         paste(qa2$project,   qa2$shortID))), "\n")
cat("Omnibus 'Scan' values that are not whole numbers (these become NA):\n")
print(omni %>% filter(is.na(suppressWarnings(as.integer(scan)))) %>% count(scan, sort = TRUE), n = 20)

# Side-by-side look at one infant who is in both sheets
ex <- intersect(omni2$shortID, qa2$shortID)[1]
cat("\nExample infant in both sheets:", ex, "\n--- omnibus rows (raw Scan column) ---\n")
print(omni %>% filter(str_to_lower(str_trim(as.character(id))) == ex) %>%
        select(id, scan, scanpma, fullid), n = 20)
cat("--- QA rows ---\n")
print(qa2 %>% filter(shortID == ex) %>% select(shortID, scanNo, scanPMA, scanDate), n = 20)

# ---- 4. One row per scan ---------------------------------------------------------
# Where a scan is in both sheets, the newer QA sheet's value wins; the omnibus
# value fills any gap.
scans <- full_join(qa2, omni2, by = keys, suffix = c("", "_omni")) %>%
  mutate(
    sex      = coalesce(sex,      sex_omni),
    birthPMA = coalesce(birthPMA, birthPMA_omni),
    scanPMA  = coalesce(scanPMA,  scanPMA_omni),
    scanner  = coalesce(scanner,  scanner_omni),
    kidokoro = coalesce(kidokoro, kidokoro_omni),
    vol      = coalesce(vol,      vol_omni),
    participant = paste(project, shortID, sep = "_")
  ) %>%
  select(-ends_with("_omni"))

# Disagreements on scan PMA between sheets (should be empty or tiny)
disagree <- inner_join(
  select(qa2,   all_of(keys), scanPMA_qa   = scanPMA),
  select(omni2, all_of(keys), scanPMA_omni = scanPMA), by = keys) %>%
  filter(abs(scanPMA_qa - scanPMA_omni) > 0.1)
cat("Scans whose PMA differs between sheets:", nrow(disagree), "\n")

# ---- 5. Assign study cohorts -------------------------------------------------------
# Looks for the study name in project, shortID, or FullID, ignoring case.
# Onesie is tested first in case its IDs also contain "preterm".
if (!"FullID" %in% names(scans)) scans$FullID <- NA_character_

scans <- scans %>%
  mutate(
    id_text = str_to_lower(paste(project, shortID, FullID, sep = " | ")),
    cohort = case_when(
      str_detect(id_text, "onesi|one_si|1z") ~ "Onesie",   # CHECK: matches onesize, onesie, 1z
      str_detect(id_text, "protect") ~ "PROTECT",
      str_detect(id_text, "preterm") ~ "Preterm",   # CHECK: add |preemie if needed
      TRUE ~ "UNASSIGNED"
    )
  )

cat("\n--- Scans per cohort (nothing should be UNASSIGNED) ---\n")
print(count(scans, cohort))
cat("\n--- Which sheet each cohort's scans come from (project label x cohort) ---\n")
print(count(scans, project, cohort), n = 50)

# What do the IDs actually look like? Shows each distinct ID "shape" (letters kept,
# digits replaced with #) with how many infants and an example, by cohort.
cat("\n--- ID shapes by cohort (use this to fix the cohort patterns) ---\n")
id_shapes <- scans %>%
  mutate(shortID_shape = str_replace_all(shortID, "\\d+", "#"),
         fullID_shape  = str_replace_all(str_to_lower(as.character(FullID)), "\\d+", "#")) %>%
  group_by(cohort, project, shortID_shape, fullID_shape) %>%
  summarize(n_infants = n_distinct(participant), example = first(shortID), .groups = "drop") %>%
  arrange(cohort, desc(n_infants))
print(id_shapes, n = 100, width = Inf)
write_csv(id_shapes, file.path(dir, "id_shapes.csv"))
if (any(scans$cohort == "UNASSIGNED")) {
  cat("Unassigned IDs:\n")
  print(distinct(filter(scans, cohort == "UNASSIGNED"), project, shortID, FullID))
}

# ---- 6. Scan-level calculations ---------------------------------------------------
# DOB is not in either sheet, so it is derived from the QA sheet as
# scanDate minus the weeks of life at the scan (scanPMA - birthPMA).
# Infants that appear only in the omnibus sheet have no scanDate, so no DOB.
scans <- scans %>%
  mutate(
    scanDate = as_date(parse_date_time(as.character(scanDate), orders = c("mdy", "ymd", "dmy"), quiet = TRUE)),
    period   = if_else(scanPMA < TERM_PMA, "preterm", "term"),
    dob_est  = scanDate - round(7 * (scanPMA - birthPMA))
  )

# ---- 7. Participant level -----------------------------------------------------------
infants <- scans %>%
  group_by(cohort, participant) %>%   # UNASSIGNED stays in so unmatched infants are visible, not dropped
  summarize(
    n_scans   = n(),
    n_preterm = sum(period == "preterm", na.rm = TRUE),
    n_term    = sum(period == "term",    na.rm = TRUE),
    first_pma = suppressWarnings(min(scanPMA, na.rm = TRUE)),
    last_pma  = suppressWarnings(max(scanPMA, na.rm = TRUE)),
    birthGA   = first(na.omit(birthPMA)),
    dob       = as_date(median(dob_est, na.rm = TRUE)),   # median across scans smooths noise
    .groups = "drop"
  ) %>%
  mutate(current_age_yrs = as.numeric(Sys.Date() - dob) / 365.25)

# ---- 8. Cohort summary (the ten requested metrics) -----------------------------------
# min/max that return NA (not Inf/-Inf) when a cohort has no values, e.g. no derivable DOB
safe_min <- function(x) if (all(is.na(x))) x[NA_integer_][1] else min(x, na.rm = TRUE)
safe_max <- function(x) if (all(is.na(x))) x[NA_integer_][1] else max(x, na.rm = TRUE)

cohort_summary <- infants %>%
  group_by(cohort) %>%
  summarize(
    `1 N infants`              = n(),
    `2 N scans`                = sum(n_scans),
    `3 Mean scans/infant`      = round(mean(n_scans), 2),
    `4 Preterm-period scans`   = sum(n_preterm),
    `5 Term-age scans`         = sum(n_term),
    `6 Mean GA at birth (wk)`  = round(mean(birthGA, na.rm = TRUE), 1),
    `7 Earliest DOB`           = safe_min(dob),
    `8 Latest DOB`             = safe_max(dob),
    `9 Current age range (yr)` = if_else(all(is.na(current_age_yrs)), NA_character_,
                                         paste0(round(safe_min(current_age_yrs), 1), "-",
                                                round(safe_max(current_age_yrs), 1))),
    `10 Mean current age (yr)` = round(mean(current_age_yrs, na.rm = TRUE), 1),
    n_missing_dob              = sum(is.na(dob)),   # infants with no derivable DOB
    .groups = "drop"
  )

cat("\n--- Cohort summary ---\n")
print(cohort_summary, width = Inf)

# Distribution of serial scans per infant
cat("\n--- Scans per infant, by cohort ---\n")
print(infants %>% count(cohort, n_scans) %>%
        pivot_wider(names_from = n_scans, values_from = n, values_fill = 0))

# ---- 9. Swimmer plot -------------------------------------------------------------------
plot_df <- scans %>%
  filter(cohort != "UNASSIGNED", !is.na(scanPMA)) %>%
  inner_join(select(infants, participant, first_pma, last_pma), by = "participant")

# Order infants within each cohort by first scan, then last scan (staircase look)
plot_order <- plot_df %>%
  distinct(cohort, participant, first_pma, last_pma) %>%
  arrange(cohort, first_pma, last_pma)
plot_df <- plot_df %>%
  mutate(participant = factor(participant, levels = unique(plot_order$participant)))

p <- ggplot(plot_df, aes(y = participant, colour = cohort, fill = cohort)) +
  geom_vline(xintercept = TERM_PMA, linetype = "dashed", colour = "#898781", linewidth = 0.4) +
  geom_segment(aes(x = first_pma, xend = last_pma, yend = participant),
               linewidth = 0.5, alpha = 0.7) +
  geom_point(aes(x = scanPMA), shape = 21, size = 1.8, stroke = 0.4, colour = "#fcfcfb") +
  facet_grid(cohort ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_manual(values = cohort_cols, name = "Study cohort") +
  scale_fill_manual(values = cohort_cols, name = "Study cohort") +
  scale_x_continuous(breaks = seq(24, 80, 4), name = "Postmenstrual age at scan (weeks)") +
  annotate("text", x = TERM_PMA + 0.3, y = Inf, label = paste0("term (", TERM_PMA, " wk)"),
           hjust = 0, vjust = 1.5, size = 3, colour = "#52514e") +
  labs(y = "Participants (one row per infant)", title = "Serial MRI scan timing by participant") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "#e1e0d9"),
    strip.text.y.left  = element_text(angle = 0, hjust = 1, colour = "#0b0b0b"),
    plot.background    = element_rect(fill = "#fcfcfb", colour = NA),
    legend.position    = "top"
  )

print(p)

# ---- 10. Save outputs -------------------------------------------------------------------
write_csv(scans,          file.path(dir, "combined_scans.csv"))
write_csv(infants,        file.path(dir, "participants.csv"))
write_csv(cohort_summary, file.path(dir, "cohort_summary.csv"))
ggsave(file.path(dir, "swimmer_plot.png"), p, width = 8, height = 10, dpi = 300)
cat("\nDone. Outputs written to:", dir, "\n")
