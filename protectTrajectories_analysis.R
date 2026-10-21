###############################################################################
# PROTECT TRAJECTORIES
# Fetal vs. preterm resting-state functional connectivity, and how connectivity
# changes with postmenstrual age (PMA) across degrees of prematurity.
#
# CLEANED, CONSOLIDATED SCRIPT. Run top to bottom (or section by section).
# Replaces protectTrajectoriesAnalyses.R and the code blocks developed in chat.
# Run it interactively (RStudio: Ctrl+Enter by section, or Source with Echo) so
# the summary numbers print; results are also kept in named objects.
# SPEED: the edge-wise loops run in parallel (N_CORES) and slow results are saved
# to figureDir/cache and reloaded on later runs (USE_CACHE); see section 0.
#
# DESIGN IN ONE PARAGRAPH
#   Fetal scans (all 1.5T) and preterm scans (1.5T or 3T; 1-4 scans per infant)
#   are harmonized with LongComBat (scanner = batch). Connectivity is Fisher z.
#   For each of the 3,003 cortical connections ("edges"):
#     (A) Fetal vs. preterm:   FC ~ PMA + scanner + group, random intercept/infant
#     (B) Preterm only:        FC ~ scanner + PMA * degree of prematurity,
#                              random intercept per infant; the PMA x prematurity
#                              interaction is tested with a 2-df Wald chi-square.
#   FDR (Benjamini-Hochberg) is applied across all 3,003 edges for every test.
#
# GROUP DEFINITIONS (fixed for the paper; see birth_cat() below)
#   EPT: birth GA <= 28 wk | VPT: > 28 to <= 32 wk | MPT: > 32 wk
#   i.e. exactly 28.0 -> EPT and exactly 32.0 -> VPT. The standard clinical
#   convention (<28 / 28-<32 / >=32) differs for 6 infants at 28.0 and 5 at 32.0;
#   it is run as a sensitivity analysis in section 9.
#
# CONTENTS
#   0  Setup: packages, switches, paths, helper functions
#   1  Read data and scan-level QA
#   2  Fetal controls + LongComBat harmonization (+ harmonization checks)
#   3  Analysis data: connection table, scan table, preterm table, Table 1
#   4  Fetal vs. preterm, all edges (Fig. 2: matrices, chords, homotopic bars; Table 2; Fig. 2 assembled; SI Fig. 2 mean matrices)
#   5  Preterm only: PMA x prematurity, all edges (Fig. 3, Table 3)
#   6  Fig. 3 Panels B and C (group trajectories)
#   7  Fig. 3 Panel A (BrainNet Viewer edge file; 7b: the same picture in R with rgl, and the combined Fig. 3 A-C)
#   8  Group slopes and pairwise contrasts for the significant edges
#   9  Sensitivity analyses and controls
#  10  Save tables
#  11  Sex differences (exploratory)
#  12  Four-group model (Fetal, EPT, VPT, MPT)
#  13  Group comparisons for the Table 1 variables, session info
#
# NOTES FROM THE ORIGINAL SCRIPT
#   * onesie_1z065_03 is removed from the connectivity data (extreme outlier).
#   * protect_092_01 is excluded before harmonization.
#   * The old (non-longitudinal) neuroCombat block, the lobe/class/module
#     analyses, GAMM/spline analyses, and the commented-out model blocks were
#     removed. They remain in protectTrajectoriesAnalyses.R if needed.
###############################################################################


###############################################################################
# 0. SETUP
###############################################################################

## ---- Packages --------------------------------------------------------------
# install.packages(c("ggplot2","nlme","lmerTest","car","circlize","ggpubr","scales","ragg","svglite","png","rgl"))
# devtools::install_github("jcbeer/longCombat")
library(ggplot2)
library(nlme)       # lme(): linear mixed models with a random intercept
library(lmerTest)   # lmer() with p-values (used for Table 3, with car)
library(car)        # Anova(): Type II Wald chi-square tests for lmer models
library(longCombat) # longitudinal ComBat
library(circlize)   # chord diagrams (Fig. 2)
library(ggpubr)     # ggarrange() for the multi-panel Fig. 3B

## ---- Switches --------------------------------------------------------------
SAVE_FIGS                 <- FALSE  # TRUE: write figure files to figureDir
FIG_FORMATS               <- c("tiff", "svg")   # formats written when SAVE_FIGS is TRUE: any of "tiff","png","svg","pdf"
SAVE_TABLES               <- FALSE  # TRUE: write .csv result tables to figureDir
RUN_QA_PLOTS              <- TRUE   # scan-QA figures (section 1)
RUN_HARMONIZATION_CHECKS  <- TRUE   # scanner effect / PCA / PMA preservation (section 2)
RUN_SENSITIVITY           <- TRUE   # section 9: five lmer scans over all edges (slow serially; fast with N_CORES, cached)
RUN_EXTRA_CONTROLS        <- TRUE   # sections 9f-9h: age-matched fetal vs. preterm, postnatal-age check, unharmonized data
RUN_DISTANCE_CONTROL      <- TRUE   # section 9i: do fetal-vs-preterm differences depend on the distance between ROIs?
RUN_QC_COVARIATE_CONTROL  <- TRUE   # section 9e: finalVol and motion (framewise displacement) covariates and similar-volume / similar-motion reruns
RUN_SEX_ANALYSES          <- TRUE   # section 11: sex differences (fetal vs. preterm and within preterm; exploratory)
RUN_TABLE1_TESTS          <- TRUE   # section 13: group comparisons for the Table 1 variables (preterm vs. fetal; EPT / VPT / MPT)
RUN_FOUR_GROUP            <- TRUE   # section 12: four-group model (Fetal, EPT, VPT, MPT)
RUN_EQUAL_VOLUME_CONTROL <- TRUE   # section 9k: fetal vs. preterm with every scan cut to the same number of volumes (needs data/rawTS or data/rawRS)
N_CORES                   <- 8      # cores for the edge-wise loops (this machine: 8 physical, 16 logical); 1 = serial
USE_CACHE                 <- TRUE   # reuse saved results of slow steps (see `cached()`); FALSE = always recompute

## ---- Paths -----------------------------------------------------------------
mainDirs <- c("C:/Users/Kevin/OneDrive - Childrens National Hospital/ComputerSyncing/Projects/Active Projects/Protect Trajectories",
              "C:/Users/kmcook/OneDrive - Childrens National Hospital/ComputerSyncing/Projects/Active Projects/Protect Trajectories")
mainDir   <- mainDirs[dir.exists(mainDirs)][1]
stopifnot(!is.na(mainDir))
dataDir   <- file.path(mainDir, "data")
figureDir <- file.path(mainDir, "figure")
dir.create(figureDir, showWarnings = FALSE)
setwd(dataDir)

## ---- Constants and helper functions ---------------------------------------
# Optimizer settings for lme(): more iterations and the "optim" optimizer help
# the occasional edge that fails to converge with the defaults.
ctrl <- lmeControl(msMaxIter = 200, tolerance = 1e-4, msTol = 1e-5, opt = "optim")

# Birth-category definition used for the paper AND inside the harmonization.
#   rule = "paper"    : EPT <= 28 | VPT > 28 to <= 32 | MPT > 32   (right = TRUE)
#   rule = "clinical" : EPT <  28 | VPT 28 to < 32   | MPT >= 32   (right = FALSE)
birth_cat <- function(ga, rule = "paper") {
  cut(ga, breaks = c(0, 28, 32, 50), labels = c("EPT", "VPT", "MPT"),
      right = (rule == "paper"), include.lowest = TRUE)
}

# Colors: fetal = red, preterm = blue. Lobe sidebar/sector colors (frontal pink,
# temporal yellow, parietal green, occipital blue) are separate and unchanged.
# Preterm blue #2a78d6 (passes the color-blind/contrast checks); the preterm heatmap is built from the same blue.
# Fig. 3B/C groups run cyan (EPT) -> blue (VPT) -> navy (MPT).
pal <- list(fetal = "#cc2936", preterm = "#2a78d6",
            cols_pt      = c(EPT = "#12b5cf", VPT = "#2a78d6", MPT = "#0d366b"),
            heat_preterm = c("#cde2fb", "#2a78d6", "#0d366b"),                     # low, mid, high
            heat_fetal   = c("#f4b1ac", "#cc2936", "#650f16"))
col_fetal   <- pal$fetal
col_preterm <- pal$preterm
cols_pt     <- pal$cols_pt
grp_labels  <- c(EPT = "Extremely Premature", VPT = "Very Premature", MPT = "Moderately Premature")

# Lobe colors (chord sectors, matrix sidebars, bar-chart strip) and the ROIs in each lobe
lobe_pal  <- c(Frontal = "#ef476f", Temporal = "#ffd166", Parietal = "#06d6a0", Occipital = "#118ab2")
lobe_rois <- list(Frontal  = c("PreCG","SFGdor","ORBsupb","MFG","ORBmid","IFGoperc","IFGtriang","ORBinf","ROL","SMA","OLF","SFGmed","ORBmed","REC","ACG"),
                  Temporal = c("INS","PHG","HES","STG","TPOsup","MTG","TPOmid","ITG"),
                  Parietal = c("MCG","PCG","PoCG","SPG","IPL","SMG","ANG","PCUN","PCL"),
                  Occipital = c("CAL","CUN","LING","SOG","MOG","IOG","FFG"))

# REGION NAMES USED IN FIGURES (chord labels and the homotopic bar charts).
# Left side = atlas code; right side = the label that appears in the figure. EDIT HERE.
# Hemisphere is dropped (the chord's left/right sides and the bilateral bars show it).
roi_display <- c(PreCG = "Precentral",               SFGdor = "Dorsal Superior Frontal", ORBsupb = "Superior Orbital",
                 MFG = "Middle Frontal",             ORBmid = "Middle Orbital",          IFGoperc = "Inferior Operculum",
                 IFGtriang = "Inferior Triangle",    ORBinf = "Inferior Orbital",        ROL = "Rolandic",
                 SMA = "Supplemental Motor",         OLF = "Olfactory",                  SFGmed = "Medial Superior Frontal",
                 ORBmed = "Medial Orbital",          REC = "Rectus",                     ACG = "Anterior Cingulate",
                 INS = "Insula",                     PHG = "Parahippocampal",            HES = "Heschl",
                 STG = "Superior Temporal",          TPOsup = "Superior Pole",           MTG = "Middle Temporal",
                 TPOmid = "Middle Pole",             ITG = "Inferior Temporal",          MCG = "Middle Cingulate",
                 PCG = "Posterior Cingulate",        PoCG = "Postcentral",               SPG = "Superior Parietal",
                 IPL = "Inferior Parietal Lobule",   SMG = "Supramarginal",              ANG = "Angular",
                 PCUN = "Precuneus",                 PCL = "Paracentral Lobule",         CAL = "Calcarine",
                 CUN = "Cuneus",                     LING = "Lingual",                   SOG = "Superior Occipital",
                 MOG = "Middle Occipital",           IOG = "Inferior Occipital",         FFG = "Fusiform")
# display_name("PreCG.L"), display_name("PreCG-R") or display_name("PreCG") -> "Precentral" (vectorized)
display_name <- function(x) {
  base <- sub("[-.][LR]$", "", x)
  out  <- unname(roi_display[base])
  ifelse(is.na(out), base, out)                      # unknown codes fall back to the code itself
}

# ggplot2 >= 3.5 places legends inside the panel with legend.position = "inside"
legend_inside <- function(x, y) {
  if (utils::packageVersion("ggplot2") >= "3.5.0") theme(legend.position = "inside", legend.position.inside = c(x, y))
  else theme(legend.position = c(x, y))
}

# Shared ggplot theme for the figures that use the lab's house style
theme_paper <- function(grid_y = TRUE, grid_x = FALSE, legend = "none") {
  theme(axis.title        = element_text(color = "#000000", size = 16),
        axis.text         = element_text(color = "#000000", size = 12),
        axis.ticks        = element_blank(),
        axis.line         = element_line(color = "#000000", linewidth = 1),
        plot.title        = element_text(color = "#000000", size = 24, hjust = 0.5),
        panel.grid.major.y = if (grid_y) element_line(color = "#aaaaaa") else element_blank(),
        panel.grid.major.x = if (grid_x) element_line(color = "#aaaaaa55") else element_blank(),
        panel.background  = element_rect(fill = "white"),
        legend.position   = legend,
        legend.text       = element_text(color = "#000000", size = 12),
        legend.title      = element_blank())
}

# Writes the figure in every format in FIG_FORMATS (when SAVE_FIGS is TRUE), or just
# draws it on screen. `file` may be given with or without an extension; the extension
# is replaced by each format. `draw` is a function that prints a ggplot / draws base graphics.
#   tiff/png : raster via ragg (or, without ragg, the cairo device). R's default tiff()
#              device leaves a notch at the top-left corner of bar outlines, so avoid it.
#   svg      : editable vector for Inkscape / Illustrator (svglite, or svg() without it).
#              Text stays live text; install the fonts or convert text to outlines on opening.
#   pdf      : vector via cairo_pdf (not written by default; add "pdf" to FIG_FORMATS).
# install.packages(c("ragg", "svglite")) once for the best output.
show_fig <- function(file, w, h, draw, res = 600, save = SAVE_FIGS) {
  if (!save) { draw(); return(invisible(NULL)) }
  base <- file.path(figureDir, tools::file_path_sans_ext(file))
  open_dev <- function(fmt) {
    path <- paste0(base, ".", fmt)
    switch(fmt,
      tiff = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_tiff(path, width = w, height = h, units = "in", res = res, compression = "lzw")
             else tiff(path, units = "in", width = w, height = h, res = res, compression = "lzw", type = "cairo"),
      png  = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png(path, width = w, height = h, units = "in", res = res)
             else png(path, units = "in", width = w, height = h, res = res, type = "cairo"),
      svg  = if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite(path, width = w, height = h)
             else svg(path, width = w, height = h),
      pdf  = cairo_pdf(path, width = w, height = h),
      stop("Unknown figure format: ", fmt))
  }
  for (fmt in FIG_FORMATS) {
    open_dev(fmt)
    tryCatch(draw(), finally = dev.off())
  }
  invisible(NULL)
}

# Trims the white border off a picture (a height x width x 3 array from png::readPNG), keeping m px of white.
# With several pictures, the same crop (the union of their drawings) is applied to all, so pieces that belong together
# keep the same scale and alignment.
trim_white <- function(im, m = 6, thr = .99) {
  one <- is.array(im)                                             # one picture (an array) or a list of pictures
  ims <- if (one) list(im) else im
  nw  <- lapply(ims, function(x) x[, , 1] < thr | x[, , 2] < thr | x[, , 3] < thr)
  r   <- range(unlist(lapply(nw, function(z) which(rowSums(z) > 0)))); cc <- range(unlist(lapply(nw, function(z) which(colSums(z) > 0))))
  out <- lapply(ims, function(x) { rr <- c(max(1, r[1] - m), min(dim(x)[1], r[2] + m)); c2 <- c(max(1, cc[1] - m), min(dim(x)[2], cc[2] + m))
                                   x[rr[1]:rr[2], c2[1]:c2[2], , drop = FALSE] })
  if (one) out[[1]] else out
}

# p-value formatting for tables: "<.001" or ".034"
fmt_p <- function(p) ifelse(p < .001, "<.001", sub("^0", "", sprintf("%.3f", p)))
# 3-decimal formatting that never prints a negative zero ("-0.000" -> "0.000")
fmt3 <- function(x) sub("^-(0\\.000)$", "\\1", sprintf("%.3f", x))

## ---- Output naming ------------------------------------------------------------------------------------------------
# Every saved file follows one pattern, so related items sort together in figureDir:
#   main text      Fig<n>_<Name>.tiff/.svg           Table<n>_<Name>.csv        (Fig2_part_<Name> = pieces of Fig. 2)
#   full results   Data_<Name>.csv                   (edge-wise results behind the main figures and tables)
#   supplement     SI<nn>_Fig_<Name>.tiff/.svg       SI<nn>_Table_<Name>.csv    (one SI number per topic; its tables and figures sort together)
# <Name> is UpperCamelCase; sub-parts are separated by underscores. The SI numbers live only in this registry: renumber here
# (e.g. SI_ID[["gt32"]] <- "SI12") and every file name follows. Fig. 1 (flowchart) is made by flowchart/make_flowchart.py.
SI_ID <- c(usable_scans = "SI01", mean_matrices = "SI02", harmonization = "SI03", posthoc = "SI04",
           clinical = "SI05", continuous_GA = "SI06", gt32 = "SI07", within_between = "SI08",
           volumes = "SI09", motion_covariate = "SI10", age_matched = "SI11", postnatal = "SI12",
           unharmonized = "SI13", distance = "SI14", motion_followup = "SI15", equal_volume = "SI16",
           sex_fetal_preterm = "SI17", sex_preterm = "SI18", four_group = "SI19", table1_tests = "SI20")
si_name  <- function(x) gsub("(^|_)([a-z])", "\\1\\U\\2", x, perl = TRUE)          # clinicalCutoffs_summary -> ClinicalCutoffs_Summary
si_table <- function(df, key, name) {
  df <- as.data.frame(df)
  if (SAVE_TABLES) write.csv(df, file.path(figureDir, paste0(SI_ID[[key]], "_Table_", si_name(name), ".csv")), row.names = FALSE)
  invisible(df)
}
si_fig <- function(key, name, w, h, draw, ...) show_fig(paste0(SI_ID[[key]], "_Fig_", si_name(name)), w, h, draw, ...)

## ---- Speed: parallel loops and cached results -------------------------------
# par_lapply(): lapply() over edges, spread across N_CORES. Linux/Mac use forked
# workers (they see all objects); Windows starts a cluster and sends it the
# objects named in `export` (the data and settings the edge functions use).
# Progress messages are only printed when N_CORES = 1.
par_lapply <- function(X, FUN, export = c("goodFCData", "ctrl", "scn_term", "scanInfo", "pt"), ...) {
  if (N_CORES <= 1) return(lapply(X, FUN, ...))
  message("Fitting ", length(X), " edges on ", N_CORES, " cores...")
  if (.Platform$OS.type == "unix") {
    res <- parallel::mclapply(X, FUN, ..., mc.cores = N_CORES)
    bad <- vapply(res, inherits, NA, what = "try-error")
    if (any(bad)) stop("Parallel worker failed: ", as.character(res[[which(bad)[1]]]))
    return(res)
  }
  cl <- parallel::makeCluster(N_CORES)
  on.exit(parallel::stopCluster(cl))
  parallel::clusterEvalQ(cl, { library(nlme); library(lmerTest); library(car) })
  parallel::clusterExport(cl, export, envir = globalenv())
  parallel::parLapply(cl, X, FUN, ...)
}

# cached(): run `expr` once and save the result to figureDir/cache. Next time the
# same name + key is found on disk and loaded instead of recomputed. The key
# fingerprints the data, so changed data gets a new file automatically. If you
# change the MODEL CODE itself, delete figureDir/cache (or set USE_CACHE <- FALSE).
cache_dir <- file.path(figureDir, "cache")
dir.create(cache_dir, showWarnings = FALSE)
cached <- function(name, key, expr) {
  h <- sprintf("%08x", as.integer(sum(utf8ToInt(key) * seq_along(utf8ToInt(key))) %% 2^31))
  f <- file.path(cache_dir, paste0(name, "_", h, ".rds"))
  if (USE_CACHE && file.exists(f)) { message("Loaded cached result: ", basename(f)); return(readRDS(f)) }
  val <- expr                                   # evaluated only when there is no cache
  if (USE_CACHE) saveRDS(val, f)
  val
}
# Fingerprint of a scan table: row count plus sums of its numeric columns and factor codes
table_key <- function(d) {
  num <- sum(vapply(d[vapply(d, is.numeric, NA)], sum, 0, na.rm = TRUE))
  fac <- sum(vapply(d[vapply(d, is.factor, NA)], function(f) sum(as.integer(f), na.rm = TRUE), 0))
  sprintf("n%d_%.8g_%.8g", nrow(d), num, fac)
}


###############################################################################
# 1. READ DATA AND SCAN-LEVEL QA
###############################################################################

# Connectivity (Fisher z), one row per scan, one column per connection.
fcData <- read.csv("correlationMatrixFull_V2.csv", stringsAsFactors = FALSE, row.names = 1)
fcData <- fcData[rownames(fcData) != "onesie_1z065_03", ]      # extreme outlier

# Participant-level demographics (Table 1)
demoData <- read.csv("demoData.csv")

# ROI labels, coordinates and lobe/class/module memberships (kept = 111 + 4 ROIs)
roiData <- read.csv(file.path(dataDir, "labelsAndCoords", "labelKey_Update.csv"), stringsAsFactors = FALSE)
roiData <- roiData[c(1:111, 120:123), ]

# Scan QA table (also holds scan-level demographics)
scanQA <- read.csv("preemieScanQAUpdate_07022025.csv", stringsAsFactors = FALSE)
scanQA$fullID <- paste0(scanQA$shortID, "_0", scanQA$scanNo)
scanQA <- subset(scanQA, !is.na(finalVol) & (kidokoro < 9 | is.na(kidokoro)))   # exclude severe injury (Kidokoro >= 9)

# A scan is usable if it has >= 4 min of data: >= 120 volumes for preterm scans (TR = 2 s); fetal scans
# need >= 80 volumes (TR = 3 s)
scanQA$good4Min <- as.numeric(scanQA$finalVol > 119)
goodData <- subset(scanQA, good4Min == 1)

## ---- Data-integrity check: birth PMA must be ONE value per infant -----------
# (A mistyped birth PMA once split one infant across birth categories.)
birth_range <- tapply(scanQA$birthPMA, scanQA$shortID, function(x) diff(range(x, na.rm = TRUE)))
bad_birth_ids <- names(birth_range)[!is.na(birth_range) & birth_range > 1e-8]
if (length(bad_birth_ids) > 0) {
  stop(length(bad_birth_ids), " infant(s) have more than one birth PMA across scans. ",
       "Inspect `bad_birth_ids`, fix the source file, and rerun.")
}

## ---- QA summary -------------------------------------------------------------
round(c(pct_usable_all   = mean(scanQA$finalVol > 119) * 100,
        pct_usable_lt37  = mean(scanQA$finalVol[which(scanQA$scanPMA <  37)] > 119) * 100,     # which(): scans with a missing scanPMA are left out
        pct_usable_gt37  = mean(scanQA$finalVol[which(scanQA$scanPMA >  37)] > 119) * 100), 1)
c(scans_missing_scanPMA = sum(is.na(scanQA$scanPMA)))                                         # these scans are dropped from the QA plots below (ggplot warns about removed rows)
c(infants_all = length(unique(scanQA$shortID)), infants_usable = length(unique(goodData$shortID)),
  scans_all = nrow(scanQA), scans_usable = nrow(goodData))
table(table(goodData$shortID))                                  # infants by number of usable scans
round(c(minutes_mean = mean(goodData$finalVol * 2 / 60),   # minutes of data per scan
        sd = sd(goodData$finalVol * 2 / 60),
        min = min(goodData$finalVol * 2 / 60), max = max(goodData$finalVol * 2 / 60)), 2)

if (RUN_QA_PLOTS) {
  # (a) every scan by PMA; usable (circle) vs. excluded (x)
  p_qa_all <- ggplot(scanQA, aes(x = shortID, y = scanPMA)) +
    scale_y_continuous(limits = c(24, 46), breaks = seq(20, 60, 2), expand = c(0, 0)) +
    scale_shape_manual(values = c(4, 21)) +
    scale_color_manual(values = c("#911717", "#095c00")) +
    geom_line(aes(group = shortID), linewidth = 1) +
    geom_point(aes(shape = as.factor(good4Min), color = as.factor(good4Min)), size = 3, stroke = 2, fill = "#7ef472") +
    xlab("Participants") + ylab("PMA at Scan") +
    theme_paper() + theme(axis.text.x = element_blank())

  # (b) salvage histogram: PMA of usable vs. excluded scans
  p_qa_hist <- ggplot(scanQA, aes(x = scanPMA)) +
    scale_x_continuous(limits = c(24, 46), breaks = seq(20, 60, 2), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 30), breaks = seq(0, 40, 2), expand = c(0, 0)) +
    geom_histogram(data = subset(scanQA, good4Min == 1), color = "#095c00", fill = "#7ef472aa", bins = 30) +
    geom_histogram(data = subset(scanQA, good4Min == 0), color = "#911717", fill = "#c75b5baa", bins = 30) +
    xlab("PMA at Scan") + theme_paper() + theme(axis.title.y = element_blank())

  # (c) volumes per scan with the 120-volume cutoff
  p_qa_vol <- ggplot(scanQA, aes(x = finalVol)) +
    geom_vline(xintercept = 120, linetype = "dashed") +
    geom_histogram(color = "#000000", fill = "#000055aa", bins = 30) +
    xlab("Final Volumes") + theme_paper() + theme(axis.title.y = element_blank())

  print(p_qa_all); print(p_qa_hist); print(p_qa_vol)
}


###############################################################################
# 2. FETAL CONTROLS + LONGITUDINAL COMBAT
#    Scanner (1.5T "450" vs 3T "750") is the batch. Fetal scans are all 1.5T.
#    Protected (fixed) effects: group (Fetal/Preterm) and scan PMA.
#    Random effects: participant intercept + a preterm-only birth-subgroup term.
#    NOTE: group x PMA interactions are NOT protected (disclose as a limitation).
###############################################################################

fetalFCData   <- read.csv(file.path(dataDir, "fetalControl", "fetalCorrelationMatrixFull_Update.csv"),
                          stringsAsFactors = FALSE, row.names = 1)
fetalDemoData <- read.csv(file.path(dataDir, "fetalControl", "fetalDemoData.csv"))
fetalDemoData$good4Min <- 1
fetalDemoData$scanType <- "Fetal"
goodData$scanType      <- "Preterm"

## ---- Supplementary figure: usable scans per participant, preterm (top) and fetal (bottom) ----
# Usable scans only, colored by scanner; both panels share the PMA axis so the groups are comparable.
# Fetal scans are all 1.5T; the legend lists both scanners.
if (RUN_QA_PLOTS) {
  qa_good_pt <- subset(scanQA, good4Min == 1)
  qa_good_ft <- fetalDemoData
  qa_ylim <- c(floor(min(qa_good_pt$scanPMA, qa_good_ft$scanPMA) / 2) * 2, ceiling(max(qa_good_pt$scanPMA, qa_good_ft$scanPMA) / 2) * 2)
  qa_panel <- function(d, ttl) {
    d$scanner_f <- factor(as.character(d$scanner), levels = c("450", "750"))
    ggplot(d, aes(x = shortID, y = scanPMA)) +
      scale_y_continuous(limits = qa_ylim, breaks = seq(20, 60, 2), expand = c(0, 0)) +
      scale_shape_manual(values = c("450" = 21, "750" = 23), labels = c("450" = "1.5T", "750" = "3T"), limits = c("450", "750")) +
      scale_fill_manual(values = c("450" = "#0096FF", "750" = "#90EE90"), labels = c("450" = "1.5T", "750" = "3T"), limits = c("450", "750")) +
      scale_color_manual(values = c("450" = "#00008B", "750" = "#095c00"), labels = c("450" = "1.5T", "750" = "3T"), limits = c("450", "750")) +
      geom_line(aes(group = shortID), linewidth = 1) +
      geom_point(aes(shape = scanner_f, color = scanner_f, fill = scanner_f), size = 3, stroke = 2) +
      labs(x = "Participants", y = "PMA at Scan", title = ttl) +
      theme_paper(legend = "bottom") + theme(axis.text.x = element_blank(), plot.title = element_text(size = 18, hjust = 0))
  }
  p_qa_pt_s <- qa_panel(qa_good_pt, sprintf("Preterm (%d infants, %d scans)", length(unique(qa_good_pt$shortID)), nrow(qa_good_pt)))
  p_qa_ft_s <- qa_panel(qa_good_ft, sprintf("Fetal (%d participants, %d scans)", length(unique(qa_good_ft$shortID)), nrow(qa_good_ft)))
  p_qa_stack <- ggpubr::ggarrange(p_qa_pt_s, p_qa_ft_s, ncol = 1, nrow = 2, align = "v", common.legend = TRUE, legend = "bottom")
  si_fig("usable_scans", "UsableScans", 18, 12, function() print(p_qa_stack))
}

combFCData   <- rbind(fcData, fetalFCData)
combDemoData <- rbind(goodData, fetalDemoData)
combDemoData$scanType <- factor(combDemoData$scanType)

## ---- Harmonization input ----------------------------------------------------
harmKeep <- combDemoData[combDemoData$fullID != "protect_092_01", ]      # excluded as in original script
harmFC   <- combFCData[rownames(combFCData) %in% harmKeep$fullID, ]
harmDemo <- harmKeep[harmKeep$fullID %in% rownames(combFCData), ]

harmDemo$group <- factor(ifelse(harmDemo$scanType == "Fetal", "Fetal", "Preterm"),
                         levels = c("Fetal", "Preterm"))                 # fixed effect
harmDemo$preterm <- as.numeric(harmDemo$group == "Preterm")              # switches the subgroup effect on for preterm only
# Preterm birth subgroup (same cutoffs as the analyses); fetal rows get a placeholder level
harmDemo$subgroup_pt <- factor(ifelse(harmDemo$group == "Preterm",
                                      as.character(birth_cat(harmDemo$birthPMA)), "none"))
stopifnot(!anyNA(harmDemo$subgroup_pt))

# Drop columns that must not enter the harmonization, attach connectivity
harmDemoTemp  <- subset(harmDemo, select = -c(finalVol, scanDate, kidokoro, birthPMA))
finalHarmData <- merge(harmDemoTemp, harmFC, by.x = "fullID", by.y = "row.names")
featCols      <- (ncol(harmDemoTemp) + 1):ncol(finalHarmData)

combat_result <- cached("longCombat", key = paste(table_key(finalHarmData[, -featCols]), sum(as.matrix(finalHarmData[, featCols]))),
                        longCombat(
  idvar    = "shortID",
  timevar  = "scanNo",
  batchvar = "scanner",
  features = featCols,
  formula  = "group + scanPMA",
  ranef    = "(1|shortID) + (0 + preterm | subgroup_pt)",
  data     = finalHarmData,
  verbose  = TRUE
))
harmonized <- combat_result$data_combat

## ---- Harmonization checks (supplement) -------------------------------------
# Goals: (1) scanner differences shrink; (2) PMA effects are preserved.
# Scanner comparisons use preterm scans only (fetal scans are all 1.5T).
if (RUN_HARMONIZATION_CHECKS) {
  # Match each pre-harmonization scan to its harmonized counterpart by ID
  key_pre  <- paste0(finalHarmData$shortID, "_0", finalHarmData$scanNo)
  key_post <- paste0(harmonized$shortID,    "_0", harmonized$scanNo)
  pre_x  <- as.matrix(finalHarmData[, featCols])
  post_x <- as.matrix(harmonized[match(key_pre, key_post), paste0(names(finalHarmData)[featCols], ".combat")])
  stopifnot(!anyNA(post_x))

  isPreterm <- finalHarmData$group == "Preterm"
  scn       <- factor(finalHarmData$scanner[isPreterm])

  # (1) Scanner effect per edge: Cohen's d (3T vs 1.5T) before and after
  scan_d <- function(x) {
    a <- x[isPreterm, ][scn == levels(scn)[2], ]
    b <- x[isPreterm, ][scn == levels(scn)[1], ]
    (colMeans(a) - colMeans(b)) / sqrt((apply(a, 2, var) + apply(b, 2, var)) / 2)
  }
  d_pre <- scan_d(pre_x); d_post <- scan_d(post_x)
  print(c(median_abs_d_pre = median(abs(d_pre)), median_abs_d_post = median(abs(d_post))))
  print(rbind(pre = range(d_pre), post = range(d_post)))

  p_scan_d <- ggplot(data.frame(d_pre, d_post)) +
    scale_y_continuous(limits = c(0, 750), breaks = seq(0, 800, 100), expand = c(0, 0)) +
    geom_histogram(aes(x = d_pre),  color = "#6ca7ed", fill = "#6ca7ed70", bins = 65, linewidth = 1) +
    geom_histogram(aes(x = d_post), color = "#e5af74", fill = "#e5af7470", bins = 65, linewidth = 1) +
    geom_vline(xintercept = 0, linewidth = 1, linetype = "dashed") +
    xlab("Scanner difference (Cohen's d)") + ylab("Count") + theme_paper(legend = "bottom")
  si_fig("harmonization", "ScannerDifference", 6, 6, function() print(p_scan_d))

  # (2) PCA of preterm scans colored by scanner, before vs. after
  par(mfrow = c(1, 2))
  for (k in 1:2) {
    pc <- prcomp(list(pre_x[isPreterm, ], post_x[isPreterm, ])[[k]], scale. = TRUE)$x
    plot(pc[, 1:2], col = as.numeric(scn) + 1, pch = 16, cex = .7,
         main = c("Before harmonization", "After harmonization")[k])
  }
  legend("topright", levels(scn), col = 2:3, pch = 16)
  par(mfrow = c(1, 1))

  # (3) PMA effect per edge before vs. after (PMA was protected -> identity line)
  pma_effect <- function(x) apply(x, 2, function(y) coef(lm(y ~ finalHarmData$scanPMA))[2])
  bioPlot <- data.frame(b_pre = pma_effect(pre_x), b_post = pma_effect(post_x))
  print(c(cor_PMA_effect_pre_post = cor(bioPlot$b_pre, bioPlot$b_post)))
  p_bio <- ggplot(bioPlot, aes(b_pre, b_post)) +
    scale_x_continuous(limits = c(-.036, .036), breaks = seq(-.04, .04, .01), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-.036, .036), breaks = seq(-.04, .04, .01), expand = c(0, 0)) +
    geom_point(shape = 21, color = "#75151E", fill = "#75151E55") +
    geom_abline(intercept = 0, slope = 1, linewidth = 1, linetype = "dashed") +
    xlab("Pre-Harmonization Beta for PMA at Scan") + ylab("Post-Harmonization Beta for PMA at Scan") +
    theme_paper(grid_x = TRUE)
  si_fig("harmonization", "PmaPreservation", 6, 6, function() print(p_bio))
}

## ---- Harmonized data replaces the raw data from here on ---------------------
goodData <- harmDemo
fcData   <- harmonized[, grep("\\.combat$", names(harmonized))]
rownames(fcData) <- paste0(harmonized$shortID, "_0", harmonized$scanNo)
colnames(fcData) <- sub("\\.combat$", "", colnames(fcData))
stopifnot(ncol(fcData) == length(featCols))
rm(harmKeep, harmFC, harmDemo, harmDemoTemp, finalHarmData, combat_result, harmonized)


###############################################################################
# 3. ANALYSIS DATA
###############################################################################

## ---- 3a. Connection table (ROI, lobe, hemisphere for every edge) -----------
goodFCData <- fcData[rownames(fcData) %in% goodData$fullID, ]
goodFCData <- goodFCData[, !grepl("Verm", colnames(goodFCData))]        # drop vermis
allEdges   <- colnames(goodFCData)

parts <- strsplit(allEdges, "_", fixed = TRUE)
roi1  <- gsub(".", "-", sapply(parts, `[`, 1), fixed = TRUE)            # "PreCG.L" -> "PreCG-L"
roi2  <- gsub(".", "-", sapply(parts, `[`, 2), fixed = TRUE)
i1 <- match(roi1, roiData$FullLabel); i2 <- match(roi2, roiData$FullLabel)
stopifnot(!anyNA(i1), !anyNA(i2))                                       # every ROI must be in the label key
columnMembership <- data.frame(Connection = allEdges, ROI1 = roi1, ROI2 = roi2,
                               Lobe1 = roiData$Lobe[i1],   Lobe2 = roiData$Lobe[i2],
                               Cortical1 = roiData$Cortical[i1], Cortical2 = roiData$Cortical[i2],
                               stringsAsFactors = FALSE)
columnMembershipCortOnly <- subset(columnMembership, Cortical1 == 1 & Cortical2 == 1)
allConn <- columnMembershipCortOnly$Connection                          # the 3,003 cortical edges
stopifnot(length(allConn) == 3003)

# Full region names for table/figure labels, e.g. "ORBmid.R" -> "Orbitofrontal cortex (middle) (R)"
full_name <- function(x) {
  r <- roiData[roiData$FullLabel == gsub(".", "-", x, fixed = TRUE), ][1, ]
  paste0(r$LongLabel, " (", r$Hemi, ")")
}

## ---- 3b. Scan table: one row per scan, aligned with goodFCData -------------
idx <- match(rownames(goodFCData), goodData$fullID)
stopifnot(!anyNA(idx))
scanInfo <- data.frame(longID   = rownames(goodFCData),
                       shortID  = factor(goodData$shortID[idx]),
                       scanType = factor(goodData$scanType[idx], levels = c("Fetal", "Preterm")),  # Fetal = reference
                       BirthPMA = goodData$birthPMA[idx],
                       ScanPMA  = goodData$scanPMA[idx],
                       scanner  = factor(goodData$scanner[idx]),        # "450" = 1.5T (reference), "750" = 3T
                       sex      = goodData$sex[idx],
                       finalVol = goodData$finalVol[idx],
                       stringsAsFactors = FALSE)
scn_term <- paste0("scanner", levels(scanInfo$scanner)[2])             # coefficient name for 3T - 1.5T
is_pt <- scanInfo$scanType == "Preterm"
scanInfo$BirthCat <- birth_cat(ifelse(is_pt, scanInfo$BirthPMA, NA))   # NA for fetal scans

## ---- 3c. Preterm scan table (used for everything in sections 5-9) ----------
pt <- droplevels(scanInfo[is_pt, ])
stopifnot(!anyNA(pt$BirthCat),
          all(tapply(pt$BirthPMA, pt$shortID, function(x) diff(range(x))) < 1e-8))   # one birth PMA per infant

# Fingerprint of the analysis data; used to name cached results
data_key <- paste(table_key(scanInfo), table_key(pt), nrow(goodFCData), ncol(goodFCData), sum(as.matrix(goodFCData)))

## ---- 3d. Table 1: sample description (infant level) ------------------------
group_counts <- function(d) rbind(scans   = table(d$BirthCat),
                                  infants = tapply(d$shortID, d$BirthCat, function(x) length(unique(x))))
group_counts(pt)                                  # scans and infants per group (EPT / VPT / MPT)
c(scans = nrow(pt), infants = nlevels(pt$shortID))
table(pt$BirthCat, pt$scanner)                    # scans by group and scanner

# Infants by scanner use: only 1.5T / only 3T / both
scn_use   <- table(pt$shortID, pt$scanner)
cls_levels <- c("1.5T only", "3T only", "both")
use_class <- factor(ifelse(scn_use[, 1] > 0 & scn_use[, 2] > 0, "both", ifelse(scn_use[, 1] > 0, "1.5T only", "3T only")),
                    levels = cls_levels)
names(use_class) <- rownames(scn_use)
grp_of_inf <- tapply(as.character(pt$BirthCat), pt$shortID, `[`, 1)[names(use_class)]
rbind(All = c(table(use_class)), t(sapply(c("EPT", "VPT", "MPT"), function(g) c(table(use_class[grp_of_inf == g])))))

# Demographics: one row per infant, birth PMA from the corrected scan table
inf <- pt[!duplicated(pt$shortID), c("shortID", "BirthPMA", "BirthCat")]
inf <- merge(inf, demoData[, setdiff(names(demoData), "BirthPMA")], by.x = "shortID", by.y = "Subject", all.x = TRUE)
matEd <- grep("^Maternal.Education", names(demoData), value = TRUE)
inf$MatEd <- inf[[matEd[1]]]
ms <- function(x) sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
md <- function(x) sprintf("%.0f [%.0f-%.0f]", median(x, na.rm = TRUE), min(x, na.rm = TRUE), max(x, na.rm = TRUE))
by_group <- function(f, v) c(All = f(inf[[v]]), sapply(c("EPT", "VPT", "MPT"), function(g) f(inf[[v]][inf$BirthCat == g])))
table1 <- rbind(
  `Birth GA, wk, mean (SD)`        = by_group(ms, "BirthPMA"),
  `Birth weight, g, mean (SD)`     = by_group(ms, "BirthWeight"),
  `Birth length, cm, mean (SD)`    = by_group(ms, "BirthLength"),
  `APGAR 1 min, median [range]`    = by_group(md, "APGAR.1min"),
  `APGAR 5 min, median [range]`    = by_group(md, "APGAR.5min"),
  `Maternal age, y, mean (SD)`     = by_group(ms, "Maternal.Age"),
  `Parity, median [range]`         = by_group(md, "Parity"),
  `Male sex, %`                    = by_group(function(x) sprintf("%.0f", mean(x == 1, na.rm = TRUE) * 100), "Gender..0.F..1.M."))
table1
# Maternal education distribution (% of infants with known education; 0 = unknown)
round(prop.table(table(inf$MatEd[inf$MatEd != 0], useNA = "no")) * 100, 1)

# Scan-level description
c(PMA_at_scan_mean = mean(pt$ScanPMA), PMA_sd = sd(pt$ScanPMA), PMA_min = min(pt$ScanPMA), PMA_max = max(pt$ScanPMA),
  scans_per_infant_median = median(table(pt$shortID)), pct_scans_3T = mean(pt$scanner == "750") * 100)

## ---- 3e. Usable volumes per scan, from the ROI time series -----------------------------------------------------
# The scan QA table's finalVol and the ROI time series disagree for about 1 scan in 5 (426 of 538 match exactly: fetal 186
# of 192, preterm 240 of 346; the rest differ by -87 to +109 volumes). The stored connectivity of EVERY scan (all 538, r =
# 1.000) is reproduced from the non-zero rows of its time series, so that count (scan_vol) is the number of volumes the
# connectivity is based on. It is used for Table 1 and the volume controls (9e, 9k); finalVol is kept only for the QA
# plots in section 1. Every scan meets the inclusion rule with either count (preterm >= 120, fetal >= 80 volumes).
# Time series: one tab-separated .1D file per scan ("<fullID>_ts.1D"; rows = volumes, columns = the 111 ROIs, censored
# volumes are all-zero rows) in data/rawTS and/or data/rawRS (subfolders fetal and preterm).
ts_dirs <- file.path(dataDir, c("rawTS", "rawRS")); ts_dirs <- ts_dirs[dir.exists(ts_dirs)]
ts_files <- if (length(ts_dirs) > 0) list.files(ts_dirs, pattern = "_ts\\.1D$", recursive = TRUE, full.names = TRUE) else character(0)
names(ts_files) <- sub("_ts\\.1D$", "", basename(ts_files))
read_ts <- function(f) {                                              # AFNI .1D text: tab-separated, first field empty
  fields <- lapply(strsplit(readLines(f), "\t", fixed = TRUE), function(x) x[nzchar(x)])
  stopifnot(length(unique(lengths(fields))) == 1)
  matrix(as.numeric(unlist(fields)), nrow = length(fields), byrow = TRUE)
}
if (length(ts_files) > 0 && all(scanInfo$longID %in% names(ts_files))) {
  stopifnot(!anyDuplicated(names(ts_files)))
  scan_vol <- cached("usableVolumes", key = paste(nrow(scanInfo), length(ts_files), sum(file.size(ts_files[scanInfo$longID]))),
                     unname(vapply(scanInfo$longID, function(id) sum(rowSums(abs(read_ts(ts_files[[id]]))) > 0), numeric(1))))
  print(c(scans_where_usable_ne_finalVol = sum(scan_vol != scanInfo$finalVol)))
} else {
  message("Time series not found for every scan: using finalVol from the QA table for the volume counts (3e, Table 1, 9e).")
  scan_vol <- scanInfo$finalVol
}

# Table 1 data-amount rows (scan level): usable volumes (3e) and minutes of data per scan, mean (SD).
# Minutes = volumes x TR / 60 (TR = 2 s preterm, 3 s fetal). Columns: fetal, all preterm, then preterm groups.
TR_SEC  <- c(Fetal = 3, Preterm = 2)
scan_min <- scan_vol * TR_SEC[as.character(scanInfo$scanType)] / 60   # kept out of scanInfo so cache keys do not change
# Motion: mean framewise displacement per scan (motionSummary.csv: project, fullID, shortID, scan, meanFD)
motion <- read.csv("motionSummary.csv", stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(motion$fullID))
scan_fd <- motion$meanFD[match(scanInfo$longID, motion$fullID)]                 # kept out of scanInfo so cache keys do not change
if (anyNA(scan_fd)) message(sum(is.na(scan_fd)), " of ", length(scan_fd), " scans have no motion value (fetal: ",
                            sum(is.na(scan_fd) & scanInfo$scanType == "Fetal"), ", preterm: ", sum(is.na(scan_fd) & scanInfo$scanType == "Preterm"), ")")
scan_cols <- list(Fetal = scanInfo$scanType == "Fetal", `All preterm` = scanInfo$scanType == "Preterm",
                  EPT = scanInfo$BirthCat %in% "EPT", VPT = scanInfo$BirthCat %in% "VPT", MPT = scanInfo$BirthCat %in% "MPT")
table1_scan <- rbind(
  `Scans, n`                          = sapply(scan_cols, sum),
  `Volumes per scan, mean (SD)`       = sapply(scan_cols, function(i) ms(scan_vol[i])),
  `Minutes of data per scan, mean (SD)` = sapply(scan_cols, function(i) sprintf("%.1f (%.1f)", mean(scan_min[i]), sd(scan_min[i]))),
  `Mean framewise displacement, mm, mean (SD)` = sapply(scan_cols, function(i) sprintf("%.2f (%.2f)", mean(scan_fd[i], na.rm = TRUE), sd(scan_fd[i], na.rm = TRUE))),
  `Scans with motion data, n`         = sapply(scan_cols, function(i) sum(!is.na(scan_fd[i]))))
table1_scan

# Table 1 in two versions with identical rows and columns (EPT, VPT, MPT, Total preterm, Fetal):
#   Table1_ByScan        every scan is a row: scan-based and infant-level variables are summarized over scans, so an
#                        infant with three scans is counted three times (as in the original Table 1). n = infants.
#   Table1_ByParticipant every infant is one row: scan-based variables (volumes, minutes, motion, PMA at scan) are first
#                        averaged over the infant's scans, and birth, sex, APGAR and clinical items are counted once.
# Birth and clinical items come from the clinical workbook (`clin_file`, edit the name below), which has two sheets:
#   PRETERM: ID, BirthWeight, APGAR 1min, APGAR 5min, Clinical Infection, Bronchopulmonary dysplasia, Patent ductus arteriosus,
#            Retinopathy of prematurity, NEC, Days of mechanical ventilation (yes/no items coded 1/0)
#   FETAL:   ID, BirthWeight, Apgar1, Apgar5
# Columns are found by name pattern, and ID must equal the infant ID used in the connectivity tables (shortID); the number of
# infants matched is printed and the script stops if too few match. These are infant-level facts: the by-scan table counts an
# infant once per scan, the by-participant table once. Fetal infants have APGAR but no clinical items ("-").
male_scan <- if (is.numeric(scanInfo$sex)) ifelse(scanInfo$sex %in% c(0, 1), scanInfo$sex, NA) else {
  s1 <- toupper(substr(trimws(as.character(scanInfo$sex)), 1, 1)); ifelse(s1 == "M", 1, ifelse(s1 == "F", 0, NA)) }     # 0 = female, 1 = male
scan_df <- data.frame(shortID = as.character(scanInfo$shortID),
                      grp = factor(ifelse(scanInfo$scanType == "Fetal", "Fetal", as.character(scanInfo$BirthCat)), c("Fetal", "EPT", "VPT", "MPT")),
                      pma = scanInfo$ScanPMA, birthPMA = scanInfo$BirthPMA, male = male_scan,
                      vol = scan_vol, mins = scan_min, fd = scan_fd, stringsAsFactors = FALSE)
if (anyDuplicated(unique(scan_df[c("shortID", "grp")])$shortID)) message("Some infants appear in more than one group; infant-level counts use each infant's first group.")
inf_t <- do.call(rbind, lapply(split(scan_df, scan_df$shortID), function(d) data.frame(
  shortID = d$shortID[1], grp = d$grp[1], n_scans = nrow(d), pma = mean(d$pma), birthPMA = d$birthPMA[1], male = d$male[1],
  vol = mean(d$vol), mins = mean(d$mins), fd = if (all(is.na(d$fd))) NA_real_ else mean(d$fd, na.rm = TRUE), stringsAsFactors = FALSE)))
clin_file <- "2026-08-03_KC_fetal_preterm.xlsx"      # clinical workbook (sheets PRETERM and FETAL); relative to the working directory (the data folder), or a full path
if (!file.exists(clin_file)) stop("Clinical workbook not found: ", clin_file, ". Set clin_file (sheets PRETERM and FETAL).")
if (!requireNamespace("readxl", quietly = TRUE)) stop("Install the readxl package: install.packages('readxl')")
read_clin <- function(sheet) { d <- as.data.frame(readxl::read_excel(clin_file, sheet = sheet)); names(d) <- trimws(gsub("\\s+", " ", names(d))); d }
clin_pick <- function(d, pat) { h <- grep(pat, names(d), ignore.case = TRUE, value = TRUE); if (length(h)) h[1] else NA_character_ }
clin_pt <- read_clin("PRETERM"); clin_ft <- read_clin("FETAL")
pt_cols <- c(`Birth weight (g)` = "birth.?weight", `APGAR 1 minute` = "apgar.*1", `APGAR 5 minute` = "apgar.*5", `Bronchopulmonary dysplasia` = "bronch", `Patent ductus arteriosus` = "ductus",
             `Retinopathy of prematurity` = "retinopath", `Clinical infection` = "infect", `Necrotizing enterocolitis` = "necrotiz|^nec", `Days on ventilation` = "ventil")
ft_cols <- c(`Birth weight (g)` = "birth.?weight", `APGAR 1 minute` = "apgar.*1", `APGAR 5 minute` = "apgar.*5")
clin_map <- function(d, cols) { m <- sapply(cols, function(p) clin_pick(d, p)); m }
pt_map <- clin_map(clin_pt, pt_cols); ft_map <- clin_map(clin_ft, ft_cols)
print(data.frame(item = names(pt_map), `PRETERM column` = unname(pt_map), check.names = FALSE), row.names = FALSE)   # check: each item should name the right column
stopifnot(!anyNA(pt_map), !anyNA(ft_map))
id_pt <- as.character(inf_t$shortID[inf_t$grp != "Fetal"]); id_ft <- as.character(inf_t$shortID[inf_t$grp == "Fetal"])
sheet_id <- function(d) trimws(as.character(d[[clin_pick(d, "^id$")]]))
m_pt <- match(id_pt, sheet_id(clin_pt)); m_ft <- match(id_ft, sheet_id(clin_ft))
message("Clinical workbook match: PRETERM ", sum(!is.na(m_pt)), " of ", length(id_pt), " infants; FETAL ", sum(!is.na(m_ft)), " of ", length(id_ft), " infants",
        "; repeated IDs: PRETERM ", sum(duplicated(sheet_id(clin_pt))), ", FETAL ", sum(duplicated(sheet_id(clin_ft))), " (the first non-empty value per infant is used)")
if (mean(!is.na(m_pt)) < 0.8) stop("Fewer than 80% of preterm infants match the PRETERM sheet's ID column: the IDs probably need reformatting (compare head(sheet_id(clin_pt)) with head(id_pt)).")
clin_items <- names(pt_cols)
# Cell coding: yes/no items may be typed as 1/0, Y/N or yes/no (blank = missing); ventilation days may carry text ("22 days").
# Entries that cannot be read are set to missing and listed below, so nothing is lost silently.
code_cell <- function(x, kind) {
  x <- trimws(as.character(x)); x[x %in% c("", "NA", "N/A", "na", "n/a", ".")] <- NA
  if (kind == "num") return(suppressWarnings(as.numeric(sub("^[^0-9.-]*(-?[0-9]*\\.?[0-9]+).*$", "\\1", x))))
  v <- tolower(x)
  ifelse(is.na(v), NA_real_, ifelse(grepl("^(1|y|yes|true)\\b", v), 1, ifelse(grepl("^(0|n|no|false)\\b", v), 0, NA_real_)))
}
item_kind <- function(k) if (grepl("^APGAR|ventilation|weight", k)) "num" else "yn"
unread <- list()
by_id <- function(d, col, ids, kind, item) {                           # first non-empty value per infant (the sheet may have several rows per infant)
  raw <- trimws(as.character(d[[col]])); x <- code_cell(raw, kind)
  bad <- unique(raw[!is.na(raw) & raw != "" & is.na(x)]); if (length(bad)) unread[[paste(item, nrow(d))]] <<- bad
  v <- sapply(split(x, sheet_id(d)), function(z) if (all(is.na(z))) NA_real_ else z[!is.na(z)][1])
  unname(v[ids])
}
clin_inf <- data.frame(shortID = c(id_pt, id_ft), stringsAsFactors = FALSE)
for (k in clin_items) clin_inf[[k]] <- c(by_id(clin_pt, pt_map[[k]], id_pt, item_kind(k), k), if (k %in% names(ft_map)) by_id(clin_ft, ft_map[[k]], id_ft, item_kind(k), k) else rep(NA_real_, length(id_ft)))
if (length(unread)) { cat("\nEntries that could not be read as 1/0 or a number (set to missing); check them:\n"); print(unread) }
print(round(100 * sapply(clin_inf[clin_items], function(x) mean(is.na(x[seq_along(id_pt)]))), 1))        # % of preterm infants with no value, per item
as_yes <- function(x) ifelse(is.na(x), NA, x > 0)
for (k in setdiff(clin_items, c("Days on ventilation", "APGAR 1 minute", "APGAR 5 minute"))) { cat("\n", k, ":\n"); print(table(clin_inf[[k]][seq_along(id_pt)], useNA = "ifany")) }   # check the coding
# Coverage check: infants with a value for each item, by group (and how many infants the group has), so low coverage is visible.
clin_cov <- do.call(rbind, lapply(clin_items, function(k) {
  v <- clin_inf[[k]][match(inf_t$shortID, clin_inf$shortID)]
  cov_grp <- list(EPT = inf_t$grp == "EPT", VPT = inf_t$grp == "VPT", MPT = inf_t$grp == "MPT", Fetal = inf_t$grp == "Fetal")
  data.frame(Item = k, t(sapply(cov_grp, function(i) sprintf("%d of %d", sum(!is.na(v[i])), sum(i)))), check.names = FALSE) }))
rownames(clin_cov) <- NULL
cat("\nInfants with a value, by item and group:\n"); print(clin_cov, row.names = FALSE)
if (SAVE_TABLES) write.csv(clin_cov, file.path(figureDir, "Table1_ClinicalCoverage.csv"), row.names = FALSE)
inf_t   <- merge(inf_t,   clin_inf, by = "shortID", all.x = TRUE)       # infant-level items; fetal infants get NA
scan_df <- merge(scan_df, clin_inf, by = "shortID", all.x = TRUE)       # each scan carries its infant's values

msd <- function(x, d = 1) sprintf(paste0("%.", d, "f (%.", d, "f)"), mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
npc <- function(x) if (all(is.na(x))) "no data" else sprintf("%d (%.0f%%)", sum(x, na.rm = TRUE), 100 * mean(x, na.rm = TRUE))
grp_sel <- function(g) list(EPT = g == "EPT", VPT = g == "VPT", MPT = g == "MPT", Total = g != "Fetal", Fetal = g == "Fetal")
make_table1 <- function(d) {                                           # d = scan_df (by scan) or inf_t (by participant)
  sel  <- grp_sel(d$grp)
  ic   <- function(f) sapply(sel, function(i) f(d[i, ]))
  pre  <- function(f) ic(function(x) if (all(x$grp == "Fetal")) "-" else f(x))
  rows <- list(
    `n`                           = sapply(sel, function(i) length(unique(d$shortID[i]))),
    `Total scans`                 = sapply(grp_sel(scan_df$grp), function(i) sum(i)),
    `Mean number of scans`        = sapply(grp_sel(inf_t$grp), function(i) sprintf("%.1f", mean(inf_t$n_scans[i]))),
    `Mean Volumes`                = ic(function(x) msd(x$vol)),
    `Mean Scan Time (min)`        = ic(function(x) msd(x$mins)),
    `Mean Framewise Displacement` = ic(function(x) msd(x$fd, 2)),
    `PMA at scan (Weeks)`         = ic(function(x) msd(x$pma, 2)),
    `PMA at birth (Weeks)`        = ic(function(x) msd(x$birthPMA, 2)),
    `Sex (Male)`                  = ic(function(x) npc(x$male == 1)))
  for (k in clin_items)                                                # APGAR and clinical items: mean (SD) for APGAR and ventilation days, n (%) for yes/no items
    rows[[k]] <- if (grepl("weight", k)) ic(function(x) msd(x[[k]], 0)) else if (grepl("^APGAR", k)) ic(function(x) msd(x[[k]], 2)) else if (grepl("ventilation", k)) pre(function(x) msd(x[[k]], 2)) else pre(function(x) npc(as_yes(x[[k]])))
  nm <- setdiff(names(rows), "Birth weight (g)")                      # birth weight goes right after PMA at birth
  rows <- rows[append(nm, "Birth weight (g)", after = match("PMA at birth (Weeks)", nm))]
  m <- do.call(rbind, lapply(rows, as.character)); dimnames(m) <- list(names(rows), names(sel)); m
}
Table1_ByScan        <- make_table1(scan_df)
Table1_ByParticipant <- make_table1(inf_t)
stopifnot(all(as.numeric(Table1_ByParticipant["n", c("EPT", "VPT", "MPT")]) == c(nlevels(droplevels(pt$shortID[pt$BirthCat == "EPT"])), nlevels(droplevels(pt$shortID[pt$BirthCat == "VPT"])), nlevels(droplevels(pt$shortID[pt$BirthCat == "MPT"])))))
cat("\nTable 1, by scan (n = infants; other rows over scans):\n");  print(Table1_ByScan, quote = FALSE)
cat("\nTable 1, by participant (each infant once):\n");             print(Table1_ByParticipant, quote = FALSE)

# Fetal vs. preterm sample comparison (PMA at scan; sex)
inf_all <- scanInfo[!duplicated(scanInfo$shortID), ]
t.test(ScanPMA ~ scanType, data = scanInfo, var.equal = TRUE)
tryCatch(chisq.test(table(inf_all$sex, inf_all$scanType)[2:3, 1:2]), error = function(e) message("Sex chi-square skipped: ", e$message))


###############################################################################
# 4. FETAL vs. PRETERM, ALL EDGES
#    Model (per edge):  FC ~ PMA + scanner + group,  random intercept per infant
#    group coefficient = Preterm - Fetal  (positive = higher FC in preterm)
###############################################################################

# Fits the model for one edge. `extra` adds numeric covariates (section 9e).
# Group means are the model-adjusted means (equivalent to emmeans): covariates at
# their sample mean, scanner averaged with equal weight.
fit_fetal_edge <- function(cn, extra = character(0), data = scanInfo) {
  d <- data; d$FC <- as.numeric(goodFCData[d$longID, cn])
  na_row <- data.frame(Connection = cn, scanPMAb = NA, scanPMAp = NA, scannerb = NA, scannerp = NA,
                       scanTypeb = NA, scanTypese = NA, scanTypep = NA,
                       emmMean.Fetal = NA, emmSE.Fetal = NA, emmMean.Preterm = NA, emmSE.Preterm = NA)
  m <- tryCatch(lme(reformulate(c("ScanPMA", extra, "scanner", "scanType"), response = "FC"),
                    random = ~ 1 | shortID, data = d, control = ctrl), error = function(e) NULL)
  if (is.null(m)) return(na_row)
  tt <- summary(m)$tTable; b <- nlme::fixef(m); V <- vcov(m)
  L <- function(group) {
    x <- setNames(numeric(length(b)), names(b)); x["(Intercept)"] <- 1
    for (v in c("ScanPMA", extra)) x[v] <- mean(d[[v]], na.rm = TRUE)
    x[scn_term] <- 0.5; x["scanTypePreterm"] <- as.numeric(group == "Preterm"); x
  }
  emm <- function(group) { x <- L(group); c(sum(x * b), sqrt(drop(t(x) %*% V %*% x))) }
  eF <- emm("Fetal"); eP <- emm("Preterm")
  data.frame(Connection = cn,
             scanPMAb = tt["ScanPMA", "Value"],  scanPMAp = tt["ScanPMA", "p-value"],
             scannerb = tt[scn_term, "Value"],   scannerp = tt[scn_term, "p-value"],
             scanTypeb = tt["scanTypePreterm", "Value"], scanTypese = tt["scanTypePreterm", "Std.Error"],
             scanTypep = tt["scanTypePreterm", "p-value"],
             emmMean.Fetal = eF[1], emmSE.Fetal = eF[2], emmMean.Preterm = eP[1], emmSE.Preterm = eP[2])
}

fetalRes <- cached("fetalRes", key = paste(data_key, "FC~ScanPMA+scanner+scanType"),
                   do.call(rbind, par_lapply(allConn, fit_fetal_edge,
                                             export = c("goodFCData", "ctrl", "scn_term", "scanInfo"))))
fetalRes$scanPMAp.Adj  <- p.adjust(fetalRes$scanPMAp,  method = "fdr")
fetalRes$scanTypep.Adj <- p.adjust(fetalRes$scanTypep, method = "fdr")

sigFetal <- fetalRes$scanTypep.Adj < .05 & !is.na(fetalRes$scanTypep.Adj)
c(failed_fits       = sum(is.na(fetalRes$scanTypep)),
  uncorrected_p05   = sum(fetalRes$scanTypep < .05, na.rm = TRUE),
  fdr05             = sum(sigFetal),
  pct_of_3003       = round(mean(sigFetal) * 100, 1),
  preterm_greater   = sum(sigFetal & fetalRes$scanTypeb > 0),
  fetal_greater     = sum(sigFetal & fetalRes$scanTypeb < 0))
c(pct_preterm_greater = round(mean(fetalRes$scanTypeb[sigFetal] > 0) * 100, 1))
sum(fetalRes$scanPMAp.Adj < .05, na.rm = TRUE)                         # edges with a PMA effect in this model

## ---- 4a. Where the differences are: hemisphere and lobe tallies ------------
tally_pairs <- function(conns) {
  cm <- columnMembership[match(conns, columnMembership$Connection), ]
  hemi <- function(x) factor(substr(x, nchar(x), nchar(x)), levels = c("R", "L"))
  lobe <- function(x) factor(substr(x, 1, 1), levels = c("F", "T", "P", "O"))
  list(hemisphere = table(ROI1 = hemi(cm$ROI1), ROI2 = hemi(cm$ROI2)),
       lobe       = table(Lobe1 = lobe(cm$Lobe1), Lobe2 = lobe(cm$Lobe2)))
}
tally_preterm_greater <- tally_pairs(fetalRes$Connection[sigFetal & fetalRes$scanTypeb > 0])
tally_fetal_greater   <- tally_pairs(fetalRes$Connection[sigFetal & fetalRes$scanTypeb < 0])
tally_preterm_greater; tally_fetal_greater

## ---- 4b. Fig. 2 matrices: triangular heatmaps of the significant edges -----
# One panel per group. Cells show |beta| of the preterm-minus-fetal difference for edges
# that are FDR-significant and larger in that group (fetal: beta < 0, preterm: beta > 0).
# ROIs are ordered by lobe; lobe boundaries are drawn as lines and each lobe gets a
# colored sidebar (left and bottom). Panel letters and titles are added by hand.
lobeOrder <- c("Frontal", "Temporal", "Parietal", "Occipital")
cp <- data.frame(Conn = fetalRes$Connection,
                 ROI1 = sapply(strsplit(fetalRes$Connection, "_", fixed = TRUE), `[`, 1),
                 ROI2 = sapply(strsplit(fetalRes$Connection, "_", fixed = TRUE), `[`, 2),
                 B = fetalRes$scanTypeb, pAdj = fetalRes$scanTypep.Adj, stringsAsFactors = FALSE)

cort       <- roiData[roiData$Cortical == 1, ]
cort       <- cort[order(factor(cort$Lobe, levels = lobeOrder)), ]            # stable within lobe
roi_levels <- gsub("-", ".", cort$FullLabel, fixed = TRUE)                      # "PreCG-L" -> "PreCG.L"
stopifnot(all(c(cp$ROI1, cp$ROI2) %in% roi_levels))
n_roi      <- length(roi_levels)
lobe_runs  <- rle(as.character(cort$Lobe))                                      # contiguous lobe blocks
lobe_end   <- cumsum(lobe_runs$lengths); lobe_start <- lobe_end - lobe_runs$lengths + 1
heat_lim   <- ceiling(max(abs(fetalRes$scanTypeb[sigFetal])) * 100) / 100       # top of the color scale

# Matrix layout settings (tune here)
tri_strip_thick <- 2      # lobe sidebar thickness, in matrix cells (the homotopic chart's strip is about 3% of its panel width)
tri_strip_gap   <- 1      # gap between the axis line and the sidebar, in cells
tri_strip_split <- 0.8    # small white gap between the lobe segments of each sidebar, in cells (0 = none)
tri_bar_len_in  <- 1.5    # colorbar length, inches
tri_bar_w_in    <- 0.17   # colorbar width, inches
tri_axis_lw     <- 0.5    # line width of the left/bottom axis lines and the hypotenuse

tri_heat <- function(keep, ramp, lim = heat_lim, margin_pt = c(8, 8, 46, 50)) {
  keep[is.na(keep)] <- FALSE
  d <- cp[keep, ]
  i <- match(d$ROI1, roi_levels); j <- match(d$ROI2, roi_levels)
  d$x <- pmin(i, j); d$y <- pmax(i, j)                                         # lower triangle: y (row) >= x (column)
  d$beta <- abs(d$B)
  bnd <- lobe_end[-length(lobe_end)] + 0.5                                     # lobe boundaries
  gap <- tri_strip_gap; thick <- tri_strip_thick                               # sidebar offset and thickness (cell units)
  mid <- (lobe_start + lobe_end) / 2
  K   <- length(lobe_start)
  sp  <- tri_strip_split / 2
  seg_lo <- lobe_start - 0.5 + c(0, rep(sp, K - 1))                           # segment start, trimmed where it meets a neighbor
  seg_hi <- lobe_end + 0.5 - c(rep(sp, K - 1), 0)                             # segment end
  # Black frame and ticks on the colorbar: ggplot2 >= 3.5 styles them through the theme (the old
  # guide arguments are ignored there); older versions through the guide.
  new_gg    <- utils::packageVersion("ggplot2") >= "3.5.0"
  bar_guide <- guides(fill = if (new_gg) guide_colourbar() else guide_colourbar(frame.colour = "#000000", ticks.colour = "#000000"))
  bar_theme <- if (new_gg) theme(legend.frame = element_rect(colour = "#000000", fill = NA, linewidth = 0.6),
                                 legend.ticks = element_line(colour = "#000000", linewidth = 0.4)) else theme()
  ggplot(d) +
    geom_tile(aes(x, y, fill = beta), width = 1, height = 1) +
    scale_fill_gradient2(low = ramp[1], mid = ramp[2], high = ramp[3], midpoint = lim / 2,
                         limits = c(0, lim), breaks = seq(0.1, lim, by = 0.1)) +
    # axis lines (left and bottom edges), hypotenuse, and lobe boundaries (inside the triangle)
    annotate("segment", x = 0.5, xend = 0.5, y = 0.5, yend = n_roi + 0.5, linewidth = tri_axis_lw) +
    annotate("segment", x = 0.5, xend = n_roi + 0.5, y = n_roi + 0.5, yend = n_roi + 0.5, linewidth = tri_axis_lw) +
    annotate("segment", x = 0.5, xend = n_roi + 0.5, y = 0.5, yend = n_roi + 0.5, linewidth = tri_axis_lw) +
    annotate("segment", x = 0.5, xend = bnd, y = bnd, yend = bnd, linewidth = 0.25) +
    annotate("segment", x = bnd, xend = bnd, y = bnd, yend = n_roi + 0.5, linewidth = 0.25) +
    # lobe sidebars + labels: left (vertical) and bottom (horizontal)
    annotate("rect", xmin = 0.5 - gap - thick, xmax = 0.5 - gap, ymin = seg_lo, ymax = seg_hi,
             fill = lobe_pal[lobe_runs$values]) +
    annotate("rect", xmin = seg_lo, xmax = seg_hi, ymin = n_roi + 0.5 + gap, ymax = n_roi + 0.5 + gap + thick,
             fill = lobe_pal[lobe_runs$values]) +
    annotate("text", x = 0.5 - gap - thick - 2.2, y = mid, label = lobe_runs$values, angle = 90, size = 3.4) +
    annotate("text", x = mid, y = n_roi + 0.5 + gap + thick + 2.6, label = lobe_runs$values, size = 3.4) +
    scale_x_continuous(expand = c(0, 0)) + scale_y_reverse(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0.5, n_roi + 0.5), ylim = c(n_roi + 0.5, 0.5), expand = FALSE, clip = "off") +
    labs(fill = "\u03b2") +
    bar_guide +
    theme_void(base_size = 10) +
    theme(aspect.ratio = 1, plot.margin = margin(margin_pt[1], margin_pt[2], margin_pt[3], margin_pt[4], "pt"),
          # a colorbar is 5 x legend.key.height long, so divide the wanted length by 5
          legend.key.height = unit(tri_bar_len_in / 5, "in"), legend.key.width = unit(tri_bar_w_in, "in"),
          legend.title = element_text(size = 12), legend.text = element_text(size = 8)) +
    bar_theme +
    legend_inside(0.88, 0.80)
}
p_tri_fetal   <- tri_heat(cp$B < 0 & cp$pAdj < .05, pal$heat_fetal)
p_tri_preterm <- tri_heat(cp$B > 0 & cp$pAdj < .05, pal$heat_preterm)
show_fig("Fig2_part_MatrixFetal",   6, 6, function() print(p_tri_fetal))
show_fig("Fig2_part_MatrixPreterm", 6, 6, function() print(p_tri_preterm))

## ---- 4c. Fig. 2 chord diagrams (one file per group) ------------------------
# ROI order and gaps: right hemisphere then left (mirrored), grouped by lobe
ord       <- order(factor(columnMembershipCortOnly$Lobe1, levels = lobeOrder))
bilatROI  <- unique(sub("-[LR]$", "", columnMembershipCortOnly$ROI1[ord]))   # bilateral ROI names in lobe order
roiList   <- c(paste0(bilatROI, "-R"), rev(paste0(bilatROI, "-L")))
roi_hemi  <- substr(roiList, nchar(roiList), nchar(roiList))
roi_lobe  <- roiData$Lobe[match(roiList, roiData$FullLabel)]
gapSize   <- sapply(seq_along(roiList), function(i) {
  if (i == length(roiList) || roi_hemi[i] != roi_hemi[i + 1]) return(10)      # hemisphere break
  if (roi_lobe[i] == roi_lobe[i + 1]) 2 else 5                                # within lobe vs. between lobes
})

# Sector colors by lobe
roiColor <- unlist(lapply(names(lobe_rois), function(l)
  setNames(rep(lobe_pal[[l]], 2 * length(lobe_rois[[l]])), c(paste0(lobe_rois[[l]], "-L"), paste0(lobe_rois[[l]], "-R")))))

chord_df <- data.frame(From = columnMembershipCortOnly$ROI1, To = columnMembershipCortOnly$ROI2,
                       Value = as.numeric(sigFetal[match(allConn, fetalRes$Connection)]), stringsAsFactors = FALSE)
b_dir    <- fetalRes$scanTypeb[match(allConn, fetalRes$Connection)] > 0
chord_col_fetal   <- ifelse(b_dir, "#00000000", col_fetal)       # fetal-greater edges only
chord_col_preterm <- ifelse(b_dir, col_preterm, "#00000000")     # preterm-greater edges only

# Chord layout settings (tune here)
chord_canvas    <- 1.65   # half-width of the drawing canvas; raise it if labels are still cut off (the circle gets smaller)
chord_label_cex <- 0.725  # label size
chord_label_gap <- 1.5    # mm between a colored sector and its label: the SAME for every label, left and right

# Sector labels use roi_display (no hemisphere). Each label starts a fixed distance (mm) beyond its
# colored sector and runs outward; with adj = c(0, .5) and niceFacing, circlize mirrors the anchor for
# the flipped (left-side) labels, so their near end sits the same distance from the sector.
draw_chord <- function(cols, label_cex = chord_label_cex) {
  circos.clear()
  circos.par(start.degree = 86, gap.after = gapSize,
             canvas.xlim = c(-chord_canvas, chord_canvas), canvas.ylim = c(-chord_canvas, chord_canvas))
  chordDiagram(chord_df, order = roiList, grid.col = roiColor, col = cols, symmetric = TRUE, scale = TRUE,
               annotationTrack = "grid", annotationTrackHeight = c(0.03, 0.03))
  circos.track(track.index = 1, bg.border = NA, panel.fun = function(x, y)
    circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + mm_y(chord_label_gap), display_name(CELL_META$sector.index),
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = label_cex))
  circos.clear()
}
show_fig("Fig2_part_ChordFetal",   9, 9, function() draw_chord(chord_col_fetal))
show_fig("Fig2_part_ChordPreterm", 9, 9, function() draw_chord(chord_col_preterm))

## ---- 4d. Table 2 + data for the homotopic chart --------------------------------
# Homotopic edges = the same region in the left and right hemispheres.
isHomotopic <- sub("-[LR]$", "", columnMembershipCortOnly$ROI1) == sub("-[LR]$", "", columnMembershipCortOnly$ROI2)
homo <- fetalRes[match(allConn[isHomotopic], fetalRes$Connection), ]
homo$ROI <- sub("-[LR]$", "", columnMembershipCortOnly$ROI1[isHomotopic])
homo <- homo[match(bilatROI[bilatROI %in% homo$ROI], homo$ROI), ]       # lobe order
table2 <- data.frame(Region = unname(sapply(homo$ROI, function(r) roiData$LongLabel[match(paste0(r, "-L"), roiData$FullLabel)])),
                     B = round(homo$scanTypeb, 3), SE = round(homo$scanTypese, 3),
                     p = fmt_p(homo$scanTypep), p_FDR = fmt_p(homo$scanTypep.Adj),
                     Direction = ifelse(homo$scanTypep.Adj >= .05, "n.s.",          # not significant after FDR
                                        ifelse(homo$scanTypeb > 0, "Preterm > Fetal", "Fetal > Preterm")),
                     row.names = NULL)
table2                                                              # B = Preterm - Fetal (Fisher z); FDR across all 3,003 edges
c(homotopic_edges = nrow(homo), significant = sum(homo$scanTypep.Adj < .05),
  preterm_greater = sum(homo$scanTypep.Adj < .05 & homo$scanTypeb > 0))

homoPlot <- rbind(data.frame(Connection = homo$ROI, Mean = homo$emmMean.Fetal,   SE = homo$emmSE.Fetal,   Group = "Fetal"),
                  data.frame(Connection = homo$ROI, Mean = homo$emmMean.Preterm, SE = homo$emmSE.Preterm, Group = "Preterm"))
homoPlot$Connection <- factor(homoPlot$Connection, levels = bilatROI)
# Bar settings for the homotopic chart below (section 4e)
bar_width   <- 0.8     # space used by the pair of bars for one ROI
bar_pad     <- 0.1     # gap between the two bars (dodge2 padding; 0.1 is the ggplot default)
errbar_frac <- 0.75    # error-bar cap width as a fraction of ONE bar's width
bar_each    <- bar_width / 2 * (1 - bar_pad)   # width of a single bar in plot units

## ---- 4e. Homotopic connectivity, horizontal version with lobe strip (Panel C) ---
# One plot: rows = bilateral ROIs in lobe order (Precentral at the top), Preterm above Fetal.
# The lobe-colored strip and the region names (single line, right-aligned) are drawn just left of the
# axis, outside the panel, so the plot margin on the left must leave room for the longest name.
homo_label_margin_in <- 1.75    # inches reserved on the left for region names (raise if a name is cut off)
homo_label_pt        <- 9       # region-name font size (points)
homo_strip_split     <- 0.4     # small white gap between lobe segments of the strip, in rows (0 = none)
homo_show_stars      <- TRUE    # TRUE: asterisk in a column at the far right for regions with an FDR-significant fetal-vs-preterm difference
homo_star_pt         <- 14      # asterisk size (points)
homo_star_vjust      <- 0.25    # nudges the asterisk up/down within its row (smaller = lower)
homo_rois  <- bilatROI[bilatROI %in% homo$ROI]
n_homo     <- length(homo_rois)
homo_lobe  <- roiData$Lobe[match(paste0(homo_rois, "-L"), roiData$FullLabel)]
row_pos    <- n_homo:1                                                   # row 1 (Precentral) sits at the top
homoH      <- homoPlot
homoH$ROI  <- factor(as.character(homoH$Connection), levels = rev(homo_rois))   # ggplot draws the first level at the bottom
homoH$Group <- factor(homoH$Group, levels = c("Fetal", "Preterm"))      # dodge order: the LAST level ends up on top, so Preterm is above Fetal

# x range: the data range rounded outward to 0.05 (tight), with a small extra pad on the left so the lowest
# error bar does not touch the axis. Set homo_x_limits <- c(lo, hi) to override.
homo_x_limits <- NULL
homo_x_pad    <- 0.02                                                  # extra room on the negative side (x units)
x_lo   <- if (is.null(homo_x_limits)) floor(min(homoH$Mean - homoH$SE) * 20) / 20 - homo_x_pad else homo_x_limits[1]
x_hi   <- if (is.null(homo_x_limits)) ceiling(max(homoH$Mean + homoH$SE) * 20) / 20 else homo_x_limits[2]
x_span <- x_hi - x_lo
strip_right <- x_lo - 0.008 * x_span                                   # strip sits just left of the axis line
strip_left  <- strip_right - 0.035 * x_span

# lobe strip as one segment per lobe, trimmed slightly where neighbors meet
homo_runs  <- rle(as.character(homo_lobe))
homo_end   <- cumsum(homo_runs$lengths); homo_start <- homo_end - homo_runs$lengths + 1
J          <- length(homo_runs$values)
strip_ymax <- row_pos[homo_start] + 0.5 - c(0, rep(homo_strip_split / 2, J - 1))
strip_ymin <- row_pos[homo_end]   - 0.5 + c(rep(homo_strip_split / 2, J - 1), 0)

# FDR-significant regions (fetal vs. preterm), in row order; the asterisk column sits just right of the panel
homo_sig <- homo$scanTypep.Adj[match(homo_rois, homo$ROI)] < .05
star_x   <- x_hi + 0.025 * x_span

p_homo_h <- ggplot(homoH, aes(x = Mean, y = ROI, group = Group)) +
  geom_col(aes(color = Group, fill = Group), width = bar_width, orientation = "y",
           position = position_dodge2(padding = bar_pad), linewidth = 0.6) +
  geom_errorbar(aes(xmin = Mean - SE, xmax = Mean + SE), orientation = "y", width = bar_each * errbar_frac,
                position = position_dodge(width = bar_width), linewidth = 0.4) +
  geom_vline(xintercept = 0, linewidth = 1) +
  # lobe strip (one segment per lobe) and single-line region names
  annotate("rect", xmin = strip_left, xmax = strip_right, ymin = strip_ymin, ymax = strip_ymax,
           fill = unname(lobe_pal[homo_runs$values])) +
  annotate("text", x = strip_left - 0.012 * x_span, y = row_pos, label = display_name(homo_rois),
           hjust = 1, size = homo_label_pt / ggplot2::.pt, color = "#000000") +
  scale_color_manual(values = c(Fetal = col_fetal, Preterm = col_preterm)) +
  scale_fill_manual(values = c(Fetal = paste0(col_fetal, "aa"), Preterm = paste0(col_preterm, "aa"))) +
  scale_x_continuous(breaks = seq(ceiling(x_lo * 4) / 4, floor(x_hi * 4) / 4, by = 0.25), labels = scales::number_format(accuracy = 0.01), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0.5)) +                                 # half a row of padding top and bottom
  coord_cartesian(xlim = c(x_lo, x_hi), expand = FALSE, clip = "off") +
  guides(color = guide_legend(reverse = TRUE), fill = guide_legend(reverse = TRUE)) +   # legend reads Preterm, Fetal (same order as the bars)
  labs(x = "Estimated marginal mean connectivity (Fisher z)", y = NULL) +
  theme_classic(base_size = 10) +
  theme(axis.line.y = element_line(color = "#000000", linewidth = 0.8),
        axis.line.x = element_line(color = "#000000", linewidth = 0.5),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(color = "#000000", size = 8),
        axis.title.x = element_text(color = "#000000", size = 10),
        panel.grid.major.x = element_line(color = "#d9d9d9", linewidth = 0.3),
        panel.grid.major.y = element_line(color = "#d9d9d9", linewidth = 0.3),   # one light-gray line per region, like the vertical ones
        legend.title = element_blank(), legend.text = element_text(size = 10),
        legend.background = element_rect(color = "#000000", fill = "white"),
        legend.key.size = unit(0.18, "in"),
        plot.margin = margin(6, if (homo_show_stars) 22 else 8, 6, homo_label_margin_in * 72, "pt")) +
  legend_inside(0.78, 0.04)
if (homo_show_stars) {       # * = FDR-corrected p < .05 (say so in the caption)
  p_homo_h <- p_homo_h +
    annotate("text", x = star_x, y = row_pos[homo_sig], label = "*", hjust = 0.5, vjust = homo_star_vjust,
             size = homo_star_pt / ggplot2::.pt, color = "#000000")
}
show_fig("Fig2_part_HomotopicBars", 6.2, 10, function() print(p_homo_h))

## ---- 4f. Figure 2 assembled: A = preterm (top), B = fetal (bottom), C = homotopic bars -------------
# Each piece (chord diagram + triangular matrix per group, and the homotopic chart) is drawn at its own design size
# (the sizes used for the separate files above), trimmed, and placed in one figure, so the proportions of the text
# stay as in the separate pieces. The chords are placed as pictures; the
# separate Fig2_chord_* files remain available as vectors for editing. Whole-figure size is set by fig2_box_h_in.
# Rounded frames around the panels (Figs. 2 and 3). Set fig2_boxes / fig3_boxes to FALSE for no frames.
fig2_boxes   <- TRUE
fig3_boxes   <- TRUE
fig3_pad_in  <- 0.15        # Fig. 3: space inside each frame above, below and beside the panels, inches
fig_box_col  <- "#2b2b2b"   # frame color (dark gray; "#808080" was the lighter version)
fig_box_lwd  <- 1.0         # frame line width
fig_box_r_in <- 0.12        # corner radius, inches
fig_box_inset_in <- 0.03    # frames are drawn this far inside the stated rectangle so the line is not clipped and neighbouring frames do not touch
round_box <- function(x, y, w, h) {                  # x, y = lower-left corner, inches; draws in the current viewport
  d <- fig_box_inset_in
  grid::grid.roundrect(x = grid::unit(x + d, "in"), y = grid::unit(y + d, "in"), width = grid::unit(w - 2 * d, "in"), height = grid::unit(h - 2 * d, "in"),
                       just = c("left", "bottom"), r = grid::unit(fig_box_r_in, "in"), gp = grid::gpar(col = fig_box_col, fill = NA, lwd = fig_box_lwd))
}
RUN_FIG2_ASSEMBLY <- TRUE
fig2_box_h_in   <- 4.2    # base height of a box, inches: sets the size of the matrices (and the first guess for the chords)
fig2_matrix_scale <- 1.1  # matrix height relative to the base height (content_h); larger = bigger matrix next to the chord (the chord is 1.15 on the same scale)
fig2_chord_scale <- 1.15  # chord diagram height relative to that (1 = same as before); the boxes, and panel C, grow with it. The figure gets wider and taller
fig2_title_h_in <- 0.5    # band at the top of each box for the title (Preterm / Fetal)
fig2_pad_in     <- 0.15   # space around the pieces
fig2_gap_in     <- 0.3    # space between the A and B rows and before C
fig2_title_pt   <- 18     # box titles (points)
fig2_letter_pt  <- 18     # panel letters (points)
fig2_dpi        <- 220    # resolution at which the pieces are drawn before they are placed (their design sizes are 6-10 in)
if (RUN_FIG2_ASSEMBLY && all(vapply(c("ragg", "png"), requireNamespace, NA, quietly = TRUE))) {
  render_piece <- function(draw, w_in, h_in) {
    f_tmp <- tempfile(fileext = ".png")
    ragg::agg_png(f_tmp, width = w_in, height = h_in, units = "in", res = fig2_dpi)
    tryCatch(draw(), finally = dev.off())
    im <- png::readPNG(f_tmp)[, , 1:3, drop = FALSE]; unlink(f_tmp); im
  }
  ch  <- trim_white(list(render_piece(function() draw_chord(chord_col_preterm), 9, 9), render_piece(function() draw_chord(chord_col_fetal), 9, 9)))
  mx  <- trim_white(list(render_piece(function() print(p_tri_preterm), 6, 6), render_piece(function() print(p_tri_fetal), 6, 6)))
  hm  <- trim_white(render_piece(function() print(p_homo_h), 6.2, 10))
  asp <- function(im) dim(im)[2] / dim(im)[1]

  content_h <- fig2_box_h_in - fig2_title_h_in - fig2_pad_in                   # height available to a matrix
  chord_h   <- content_h * fig2_chord_scale                                    # the chord can be larger than the matrix; the boxes grow to fit it
  box_h     <- fig2_title_h_in + fig2_pad_in + chord_h
  ch_w  <- chord_h * asp(ch[[1]]); mx_h <- content_h * fig2_matrix_scale; mx_w <- mx_h * asp(mx[[1]])
  left_w  <- fig2_pad_in + ch_w + 0.1 + mx_w + fig2_pad_in
  homo_h  <- 2 * box_h + fig2_gap_in - 2 * fig2_pad_in - 0.3
  homo_w  <- homo_h * asp(hm)
  right_w <- homo_w + 2 * fig2_pad_in
  W2 <- left_w + fig2_gap_in + right_w + 0.1; H2 <- 2 * box_h + fig2_gap_in + 0.1

  draw_fig2 <- function() {
    grid::grid.newpage()
    u <- function(v) grid::unit(v, "in")
    img <- function(im, x, y, w, h) grid::grid.raster(im, x = u(x), y = u(y), width = u(w), height = u(h), just = c("left", "bottom"), interpolate = TRUE)
    txt <- function(label, x, y, pt, just = "centre") grid::grid.text(label, x = u(x), y = u(y), just = just, gp = grid::gpar(fontsize = pt, fontface = "bold"))
    x0 <- 0.05
    for (k in 1:2) {                                              # k = 1: preterm (top), 2: fetal (bottom)
      y_bot <- 0.05 + (2 - k) * (box_h + fig2_gap_in)
      if (fig2_boxes) round_box(x0, y_bot, left_w, box_h)
      img(ch[[k]], x0 + fig2_pad_in, y_bot + fig2_pad_in * 0.5, ch_w, chord_h)
      img(mx[[k]], x0 + fig2_pad_in + ch_w + 0.1, y_bot + fig2_pad_in * 0.5 + (chord_h - mx_h) / 2, mx_w, mx_h)
      txt(c("Preterm", "Fetal")[k], x0 + left_w / 2, y_bot + box_h - fig2_title_h_in / 2, fig2_title_pt)
      txt(c("A", "B")[k], x0 + 0.1, y_bot + box_h - 0.08, fig2_letter_pt, just = c("left", "top"))
    }
    xc <- x0 + left_w + fig2_gap_in
    if (fig2_boxes) round_box(xc, 0.05, right_w, 2 * box_h + fig2_gap_in)
    img(hm, xc + fig2_pad_in, 0.05 + fig2_pad_in, homo_w, homo_h)
    txt("C", xc + 0.1, 0.05 + 2 * box_h + fig2_gap_in - 0.08, fig2_letter_pt, just = c("left", "top"))
  }
  show_fig("Fig2_FetalVsPreterm", W2, H2, draw_fig2, save = TRUE)   # always written to file: drawing this many large pictures in the RStudio plot pane times out
}

## ---- 4g. SI Fig. 2: mean connectivity matrices, one per group (Fetal, EPT, VPT, MPT) ----
# Full symmetric matrices of mean connectivity (r = tanh of the mean Fisher z per edge, averaged over the
# group's scans). ROIs are ordered by lobe (same order as the Fig. 2 matrices); the diagonal runs bottom-left
# to top-right and is drawn gray. Lobe boundaries are black lines; lobe color bars sit on top and on the right.
# Only the combined figure is saved: preterm groups on top (EPT, VPT, MPT), fetal on the bottom with the colorbar.
# Panel letters and titles are added by hand (Illustrator).
RUN_SI_MATRICES <- TRUE
MAT_SOURCE      <- "harmonized"   # "harmonized" (LongComBat, as in all analyses) or "unharmonized" (data before ComBat)
mat_lim         <- 1              # color scale runs from -mat_lim to +mat_lim
mat_ramp        <- c(low = "#2b2d9e", mid = "#ffffff", high = "#e4281f")
mat_diag_col    <- "#7f7f7f"      # diagonal cells
mat_label_pt_x  <- 2.3            # x-axis ROI tick-label size (points); 78 labels per axis, so about 2.3 pt is the largest that does not overlap
mat_label_pt_y  <- 2.3            # y-axis ROI tick-label size (points)
mat_tick_len_pt <- 1.5            # length of the tick marks, points
mat_tick_gap_pt <- 1.5            # space between the tick marks and the labels, points
mat_label_hemi  <- TRUE           # TRUE: "Precentral (L)"; FALSE: "Precentral" (rows alternate L, R, so the hemisphere is implied)
mat_line_w      <- 0.25           # lobe boundary lines
mat_letter_pt   <- 16             # panel letters A-D (points)
mat_letter_x_in <- 0.05           # letter position from the left edge of its cell, inches
mat_letter_y_in <- 0.02           # letter position from the top edge of its cell, inches
mat_title_pt    <- 11             # group titles above the matrices (points)
mat_titles      <- c(EPT = "EPT", VPT = "VPT", MPT = "MPT", Fetal = "Fetal")   # edit the wording here
mat_bar_len_in  <- 1.5            # colorbar length, inches
mat_bar_w_in    <- 0.12           # colorbar width, inches
mat_margin_top_pt <- 24           # top margin of every matrix, points (room for the group title above the lobe bar)
mat_margin_right_pt <- 14         # right margin of every matrix, points: must exceed the lobe bar (gap + thickness) or it is cut off at the cell edge; the fetal lobe names spill into the empty cell beside it
si2_cell_w_in   <- 3.0            # width of each of the 3 columns, inches (figure width = 3 x this)
si2_cell_h_in   <- 3.1            # height of each of the 2 rows, inches (figure height = 2 x this); lower it to remove vertical white space

if (RUN_SI_MATRICES) {
  fc_mat_src <- if (MAT_SOURCE == "unharmonized") combFCData else goodFCData
  mat_groups <- list(Fetal = scanInfo$longID[scanInfo$scanType == "Fetal"],
                     EPT   = scanInfo$longID[scanInfo$BirthCat %in% "EPT"],
                     VPT   = scanInfo$longID[scanInfo$BirthCat %in% "VPT"],
                     MPT   = scanInfo$longID[scanInfo$BirthCat %in% "MPT"])
  stopifnot(all(unlist(mat_groups) %in% rownames(fc_mat_src)))
  edge_i <- match(sapply(strsplit(allConn, "_", fixed = TRUE), `[`, 1), roi_levels)
  edge_j <- match(sapply(strsplit(allConn, "_", fixed = TRUE), `[`, 2), roi_levels)
  stopifnot(!anyNA(edge_i), !anyNA(edge_j))
  roi_tick <- if (mat_label_hemi) paste0(display_name(roi_levels), " (", substr(roi_levels, nchar(roi_levels), nchar(roi_levels)), ")")   # "PreCG.L" -> "Precentral (L)"
              else display_name(roi_levels)

  mat_long <- function(ids) {
    r <- tanh(colMeans(as.matrix(fc_mat_src[ids, allConn])))
    data.frame(x = c(edge_i, edge_j), y = c(edge_j, edge_i), r = c(r, r))
  }
  mat_plot <- function(ids, title = NULL, lobe_labels = FALSE, margin_pt = c(6, 6, 2, 2)) {
    d    <- mat_long(ids)
    bnd  <- lobe_end[-length(lobe_end)] + 0.5
    gap  <- tri_strip_gap; thick <- tri_strip_thick; K <- length(lobe_start); sp <- tri_strip_split / 2
    seg_lo <- lobe_start - 0.5 + c(0, rep(sp, K - 1)); seg_hi <- lobe_end + 0.5 - c(rep(sp, K - 1), 0)
    mid  <- (lobe_start + lobe_end) / 2
    edge <- n_roi + 0.5 + gap + thick                                               # outer edge of the lobe bars
    p <- ggplot() +
      geom_tile(data = d, aes(x, y, fill = r), width = 1, height = 1) +
      geom_tile(data = data.frame(x = seq_len(n_roi)), aes(x, x), fill = mat_diag_col, width = 1, height = 1) +
      scale_fill_gradient2(low = mat_ramp[["low"]], mid = mat_ramp[["mid"]], high = mat_ramp[["high"]],
                           midpoint = 0, limits = c(-mat_lim, mat_lim), guide = "none") +
      geom_hline(yintercept = bnd, linewidth = mat_line_w) + geom_vline(xintercept = bnd, linewidth = mat_line_w) +
      annotate("rect", xmin = seg_lo, xmax = seg_hi, ymin = n_roi + 0.5 + gap, ymax = edge, fill = lobe_pal[lobe_runs$values]) +   # top bar
      annotate("rect", xmin = n_roi + 0.5 + gap, xmax = edge, ymin = seg_lo, ymax = seg_hi, fill = lobe_pal[lobe_runs$values]) +   # right bar
      scale_x_continuous(breaks = seq_len(n_roi), labels = roi_tick, expand = c(0, 0)) +
      scale_y_continuous(breaks = seq_len(n_roi), labels = roi_tick, expand = c(0, 0)) +
      coord_cartesian(xlim = c(0.5, n_roi + 0.5), ylim = c(0.5, n_roi + 0.5), expand = FALSE, clip = "off") +
      labs(x = NULL, y = NULL) +
      theme_classic(base_size = 8) +
      theme(aspect.ratio = 1, plot.margin = margin(margin_pt[1], margin_pt[2], margin_pt[3], margin_pt[4], "pt"),
            axis.line = element_blank(), panel.border = element_rect(color = "#000000", fill = NA, linewidth = 0.5),
            axis.ticks = element_line(color = "#000000", linewidth = 0.2), axis.ticks.length = unit(mat_tick_len_pt, "pt"),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = mat_label_pt_x, color = "#000000", margin = margin(t = mat_tick_gap_pt)),
            axis.text.y = element_text(size = mat_label_pt_y, color = "#000000", margin = margin(r = mat_tick_gap_pt)),
            plot.background = element_blank())                                       # transparent: neighbouring panels never paint over the lobe bars
    if (!is.null(title)) p <- p +                                                  # group title centered above the top lobe bar
      annotate("text", x = n_roi / 2 + 0.5, y = edge + 2.5, label = title, vjust = 0, size = mat_title_pt / ggplot2::.pt, fontface = "bold")
    if (lobe_labels) p <- p +                                                      # lobe names with leader lines (right of the right bar)
      annotate("segment", x = edge + 0.5, xend = edge + 6, y = mid, yend = mid, linewidth = 0.3) +
      annotate("text", x = edge + 7, y = mid, label = lobe_runs$values, hjust = 0, size = 3.4)
    p
  }
  # Colorbar drawn in inches: the plot fills its viewport exactly, so the bar is exactly mat_bar_w_in wide
  # and mat_bar_len_in long; tick labels and "r" sit to its right.
  bar_pad_in <- 0.08                                                                # room above/below for the end labels
  bar_vp_w   <- mat_bar_w_in + 0.42; bar_vp_h <- mat_bar_len_in + 2 * bar_pad_in
  per_in     <- 2 * mat_lim / mat_bar_len_in                                        # data units per inch
  bar_breaks <- seq(-mat_lim, mat_lim, length.out = 5)
  mat_bar <- ggplot() +
    geom_tile(data = data.frame(y = seq(-mat_lim, mat_lim, length.out = 200)), aes(x = mat_bar_w_in / 2, y = y, fill = y),
              width = mat_bar_w_in, height = 2 * mat_lim / 199) +
    scale_fill_gradient2(low = mat_ramp[["low"]], mid = mat_ramp[["mid"]], high = mat_ramp[["high"]],
                         midpoint = 0, limits = c(-mat_lim, mat_lim), guide = "none") +
    annotate("rect", xmin = 0, xmax = mat_bar_w_in, ymin = -mat_lim, ymax = mat_lim, fill = NA, color = "#000000", linewidth = 0.5) +
    annotate("segment", x = mat_bar_w_in, xend = mat_bar_w_in + 0.04, y = bar_breaks, yend = bar_breaks, linewidth = 0.3) +
    annotate("text", x = mat_bar_w_in + 0.06, y = bar_breaks, label = format(bar_breaks, trim = TRUE), hjust = 0, size = 8 / ggplot2::.pt) +
    annotate("text", x = mat_bar_w_in + 0.30, y = 0, label = "r", hjust = 0, size = 9 / ggplot2::.pt) +
    coord_cartesian(xlim = c(0, bar_vp_w), ylim = c(-mat_lim - bar_pad_in * per_in, mat_lim + bar_pad_in * per_in), expand = FALSE, clip = "off") +
    theme_void() + theme(plot.margin = margin(0, 0, 0, 0))

  print(sapply(mat_groups, length))                                                 # scans per group
  # Combined figure: EPT, VPT, MPT on top; empty cell, fetal (with lobe names), colorbar on the bottom.
  # Every matrix gets the same margins so all four panels come out the same size.
  m_pt <- c(mat_margin_top_pt, mat_margin_right_pt, 2, 2)
  p_si_ept <- mat_plot(mat_groups$EPT, mat_titles[["EPT"]], margin_pt = m_pt)
  p_si_vpt <- mat_plot(mat_groups$VPT, mat_titles[["VPT"]], margin_pt = m_pt)
  p_si_mpt <- mat_plot(mat_groups$MPT, mat_titles[["MPT"]], margin_pt = m_pt)
  p_si_ft  <- mat_plot(mat_groups$Fetal, mat_titles[["Fetal"]], lobe_labels = TRUE, margin_pt = m_pt)      # lobe names overflow to the right (no clipping)
  # Fixed-size layout in inches (no automatic stretching): 3 columns x 2 rows of equal cells
  draw_si2 <- function() {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 3, widths = grid::unit(rep(si2_cell_w_in, 3), "in"),
                                                                  heights = grid::unit(rep(si2_cell_h_in, 2), "in"))))
    cell <- function(r, c) grid::viewport(layout.pos.row = r, layout.pos.col = c)
    print(p_si_ept, vp = cell(1, 1)); print(p_si_vpt, vp = cell(1, 2)); print(p_si_mpt, vp = cell(1, 3))
    print(p_si_ft,  vp = cell(2, 2))
    grid::pushViewport(cell(2, 3))                                                  # colorbar: exact size, centered in the bottom-right cell
    print(mat_bar, vp = grid::viewport(width = grid::unit(bar_vp_w, "in"), height = grid::unit(bar_vp_h, "in")))
    grid::popViewport()
    # Panel letters at the top-left corner of each matrix's cell: A = EPT, B = VPT, C = MPT, D = Fetal
    letter <- function(txt, r, c) { grid::pushViewport(cell(r, c))
      grid::grid.text(txt, x = grid::unit(mat_letter_x_in, "in"), y = grid::unit(1, "npc") - grid::unit(mat_letter_y_in, "in"),
                      just = c("left", "top"), gp = grid::gpar(fontsize = mat_letter_pt, fontface = "bold")); grid::popViewport() }
    letter("A", 1, 1); letter("B", 1, 2); letter("C", 1, 3); letter("D", 2, 2)
    grid::popViewport()
  }
  si_fig("mean_matrices", "MeanConnectivity", 3 * si2_cell_w_in, 2 * si2_cell_h_in, draw_si2)   # the only file written for this figure
}


###############################################################################
# 5. PRETERM ONLY: PMA x DEGREE OF PREMATURITY, ALL EDGES
#    Model (per edge):  FC ~ scanner + PMA * BirthCat,  random intercept/infant
#    Fit with lmer; the interaction is tested with a Type II Wald chi-square
#    (car::Anova, 2 df), then FDR across all 3,003 edges. The same test supplies
#    the Table 3 statistics, so selection and reporting use one method. In this
#    parameterization the PMA coefficient is the EPT slope and the interaction
#    coefficients are VPT/MPT slope minus EPT slope.
###############################################################################

# Generic scan over all edges: Type II Wald chi-square for one term of an lmer model
# (REML fit, random intercept per infant added to the fixed-effects formula `fml`).
# Used for the sensitivity analyses in section 9; the primary scan below is
# the same test computed together with the Table 3 statistics.
wald_scan <- function(d0, fml, term) {
  key <- paste(data_key, paste(deparse(fml), collapse = ""), term, table_key(d0), "lmer+car")
  cached("waldScan", key = key, wald_scan_run(d0, fml, term))
}
wald_scan_run <- function(d0, fml, term) {
  fml_re <- update(fml, . ~ . + (1 | shortID))
  out <- par_lapply(allConn, function(cn) {
    d <- d0; d$FCz <- as.numeric(goodFCData[d$longID, cn])
    tryCatch({
      m  <- suppressMessages(suppressWarnings(lmer(fml_re, data = d)))
      an <- car::Anova(m, type = "II", test = "Chisq")
      c(chi = an[term, "Chisq"], df = an[term, "Df"], p = an[term, "Pr(>Chisq)"])
    }, error = function(e) c(chi = NA, df = NA, p = NA))
  })
  res <- data.frame(Connection = allConn, do.call(rbind, out), row.names = NULL)
  res$pAdj <- p.adjust(res$p, "fdr")
  res
}

fml_main <- FCz ~ scanner + ScanPMA * BirthCat       # also used for Panel B and the slope tables

## ---- 5a. Statistics for every edge (primary scan + Table 3 columns) ---------
#   PMA_B / Scanner_B : fixed-effect coefficients (PMA = EPT slope; scanner = 3T - 1.5T)
#   Prem_chi2         : main effect of prematurity group, chi2(2)  (Type II)
#   Intx_chi2         : PMA x prematurity, chi2(2)                 (Type II)
tab3_row <- function(cn) {
  na_row <- data.frame(Connection = cn, PMA_B = NA, PMA_p = NA, Scanner_B = NA, Scanner_p = NA,
                       Prem_chi2 = NA, Prem_p = NA, Intx_chi2 = NA, Intx_p = NA)
  d <- pt; d$FC <- as.numeric(goodFCData[d$longID, cn])
  tryCatch({
    m  <- suppressMessages(suppressWarnings(lmer(FC ~ ScanPMA + scanner + BirthCat + ScanPMA:BirthCat + (1 | shortID), data = d)))
    cf <- summary(m)$coefficients
    an <- car::Anova(m, type = "II", test = "Chisq")
    data.frame(Connection = cn,
               PMA_B = cf["ScanPMA", "Estimate"],   PMA_p = cf["ScanPMA", "Pr(>|t|)"],
               Scanner_B = cf[scn_term, "Estimate"], Scanner_p = cf[scn_term, "Pr(>|t|)"],
               Prem_chi2 = an["BirthCat", "Chisq"], Prem_p = an["BirthCat", "Pr(>Chisq)"],
               Intx_chi2 = an["ScanPMA:BirthCat", "Chisq"], Intx_p = an["ScanPMA:BirthCat", "Pr(>Chisq)"])
  }, error = function(e) na_row)
}

# ~10-20 min on one core; parallel + cached after the first run
tab3_all <- cached("tab3_all", key = paste(data_key, "tab3_row"),
                   do.call(rbind, par_lapply(allConn, tab3_row, export = c("goodFCData", "scn_term", "pt"))))
for (v in c("PMA", "Scanner", "Prem", "Intx")) tab3_all[[paste0(v, "_pAdj")]] <- p.adjust(tab3_all[[paste0(v, "_p")]], "fdr")

## ---- 5b. Primary result: PMA x prematurity, FDR across all edges ------------
intxRes  <- data.frame(Connection = tab3_all$Connection, chi = tab3_all$Intx_chi2, df = 2,
                       p = tab3_all$Intx_p, pAdj = tab3_all$Intx_pAdj)
sigEdges <- intxRes$Connection[intxRes$pAdj < .05 & !is.na(intxRes$pAdj)]
c(failed_fits = sum(is.na(intxRes$p)), uncorrected_p05 = sum(intxRes$p < .05, na.rm = TRUE), fdr05 = length(sigEdges))

# ORDER OF THE SIGNIFICANT EDGES for Table 3, Fig. 3B (left to right, top to bottom),
# the BrainNet file and the slope tables. EDIT THIS LIST to reorder: keep one edge
# name per line, in the order you want. It must contain exactly the significant
# edges; if the data change and the set changes, the script stops and tells you
# which edges to add or remove.
edge_order <- c(
  # other frontal
  "ORBinf.L_OLF.R",
  # orbitofrontal / rectus
  "ORBsupb.L_ORBsupb.R",
  "ORBsupb.L_ORBmid.R",
  "ORBsupb.L_REC.L",
  "ORBsupb.L_REC.R",
  "ORBsupb.R_ORBmid.R",
  "ORBsupb.R_REC.L",
  "ORBmid.L_ORBmid.R",
  "ORBmid.L_REC.L",
  "ORBmid.L_REC.R",
  "ORBmid.R_REC.L",
  "ORBmid.R_REC.R",
  "REC.L_REC.R",
  # oddballs
  "PreCG.L_ORBmid.R",
  "SPG.L_HES.R",
  "SMG.L_SMG.R"
)
if (!setequal(edge_order, sigEdges)) {
  stop("edge_order does not match the significant edges.\n  Missing from edge_order: ",
       paste(setdiff(sigEdges, edge_order), collapse = ", "),
       "\n  In edge_order but not significant: ", paste(setdiff(edge_order, sigEdges), collapse = ", "))
}
stopifnot(!anyDuplicated(edge_order))
sigEdges <- edge_order
intxRes[match(sigEdges, intxRes$Connection), ]

## ---- 5c. Table 3 -------------------------------------------------------------
tab3 <- tab3_all[match(sigEdges, tab3_all$Connection), ]
# Footnote numbers: smallest FDR-adjusted p among the significant edges, for each effect
adj_cols <- c("PMA_pAdj", "Scanner_pAdj", "Prem_pAdj")
sapply(tab3[adj_cols], min, na.rm = TRUE)

tab3_display <- data.frame(
  `Region 1` = sapply(strsplit(tab3$Connection, "_", fixed = TRUE), function(p) full_name(p[1])),
  `Region 2` = sapply(strsplit(tab3$Connection, "_", fixed = TRUE), function(p) full_name(p[2])),
  `PMA B` = fmt3(tab3$PMA_B), `PMA p` = fmt_p(tab3$PMA_p),
  `Scanner B` = fmt3(tab3$Scanner_B), `Scanner p` = fmt_p(tab3$Scanner_p),
  `Prematurity chi2(2)` = sprintf("%.2f", tab3$Prem_chi2), `Prematurity p` = fmt_p(tab3$Prem_p),
  `Interaction chi2(2)` = sprintf("%.2f", tab3$Intx_chi2),
  `Interaction p` = paste0(fmt_p(tab3$Intx_p), ifelse(!is.na(tab3$Intx_pAdj) & tab3$Intx_pAdj < .05, "*", "")),   # * = FDR < .05
  check.names = FALSE)
tab3_display
# Footnotes: PMA B = slope in the EPT group (weeks); Scanner B > 0 means 3T > 1.5T;
# * = FDR-corrected p < .05 across 3,003 connections.


###############################################################################
# 6. FIG. 3 PANEL B: GROUP TRAJECTORIES FOR THE SIGNIFICANT EDGES (Fisher z)
#    Points: scans, adjusted to the average of the two scanners.
#    Lines: model-estimated trends over each group's observed ages (no
#    extrapolation). Layout: Panel B (the 13 consistent connections) and Panel C (the 3 oddballs), panel_ncol_split columns wide.
###############################################################################

# Appearance (tune here)
show_ci         <- FALSE   # TRUE: 95% bands
y_span          <- 2.4     # same vertical range (Fisher z) for every panel, centered on each panel's data
pt_size         <- 0.75    # point size
pt_alpha        <- 0.45    # point opacity
line_w          <- 0.9     # colored line width
halo_w          <- 0.2     # extra width of the white halo under each line
panel_ncol_split <- 5      # panels per row (Panel B and Panel C are both this many columns wide)
# Size: design the figure at the size it will be printed, so the text sizes below are the printed sizes.
panelB_w_in      <- 7.5    # export width (inches); a bit wider than tall leaves room for the brain image above
panelB_row_h_in  <- 1.3    # height of each row of panels (inches); the plots are square, so this only needs to fit the square plus its x-axis label
axis_text_size   <- 5.5    # tick labels (points)
axis_title_size  <- 6      # axis titles (points), both axes
legend_text_size <- 7      # legend (points)
panelB_margin_r_pt <- 2    # white space to the right of each plot (points); smaller = more room for the plot itself
panelB_labels    <- "short"   # y-axis names: "short" = display names with hemisphere, e.g. "Superior Orbital (L)"; "long" = atlas long names

panel_data <- function(cn, d0 = pt) {
  d <- d0; d$FCz <- as.numeric(goodFCData[d$longID, cn])
  m <- lme(fml_main, random = ~ 1 | shortID, data = d, method = "REML", control = ctrl)   # same model as the Wald test
  b <- nlme::fixef(m); V <- vcov(m)
  sc <- scn_term; s2 <- levels(d$scanner)[2]

  d$val <- d$FCz - b[sc] * ((d$scanner == s2) - 0.5)                  # scanner-adjusted points

  nd <- do.call(rbind, lapply(levels(d$BirthCat), function(g) {
    a <- d$ScanPMA[d$BirthCat == g]
    data.frame(ScanPMA = seq(min(a), max(a), length.out = 50),
               BirthCat = factor(g, levels = levels(d$BirthCat)),
               scanner  = factor(levels(d$scanner)[1], levels = levels(d$scanner)))
  }))
  X <- model.matrix(~ scanner + ScanPMA * BirthCat, nd)[, names(b)]
  X[, sc] <- 0.5                                                      # average of the two scanners
  fit <- drop(X %*% b); se <- sqrt(rowSums((X %*% V) * X))
  nd$val <- fit; nd$lo <- fit - 1.96 * se; nd$hi <- fit + 1.96 * se

  lab <- sub("_", " \u2013 ", gsub("\\.", "-", cn))
  d$label <- lab; nd$label <- lab
  list(pts = d[, c("ScanPMA", "BirthCat", "val", "label")], fit_lines = nd[, c("ScanPMA", "BirthCat", "val", "lo", "hi", "label")])
}

short_name <- function(x) paste0(display_name(x), " (", substr(x, nchar(x), nchar(x)), ")")   # "ORBsupb.L" -> "Superior Orbital (L)"
ylab_for <- function(cn) {                                            # two lines: ROI 1 -, ROI 2
  p  <- strsplit(cn, "_", fixed = TRUE)[[1]]
  nm <- if (panelB_labels == "short") short_name else full_name
  paste0(nm(p[1]), " \u2013\n", nm(p[2]))
}

# data: scan table to plot/fit (default: all preterm scans); xlim/xbreaks: PMA axis
make_panel_b <- function(edges, ncol = panel_ncol_split, data = pt, xlim = c(24, 49), xbreaks = seq(25, 45, 5)) {
  pd        <- lapply(edges, panel_data, d0 = data)
  lab_order <- sub("_", " \u2013 ", gsub("\\.", "-", edges))
  pts       <- do.call(rbind, lapply(pd, `[[`, "pts"));       pts$label       <- factor(pts$label, levels = lab_order)
  fit_lines <- do.call(rbind, lapply(pd, `[[`, "fit_lines")); fit_lines$label <- factor(fit_lines$label, levels = lab_order)

  edge_plot <- function(cn, lab) {
    dp   <- pts[pts$label == lab, ]; dp <- dp[order(dp$BirthCat), ]   # EPT first, MPT last (on top)
    dl   <- fit_lines[fit_lines$label == lab, ]
    ylim <- median(dp$val) + c(-1, 1) * y_span / 2
    p <- ggplot() +
      geom_point(data = dp, aes(ScanPMA, val, fill = BirthCat), shape = 21, size = pt_size,
                 color = "transparent", stroke = 0, alpha = pt_alpha)
    if (show_ci) p <- p + geom_ribbon(data = dl, aes(ScanPMA, ymin = lo, ymax = hi, fill = BirthCat), alpha = .15)
    p +
      geom_line(data = dl, aes(ScanPMA, val, group = BirthCat), color = "#ffffff", linewidth = line_w + halo_w) +
      geom_line(data = dl, aes(ScanPMA, val, color = BirthCat), linewidth = line_w) +
      scale_color_manual(values = cols_pt, labels = grp_labels) +
      scale_fill_manual(values = cols_pt, labels = grp_labels) +
      scale_x_continuous(breaks = xbreaks) +
      scale_y_continuous(breaks = scales::breaks_width(0.5)) +
      coord_cartesian(xlim = xlim, ylim = ylim) +
      labs(x = "PMA at scan (weeks)", y = ylab_for(cn), color = NULL, fill = NULL) +
      theme_classic(base_size = 8) +
      theme(aspect.ratio = 1,                                            # square plotting area
            axis.line = element_line(color = "#000000", linewidth = .5),
            axis.ticks = element_line(color = "#000000", linewidth = .4),
            axis.text = element_text(size = axis_text_size, color = "#000000"),
            axis.title.x = element_text(size = axis_title_size, color = "#000000"),
            axis.title.y = element_text(size = axis_title_size, color = "#000000"),
            plot.margin = margin(3, panelB_margin_r_pt, 1, 1),
            panel.grid.major.y = element_line(color = "#e6e6e6", linewidth = .3),
            legend.position = "none")                                        # the legend is drawn separately (pt_legend_plot)
  }
  plots <- lapply(seq_along(edges), function(j) edge_plot(edges[j], lab_order[j]))
  nrow_fig <- ceiling(length(plots) / ncol)
  list(fig = ggpubr::ggarrange(plotlist = plots, ncol = ncol, nrow = nrow_fig, align = "hv",
                               legend = "none"),
       nrow = nrow_fig)
}

# Legend-only plot for the empty cells beside Panel C
pt_legend_plot <- function() {
  d <- data.frame(g = factor(rep(names(cols_pt), each = 2), levels = names(cols_pt)), x = rep(1:2, 3), y = 1)
  p <- ggplot(d, aes(x, y, color = g, group = g)) + geom_line(linewidth = line_w) + geom_point(size = 2.2) +   # key = line with a dot (dots are scans, lines are model fits)
    scale_color_manual(values = cols_pt, labels = grp_labels, name = NULL) +
    guides(color = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(linewidth = 1.4, size = 2.2))) +
    theme_void() + theme(legend.position = "left", legend.text = element_text(size = legend_text_size, color = "#000000"),
                         legend.key.width = unit(0.45, "in"), legend.key.height = unit(0.24, "in"),
                         legend.background = element_rect(color = "#000000", fill = "#ffffff", linewidth = 0.3),   # box around the legend
                         legend.margin = margin(7, 10, 7, 10))
  ggpubr::as_ggplot(ggpubr::get_legend(p))
}
# draws the legend in the right-hand part of the Panel C area of the given size (columns beyond the three plots)
draw_pt_legend <- function(w_in, h_in, x0_in = 0, y0_in = 0, first_col = 3, ncols = panel_ncol_split) {
  print(pt_legend_plot(), vp = grid::viewport(x = grid::unit(x0_in + w_in * first_col / ncols, "in"), y = grid::unit(y0_in, "in"),
                                              width = grid::unit(w_in * (ncols - first_col) / ncols, "in"), height = grid::unit(h_in, "in"), just = c("left", "bottom")))
}

# Draws Panel B above Panel C (and the group legend beside C) in the current device. letters: labels for the two panels.
# pad_in > 0 (used with the rounded frames of Fig. 3): each of B and C gets pad_in of space above and below inside its own
# frame (B and C stay stacked, C below B) and the plots are shifted right by pad_in; the legend then sits inside B's frame,
# centered in the empty block to the right of the last row of plots. With pad_in = 0 the layout is as before (legend centered
# on the boundary between B's last row and C).
draw_panels_bc <- function(pB, pC, w_in, letters = c("B", "C"), x0_in = bn_panel_left_in, y0_in = 0, pad_in = 0) {
  b_h <- panelB_row_h_in * pB$nrow; c_h <- panelB_row_h_in * pC$nrow; bw <- w_in - x0_in; p <- pad_in
  yc <- y0_in + p; yb <- y0_in + c_h + 3 * p                       # bottoms of the plot areas of C and B
  vp  <- function(y, h) grid::viewport(x = grid::unit(x0_in + p, "in"), y = grid::unit(y, "in"), width = grid::unit(bw, "in"),
                                       height = grid::unit(h, "in"), just = c("left", "bottom"))
  print(pB$fig, vp = vp(yb, b_h)); print(pC$fig, vp = vp(yc, c_h))
  if (p > 0) draw_pt_legend(bw, panelB_row_h_in, x0_in = x0_in + p, y0_in = yb)
  else       draw_pt_legend(bw, 2 * panelB_row_h_in, x0_in = x0_in, y0_in = y0_in + c_h - panelB_row_h_in)
  lab <- function(txt, y_top) grid::grid.text(txt, x = grid::unit(bn_letter_in + 0.6 * p, "in"), y = grid::unit(y_top - 0.4 * p, "in") - grid::unit(0.02, "in"),
                                              just = c("left", "top"), gp = grid::gpar(fontsize = bn_letter_pt, fontface = "bold"))
  lab(letters[1], y0_in + c_h + b_h + 4 * p); lab(letters[2], y0_in + c_h + 2 * p)      # tops of the B and C frames
}

# Panel B = the connections where MPT differs from VPT and EPT (which look alike); Panel C = the three oddballs (listed
# below). Both are panel_ncol_split columns wide; the legend sits in the empty cells to the right of C.
panelC_edges     <- c("PreCG.L_ORBmid.R", "SPG.L_HES.R", "SMG.L_SMG.R")
stopifnot(all(panelC_edges %in% sigEdges))
edgesB <- setdiff(sigEdges, panelC_edges)                                        # keeps the order set in edge_order
panelB <- make_panel_b(edgesB)
panelC <- make_panel_b(panelC_edges)


###############################################################################
# 7. FIG. 3 PANEL A: EDGE FILE FOR BRAINNET VIEWER
#    Writes a symmetric binary adjacency matrix (rows/columns in node-file order)
#    with 1 for each significant PMA x prematurity edge.
###############################################################################

write_brainnet_edges <- function(edges, nodeFile, outFile) {
  nodes <- read.table(nodeFile)                                       # x y z color size label
  labs  <- as.character(nodes$V6)
  A <- matrix(0, length(labs), length(labs), dimnames = list(labs, labs))
  for (cn in edges) {
    p <- strsplit(cn, "_", fixed = TRUE)[[1]]
    r1 <- match(gsub(".", "-", p[1], fixed = TRUE), labs); r2 <- match(gsub(".", "-", p[2], fixed = TRUE), labs)
    if (anyNA(c(r1, r2))) stop("Edge ", cn, " not found in the node-file labels (expected e.g. 'PreCG-L').")
    A[r1, r2] <- 1; A[r2, r1] <- 1
  }
  write.table(A, outFile, sep = " ", col.names = FALSE, row.names = FALSE)
  invisible(A)
}
nodeFile <- file.path(figureDir, "bainNetPlotting", "new_AALNode_80_cortical.node")
if (file.exists(nodeFile)) {
  edgeMat <- write_brainnet_edges(sigEdges, nodeFile, file.path(figureDir, "bainNetPlotting", "Fig3_BrainNet_Edges.edge"))
  print(c(edges_written = sum(edgeMat) / 2, n_sigEdges = length(sigEdges)))   # should match
} else message("Node file not found; skipping BrainNet edge export.")

## ---- 7b. Fig. 3 Panel A in R: nodes and significant edges on the neonatal surface (rgl) ----
# Reads the neonatal surface exported from Surf Ice (preterm_mesh.obj) and the same node file used for
# BrainNet Viewer, then draws the nodes (colored by lobe) and the significant edges on a semi-transparent
# surface, then combines left lateral, right lateral and a dorsal view (with a box around the densest cluster) plus a
# magnified panel of that box into ONE figure with white space around each image for labeling by hand.
# Needs the rgl and png packages and an OpenGL display (a normal desktop session); the rgl windows must stay visible.
# Orientation: node and mesh coordinates are assumed to share one space with x = right, y = anterior,
# z = superior (as in BrainNet Viewer); they are rotated into rgl's axes (x, z, -y), a pure rotation.
# A warning is given if the nodes fall outside the mesh (different coordinate spaces).
RUN_BRAIN_RENDER <- TRUE
bn_mesh_name     <- "preterm_mesh.obj"
bn_nodes_shown   <- "edges_only"      # "edges_only": only ROIs touched by a significant edge; "all": all 78 ROIs
bn_node_scale    <- 0.022             # node radius as a fraction of the brain's width
bn_edge_scale    <- 0.006             # edge (cylinder) radius as a fraction of the brain's width
bn_edge_col      <- "#333333"
bn_surface_col   <- "#c8c8c8"
bn_surface_alpha <- 0.12              # surface opacity (lower = see the nodes through the brain more easily)
bn_fov           <- 10                # 0 = orthographic (nodes the same size at any depth, like a glass-brain figure); 30-60 = perspective
bn_zoom          <- c(right_lateral = 0.6, left_lateral = 0.6, dorsal = 0.7, ventral = 0.7, anterior = 0.6, posterior = 0.6)
                                      # per view; SMALLER = brain fills more of the image. The dorsal view is taller than the lateral ones, so it needs a larger value
bn_png_px        <- c(2000, 1600)     # snapshot width x height, pixels (= the rgl window size; keep it within your screen: a window larger than the screen is cut off)
# Views to save: theta (rotation about the vertical axis, degrees) and phi (elevation, degrees).
# If left and right come out swapped, swap the signs of theta in the two lateral views.
bn_views <- list(right_lateral = c(theta = 90,  phi = 0),  left_lateral = c(theta = -90, phi = 0),
                 dorsal        = c(theta = 0,   phi = 90), ventral      = c(theta = 0,   phi = -90),
                 anterior      = c(theta = 180, phi = 0),  posterior    = c(theta = 0,   phi = 0))
bn_panels         <- c("right_lateral", "anterior", "dorsal")   # views in the figure, left to right: any of right_lateral, left_lateral, anterior (coronal), posterior, dorsal, ventral
bn_inset          <- TRUE             # add a magnified panel of the most medial frontal connections (the dorsal view gets a rectangle around them)
bn_inset_view     <- "dorsal"         # view used for the main picture with the box and for the magnification
bn_inset_n_edges  <- 6                # the magnified region frames the n most medial frontal-frontal connections (rectangle around their nodes)
bn_inset_margin_mm <- 6               # extra space around those nodes inside the rectangle (mm)
bn_inset_shift_mm <- c(x = -3, y = 0)  # nudge the rectangle and the magnified view (mm): negative x = left, positive y = toward the front
bn_inset_scale    <- 0.85             # height of the magnified panel relative to the other images (it is wide, so a smaller height saves width)
bn_inset_zoom_trim <- 0.955          # fine correction of the automatic magnification: the scale is measured at the brain's widest depth, the nodes sit deeper, so the view came out ~4.5% too wide (smaller = more magnified)
bn_inset_zoom     <- NULL             # NULL = estimated from the box size; set a number (smaller = more magnified) to tune by hand
bn_box_col        <- "#000000"        # color of the box and the magnification frame
bn_panel_h_in     <- 2.3              # height of every image in the combined figure, inches
bn_pad_x_in       <- 0.2              # white space between images, inches (room for ROI labels; raise it if you need more)
bn_pad_y_in       <- 0.1              # white space above and below the images, inches
bn_key_h_in      <- 0.24             # height of the lobe color key under Panel A, inches
bn_key_pt        <- 7                # lobe key text size (points)
bn_gap_ab_in      <- 0.1              # space between Panel A and Panel B in the combined figure, inches
bn_letter_pt      <- 14               # panel letters A and B (points)
bn_panel_left_in  <- 0.22             # Panels B and C are indented by this much (inches) so the letters do not sit on the y-axis labels
bn_letter_in      <- 0.03             # distance of the letters from the top-left corner of each panel, inches

bn_dir <- c(file.path(figureDir, "brainNetPlotting"), file.path(figureDir, "bainNetPlotting"))   # either spelling of the folder
bn_dir <- bn_dir[dir.exists(bn_dir)][1]
if (RUN_BRAIN_RENDER && !all(c(requireNamespace("rgl", quietly = TRUE), requireNamespace("png", quietly = TRUE)))) {
  message("Packages rgl and png are needed for the brain render; skipping (install.packages(c('rgl', 'png'))).")
} else if (RUN_BRAIN_RENDER && (is.na(bn_dir) || !file.exists(file.path(bn_dir, bn_mesh_name)) || !file.exists(file.path(bn_dir, basename(nodeFile))))) {
  message("Mesh or node file not found in the BrainNet folder; skipping the brain render.")
} else if (RUN_BRAIN_RENDER) {
  to_rgl <- function(m) cbind(m[, 1], m[, 3], -m[, 2])                  # (x, y, z) -> rgl (x, z, -y)
  mesh <- rgl::readOBJ(file.path(bn_dir, bn_mesh_name))
  xyz  <- t(mesh$vb[1:3, ]); mesh$vb[1:3, ] <- t(to_rgl(xyz))
  mesh <- rgl::addNormals(mesh)
  bn_nodes <- read.table(file.path(bn_dir, basename(nodeFile)), stringsAsFactors = FALSE)
  names(bn_nodes) <- c("x", "y", "z", "color", "size", "label")
  bn_nodes$lobe <- roiData$Lobe[match(bn_nodes$label, roiData$FullLabel)]
  stopifnot(!anyNA(bn_nodes$lobe))
  nd_rng <- apply(bn_nodes[, c("x", "y", "z")], 2, range); ms_rng <- apply(xyz, 2, range)
  if (any(nd_rng[1, ] < ms_rng[1, ] - 0.1 * (ms_rng[2, ] - ms_rng[1, ]) | nd_rng[2, ] > ms_rng[2, ] + 0.1 * (ms_rng[2, ] - ms_rng[1, ])))
    warning("Some nodes lie well outside the mesh: the node file and the mesh may not share a coordinate space.")
  width_mm <- diff(range(xyz[, 1]))

  ends  <- do.call(rbind, strsplit(gsub(".", "-", sigEdges, fixed = TRUE), "_", fixed = TRUE))
  p1    <- to_rgl(as.matrix(bn_nodes[match(ends[, 1], bn_nodes$label), c("x", "y", "z")]))
  p2    <- to_rgl(as.matrix(bn_nodes[match(ends[, 2], bn_nodes$label), c("x", "y", "z")]))
  stopifnot(!anyNA(p1), !anyNA(p2))
  shown <- if (bn_nodes_shown == "all") seq_len(nrow(bn_nodes)) else which(bn_nodes$label %in% c(ends))
  node_xyz <- to_rgl(as.matrix(bn_nodes[shown, c("x", "y", "z")]))

  draw_brain <- function(box = FALSE) {
    rgl::bg3d("white")
    rgl::shade3d(mesh, color = bn_surface_col, alpha = bn_surface_alpha, specular = "black")
    for (i in seq_len(nrow(p1)))
      rgl::shade3d(rgl::cylinder3d(rbind(p1[i, ], p2[i, ]), radius = bn_edge_scale * width_mm, sides = 10, closed = -2),
                   color = bn_edge_col, specular = "black")
    rgl::spheres3d(node_xyz, radius = bn_node_scale * width_mm, color = lobe_pal[bn_nodes$lobe[shown]], specular = "black")
    if (box) for (i in 1:4)                                                          # square marking the magnified region (dorsal view)
      rgl::shade3d(rgl::cylinder3d(rbind(box_pts[i, ], box_pts[i %% 4 + 1, ]), radius = 0.25 * bn_edge_scale * width_mm, sides = 6, closed = -2),
                   color = bn_box_col, specular = "black")
  }

  # ---- magnified region: the bn_inset_n_edges most medial frontal-frontal connections ----
  # "Medial" = smallest mean distance of the two end nodes from the midline (x centre of the surface).
  nd_xyz <- as.matrix(bn_nodes[, c("x", "y", "z")]); rownames(nd_xyz) <- bn_nodes$label
  mid_x  <- mean(range(xyz[, 1]))
  is_fr  <- bn_nodes$lobe[match(ends[, 1], bn_nodes$label)] == "Frontal" & bn_nodes$lobe[match(ends[, 2], bn_nodes$label)] == "Frontal"
  medial <- (abs(nd_xyz[ends[, 1], 1] - mid_x) + abs(nd_xyz[ends[, 2], 1] - mid_x)) / 2
  pick   <- head(order(ifelse(is_fr, medial, Inf)), min(bn_inset_n_edges, sum(is_fr)))
  cl     <- unique(c(ends[pick, ]))                                                 # nodes of the picked connections
  box_lo <- apply(nd_xyz[cl, , drop = FALSE], 2, min); box_hi <- apply(nd_xyz[cl, , drop = FALSE], 2, max)
  hx <- (box_hi[1] - box_lo[1]) / 2 + bn_inset_margin_mm; hy <- (box_hi[2] - box_lo[2]) / 2 + bn_inset_margin_mm   # half-width / half-height of the rectangle (mm)
  clampc <- function(v, r, h) if (2 * h + 2 >= diff(r)) mean(r) else min(max(v, r[1] + h + 1), r[2] - h - 1)
  cl_ctr  <- (box_lo + box_hi) / 2 + c(bn_inset_shift_mm, 0)                         # centre of the rectangle, nudged by bn_inset_shift_mm (x = left/right, y = back/front)
  box_ctr <- c(clampc(cl_ctr[1], range(xyz[, 1]), hx), clampc(cl_ctr[2], range(xyz[, 2]), hy), cl_ctr[3])   # centre, kept inside the brain
  box_pts <- to_rgl(cbind(box_ctr[1] + c(-1, 1, 1, -1) * hx, box_ctr[2] + c(-1, -1, 1, 1) * hy, box_ctr[3]))   # at the depth of the nodes, so it lines up with the magnification
  message("Magnified region: connections ", paste(sigEdges[pick], collapse = ", "), "; box ", round(2 * hx, 1), " x ", round(2 * hy, 1), " mm")

  # ---- snapshots ----
  snap <- function(v, file, zoom, center_on = NULL, box = FALSE) {
    rgl::clear3d(); draw_brain(box)
    rgl::view3d(theta = bn_views[[v]][["theta"]], phi = bn_views[[v]][["phi"]], fov = bn_fov, zoom = zoom)
    if (!is.null(center_on)) {                                                      # shift the view so `center_on` (rgl coordinates) sits in the middle
      bb <- rgl::par3d("bbox"); c0 <- c(mean(bb[1:2]), mean(bb[3:4]), mean(bb[5:6]))
      um <- rgl::par3d("userMatrix"); um[1:3, 4] <- um[1:3, 4] - drop(um[1:3, 1:3] %*% (center_on - c0))
      rgl::par3d(userMatrix = um)
    }
    rgl::snapshot3d(file, webshot = FALSE)                                          # plain screen capture; no webshot2/Chrome needed
  }
  f_png <- function(nm) file.path(tempdir(), paste0("bn_", nm, ".png"))                 # intermediate snapshots are not kept
  rgl::open3d(windowRect = c(20, 20, 20 + bn_png_px[1], 20 + bn_png_px[2]))      # the saved PNG has the size of this window
  for (v in bn_panels) snap(v, f_png(v), bn_zoom[[v]], box = bn_inset && v == bn_inset_view)
  if (bn_inset) {
    # Calibrate the magnification from the dorsal picture itself: pixels per mm there (brain width in pixels / width in mm),
    # then choose the zoom at which the rectangle just fits the window. Zoom is inversely proportional to scale.
    im_d  <- png::readPNG(f_png(bn_inset_view))[, , 1:3, drop = FALSE]
    nw_d  <- im_d[, , 1] < .985 | im_d[, , 2] < .985 | im_d[, , 3] < .985
    px_mm <- diff(range(which(colSums(nw_d) > 0))) / diff(range(xyz[, 1]))
    px_mm_t <- min(dim(im_d)[2] / (2 * hx), dim(im_d)[1] / (2 * hy))                 # target pixels per mm in the magnified picture
    inset_zoom <- if (is.null(bn_inset_zoom)) bn_zoom[[bn_inset_view]] * px_mm / px_mm_t * bn_inset_zoom_trim else bn_inset_zoom
    message("Magnification zoom: ", round(inset_zoom, 3), " (set bn_inset_zoom to override)")
    snap(bn_inset_view, f_png("inset"), inset_zoom, center_on = to_rgl(matrix(box_ctr, 1))[1, ])
  }

  # ---- combine into one figure: the views in bn_panels, then the magnification (all in one row) ----
  # Each image is trimmed to the brain and surrounded by bn_pad_x_in and bn_pad_y_in of white space for labeling the ROIs by hand.
  # The magnification is calibrated from the dorsal image so that the box in the dorsal view and the magnified panel cover the same region.
  {
    rd <- function(f) png::readPNG(f)[, , 1:3, drop = FALSE]
    centre_rect <- function(im, w, h) { h <- min(h, dim(im)[1]); w <- min(w, dim(im)[2])                    # w x h pixels from the middle of the window
                                        r0 <- floor((dim(im)[1] - h) / 2); c0 <- floor((dim(im)[2] - w) / 2)
                                        im[r0 + seq_len(h), c0 + seq_len(w), , drop = FALSE] }
    ims   <- setNames(lapply(bn_panels, function(v) trim_white(rd(f_png(v)))), bn_panels)
    if (bn_inset) ims$inset <- centre_rect(rd(f_png("inset")), round(2 * hx * px_mm_t), round(2 * hy * px_mm_t))
    H     <- bn_panel_h_in; padx <- bn_pad_x_in; pady <- bn_pad_y_in
    hgt   <- setNames(rep(H, length(ims)), names(ims)); if (bn_inset) hgt[["inset"]] <- H * bn_inset_scale
    wid   <- sapply(names(ims), function(nm) hgt[[nm]] * dim(ims[[nm]])[2] / dim(ims[[nm]])[1])
    rows  <- list(c(bn_panels, if (bn_inset) "inset"))                                          # all images in one row
    row_w <- function(r) sum(wid[r]) + (length(r) + 1) * padx
    W     <- max(sapply(rows, row_w)); Hfig <- length(rows) * H + (length(rows) + 1) * pady
    # Draws the image strip (Panel A) in the current viewport; sc scales every size (1 = natural size in inches).
    draw_strip <- function(sc = 1, lwd = 1) {
      for (k in seq_along(rows)) {
        r <- rows[[k]]; x <- (W - row_w(r)) / 2 + padx; y0 <- Hfig - k * (H + pady)               # row k from the top
        for (nm in r) {
          yy <- y0 + (H - hgt[[nm]]) / 2                                                          # centred vertically in the row
          grid::grid.raster(ims[[nm]], x = grid::unit(x * sc, "in"), y = grid::unit(yy * sc, "in"), width = grid::unit(wid[[nm]] * sc, "in"),
                            height = grid::unit(hgt[[nm]] * sc, "in"), just = c("left", "bottom"), interpolate = TRUE)
          if (nm == "inset") grid::grid.rect(x = grid::unit(x * sc, "in"), y = grid::unit(yy * sc, "in"), width = grid::unit(wid[[nm]] * sc, "in"),
                                             height = grid::unit(hgt[[nm]] * sc, "in"), just = c("left", "bottom"), gp = grid::gpar(col = bn_box_col, fill = NA, lwd = lwd))
          x <- x + wid[[nm]] + padx
        }
      }
    }
    # Figure 3, Panels A + B: the image strip on top (scaled to the width of Panel B), the interaction plots below.
    # Panel B is drawn at its own design size (section 6), so its text sizes are unchanged.
    sc_a    <- panelB_w_in / W
    a_h_in  <- Hfig * sc_a
    b_h_in  <- panelB_row_h_in * panelB$nrow                                                                # Panel B
    c_h_in  <- panelB_row_h_in * panelC$nrow                                                                # Panel C (its legend sits beside the plots)
    p3      <- if (fig3_boxes) fig3_pad_in else 0                                                           # space inside each frame
    fig3_w  <- panelB_w_in + 2 * p3
    fig3_h  <- a_h_in + bn_key_h_in + b_h_in + c_h_in + 6 * p3 + (if (fig3_boxes) 0 else bn_gap_ab_in)
    draw_fig3 <- function() {
      grid::grid.newpage()
      grid::pushViewport(grid::viewport(x = grid::unit(p3, "in"), y = grid::unit(fig3_h - p3, "in"), width = grid::unit(panelB_w_in, "in"), height = grid::unit(a_h_in, "in"), just = c("left", "top")))
      draw_strip(sc_a, lwd = 1)                                                                   # Panel A (same drawing, scaled)
      grid::popViewport()
      # Lobe color key under Panel A: a dot per lobe, same colors as the nodes
      key_y  <- fig3_h - p3 - a_h_in - bn_key_h_in / 2
      key_w  <- c(0.62, 0.7, 0.62, 0.55)                                                          # width of each dot + label (inches); the last is the label's own width
      key_x0 <- p3 + (panelB_w_in - sum(key_w)) / 2
      grid::grid.rect(x = grid::unit(p3 + panelB_w_in / 2, "in"), y = grid::unit(key_y, "in"), width = grid::unit(sum(key_w) + 0.2, "in"),
                      height = grid::unit(bn_key_h_in - 0.04, "in"), gp = grid::gpar(col = "#000000", fill = NA, lwd = 0.85))   # box around the lobe key (same weight as the group legend box)
      for (k in seq_along(lobe_pal)) {
        xk <- key_x0 + sum(key_w[seq_len(k - 1)])
        grid::grid.circle(x = grid::unit(xk + 0.07, "in"), y = grid::unit(key_y, "in"), r = grid::unit(0.065, "in"),
                          gp = grid::gpar(fill = lobe_pal[[k]], col = "#00000066", lwd = 0.6))
        grid::grid.text(names(lobe_pal)[k], x = grid::unit(xk + 0.17, "in"), y = grid::unit(key_y, "in"), just = "left",
                        gp = grid::gpar(fontsize = bn_key_pt))
      }
      draw_panels_bc(panelB, panelC, panelB_w_in, pad_in = p3)                                    # Panels B and C below, indented so the letters fit
      grid::grid.text("A", x = grid::unit(bn_letter_in + 0.6 * p3, "in"), y = grid::unit(fig3_h - 0.4 * p3, "in") - grid::unit(0.02, "in"), just = c("left", "top"),
                      gp = grid::gpar(fontsize = bn_letter_pt, fontface = "bold"))
      if (fig3_boxes) {                                                                           # rounded frames: A (strip + lobe key), B, C
        round_box(0, c_h_in + b_h_in + 4 * p3, fig3_w, a_h_in + bn_key_h_in + 2 * p3)
        round_box(0, c_h_in + 2 * p3, fig3_w, b_h_in + 2 * p3)
        round_box(0, 0, fig3_w, c_h_in + 2 * p3)
      }
    }
    show_fig("Fig3_PmaByPrematurity", fig3_w, fig3_h, draw_fig3)
  }
  # The rgl window stays open: rotate it by hand, or use rgl::view3d(theta, phi) to try other angles.
}


###############################################################################
# 8. GROUP SLOPES AND PAIRWISE CONTRASTS FOR THE SIGNIFICANT EDGES
#    For each edge: PMA slope (Fisher z per week) in EPT, VPT and MPT, and the
#    pairwise slope differences, each with a 95% CI (conservative: the smallest
#    lme denominator df of the terms involved).
###############################################################################

slope_table <- function(d0, edges, fml, pre = "ScanPMA") {
  nV <- paste0(pre, ":BirthCatVPT"); nM <- paste0(pre, ":BirthCatMPT")
  do.call(rbind, lapply(edges, function(cn) {
    d <- d0; d$FCz <- as.numeric(goodFCData[d$longID, cn])
    m <- lme(fml, random = ~ 1 | shortID, data = d, method = "REML", control = ctrl)
    b <- nlme::fixef(m); V <- vcov(m)
    stopifnot(all(c(pre, nV, nM) %in% names(b)))
    w <- function(a) { x <- setNames(numeric(length(b)), names(b)); x[names(a)] <- a; x }
    L <- list(slope_EPT  = w(setNames(1, pre)),
              slope_VPT  = w(setNames(c(1, 1), c(pre, nV))),
              slope_MPT  = w(setNames(c(1, 1), c(pre, nM))),
              VPT_vs_EPT = w(setNames(1, nV)),
              MPT_vs_EPT = w(setNames(1, nM)),
              MPT_vs_VPT = w(setNames(c(1, -1), c(nM, nV))))
    dfd <- min(m$fixDF$X[c(pre, nV, nM)]); tq <- qt(.975, dfd)
    do.call(rbind, lapply(names(L), function(k) {
      x <- L[[k]]; est <- sum(x * b); se <- sqrt(drop(t(x) %*% V %*% x))
      data.frame(Connection = cn, term = k, est = est, lo = est - tq * se, hi = est + tq * se,
                 p = 2 * stats::pt(-abs(est / se), dfd))
    }))
  }))
}

# Wide display table: "slope [95% CI]" per group and contrast, plus contrast p-values
slope_wide <- function(tl) {
  fmt <- function(r) sprintf("%s [%s, %s]", fmt3(r$est), fmt3(r$lo), fmt3(r$hi))
  g   <- function(term) tl[tl$term == term, ]
  data.frame(Connection = g("slope_EPT")$Connection,
             EPT = fmt(g("slope_EPT")), VPT = fmt(g("slope_VPT")), MPT = fmt(g("slope_MPT")),
             `MPT-EPT` = fmt(g("MPT_vs_EPT")), `p (MPT-EPT)` = fmt_p(g("MPT_vs_EPT")$p),
             `MPT-VPT` = fmt(g("MPT_vs_VPT")), `p (MPT-VPT)` = fmt_p(g("MPT_vs_VPT")$p),
             `VPT-EPT` = fmt(g("VPT_vs_EPT")), `p (VPT-EPT)` = fmt_p(g("VPT_vs_EPT")$p),
             check.names = FALSE)
}
slopesMain <- slope_table(pt, sigEdges, fml_main)
slopesMain_wide <- slope_wide(slopesMain)
slopesMain_wide
# Uncorrected contrast p-values: these are follow-ups to the FDR-significant omnibus interaction.

# Supplemental table of the post hoc results: group slopes and pairwise slope contrasts on the 16 edges, with full region
# names, the Fig. 3 panel, uncorrected p and BH-FDR p (across the 16 edges, per contrast). Written under the "posthoc" SI number (see SI_ID).
si_posthoc <- local({
  tl <- slopesMain; g <- function(k) tl[tl$term == k, ]
  cell <- function(r) sprintf("%s [%s, %s]", fmt3(r$est), fmt3(r$lo), fmt3(r$hi))
  ctr  <- function(k, lab) { r <- g(k); setNames(data.frame(cell(r), fmt_p(r$p), fmt_p(p.adjust(r$p, "fdr")), stringsAsFactors = FALSE),
                                                 paste0(lab, c("", ", p", ", FDR p"))) }
  e <- g("slope_EPT")$Connection
  cbind(`Region 1` = sapply(strsplit(e, "_", fixed = TRUE), function(p) full_name(p[1])),
        `Region 2` = sapply(strsplit(e, "_", fixed = TRUE), function(p) full_name(p[2])),
        `Fig. 3 panel` = ifelse(e %in% panelC_edges, "C", "B"),
        `EPT slope` = cell(g("slope_EPT")), `VPT slope` = cell(g("slope_VPT")), `MPT slope` = cell(g("slope_MPT")),
        ctr("MPT_vs_EPT", "MPT - EPT"), ctr("MPT_vs_VPT", "MPT - VPT"), ctr("VPT_vs_EPT", "VPT - EPT"))
})
si_posthoc_counts <- sapply(c(`MPT - EPT` = "MPT_vs_EPT", `MPT - VPT` = "MPT_vs_VPT", `VPT - EPT` = "VPT_vs_EPT"), function(k) {
  p <- slopesMain$p[slopesMain$term == k]; c(`p < .05` = sum(p < .05), `FDR p < .05 (across 16)` = sum(p.adjust(p, "fdr") < .05), Edges = length(p)) })
si_posthoc_counts


###############################################################################
# 9. SENSITIVITY ANALYSES AND CONTROLS
#    Each repeats the PMA x prematurity scan under a different specification and
#    compares the FDR-significant set with the primary result (`sigEdges`).
#
#    SUPPLEMENTAL OUTPUTS: every control writes >= 1 table (SI<nn>_Table_*.csv, needs SAVE_TABLES) and >= 1 figure
#    (SI<nn>_Fig_*, needs SAVE_FIGS) to figureDir. Numbers come from the SI_ID registry (defined just before 9a), so
#    renumbering is one edit there; outputs of one control share a number and sort together:
#      SI01 usable scans (QA)       SI02 mean connectivity matrices   SI03 harmonization checks   SI04 post hoc slope contrasts
#      SI05 clinical cutoffs (9a)   SI06 continuous birth GA (9b)     SI07 scans > 32 wk (9c)    SI08 within/between (9d)
#      SI09 usable volumes (9e A-C) SI10 motion as covariate (9e D-F) SI11 age-matched (9f)      SI12 postnatal age (9g)
#      SI13 unharmonized (9h)       SI14 distance (9i)                SI15 motion follow-up (9j) SI16 equal volume (9k)
#      SI17 sex, fetal vs. preterm (11b)   SI18 sex, preterm (11c-d)   SI19 four-group model (12)   SI20 Table 1 group comparisons (13)
#
#    OVERVIEW (what each check tests and what it showed; details above each block)
#      9a  clinical cutoffs ........ 9 of the edges FDR-significant; 8 of the 16 hold (cutoff-dependent set)
#      9b  continuous birth GA ..... 0 significant (the effect is a step, not a linear gradient)
#      9c  scans after 32 wk PMA ... 6 of 16 hold; same slope pattern (less power)
#      9d  within/between infant ... omnibus not significant (underpowered); within-infant slopes match
#      9e  data amount (usable volumes; rerun pending) .. robust (direction kept for >= 97% of edges); motion: group and FD do not overlap, within-group adjustment leaves results unchanged; tSNR not planned (add on request at review)
#      9f  age-matched fetal/preterm  robust (1,174-1,354 of 1,404 edges stay significant; all same direction)
#      9g  postnatal age ........... robust (14 of 16 edges hold when adjusted); (C) PMA vs. postnatal clock: see 9g
#      9h  unharmonized data ....... fetal-vs-preterm fully robust; PMA x prematurity 9 of 16 hold
#      9i  distance between ROIs ... group differences unrelated to distance
#      9j  motion (mean FD) ........ r = -0.93 with group (VIF 7.2); in fetal scans motion is unrelated to FC (1 edge)
#      9k  strict equal-volume ..... fetal-vs-preterm unchanged at 80 volumes per scan (1,404 same direction, r = 0.99)
###############################################################################

compare_sets <- function(res_alt) {
  s <- res_alt$Connection[res_alt$pAdj < .05 & !is.na(res_alt$pAdj)]
  summary_vec <- c(raw_p05 = sum(res_alt$p < .05, na.rm = TRUE), fdr05 = length(s),
                   of_primary = sum(sigEdges %in% s), n_primary = length(sigEdges), newly_sig = sum(!s %in% sigEdges))
  list(summary = summary_vec, dropped = setdiff(sigEdges, s), entered = setdiff(s, sigEdges),
       primary_edges = res_alt[match(sigEdges, res_alt$Connection), c("Connection", "chi", "p", "pAdj")])
}

# Temporarily swaps the connectivity matrix (and the cache key) used by the scan functions (wald_scan, run_fetal, slope_table),
# so any analysis can be repeated on another version of the connectivity values (unharmonized in 9h, truncated in 9k).
with_fc <- function(fc, code, tag = "alt") {
  old <- list(goodFCData = goodFCData, data_key = data_key)
  on.exit({ goodFCData <<- old$goodFCData; data_key <<- old$data_key }, add = TRUE)
  goodFCData <<- fc
  data_key   <<- paste(old$data_key, tag, sum(as.matrix(fc)))
  force(code)
}

## ---- Supplemental (SI) outputs: helpers (the numbering registry SI_ID is in section 0) ----------------------------
si_theme <- function() theme_classic(base_size = 9) + theme(axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
                                                            strip.background = element_blank(), legend.position = "bottom")
facet_dims <- function(n, w = 2.5, h = 2.6) { nc <- min(n, 3); c(w = w * nc + 0.4, h = h * ceiling(n / nc) + 0.7) }
si_group_fill <- scale_fill_manual(values = c(Fetal = col_fetal, Preterm = col_preterm))
dist_plot <- function(d, var, xlab) {                       # distribution of a scan-level variable by group
  ggplot(d, aes(.data[[var]], fill = scanType)) + geom_density(alpha = 0.5, linewidth = 0.3, color = NA) +
    si_group_fill + labs(x = xlab, y = "Density", fill = NULL) + si_theme()
}

# PMA x prematurity reruns (wald_scan results). res_list and cmp_list are named lists (same names), one entry per control.
intx_summary <- function(cmp_list) {
  prim <- c(raw_p05 = sum(intxRes$p < .05, na.rm = TRUE), fdr05 = length(sigEdges), of_primary = length(sigEdges),
            n_primary = length(sigEdges), newly_sig = 0)
  rbind(data.frame(Analysis = "Primary", t(prim), check.names = FALSE),
        do.call(rbind, lapply(names(cmp_list), function(n) data.frame(Analysis = n, t(cmp_list[[n]]$summary), check.names = FALSE))))
}
intx_edges <- function(res_list) {                         # chi-square, p and FDR p of the 16 primary edges in each analysis
  all_res <- c(list(Primary = intxRes), res_list)
  do.call(rbind, lapply(names(all_res), function(n) {
    r <- all_res[[n]][match(sigEdges, all_res[[n]]$Connection), ]
    data.frame(Analysis = n, Connection = sigEdges, chi = round(r$chi, 2), p = signif(r$p, 3), pAdj = signif(r$pAdj, 3))
  }))
}
intx_scatter <- function(res_list) {                       # -log10 p of every edge: control vs. primary; the 16 primary edges are colored
  d <- do.call(rbind, lapply(names(res_list), function(n) {
    r <- res_list[[n]]; i <- match(r$Connection, intxRes$Connection)
    data.frame(Analysis = n, x = -log10(intxRes$p[i]), y = -log10(r$p), primary = r$Connection %in% sigEdges, fdr = !is.na(r$pAdj) & r$pAdj < .05)
  }))
  d$Analysis <- factor(d$Analysis, names(res_list)); d$fdr <- factor(d$fdr, c(TRUE, FALSE), c("FDR p < .05 in control", "Not significant in control"))
  ggplot(d, aes(x, y)) + geom_abline(linetype = "dashed", linewidth = 0.3) +
    geom_point(data = d[!d$primary, ], color = "#8a8a8a", alpha = 0.25, size = 0.5, stroke = 0) +
    geom_point(data = d[d$primary, ], aes(color = fdr), size = 1.5) +
    scale_color_manual(values = c("FDR p < .05 in control" = col_preterm, "Not significant in control" = "#e08a1e"), name = NULL) +
    facet_wrap(~Analysis, ncol = min(length(res_list), 3)) +
    labs(x = "-log10 p, primary analysis", y = "-log10 p, control") + si_theme()
}
# Writes the three standard outputs for a PMA x prematurity control: summary table, table of the 16 primary edges, p-value scatter.
si_intx <- function(key, name, res_list, cmp_list) {
  si_table(intx_summary(cmp_list), key, paste0(name, "_summary"))
  si_table(intx_edges(res_list),   key, paste0(name, "_primaryEdges"))
  dm <- facet_dims(length(res_list))
  si_fig(key, paste0(name, "_pvalues"), dm[["w"]], dm[["h"]], function() print(intx_scatter(res_list)))
}

si_table(si_posthoc, "posthoc", "pmaByPrematurity_posthocSlopeContrasts")
si_table(data.frame(Contrast = colnames(si_posthoc_counts), t(si_posthoc_counts), check.names = FALSE), "posthoc", "pmaByPrematurity_posthocCounts")

if (RUN_SENSITIVITY) {

  ## ---- 9a. Clinical cutoff convention (EPT <28 | VPT 28-<32 | MPT >=32) -----
  # WHAT IT TESTS: whether the PMA x prematurity result depends on where the group boundaries fall.
  #   The paper's rule puts infants born at exactly 28.0 wk in EPT and exactly 32.0 wk in VPT; the usual
  #   clinical convention puts them one group up (6 infants at 28.0 wk, 5 at 32.0 wk). Only group
  #   membership changes; harmonization is not rerun.
  # RESULT: 9 edges FDR-significant (primary: 16). 8 of the 16 hold (ORBsupb.L_ORBsupb.R, ORBsupb.L_REC.L,
  #   ORBsupb.L_REC.R, ORBsupb.R_REC.L, ORBmid.R_REC.R, REC.L_REC.R, SPG.L_HES.R, SMG.L_SMG.R); 1 new edge
  #   (SFGmed.R_REC.R). The 8 that drop have adjusted p of .055-.234 (two sit just above .05).
  # INTERPRETATION: the pattern is stable but the exact FDR set depends on the cutoff convention. Report
  #   the 8 edges that survive both rules as the cutoff-robust core, and state the boundary rule in Methods.
  # Changes group membership of the infants born at exactly 28.0 and 32.0 weeks.
  # Harmonization is NOT rerun (it used the paper's cutoffs).
  pt_clin <- pt; pt_clin$BirthCat <- birth_cat(pt_clin$BirthPMA, rule = "clinical")
  print(group_counts(pt_clin))
  res_clin <- wald_scan(pt_clin, fml_main, "ScanPMA:BirthCat")
  cmp_clin <- compare_sets(res_clin); print(cmp_clin)
  si_intx("clinical", "clinicalCutoffs", list(`Clinical cutoffs` = res_clin), list(`Clinical cutoffs` = cmp_clin))

  ## ---- 9b. Birth GA as a continuous moderator (no cutoffs; 1 df) ------------
  # WHAT IT TESTS: the same question with no cutoffs at all, using birth GA as a continuous moderator
  #   of the PMA slope (1 df).
  # RESULT: 0 edges FDR-significant (smallest adjusted p among the 16 edges: .238); 292 edges at raw p < .05.
  # INTERPRETATION: the effect is a step (MPT differs from EPT and VPT), not a linear gradient across
  #   birth GA, so a 1-df linear test is the wrong shape and has little power for it. This does not argue
  #   against the step; report it as a limitation (no linear dose-response with gestational age).
  pt_cont <- pt; pt_cont$GA <- pt_cont$BirthPMA - mean(pt_cont$BirthPMA)
  res_cont <- wald_scan(pt_cont, FCz ~ scanner + ScanPMA * GA, "ScanPMA:GA")
  cmp_cont <- compare_sets(res_cont); print(cmp_cont)
  si_intx("continuous_GA", "continuousBirthGA", list(`Birth GA, continuous (1 df)` = res_cont), list(`Birth GA, continuous (1 df)` = cmp_cont))

  ## ---- 9c. Scans after 32 weeks PMA only ------------------------------------
  # WHAT IT TESTS: whether the effect depends on the earliest scans (before 32 wk PMA), which exist only in
  #   EPT and VPT infants and could drive the group differences in slope.
  # RESULT: 301 scans from 176 infants remain (45 scans dropped). 6 of the 16 edges stay FDR-significant;
  #   several of the 10 that drop are near misses (adjusted p .052-.081). The group slopes (EPT and VPT near
  #   zero, MPT steep) and the MPT-EPT contrasts are essentially unchanged.
  # INTERPRETATION: the pattern does not depend on the earliest scans; the lower count of significant edges
  #   reflects fewer scans (less power).
  # Drops the earliest EPT/VPT scans (all MPT scans are already later).
  pt32 <- droplevels(pt[pt$ScanPMA > 32, ])
  print(c(scans = nrow(pt32), infants = nlevels(pt32$shortID), dropped = nrow(pt) - nrow(pt32)))
  print(table(pt32$BirthCat, pt32$scanner))                            # every group should still have both scanners
  res_gt32 <- wald_scan(pt32, fml_main, "ScanPMA:BirthCat")
  cmp_gt32 <- compare_sets(res_gt32); print(cmp_gt32)
  slopes_gt32_wide <- slope_wide(slope_table(pt32, sigEdges, fml_main)); print(slopes_gt32_wide)
  si_intx("gt32", "after32wk", list(`Scans after 32 wk PMA` = res_gt32), list(`Scans after 32 wk PMA` = cmp_gt32))
  si_table(slopes_gt32_wide, "gt32", "after32wk_slopes")
  # Supplemental figure: the Fig. 3 B/C trajectories refit on these scans only (same connections, order and layout).
  pB32 <- make_panel_b(edgesB,       data = pt32, xlim = c(32, 49), xbreaks = seq(35, 45, 5))
  pC32 <- make_panel_b(panelC_edges, data = pt32, xlim = c(32, 49), xbreaks = seq(35, 45, 5))
  si_fig("gt32", "After32wk_Trajectories", panelB_w_in, panelB_row_h_in * (pB32$nrow + pC32$nrow),
           function() { grid::grid.newpage(); draw_panels_bc(pB32, pC32, panelB_w_in, letters = c("A", "B")) },
           save = TRUE)                                                   # always written to figureDir (the plot pane shows blank space above it: the pane is taller than the figure)

  ## ---- 9d. Within- vs. between-infant PMA effects ---------------------------
  # WHAT IT TESTS: whether the MPT steepening reflects change WITHIN infants over time or differences
  #   BETWEEN infants scanned at different ages. PMA is split into an infant's mean (between) and each
  #   scan's deviation from that mean (within).
  # RESULT: infants with 2+ scans: EPT 47, VPT 53, MPT 23. Within-infant omnibus interaction: 0 edges
  #   FDR-significant (smallest adjusted p .095; 215 edges at raw p < .05). Between-infant: 0 (smallest .61).
  #   The within-infant slopes show the same pattern as the full sample (MPT about +0.04 to +0.10 per week,
  #   EPT and VPT near zero), e.g. ORBsupb.L_ORBsupb.R MPT 0.094 within vs. 0.091 overall.
  # INTERPRETATION: the MPT steepening appears within the same infants, not only between infants, but only
  #   23 MPT infants have repeat scans, so the omnibus test is underpowered. Report as supportive, not as
  #   an FDR-significant replication.
  # Splits PMA into an infant's mean (between) and the deviation of each scan
  # from that mean (within). The within term tests change inside the same infants.
  ptw <- pt
  ptw$n_scans  <- ave(ptw$ScanPMA, ptw$shortID, FUN = length)
  ptw$pma_mean <- ave(ptw$ScanPMA, ptw$shortID)
  ptw$pma_b    <- ptw$pma_mean - mean(ptw$ScanPMA)                    # between-infant
  ptw$pma_w    <- ptw$ScanPMA - ptw$pma_mean                          # within-infant
  # Infants with 2+ scans (information available for the within-infant terms)
  print(sapply(c("EPT", "VPT", "MPT"), function(g) length(unique(ptw$shortID[ptw$BirthCat == g & ptw$n_scans >= 2]))))
  fml_w <- FCz ~ scanner + pma_b + pma_w + BirthCat + pma_b:BirthCat + pma_w:BirthCat
  res_within  <- wald_scan(ptw, fml_w, "pma_w:BirthCat")
  res_between <- wald_scan(ptw, fml_w, "pma_b:BirthCat")
  cmp_within  <- compare_sets(res_within);  print(cmp_within)
  cmp_between <- compare_sets(res_between); print(cmp_between)
  slopes_within_wide  <- slope_wide(slope_table(ptw, sigEdges, fml_w, pre = "pma_w"))
  slopes_between_wide <- slope_wide(slope_table(ptw, sigEdges, fml_w, pre = "pma_b"))
  print(slopes_within_wide); print(slopes_between_wide)
  si_intx("within_between", "withinBetween", list(`Within-infant (pma_w x group)` = res_within, `Between-infant (pma_b x group)` = res_between),
          list(`Within-infant (pma_w x group)` = cmp_within, `Between-infant (pma_b x group)` = cmp_between))
  si_table(slopes_within_wide,  "within_between", "withinBetween_slopesWithin")
  si_table(slopes_between_wide, "within_between", "withinBetween_slopesBetween")
}

## ---- Shared helpers for 9e-9h (fetal vs. preterm reruns) -------------------------
# Summary of a fetal-vs-preterm rerun against the primary result (fetalRes / sigFetal)
cmp_fetal <- function(res_alt, n_scans = NA) {
  s_alt <- res_alt$scanTypep.Adj < .05 & !is.na(res_alt$scanTypep.Adj)
  c(scans = n_scans, fdr05 = sum(s_alt), preterm_greater = sum(s_alt & res_alt$scanTypeb > 0),
    fetal_greater = sum(s_alt & res_alt$scanTypeb < 0),
    primary_n = sum(sigFetal), still_sig_of_primary = sum(s_alt & sigFetal),
    same_direction_of_primary = sum(sigFetal & sign(res_alt$scanTypeb) == sign(fetalRes$scanTypeb), na.rm = TRUE),
    r_beta_with_primary = round(cor(fetalRes$scanTypeb, res_alt$scanTypeb, use = "complete.obs"), 3))
}
# Fetal-vs-preterm scan over all edges for a given scan table; `extra` adds covariates (parallel, cached)
run_fetal <- function(d, tag, extra = character(0)) {
  res <- cached(paste0("fetalRes_", tag), key = paste(data_key, "fetal", tag, paste(extra, collapse = "+"), table_key(d)),
                do.call(rbind, par_lapply(allConn, fit_fetal_edge, export = c("goodFCData", "ctrl", "scn_term", "scanInfo"),
                                          extra = extra, data = d)))
  res$scanPMAp.Adj  <- p.adjust(res$scanPMAp,  method = "fdr")
  res$scanTypep.Adj <- p.adjust(res$scanTypep, method = "fdr")
  res
}
# 1:1 nearest-neighbor matching of fetal to preterm scans on a numeric variable within a caliper
# (greedy, without replacement; repeated scans of one infant are kept)
match_scans <- function(d, var, caliper) {
  fet <- d[d$scanType == "Fetal", ]; pre <- d[d$scanType == "Preterm", ]
  free <- rep(TRUE, nrow(pre)); ki <- integer(0); kj <- integer(0)
  for (i in order(fet[[var]])) {
    dist <- abs(pre[[var]] - fet[[var]][i]); dist[!free] <- Inf
    j <- which.min(dist)
    if (is.finite(dist[j]) && dist[j] <= caliper) { ki <- c(ki, i); kj <- c(kj, j); free[j] <- FALSE }
  }
  rbind(fet[ki, ], pre[kj, ])
}
# One table row for a fetal-vs-preterm rerun. beta ratio = median over the primary significant edges of
# (beta in this analysis / primary beta); slope = regression slope of this analysis' betas on the primary betas.
ctrl_row <- function(res, label, n_scans = nrow(scanInfo)) {
  cf <- cmp_fetal(res, n_scans = n_scans)
  data.frame(Analysis = label, Scans = n_scans, `FDR-significant edges` = cf[["fdr05"]],
             `Preterm > fetal` = cf[["preterm_greater"]], `Fetal > preterm` = cf[["fetal_greater"]],
             `Primary edges still significant` = sprintf("%d (%.1f%%)", cf[["still_sig_of_primary"]], 100 * cf[["still_sig_of_primary"]] / cf[["primary_n"]]),
             `Primary edges, same direction` = sprintf("%d (%.1f%%)", cf[["same_direction_of_primary"]], 100 * cf[["same_direction_of_primary"]] / cf[["primary_n"]]),
             `Pearson r with primary betas` = round(cf[["r_beta_with_primary"]], 3),
             `Spearman rho with primary betas` = round(cor(fetalRes$scanTypeb, res$scanTypeb, method = "spearman", use = "complete.obs"), 3),
             `Slope on primary betas` = round(unname(coef(lm(res$scanTypeb ~ fetalRes$scanTypeb))[2]), 3),
             `Beta ratio (median, primary sig. edges)` = round(median(res$scanTypeb[sigFetal] / fetalRes$scanTypeb[sigFetal], na.rm = TRUE), 3),
             check.names = FALSE)
}
# Table of fetal-vs-preterm reruns against the primary result. res_list: named list of run_fetal() results; n_list: scans per analysis.
fetal_table <- function(res_list, n_list) {
  rbind(ctrl_row(fetalRes, "Primary", nrow(scanInfo)),
        do.call(rbind, lapply(names(res_list), function(n) ctrl_row(res_list[[n]], n, n_list[[n]]))))
}
# Group-difference betas of every edge: control vs. primary (the 1,404 primary significant edges are colored by direction)
beta_scatter <- function(res_list) {
  d <- do.call(rbind, lapply(names(res_list), function(n)
    data.frame(Analysis = n, x = fetalRes$scanTypeb, y = res_list[[n]]$scanTypeb, sig = sigFetal)))
  d$Analysis <- factor(d$Analysis, names(res_list))
  ann <- do.call(rbind, lapply(split(d, d$Analysis), function(g)
    data.frame(Analysis = g$Analysis[1], x = min(d$x, na.rm = TRUE), y = max(d$y, na.rm = TRUE),
               lab = sprintf("r = %.3f\nslope = %.2f", cor(g$x, g$y, use = "complete.obs"), unname(coef(lm(y ~ x, g))[2])))))
  ggplot(d, aes(x, y)) + geom_abline(linetype = "dashed", linewidth = 0.3) +
    geom_point(data = d[!d$sig, ], color = "#8a8a8a", alpha = 0.2, size = 0.4, stroke = 0) +
    geom_point(data = d[d$sig, ], color = col_preterm, alpha = 0.35, size = 0.5, stroke = 0) +
    geom_text(data = ann, aes(label = lab), hjust = 0, vjust = 1, size = 2.6) +
    facet_wrap(~Analysis, ncol = min(length(res_list), 3)) +
    labs(x = "Primary group difference (beta, preterm - fetal)", y = "Control group difference (beta)") + si_theme()
}
# Standard outputs for a fetal-vs-preterm control: results table and beta scatter (optionally with a second panel `extra`, e.g. a distribution).
si_fetal <- function(key, name, res_list, n_list, extra = NULL) {
  si_table(fetal_table(res_list, n_list), key, paste0(name, "_results"))
  dm <- facet_dims(length(res_list)); p <- beta_scatter(res_list)
  if (is.null(extra)) si_fig(key, paste0(name, "_betas"), dm[["w"]], dm[["h"]], function() print(p))
  else si_fig(key, paste0(name, "_betas"), dm[["w"]] + 3, dm[["h"]],
              function() print(ggpubr::ggarrange(p, extra, nrow = 1, widths = c(dm[["w"]], 3), labels = c("A", "B"))))
}

## ---- 9e. Fetal vs. preterm with data-quality covariates and similar data amounts -----
# WHAT IT TESTS: whether the fetal-vs-preterm group difference is explained by differences in the amount
#   or quality of data between fetal and preterm scans. Volumes = number of usable (non-censored) volumes per scan from the
#   ROI time series (section 3e; the QA table's finalVol differs from it in about 1 scan in 5).
#   (A) volumes as a covariate (raw); (B) volumes centered within group (so the covariate is not
#   collinear with group); (C) similar volumes: restrict to the volume range both groups share, and
#   1:1 match fetal to preterm scans on volumes; (D) motion (mean framewise displacement): as a raw and a
#   group-centered covariate, restricted to the range both groups share, and 1:1 matched on motion;
#   (E) volumes and motion together; (F) preterm only: the PMA x prematurity
#   test (Fig. 3) with motion as a covariate, since motion does overlap within the preterm group.
# STATUS: motion RUN (D-F). Volume controls A-C and E were first run with finalVol and then switched to the usable
#   volume count (scan_vol); the numbers marked [finalVol] below are from the first version and are to be replaced by
#   the rerun. tSNR is not planned; if a reviewer asks, add a tSNR column to si_fd (centered
#   within group) and list it in qc_more in part E.
#   Result summaries use the same columns as 9f:
#   fdr05, still_sig_of_primary (of the 1,404 primary edges), same_direction_of_primary, r_beta_with_primary.
# NOTE: this is a proxy for "equal amounts of data". The strict version (truncate every scan to the same number
#   of volumes and recompute FC) needs the ROI time series.
# DATA: finalVol differs by group: fetal 80-281 (median 117, 192 scans), preterm 120-296 (median 220, 346
#   scans); correlation with group r = 0.63. (An earlier minimum of 32 was a data-entry error, corrected to 132.)
#   102 of the 192 fetal scans (53%) have < 120 volumes while every preterm scan has >= 120, because the
#   criterion is 4 min of data in both groups: >= 120 volumes at TR = 2 s (preterm) and >= 80 volumes at
#   TR = 3 s (fetal). Fetal scans therefore have fewer volumes by design, and both groups have the same
#   duration of data; the volume controls below address the difference in volume count.
# RESULT [finalVol version, to be replaced by the usable-volume rerun] (each vs. the primary 1,404 significant edges;
#   'same direction' = of those 1,404):
#   A raw finalVol covariate ........ 1,113 significant; 1,029 still significant; 1,403 same direction; r = 0.962
#   B finalVol centered in group .... 1,426 significant; 1,387 still significant; 1,404 same direction; r = 0.999
#   C1 shared range (120-281) ....... 375 scans (90 fetal); 1,164 significant; 1,083 still; 1,403 same; r = 0.962
#   C2 1:1 matched on finalVol ...... 92 pairs (165.0 vs 167.7 volumes, SD 48.5 vs 46.5); 603 significant;
#                                      581 still; 1,373 same direction; r = 0.910
# MOTION (mean framewise displacement, mm; D, E): fetal 0.753 (SD 0.219, range 0.265-1.405), preterm 0.076 (SD 0.032,
#   range 0.023-0.225): the ranges do NOT overlap, and meanFD correlates -0.93 with group.
#   D1 raw meanFD covariate ......... 355 significant; 310 still; 1,257 same direction; r = 0.749 (group and motion
#                                      are nearly the same variable, so this model cannot separate them: do not interpret;
#                                      see 9j: VIF 7.2, group SEs inflated 2.32-fold)
#   D2 meanFD centered in group ..... 1,402 significant; 1,391 still; 1,404 same direction; r = 1.000
#   D3 shared range / 1:1 matching .. not possible (no overlap)
#   E  finalVol + meanFD centered ... 1,431 significant; 1,389 still; 1,404 same direction; r = 0.999
# INTERPRETATION (motion): between-group motion cannot be controlled: every fetal scan has more motion than every
#   preterm scan, so group and motion are confounded by design (and fetal and infant FD are probably not measured on
#   the same footing: the fetal head moves inside the mother and is registered differently). This is a limitation to
#   state. What can be shown is that motion variation WITHIN each group does not drive the result: centering within
#   group leaves every beta unchanged (r = 1.000, all 1,404 in the same direction), also with finalVol included.
# REPORT (motion): "Mean framewise displacement was higher in fetal than preterm scans, with no overlap between the
#   groups (0.75 vs 0.08 mm), so group differences cannot be separated from motion differences. Adjusting for
#   within-group variation in motion, alone or together with the number of volumes, left the results unchanged (1,402
#   and 1,431 significant edges; all 1,404 primary edges in the same direction)."
# MOTION, PRETERM ONLY (F): PMA x prematurity with meanFD as a covariate (all 346 preterm scans have a value).
#   Motion (mean FD, mm): EPT 0.082 (SD 0.034), VPT 0.074 (0.031), MPT 0.067 (0.024); meanFD falls with PMA
#   (r = -0.29 overall; -0.30 EPT, -0.27 VPT, -0.18 MPT), so motion could in principle reach the interaction.
#   RESULT: 15 edges FDR-significant; 14 of the 16 retained, 1 new (SFGmed.R_REC.R). The 2 that drop are near misses
#   (ORBinf.L_OLF.R adjusted p .051; ORBsupb.L_ORBmid.R .056). Group slopes are essentially unchanged: EPT and VPT near
#   zero, MPT 0.05-0.09 per week in the 13 consistent edges (MPT-EPT contrasts p < .001 in 12 of them), and the three
#   oddballs keep their different pattern.
#   INTERPRETATION: the PMA x prematurity result is not a motion effect. Motion is slightly higher in the youngest
#   and most premature infants, but adjusting for it leaves the interaction and the slopes intact.
#   REPORT: "The PMA x prematurity interaction was robust to adjustment for mean framewise displacement (14 of the 16
#   edges remained FDR-significant; the other two had adjusted p = .051 and .056), and the group slopes were
#   unchanged."
# INTERPRETATION: the group difference is not explained by the amount of data: direction is preserved for at
#   least 97.8% of edges and the betas correlate 0.91-0.999. Fewer significant edges in A (collinearity inflates
#   standard errors by about 30%) and C (90-92 fetal scans) reflect precision, not a changed effect. B equals the
#   primary result by construction (centering leaves the group coefficient unchanged); it only shows that
#   volumes matter little within groups. This covers data AMOUNT; for data QUALITY see MOTION above (tSNR is not planned; see STATUS).
# REPORT: "Results were unchanged when including finalVol as a covariate (1,113 significant edges; 1,403 of
#   1,404 primary edges in the same direction) or restricting to scans with >= 120 volumes in both groups
#   (1,164 significant edges; 1,403 of 1,404 in the same direction)."
if (RUN_QC_COVARIATE_CONTROL) {
  # Volume measure: usable volumes per scan from the time series (scan_vol, section 3e). Copies of the scan table carry it,
  # so scanInfo (and every cache keyed on it) stays as it is.
  si_vol <- scanInfo; si_vol$vol <- scan_vol
  stopifnot(!anyNA(si_vol$vol))                                        # every scan needs a value
  print(tapply(si_vol$vol, si_vol$scanType, function(x) c(n = length(x), min = min(x), median = median(x), max = max(x))))
  print(c(r_vol_with_group = cor(si_vol$vol, as.numeric(si_vol$scanType == "Preterm"))))                  # collinearity with group
  si_vol$vol_c <- si_vol$vol - ave(si_vol$vol, si_vol$scanType)                                          # group-mean-centered

  ## A. raw volumes as a covariate
  res_qc_raw <- run_fetal(si_vol, "usableVol", extra = "vol")
  print(cmp_fetal(res_qc_raw, n_scans = nrow(si_vol)))
  ## B. volumes centered within group
  res_qc_cen <- run_fetal(si_vol, "usableVolC", extra = "vol_c")
  print(cmp_fetal(res_qc_cen, n_scans = nrow(si_vol)))

  ## C. similar volumes between fetal and preterm scans
  vol_rng <- lapply(split(si_vol$vol, si_vol$scanType), range)
  vol_lo <- max(sapply(vol_rng, `[`, 1)); vol_hi <- min(sapply(vol_rng, `[`, 2))
  si_vol_range <- droplevels(si_vol[si_vol$vol >= vol_lo & si_vol$vol <= vol_hi, ])
  print(c(vol_lo = vol_lo, vol_hi = vol_hi)); print(table(si_vol_range$scanType))
  res_vol_range <- run_fetal(si_vol_range, "usableVolRange")
  print(cmp_fetal(res_vol_range, n_scans = nrow(si_vol_range)))
  si_vol_match <- droplevels(match_scans(si_vol, "vol", caliper = 0.1 * sd(si_vol$vol)))
  print(c(matched_pairs = sum(si_vol_match$scanType == "Fetal")))
  print(tapply(si_vol_match$vol, si_vol_match$scanType, function(x) c(mean = mean(x), sd = sd(x))))   # balance on volumes
  res_vol_match <- run_fetal(si_vol_match, "usableVolMatched")
  print(cmp_fetal(res_vol_match, n_scans = nrow(si_vol_match)))
  si_table(do.call(rbind, lapply(split(si_vol, si_vol$scanType), function(g)
    data.frame(Group = g$scanType[1], Scans = nrow(g), `Volumes, mean (SD)` = sprintf("%.1f (%.1f)", mean(g$vol), sd(g$vol)),
               `Volumes, median` = median(g$vol), `Volumes, range` = paste(range(g$vol), collapse = "-"), check.names = FALSE))),
    "volumes", "usableVolumes_descriptives")
  si_fetal("volumes", "usableVolumes",
           list(`Volumes as covariate` = res_qc_raw, `Volumes centered within group` = res_qc_cen,
                `Shared volume range` = res_vol_range, `Matched on volumes` = res_vol_match),
           list(`Volumes as covariate` = nrow(si_vol), `Volumes centered within group` = nrow(si_vol),
                `Shared volume range` = nrow(si_vol_range), `Matched on volumes` = nrow(si_vol_match)),
           extra = dist_plot(si_vol, "vol", "Usable volumes per scan"))

  ## D. motion (mean framewise displacement, from Table 1): as a covariate, and similar motion between groups.
  # Uses a copy of the scan table (si_fd) restricted to scans with a motion value, so scanInfo and the caches stay as they are.
  si_fd <- scanInfo; si_fd$meanFD <- scan_fd; si_fd$vol_c <- si_vol$vol_c[match(si_fd$longID, si_vol$longID)]
  si_fd <- droplevels(si_fd[!is.na(si_fd$meanFD), ])
  if (nlevels(si_fd$scanType) < 2) stop("meanFD is missing for every scan of one group, so a fetal-vs-preterm comparison is not possible.")
  si_fd$meanFD_c <- si_fd$meanFD - ave(si_fd$meanFD, si_fd$scanType)                                    # group-mean-centered
  print(tapply(si_fd$meanFD, si_fd$scanType, function(x) c(n = length(x), mean = mean(x), sd = sd(x), min = min(x), median = median(x), max = max(x))))
  print(c(r_meanFD_with_group = cor(si_fd$meanFD, as.numeric(si_fd$scanType == "Preterm"))))             # collinearity with group
  # D1. raw meanFD as a covariate; D2. meanFD centered within group
  res_fd_raw <- run_fetal(si_fd, "meanFD", extra = "meanFD")
  print(cmp_fetal(res_fd_raw, n_scans = nrow(si_fd)))
  res_fd_cen <- run_fetal(si_fd, "meanFDC", extra = "meanFD_c")
  print(cmp_fetal(res_fd_cen, n_scans = nrow(si_fd)))
  # D3. similar motion: the meanFD range both groups share, and 1:1 matching on meanFD. Skipped when the groups do not overlap
  # (the case here: every preterm scan is below every fetal scan, so there is no shared range to restrict to or match within).
  fd_rng <- lapply(split(si_fd$meanFD, si_fd$scanType), range)
  fd_lo <- max(sapply(fd_rng, `[`, 1)); fd_hi <- min(sapply(fd_rng, `[`, 2))
  if (fd_lo >= fd_hi) {
    message("meanFD ranges of the two groups do not overlap (fetal ", round(fd_rng$Fetal[1], 3), "-", round(fd_rng$Fetal[2], 3),
            "; preterm ", round(fd_rng$Preterm[1], 3), "-", round(fd_rng$Preterm[2], 3), "): similar-motion restriction and matching skipped.")
  } else {
    si_fd_range <- droplevels(si_fd[si_fd$meanFD >= fd_lo & si_fd$meanFD <= fd_hi, ])
    print(c(fd_lo = fd_lo, fd_hi = fd_hi)); print(table(si_fd_range$scanType))
    res_fd_range <- run_fetal(si_fd_range, "fdRange")
    print(cmp_fetal(res_fd_range, n_scans = nrow(si_fd_range)))
    si_fd_match <- droplevels(match_scans(si_fd, "meanFD", caliper = 0.1 * sd(si_fd$meanFD)))
    print(c(matched_pairs = sum(si_fd_match$scanType == "Fetal")))
    print(tapply(si_fd_match$meanFD, si_fd_match$scanType, function(x) c(mean = mean(x), sd = sd(x))))      # balance on motion
    res_fd_match <- run_fetal(si_fd_match, "fdMatched")
    print(cmp_fetal(res_fd_match, n_scans = nrow(si_fd_match)))
  }

  ## E. every QC covariate together: volumes (centered) + meanFD (centered) + anything in qc_more
  ## (tSNR: add e.g. si_fd$tSNR_c <- ... above and set qc_more <- c("tSNR_c"))
  qc_more <- "meanFD_c"
  stopifnot(all(qc_more %in% names(si_fd)), !anyNA(si_fd[qc_more]))
  res_qc_all <- run_fetal(si_fd, "allQC", extra = c("vol_c", qc_more))
  print(cmp_fetal(res_qc_all, n_scans = nrow(si_fd)))

  ## F. preterm only: does the PMA x prematurity result (Fig. 3) survive adjusting for motion? Motion overlaps within the
  ## preterm group, so unlike D this adjustment is meaningful. Same test as the primary scan plus meanFD as a covariate.
  pt_fd <- pt; pt_fd$meanFD <- scan_fd[match(pt$longID, scanInfo$longID)]
  pt_fd <- droplevels(pt_fd[!is.na(pt_fd$meanFD), ])
  print(c(scans = nrow(pt_fd), dropped_no_motion = nrow(pt) - nrow(pt_fd)))
  print(round(do.call(rbind, tapply(pt_fd$meanFD, pt_fd$BirthCat, function(x) c(mean = mean(x), sd = sd(x)))), 3))          # motion by birth group
  print(c(r_meanFD_PMA = round(cor(pt_fd$meanFD, pt_fd$ScanPMA), 3),                                                       # motion vs. age at scan
          sapply(split(pt_fd, pt_fd$BirthCat), function(d) c(round(cor(d$meanFD, d$ScanPMA), 3)))))
  fml_fd <- FCz ~ scanner + ScanPMA * BirthCat + meanFD
  res_pt_fd <- wald_scan(pt_fd, fml_fd, "ScanPMA:BirthCat")
  cmp_pt_fd <- compare_sets(res_pt_fd); print(cmp_pt_fd)
  slopes_pt_fd_wide <- slope_wide(slope_table(pt_fd, sigEdges, fml_fd)); print(slopes_pt_fd_wide)

  ## SI outputs for D-F: fetal vs. preterm with motion (table + beta scatter with the meanFD distribution), and the preterm-only PMA x group check
  fd_res <- list(`meanFD as covariate` = res_fd_raw, `meanFD centered within group` = res_fd_cen)
  fd_n   <- list(`meanFD as covariate` = nrow(si_fd), `meanFD centered within group` = nrow(si_fd))
  if (fd_lo < fd_hi) {
    fd_res[["Shared meanFD range"]] <- res_fd_range; fd_n[["Shared meanFD range"]] <- nrow(si_fd_range)
    fd_res[["Matched on meanFD"]]   <- res_fd_match; fd_n[["Matched on meanFD"]]   <- nrow(si_fd_match)
  }
  fd_res[["Volumes + meanFD (centered)"]] <- res_qc_all; fd_n[["Volumes + meanFD (centered)"]] <- nrow(si_fd)
  si_table(do.call(rbind, lapply(split(si_fd, si_fd$scanType), function(g)
    data.frame(Group = g$scanType[1], Scans = nrow(g), `meanFD (mm), mean (SD)` = sprintf("%.3f (%.3f)", mean(g$meanFD), sd(g$meanFD)),
               `meanFD, median` = round(median(g$meanFD), 3), `meanFD, range` = paste(round(range(g$meanFD), 3), collapse = "-"), check.names = FALSE))),
    "motion_covariate", "meanFD_descriptives")
  si_fetal("motion_covariate", "meanFD_fetalVsPreterm", fd_res, fd_n, extra = dist_plot(si_fd, "meanFD", "Mean framewise displacement (mm)"))
  si_intx("motion_covariate", "meanFD_pmaByPrematurity", list(`Preterm: PMA x group, adjusted for meanFD` = res_pt_fd),
          list(`Preterm: PMA x group, adjusted for meanFD` = cmp_pt_fd))
  si_table(slopes_pt_fd_wide, "motion_covariate", "meanFD_pmaByPrematurity_slopes")
}

## ---- 9f-9h. Controls that use only data already in this script -------------------
if (RUN_EXTRA_CONTROLS) {

  ## ---- 9f. Age-matched fetal vs. preterm (control #3) ---------------------------
  # WHAT IT TESTS: whether the fetal-vs-preterm difference is an artifact of the groups covering different
  #   ages (fetal 18.7-40.0 wk, preterm 25.7-48.1 wk). The primary model adjusts for PMA, but this removes
  #   the age mismatch by design.
  # RESULT: (A) shared PMA range, 25.7-40.0 wk (171 fetal / 273 preterm scans): 1,416 FDR-significant edges;
  #   1,354 of the primary 1,404 remain; all 1,404 keep the same direction; correlation of betas 0.997.
  #   (B) 1:1 matching within 1 wk (166 pairs, PMA 34.2 vs 34.3 wk, SD about 3): 1,196 significant; 1,174 of
  #   the primary 1,404 remain; all same direction; r = 0.988.
  # INTERPRETATION: the group difference is not explained by age differences between the groups. Fewer
  #   significant edges in (B) reflect fewer scans, not a change in the effect.
  # The primary model adjusts for PMA, but the groups do not cover the same ages. Two checks:
  #  A. restrict both groups to the PMA range they share;
  #  B. 1:1 nearest-neighbor matching on PMA within a caliper (greedy, without replacement).
  # Repeated scans of one infant are kept (the model has a random intercept per infant).
  pma_rng <- lapply(split(scanInfo$ScanPMA, scanInfo$scanType), range)
  am_lo <- max(sapply(pma_rng, `[`, 1)); am_hi <- min(sapply(pma_rng, `[`, 2))
  pma_tab <- rbind(Fetal = pma_rng$Fetal, Preterm = pma_rng$Preterm, Overlap = c(am_lo, am_hi))
  colnames(pma_tab) <- c("min", "max"); print(pma_tab)

  si_range <- droplevels(scanInfo[scanInfo$ScanPMA >= am_lo & scanInfo$ScanPMA <= am_hi, ])
  print(table(si_range$scanType))
  res_am_range <- run_fetal(si_range, "ageRange")
  print(cmp_fetal(res_am_range, n_scans = nrow(si_range)))

  si_match <- droplevels(match_scans(scanInfo, "ScanPMA", caliper = 1))
  print(c(matched_pairs = sum(si_match$scanType == "Fetal"), fetal_infants = nlevels(droplevels(si_match$shortID[si_match$scanType == "Fetal"])),
          preterm_infants = nlevels(droplevels(si_match$shortID[si_match$scanType == "Preterm"]))))
  print(tapply(si_match$ScanPMA, si_match$scanType, function(x) c(mean = mean(x), sd = sd(x))))   # balance on PMA after matching
  res_am_match <- run_fetal(si_match, "ageMatched")
  print(cmp_fetal(res_am_match, n_scans = nrow(si_match)))
  si_table(data.frame(Set = rownames(pma_tab), pma_tab, check.names = FALSE, row.names = NULL), "age_matched", "ageMatched_pmaRanges")
  si_table(do.call(rbind, lapply(split(si_match, si_match$scanType), function(g)                      # balance on PMA after matching
    data.frame(Group = g$scanType[1], Scans = nrow(g), Infants = nlevels(droplevels(g$shortID)),
               `PMA, mean (SD)` = sprintf("%.2f (%.2f)", mean(g$ScanPMA), sd(g$ScanPMA)), check.names = FALSE))),
    "age_matched", "ageMatched_balance")
  si_fetal("age_matched", "ageMatched",
           list(`Overlapping PMA range` = res_am_range, `Matched within 1 wk PMA` = res_am_match),
           list(`Overlapping PMA range` = nrow(si_range), `Matched within 1 wk PMA` = nrow(si_match)),
           extra = dist_plot(scanInfo, "ScanPMA", "PMA at scan (weeks)"))

  ## ---- 9g. Postnatal-age check for PMA x prematurity (control #5) -----------------
  # WHAT IT TESTS: whether the PMA x prematurity result could be a postnatal-age effect. PMA = birth GA +
  #   postnatal age, so groups differ in postnatal age at the same PMA (median postnatal age: EPT 10.3,
  #   VPT 6.3, MPT 3.9 wk; ranges 0.6-21 wk). Postnatal age is 0.92-0.96 correlated with PMA within group,
  #   so this is a demanding test with wider confidence intervals.
  # RESULT: (A) adding postnatal age as a covariate: 14 edges FDR-significant, 14 of the 16 retained, 0 new;
  #   the 2 that drop have adjusted p = .0536. Slopes are unchanged: EPT/VPT near zero, MPT +0.045 to +0.093.
  #   (B) postnatal age as the time axis instead of PMA (a different hypothesis): 12 significant, 10 of the 16
  #   retained, 2 new (SFGmed.R_REC.R, PCUN.L_HES.L); the dropped edges are near misses (adjusted p .11-.19).
  # INTERPRETATION: the effect is not a postnatal-age artifact. (A) is the direct test; (B) answers a
  #   different question and is expected to differ somewhat.
  # PMA = gestational age at birth + postnatal age, so differences in postnatal age could
  # masquerade as differences in PMA slopes. (A) Add postnatal age as a covariate and retest the
  # interaction; (B) use postnatal age as the time axis instead of PMA.
  pt_pna <- pt; pt_pna$postnatal <- pt_pna$ScanPMA - pt_pna$BirthPMA          # weeks since birth
  print(tapply(pt_pna$postnatal, pt_pna$BirthCat, function(x) round(c(min = min(x), median = median(x), max = max(x)), 1)))
  print(round(tapply(seq_len(nrow(pt_pna)), pt_pna$BirthCat, function(i) cor(pt_pna$ScanPMA[i], pt_pna$postnatal[i])), 2))   # collinearity within group
  fml_pna <- FCz ~ scanner + ScanPMA * BirthCat + postnatal
  res_pna <- wald_scan(pt_pna, fml_pna, "ScanPMA:BirthCat")                    # (A) PMA x group, adjusted for postnatal age
  cmp_pna <- compare_sets(res_pna); print(cmp_pna)
  slopes_pna_wide <- slope_wide(slope_table(pt_pna, sigEdges, fml_pna)); print(slopes_pna_wide)
  res_pna_axis <- wald_scan(pt_pna, FCz ~ scanner + postnatal * BirthCat, "postnatal:BirthCat")   # (B) postnatal age as the time axis
  cmp_pna_axis <- compare_sets(res_pna_axis); print(cmp_pna_axis)

  ## (C) Which clock organizes connectivity better: PMA or time since birth? ----------------------
  # WHAT IT TESTS: the "not how premature, but how many weeks postnatal" hypothesis. PMA = birth GA +
  #   postnatal age, so any model containing PMA and postnatal age already contains birth GA; the two
  #   ideas can only be separated by comparing single-clock models and by testing what each adds.
  #   Six models are fitted to every edge with ML (needed to compare fixed effects), preterm scans only:
  #     pma        FCz ~ scanner + PMA                         (maturation clock = PMA)
  #     post       FCz ~ scanner + postnatal                   (time since birth only)
  #     both       FCz ~ scanner + PMA + postnatal             (= PMA + linear birth GA)
  #     pma_cat    FCz ~ scanner + PMA * BirthCat              (the primary model)
  #     post_cat   FCz ~ scanner + postnatal * BirthCat        (same, on the postnatal clock)
  #     pma_x_post FCz ~ scanner + PMA * postnatal             (continuous interaction)
  #   Output: per model, the median BIC difference from the PMA-only model (negative = better) and the number
  #   of edges where that model has the lowest BIC; then likelihood-ratio tests (1 df, FDR over edges) of
  #   whether postnatal age adds to PMA and whether PMA adds to postnatal age. "Postnatal adds to PMA" is the
  #   same as "birth GA matters at a given PMA". Run for all edges and for the 16 primary edges.
  # RESULT: lowest BIC over all 3,003 edges: PMA only 2,046 (68%), postnatal only 936 (31%), the other four
  #   models 21 combined (BIC charges about 5.9 per added parameter). Median BIC vs. PMA only: postnatal only
  #   +0.8 (near ties), PMA + postnatal +5.4, PMA x BirthCat +19.6, postnatal x BirthCat +19.9, PMA x
  #   postnatal +10.1. Likelihood-ratio tests (FDR over edges): postnatal age adds to PMA in 0 of 3,003 edges
  #   (0 of the 16 primary edges); PMA adds to postnatal age in 366 edges (5 of the 16). Primary 16 edges,
  #   lowest BIC: PMA x BirthCat 8, PMA only 6, postnatal only 1, postnatal x BirthCat 1; postnatal x BirthCat
  #   is 2.6 BIC worse than PMA x BirthCat (median +2.7 vs +0.1 relative to PMA only).
  # INTERPRETATION: connectivity tracks maturational age (PMA), not time since birth: PMA fits as well or
  #   better than postnatal age, and postnatal age (equivalently, linear birth GA at a given PMA) adds
  #   nothing. BIC is conservative, so the categorical model only ties PMA-only on median BIC even for the 16
  #   edges whose interaction is FDR-significant. PMA and postnatal age are 0.92-0.96 correlated, so this is
  #   the absence of a detectable postnatal-age effect, not proof that none exists.
  # REPORT: "Connectivity was better described by PMA than by postnatal age (PMA-only model had the lowest
  #   BIC in 68% of edges vs. 31% for postnatal age only). Postnatal age added nothing beyond PMA (0 of 3,003
  #   edges, FDR-corrected likelihood-ratio tests), whereas PMA added beyond postnatal age in 366 edges."
  ic_models <- list(pma        = FCz ~ scanner + ScanPMA,
                    post       = FCz ~ scanner + postnatal,
                    both       = FCz ~ scanner + ScanPMA + postnatal,
                    pma_cat    = FCz ~ scanner + ScanPMA * BirthCat,
                    post_cat   = FCz ~ scanner + postnatal * BirthCat,
                    pma_x_post = FCz ~ scanner + ScanPMA * postnatal)
  ic_edge <- function(cn) {
    d <- pt_pna; d$FCz <- as.numeric(goodFCData[d$longID, cn])
    sapply(ic_models, function(f) {
      m <- tryCatch(lme(f, random = ~ 1 | shortID, data = d, method = "ML", control = ctrl), error = function(e) NULL)
      if (is.null(m)) c(AIC = NA, BIC = NA, logLik = NA) else c(AIC = AIC(m), BIC = BIC(m), logLik = as.numeric(logLik(m)))
    })
  }
  res_ic <- cached("icClocks", key = paste(data_key, "icClocks", table_key(pt_pna)),
                   par_lapply(allConn, ic_edge, export = c("goodFCData", "ctrl", "pt_pna", "ic_models")))
  bic <- t(sapply(res_ic, function(m) m["BIC", ])); ll <- t(sapply(res_ic, function(m) m["logLik", ]))
  ic_summary <- function(idx) {
    b <- bic[idx, , drop = FALSE]; b <- b[complete.cases(b), , drop = FALSE]
    data.frame(model = colnames(b), median_dBIC_vs_pma = round(apply(b - b[, "pma"], 2, median), 1),
               edges_lowest_BIC = as.integer(table(factor(colnames(b)[max.col(-b, ties.method = "first")], levels = colnames(b)))),
               row.names = NULL)
  }
  lrt_fdr <- function(big, small) p.adjust(pchisq(pmax(2 * (ll[, big] - ll[, small]), 0), 1, lower.tail = FALSE), "fdr")
  sig_idx <- match(sigEdges, allConn)
  cat("\nAll edges:\n");   print(ic_summary(seq_along(allConn)))
  cat("\nPrimary 16 edges:\n"); print(ic_summary(sig_idx))
  p_post_adds_pma <- lrt_fdr("both", "pma"); p_pma_adds_post <- lrt_fdr("both", "post")
  print(c(all_edges_postnatal_adds_to_PMA   = sum(p_post_adds_pma < .05, na.rm = TRUE),
          all_edges_PMA_adds_to_postnatal   = sum(p_pma_adds_post < .05, na.rm = TRUE),
          primary16_postnatal_adds_to_PMA   = sum(p_post_adds_pma[sig_idx] < .05, na.rm = TRUE),
          primary16_PMA_adds_to_postnatal   = sum(p_pma_adds_post[sig_idx] < .05, na.rm = TRUE)))

  ## SI outputs: (A, B) PMA x group with postnatal age as a covariate / as the time axis, slopes, and the model comparison of the two clocks
  pna_res <- list(`PMA x group, adjusted for postnatal age` = res_pna, `Postnatal age x group (replaces PMA)` = res_pna_axis)
  pna_cmp <- list(`PMA x group, adjusted for postnatal age` = cmp_pna, `Postnatal age x group (replaces PMA)` = cmp_pna_axis)
  si_intx("postnatal", "postnatalAge", pna_res, pna_cmp)
  si_table(slopes_pna_wide, "postnatal", "postnatalAge_slopes")
  ic_tab <- rbind(data.frame(Edges = "All 3,003 edges", ic_summary(seq_along(allConn))), data.frame(Edges = "Primary 16 edges", ic_summary(sig_idx)))
  si_table(ic_tab, "postnatal", "clocks_BICsummary")
  si_table(data.frame(Comparison = c("Postnatal age adds to PMA", "PMA adds to postnatal age"), check.names = FALSE,
                      `All edges, FDR p < .05` = c(sum(p_post_adds_pma < .05, na.rm = TRUE), sum(p_pma_adds_post < .05, na.rm = TRUE)),
                      `Primary 16 edges, FDR p < .05` = c(sum(p_post_adds_pma[sig_idx] < .05, na.rm = TRUE), sum(p_pma_adds_post[sig_idx] < .05, na.rm = TRUE))),
           "postnatal", "clocks_LRT")
  p_ic <- ggplot(ic_tab, aes(model, median_dBIC_vs_pma, fill = Edges)) + geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    scale_fill_manual(values = c("All 3,003 edges" = "#8a8a8a", "Primary 16 edges" = col_preterm)) +
    labs(x = NULL, y = "Median BIC change vs. the PMA-only model\n(lower = better fit)", fill = NULL) + si_theme() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  p_pna_scatter <- intx_scatter(pna_res)
  si_fig("postnatal", "clocks_and_pvalues", 7.5, 3.4,
         function() print(ggpubr::ggarrange(p_pna_scatter, p_ic, nrow = 1, widths = c(2, 1), labels = c("A", "B"))))

  ## ---- 9h. Unharmonized data with scanner as a covariate (control #11) ------------
  # WHAT IT TESTS: whether the findings depend on LongComBat. All models include scanner as a fixed effect,
  #   so this repeats them on the connectivity values BEFORE harmonization. It also addresses the limitation
  #   that LongComBat does not protect interaction terms.
  # RESULT: fetal vs. preterm: 1,439 edges significant, all 1,404 primary edges remain, all the same
  #   direction, correlation of betas 1.000 (scanner is already a covariate and fetal scans are all 1.5T, so
  #   harmonization changes nothing here). PMA x prematurity: 9 of the 16 edges significant, none new; the 7
  #   that drop have adjusted p of .058-.147 (most just above .05) and chi-square of 14-17.
  # INTERPRETATION: the group difference is fully robust to harmonization; the interaction pattern persists
  #   without it, but harmonization improves power (16 vs. 9 edges). Report both numbers.
  # Repeats the two main analyses on the connectivity values BEFORE LongComBat. Every model already
  # includes scanner as a fixed effect, so this shows whether the findings depend on harmonization.
  rawFC <- combFCData[scanInfo$longID, colnames(goodFCData)]
  stopifnot(!anyNA(rawFC), identical(dim(rawFC), dim(goodFCData)))
  fetalRes_raw <- with_fc(rawFC, run_fetal(scanInfo, "unharmonized"), "unharmonized")
  print(cmp_fetal(fetalRes_raw, n_scans = nrow(scanInfo)))
  res_raw <- with_fc(rawFC, wald_scan(pt, fml_main, "ScanPMA:BirthCat"), "unharmonized")
  cmp_raw <- compare_sets(res_raw); print(cmp_raw)

  ## SI outputs: tables for both analyses; one figure with the fetal-vs-preterm betas (A) and the PMA x group p-values (B)
  si_table(fetal_table(list(`Unharmonized` = fetalRes_raw), list(`Unharmonized` = nrow(scanInfo))), "unharmonized", "fetalVsPreterm_results")
  si_table(intx_summary(list(`Unharmonized` = cmp_raw)), "unharmonized", "pmaByPrematurity_summary")
  si_table(intx_edges(list(`Unharmonized` = res_raw)),   "unharmonized", "pmaByPrematurity_primaryEdges")
  si_fig("unharmonized", "betas_and_pvalues", 5.6, 3.2,
         function() print(ggpubr::ggarrange(beta_scatter(list(`Fetal vs. preterm` = fetalRes_raw)),
                                            intx_scatter(list(`PMA x group` = res_raw)), nrow = 1, labels = c("A", "B"))))
}

## ---- 9i. Do the fetal-vs-preterm differences depend on distance between ROIs? ----
# Connectivity usually falls off with distance, so a group difference could simply be larger for
# short or long connections. Distances are Euclidean distances between ROI centroids from the
# age-specific coordinate templates in labelsAndCoords (GA 23-37 and 40 weeks).
#   (1) Edge level: relate each edge's fetal-vs-preterm difference (beta) and each group's adjusted
#       mean FC to distance (Spearman; Mantel-type permutation p-values that shuffle ROI labels,
#       because edges sharing an ROI are not independent).
#   (2) Sensitivity: the same correlation with distances from every template age.
#   (3) Scan level: for each scan, the Spearman correlation between its FC and edge distance (using
#       the template nearest its own PMA) = how strongly its FC falls off with distance. The two
#       groups are compared with a mixed model adjusted for PMA and scanner.
# ---- CONCLUSION FOR THE PAPER (from the run on the final data) ---------------------------------
# Part (1) is the control that answers the question: fetal-vs-preterm differences are NOT associated
# with ROI distance (GA 36 template: rho = -0.002, permutation p = .94; within hemisphere rho = 0.025,
# p = .40; between hemispheres rho = -0.004, p = .91), and |rho| <= 0.03 with the template of every age
# (part 2). Connectivity itself falls with distance in both groups (fetal rho = -0.33, preterm
# rho = -0.44, permutation p < .001; 0.0002 is the smallest p 5,000 permutations can give), as expected
# for an early-development brain with predominantly local organization.
#
# Suggested wording.
#   Methods: "To test whether fetal-preterm differences were confounded with connection length, we
#   computed the Euclidean distance between ROI centroids using age-specific coordinate templates and
#   correlated it with the edge-wise group differences (Spearman; significance by permutation of ROI
#   labels, 5,000 permutations)."
#   Results: "Group differences were unrelated to distance (rho = -0.002, permutation p = .94; within
#   hemisphere rho = 0.03, p = .40; between hemispheres rho = -0.004, p = .91), and this held for
#   templates of every age (|rho| <= 0.03). Connectivity itself fell with distance in both groups
#   (fetal rho = -0.33; preterm rho = -0.44; both p < .001), as expected for a developing brain with
#   predominantly local organization."
#
# Part (3) (scan level) is EXPLORATORY and answers a different question: how strongly each scan's
# overall FC falls off with distance. It differs by group and by PMA (the preterm-fetal difference is
# about +0.045 rho at 28 wk, ~0 at 30 wk, and -0.03 to -0.19 at 32-40 wk; interaction p < .001, still
# present in 28-38 wk), and it differs by scanner within preterm scans (3T - 1.5T = +0.10, p < .001).
# Do NOT state that overall distance dependence is equivalent between groups. It does not change the
# result above, and it is not needed for the paper; if a reviewer asks, this is the answer: the groups
# differ in overall distance dependence in an age- and scanner-dependent way, but this does not drive
# the edge-wise group differences. The Wilcoxon test of edge distance by direction is descriptive only
# (it treats edges as independent).
# -----------------------------------------------------------------------------------------------
if (RUN_DISTANCE_CONTROL) {
  ga_avail <- c(23:37, 40)       # templates that exist (no 38 or 39)
  dist_ga  <- NA                 # template for part (1); NA = the template nearest the median scan PMA
  n_perm   <- 5000               # permutations
  nearest_ga <- function(pma) ga_avail[apply(abs(outer(pma, ga_avail, "-")), 1, which.min)]
  if (is.na(dist_ga)) dist_ga <- nearest_ga(median(scanInfo$ScanPMA))

  read_xyz <- function(ga) {
    cd <- read.csv(file.path(dataDir, "labelsAndCoords", sprintf("ga%dCoords_reformat.csv", ga)), stringsAsFactors = FALSE)
    cc <- intersect(c("X", "Y", "Z"), names(cd))
    stopifnot("ROI" %in% names(cd), all(c("X", "Y") %in% cc))
    m <- as.matrix(cd[, cc, drop = FALSE]); rownames(m) <- cd$ROI; m
  }
  e1 <- columnMembershipCortOnly$ROI1; e2 <- columnMembershipCortOnly$ROI2          # "PreCG-L" format, same order as allConn
  stopifnot(identical(fetalRes$Connection, allConn))
  nodes <- sort(unique(c(e1, e2))); i1 <- match(e1, nodes); i2 <- match(e2, nodes)
  Dnode <- lapply(ga_avail, function(ga) {
    xyz <- read_xyz(ga); stopifnot(all(nodes %in% rownames(xyz)))
    as.matrix(dist(xyz[nodes, , drop = FALSE]))                                     # node-by-node distances
  })
  names(Dnode) <- ga_avail
  Dedge <- sapply(Dnode, function(D) D[cbind(i1, i2)])                               # edges x templates
  message("Distances use coordinate columns: ", paste(colnames(read_xyz(ga_avail[1])), collapse = ", "),
          "; part (1) uses the GA ", dist_ga, " template")
  same_hemi <- substr(e1, nchar(e1), nchar(e1)) == substr(e2, nchar(e2), nchar(e2))
  subsets   <- list(`All edges` = seq_along(allConn), `Within hemisphere` = which(same_hemi), `Between hemispheres` = which(!same_hemi))

  ## ---- (1) edge level, one template ----
  Dmat      <- Dnode[[as.character(dist_ga)]]
  edge_dist <- Dedge[, as.character(dist_ga)]
  outcomes  <- list(`Group difference (beta, preterm - fetal)` = fetalRes$scanTypeb,
                    `Fetal mean FC (adjusted)`   = fetalRes$emmMean.Fetal,
                    `Preterm mean FC (adjusted)` = fetalRes$emmMean.Preterm)
  mantel_rho <- function(y, idx, nperm = n_perm) {
    ry <- rank(y[idx]); obs <- cor(rank(edge_dist[idx]), ry)
    set.seed(2026)
    perm <- replicate(nperm, { p <- sample(length(nodes)); cor(rank(Dmat[cbind(p[i1[idx]], p[i2[idx]])]), ry) })
    c(rho = obs, p_perm = (1 + sum(abs(perm) >= abs(obs))) / (1 + nperm), n_edges = length(idx))
  }
  dist_tab <- do.call(rbind, lapply(names(outcomes), function(o) do.call(rbind, lapply(names(subsets), function(g)
    data.frame(Outcome = o, Edges = g, t(mantel_rho(outcomes[[o]], subsets[[g]])), check.names = FALSE)))))
  dist_tab$rho <- round(dist_tab$rho, 3); dist_tab$p_perm <- round(dist_tab$p_perm, 4)
  print(dist_tab)
  si_table(dist_tab, "distance", "distance_mantel")

  dir_sig <- ifelse(fetalRes$scanTypeb[sigFetal] > 0, "Preterm greater", "Fetal greater")      # distance of the significant edges, by direction
  print(tapply(edge_dist[sigFetal], dir_sig, function(x) round(c(n = length(x), median = median(x), q25 = quantile(x, .25, names = FALSE), q75 = quantile(x, .75, names = FALSE)), 1)))
  print(wilcox.test(edge_dist[sigFetal] ~ dir_sig))

  ## ---- (2) sensitivity: distances from every template age ----
  rho_by_ga <- data.frame(GA = ga_avail,
                          rho_beta = sapply(as.character(ga_avail), function(g) cor(Dedge[, g], fetalRes$scanTypeb, method = "spearman")))
  rho_by_ga$rho_beta <- round(rho_by_ga$rho_beta, 3); rownames(rho_by_ga) <- NULL
  print(rho_by_ga)
  si_table(rho_by_ga, "distance", "distance_rhoByTemplate")

  ## ---- (3) scan level: how strongly does each scan's FC fall off with distance? ----
  FCm      <- as.matrix(goodFCData[scanInfo$longID, allConn])
  scan_ga  <- nearest_ga(scanInfo$ScanPMA)                                           # template for each scan
  coupling <- function(idx) {                                                        # Spearman rho of FC vs distance, per scan
    R <- t(apply(FCm[, idx, drop = FALSE], 1, rank)); out <- numeric(nrow(scanInfo))
    for (g in unique(scan_ga)) { sc <- which(scan_ga == g)
      out[sc] <- cor(t(R[sc, , drop = FALSE]), rank(Dedge[idx, as.character(g)]))[, 1] }
    out
  }
  coup_tab <- do.call(rbind, lapply(names(subsets), function(g) {
    d <- scanInfo; d$coupling <- coupling(subsets[[g]])
    m  <- lme(coupling ~ ScanPMA + scanner + scanType, random = ~ 1 | shortID, data = d, control = ctrl)
    tt <- summary(m)$tTable["scanTypePreterm", ]
    data.frame(Edges = g, fetal_mean_rho = mean(d$coupling[d$scanType == "Fetal"]),
               preterm_mean_rho = mean(d$coupling[d$scanType == "Preterm"]),
               b_preterm_minus_fetal = unname(tt["Value"]), se = unname(tt["Std.Error"]), p = unname(tt["p-value"]))
  }))
  coup_tab[-1] <- lapply(coup_tab[-1], function(x) signif(x, 3))
  si_table(coup_tab, "distance", "distance_coupling_groupDifference")
  print(coup_tab)                       # negative rho = FC falls with distance; b > 0 means a weaker fall-off in preterm scans

  # Does the group difference depend on PMA? (the fetal and preterm lines can cross, so the single
  # group coefficient above averages over different differences at different ages)
  d_all  <- scanInfo; d_all$coupling <- coupling(subsets[["All edges"]])
  m_int  <- lme(coupling ~ ScanPMA * scanType + scanner, random = ~ 1 | shortID, data = d_all, control = ctrl)
  print(round(summary(m_int)$tTable[c("ScanPMA", "scanTypePreterm", "ScanPMA:scanTypePreterm"), c("Value", "Std.Error", "p-value")], 4))
  # rows: slope with PMA in fetal scans | preterm-minus-fetal at PMA = 0 (not meaningful) | difference in slopes
  b_int <- nlme::fixef(m_int); V_int <- vcov(m_int)
  dfi   <- min(m_int$fixDF$X[c("scanTypePreterm", "ScanPMA:scanTypePreterm")])
  diff_at_age <- t(sapply(c(28, 30, 32, 35, 38, 40), function(a) {                   # preterm minus fetal at a given PMA
    w <- setNames(numeric(length(b_int)), names(b_int)); w["scanTypePreterm"] <- 1; w["ScanPMA:scanTypePreterm"] <- a
    est <- sum(w * b_int); se <- sqrt(drop(t(w) %*% V_int %*% w))
    c(PMA = a, diff = est, se = se, p = 2 * pt(-abs(est / se), dfi))
  }))
  print(round(diff_at_age, 4))
  si_table(round(diff_at_age, 4), "distance", "distance_coupling_differenceByPMA")

  # Scanner effect on distance coupling within preterm scans: harmonization removes per-edge scanner
  # means and variances, but not necessarily a difference in the PATTERN of FC across edges.
  m_sc <- lme(coupling ~ ScanPMA + scanner, random = ~ 1 | shortID, data = d_all[d_all$scanType == "Preterm", ], control = ctrl)
  print(round(summary(m_sc)$tTable[scn_term, c("Value", "Std.Error", "p-value")], 4))       # 3T minus 1.5T
  # Interaction restricted to the ages where both groups are well sampled (avoids extrapolating the fetal line)
  rng_both <- c(28, 38)                                                                     # edit as needed
  d_rng    <- d_all[d_all$ScanPMA >= rng_both[1] & d_all$ScanPMA <= rng_both[2], ]
  print(table(d_rng$scanType))
  m_rng <- lme(coupling ~ ScanPMA * scanType + scanner, random = ~ 1 | shortID, data = d_rng, control = ctrl)
  print(round(summary(m_rng)$tTable[c("ScanPMA", "scanTypePreterm", "ScanPMA:scanTypePreterm"), c("Value", "Std.Error", "p-value")], 4))

  ## ---- Figures ----
  dd <- data.frame(dist = edge_dist, beta = fetalRes$scanTypeb)
  dd$bin <- cut(dd$dist, breaks = unique(quantile(dd$dist, seq(0, 1, 0.1))), include.lowest = TRUE)
  bin_df <- do.call(rbind, lapply(split(dd, dd$bin), function(g) data.frame(dist = mean(g$dist), mean = mean(g$beta), se = sd(g$beta) / sqrt(nrow(g)))))
  p_dist <- ggplot(dd, aes(dist, beta)) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    geom_point(color = "#8a8a8a", alpha = 0.18, size = 0.5, stroke = 0) +
    geom_errorbar(data = bin_df, aes(x = dist, y = mean, ymin = mean - 1.96 * se, ymax = mean + 1.96 * se), width = 0, color = col_preterm, linewidth = 0.6) +
    geom_line(data = bin_df, aes(x = dist, y = mean), color = col_preterm, linewidth = 0.6) +
    geom_point(data = bin_df, aes(x = dist, y = mean), color = col_preterm, size = 1.6) +
    labs(x = paste0("Distance between ROIs (GA ", dist_ga, " template units)"), y = "Preterm - fetal difference (Fisher z)") +
    theme_classic(base_size = 9) + theme(axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4))
  si_fig("distance", "groupDifferenceByDistance", 3.4, 3.0, function() print(p_dist))

  dc <- scanInfo; dc$coupling <- coupling(subsets[["All edges"]])
  p_coup <- ggplot(dc, aes(ScanPMA, coupling, color = scanType)) +
    geom_point(alpha = 0.5, size = 0.9, stroke = 0) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.8) +
    scale_color_manual(values = c(Fetal = col_fetal, Preterm = col_preterm)) +
    labs(x = "PMA at scan (weeks)", y = "Distance dependence of FC (Spearman rho)", color = NULL) +
    theme_classic(base_size = 9) + theme(axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4), legend.position = "bottom")
  si_fig("distance", "distanceCouplingByPMA", 3.4, 3.2, function() print(p_coup))
}

## ---- 9j. Motion (mean framewise displacement): collinearity with group, and motion vs. FC within each group ----
# WHAT IT TESTS: two follow-ups to the motion control in 9e. (A) How collinear is meanFD with group in the fetal-vs-
#   preterm model (variance inflation, and how much the group standard errors grow when meanFD is added: this is why
#   the raw-meanFD run, D1 in 9e, cannot be interpreted). (B) Does motion relate to connectivity WITHIN each group,
#   where motion does vary? Per edge: FCz ~ PMA (+ scanner for preterm) + meanFD with a random intercept per infant,
#   fetal and preterm separately, FDR across the 3,003 edges. Then: do the edge-wise motion effects resemble the
#   edge-wise fetal-vs-preterm differences? (If the group difference were a motion artifact, edges that motion
#   changes most should also be the edges that differ most between groups, and in the matching direction.)
# HOW TO READ THE NUMBERS (for later reference):
#   * Variance inflation factor (VIF) = 1 / (1 - R2), where R2 is how well the other predictors explain this one. Above
#     about 5 is serious: the model cannot tell the two variables apart. Standard errors grow by sqrt(VIF with / VIF
#     without) = sqrt(7.16 / 1.33) = 2.32, which is what was observed, so the lost significance of the raw-meanFD run
#     (D1 in 9e) is a precision problem, not a changed effect. Centering meanFD within group (D2 in 9e) removes the
#     collinearity and leaves the group beta unchanged, so D2 only shows that motion variation WITHIN groups does not
#     drive the result; it cannot say anything about the between-group motion difference.
#   * "Raw p < .05" in the fetal scan is 203 of 3,003 edges; chance alone gives about 150 (5%), so there is a small
#     excess but only 1 edge survives FDR.
#   * "Per SD" effects are FC change (Fisher z) for a 1-SD higher motion within that group (fetal SD 0.219 mm, preterm
#     SD 0.032 mm), which makes the two groups comparable despite their different motion scales; per-mm slopes are
#     the same effects divided by the group's SD (preterm -0.38 z/mm, fetal -0.002 z/mm).
#   * The extrapolation argument (see INTERPRETATION) is a plausibility check, not a test: the preterm slope is
#     estimated over 0.02-0.22 mm and cannot be extended to the fetal range; the fetal data cover that range and show
#     no slope.
#   * LIMITATION to state: group and motion are confounded by design (every fetal scan has more motion than every
#     preterm scan), so no model can separate them; fetal and infant FD may also not be measured on the same footing.
# RESULT (A): meanFD is almost the same variable as group (r = -0.925, R2 = 0.86). Variance inflation of the group term
#   rises from 1.33 to 7.16 when meanFD is added (meanFD itself 6.97; PMA 1.70 and scanner 1.86 unchanged). The group
#   standard errors grow by a median factor of 2.32 across edges, exactly what the VIF predicts (sqrt(7.16/1.33)).
#   (B) fetal scans only (192 scans, 161 infants): 203 edges with raw p < .05 (about 150 expected by chance), 1 FDR-
#   significant; median effect -0.0005 z per SD of FD (51% of edges negative; largest |effect| 0.148 z). Preterm scans
#   only (346 scans, 179 infants): 843 raw, 392 FDR-significant; median -0.0121 z per SD of FD (70% negative; largest
#   0.111). Edge-wise Spearman correlations: fetal vs. preterm motion effects 0.034; fetal motion effect vs. group
#   difference -0.112; preterm motion effect vs. group difference 0.057.
# INTERPRETATION: (A) with VIF > 7 the raw-meanFD model (D1 in 9e) cannot separate group from motion, which explains its
#   loss of significant edges (inflated standard errors), not a change in the effect. (B) In the fetal data, where
#   motion varies most (SD 0.22 mm), motion is essentially unrelated to connectivity (1 of 3,003 edges). In preterm
#   infants motion is modestly associated with lower connectivity (about -0.012 z per SD of 0.032 mm, or -0.38 z per
#   mm, over a 0.02-0.22 mm range). Extrapolating that preterm slope across the 0.68 mm fetal-preterm gap would give
#   about 0.26 z, on the order of the group differences, but that is far outside the preterm range and the fetal data
#   (-0.002 z per mm over 0.26-1.4 mm) do not show it. The edge patterns of the motion effects also do not resemble the
#   group differences (|rho| <= 0.11). Together this argues against motion as the explanation of the fetal-preterm
#   difference, though group and motion cannot be separated statistically (state as a limitation).
# REPORT: "Mean framewise displacement was almost completely collinear with group (r = -0.93; variance inflation 7.2),
#   so it could not be included as a covariate. Within the fetal scans, where motion varied most, it was unrelated to
#   connectivity (1 of 3,003 edges significant after FDR), and the edge-wise motion effects did not resemble the
#   fetal-preterm differences (Spearman rho = -0.11 and 0.06 for fetal and preterm motion effects)."
if (RUN_QC_COVARIATE_CONTROL) {
  ## A. collinearity of meanFD with group
  dv <- data.frame(group = as.numeric(si_fd$scanType == "Preterm"), PMA = si_fd$ScanPMA,
                   scanner = as.numeric(si_fd$scanner == levels(si_fd$scanner)[2]), meanFD = si_fd$meanFD)
  vif_of <- function(d, vars) sapply(setNames(vars, vars), function(v) 1 / (1 - summary(lm(reformulate(setdiff(vars, v), v), data = d))$r.squared))
  vif_tab <- rbind(`without meanFD` = c(vif_of(dv, c("group", "PMA", "scanner")), meanFD = NA),
                   `with meanFD`    = vif_of(dv, c("group", "PMA", "scanner", "meanFD")))
  print(round(vif_tab, 2))                                                   # variance inflation factors (above ~5 = serious)
  print(c(r_group_meanFD = round(cor(dv$group, dv$meanFD), 3), R2_group_on_meanFD = round(cor(dv$group, dv$meanFD)^2, 3)))
  se_ratio <- res_fd_raw$scanTypese / fetalRes$scanTypese                    # group SE with meanFD / group SE without, per edge
  print(c(median_SE_ratio = round(median(se_ratio, na.rm = TRUE), 2),
          expected_from_VIF = round(sqrt(vif_tab["with meanFD", "group"] / vif_tab["without meanFD", "group"]), 2)))

  ## B. motion vs. connectivity within each group
  fd_edge <- function(cn, d, terms) {
    d$FCz <- as.numeric(goodFCData[d$longID, cn])
    tryCatch({
      m  <- lme(reformulate(terms, "FCz"), random = ~ 1 | shortID, data = d, control = ctrl)
      tt <- summary(m)$tTable["meanFD", ]
      c(beta = unname(tt["Value"]), se = unname(tt["Std.Error"]), p = unname(tt["p-value"]))
    }, error = function(e) c(beta = NA, se = NA, p = NA))
  }
  fd_scan <- function(d, label, terms) {
    res <- cached(paste0("fdEdge_", label), key = paste(data_key, "fdEdge", label, paste(terms, collapse = "+"), table_key(d)),
                  data.frame(Connection = allConn, do.call(rbind, par_lapply(allConn, fd_edge, d = d, terms = terms, export = c("goodFCData", "ctrl")))))
    res$pAdj      <- p.adjust(res$p, "fdr")
    res$beta_perSD <- res$beta * sd(d$meanFD)                                # FC change (Fisher z) per 1 SD of motion within this group
    res
  }
  fd_summary <- function(res, d) c(scans = nrow(d), infants = nlevels(droplevels(d$shortID)), raw_p05 = sum(res$p < .05, na.rm = TRUE),
                                   fdr05 = sum(res$pAdj < .05, na.rm = TRUE), median_beta_perSD = round(median(res$beta_perSD, na.rm = TRUE), 4),
                                   pct_beta_negative = round(100 * mean(res$beta < 0, na.rm = TRUE), 1),
                                   max_abs_beta_perSD = round(max(abs(res$beta_perSD), na.rm = TRUE), 3))
  d_fet <- droplevels(si_fd[si_fd$scanType == "Fetal", ]);   d_pre <- droplevels(si_fd[si_fd$scanType == "Preterm", ])
  res_fd_fet <- fd_scan(d_fet, "fetal",   c("ScanPMA", "meanFD"))                  # fetal scans are all 1.5T: no scanner term
  res_fd_pre <- fd_scan(d_pre, "preterm", c("ScanPMA", "scanner", "meanFD"))
  print(rbind(Fetal = fd_summary(res_fd_fet, d_fet), Preterm = fd_summary(res_fd_pre, d_pre)))

  ## Do the edge-wise motion effects resemble the group differences? (Spearman correlations across the 3,003 edges;
  ## group difference = preterm minus fetal, so an artifact driven by fetal motion would give a NEGATIVE correlation
  ## with the fetal motion effect if motion lowers FC, and a positive one if it raises it.)
  rho <- function(a, b) round(cor(a, b, method = "spearman", use = "complete.obs"), 3)
  print(c(rho_fetal_vs_preterm_motion_effect  = rho(res_fd_fet$beta_perSD, res_fd_pre$beta_perSD),
          rho_fetal_motion_vs_group_difference   = rho(res_fd_fet$beta_perSD, fetalRes$scanTypeb),
          rho_preterm_motion_vs_group_difference = rho(res_fd_pre$beta_perSD, fetalRes$scanTypeb)))

  ## SI outputs: collinearity, within-group motion effects, and the edge-wise motion-vs-group-difference comparison
  si_table(data.frame(Model = rownames(vif_tab), vif_tab, check.names = FALSE, row.names = NULL), "motion_followup", "meanFD_collinearity_VIF")
  si_table(data.frame(Check = c("r (group, meanFD)", "R2 (group on meanFD)", "Median SE ratio, group term (with / without meanFD)", "Expected from VIF"),
                      Value = c(round(cor(dv$group, dv$meanFD), 3), round(cor(dv$group, dv$meanFD)^2, 3),
                                round(median(se_ratio, na.rm = TRUE), 2), round(sqrt(vif_tab["with meanFD", "group"] / vif_tab["without meanFD", "group"]), 2))),
           "motion_followup", "meanFD_collinearity_summary")
  si_table(data.frame(Group = c("Fetal", "Preterm"), rbind(fd_summary(res_fd_fet, d_fet), fd_summary(res_fd_pre, d_pre))), "motion_followup", "meanFD_withinGroup_effects")
  si_table(data.frame(Comparison = c("Fetal vs. preterm motion effect", "Fetal motion effect vs. group difference", "Preterm motion effect vs. group difference"),
                      `Spearman rho` = c(rho(res_fd_fet$beta_perSD, res_fd_pre$beta_perSD), rho(res_fd_fet$beta_perSD, fetalRes$scanTypeb),
                                         rho(res_fd_pre$beta_perSD, fetalRes$scanTypeb)), check.names = FALSE),
           "motion_followup", "meanFD_edgewise_correlations")
  d_mot <- rbind(data.frame(Panel = "Fetal scans", x = res_fd_fet$beta_perSD, y = fetalRes$scanTypeb),
                 data.frame(Panel = "Preterm scans", x = res_fd_pre$beta_perSD, y = fetalRes$scanTypeb))
  ann_mot <- data.frame(Panel = c("Fetal scans", "Preterm scans"), x = -Inf, y = Inf,
                        lab = paste0("rho = ", c(rho(res_fd_fet$beta_perSD, fetalRes$scanTypeb), rho(res_fd_pre$beta_perSD, fetalRes$scanTypeb))))
  p_mot <- ggplot(d_mot, aes(x, y)) + geom_hline(yintercept = 0, linewidth = 0.3) + geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_point(color = "#8a8a8a", alpha = 0.2, size = 0.4, stroke = 0) +
    geom_text(data = ann_mot, aes(label = lab), hjust = -0.1, vjust = 1.5, size = 2.8) +
    facet_wrap(~Panel) +
    labs(x = "Within-group effect of motion on FC (Fisher z per SD of meanFD)", y = "Group difference (beta, preterm - fetal)") + si_theme()
  si_fig("motion_followup", "motionEffects_vs_groupDifference", 5.4, 3.0, function() print(p_mot))
}

## ---- 9k. Strict equal-volume control: ROI time series cut to the same number of volumes ----------------------
# WHAT IT TESTS: the strict version of the volume controls in 9e, for the fetal-vs-preterm comparison only. Instead of
#   adjusting for, or restricting on, the number of volumes, every scan's connectivity is RECOMPUTED from the first
#   vol_equal_n usable (non-censored) volumes of its ROI time series, so fetal and preterm scans contribute exactly the
#   same number of volumes, and the main fetal-vs-preterm model (PMA + scanner + group) is repeated. (The preterm-only
#   PMA x prematurity analysis is not repeated here: volumes are not capped within the preterm group.)
# DATA: one tab-separated .1D file per scan ("<fullID>_ts.1D"; rows = volumes, columns = the 111 ROIs in label-key
#   order, censored volumes are all-zero rows) in data/rawTS and/or data/rawRS (subfolders fetal and preterm).
#   Connectivity is the Fisher z of the correlation across the usable volumes. The truncated values are NOT
#   re-harmonized: scanner is a covariate in the model (as in 9h, where harmonization changed nothing for the
#   fetal-vs-preterm contrast).
# CHECKS built in: the time series are first used in full to reproduce the stored (unharmonized) connectivity of a few
#   scans, which would reveal a wrong column order or different preprocessing; usable-volume counts are compared with
#   finalVol; for every scan, the truncated connectivity is correlated with the full-data connectivity (fidelity).
# NOTE: all scans at vol_equal_n volumes means fetal scans keep all their data at the minimum (80 volumes = 4 min at TR
#   3 s) while preterm scans (>= 120 volumes at TR 2 s) are cut to 80 volumes = 2.7 min; duration and volume count
#   cannot both be equated.
# CHECKS (run): the full time series reproduce the stored connectivity (r = 1.000 for 6 test scans, so the column order
#   is right). Every scan has >= 80 usable (non-zero) volumes. finalVol differs from the usable count in 112 of 538 scans
#   (finalVol minus usable volumes: fetal, 6 scans, -1 to +100; preterm, 106 scans, -87 to +109; median 0). The stored
#   connectivity matches the usable volumes, so the usable count is the number of volumes the connectivity is based on.
# DESCRIPTIVES (Table S, equal-volume control), usable volumes before truncation, mean (SD) [range]:
#   fetal 192 scans / 161 infants, TR 3 s: 129.9 (47.8) [80-281] volumes = 6.5 (2.4) min; preterm 346 scans / 179
#   infants, TR 2 s: 216.7 (53.7) [120-296] volumes = 7.2 (1.8) min. After truncation: 80 volumes = 4.00 min (fetal),
#   2.67 min (preterm). Fidelity of the 80-volume connectivity to the full-data connectivity (per-scan correlation
#   across the 3,003 edges): fetal 0.93 (SD 0.07, min 0.67); preterm 0.84 (SD 0.08, min 0.65).
# RESULT (fetal vs. preterm; each row vs. the primary result):
#   primary (harmonized, all volumes) ...... 1,404 significant (1,118 preterm > fetal, 286 fetal > preterm)
#   unharmonized, all volumes (9h) ......... 1,439 (1,141 / 298); 1,404 (100%) primary edges still significant, same
#                                           direction 100%; r = 1.00, rho = 1.000, slope 1.000, beta ratio 1.001
#   unharmonized, 80 volumes per scan ...... 1,248 (1,002 / 246); 1,206 (85.9%) primary edges still significant;
#                                           same direction 1,404 (100%); r = 0.99, rho = 0.984, slope 1.000,
#                                           median beta ratio 0.996
# INTERPRETATION: the fetal-vs-preterm difference does not depend on the number of volumes. With every scan cut to 80
#   volumes the direction is preserved for all 1,404 primary edges, the betas keep their rank order (rho = 0.984) and
#   their size (regression slope 1.00, median ratio 0.996), so effects are not attenuated by using less data. The 14% of
#   primary edges that lose significance reflect larger standard errors (noisier connectivity from 80 volumes,
#   especially in preterm scans, fidelity 0.84), not smaller effects.
# OUTPUT FOR THE SUPPLEMENT (printed, and kept as objects): table_s_equal_desc = descriptives by group (scans, infants,
#   TR, usable volumes and minutes before truncation, volumes and minutes after, fidelity of the truncated connectivity
#   to the full-data connectivity); table_s_equal = results for the primary analysis, the unharmonized full-data
#   analysis (9h) and the equal-volume analysis (FDR-significant edges and direction split, how many primary edges
#   stay significant / keep their direction, Pearson and Spearman correlation of the betas with the primary betas,
#   the regression slope of the betas on the primary betas, and the median ratio of the betas on the primary
#   significant edges, which shows whether effect sizes shrink when scans are shortened). Both are written to csv when
#   SAVE_TABLES is TRUE (section 10).
# REPORT: "When every scan was reduced to the same number of volumes (80; 4.0 min for fetal scans at TR 3 s, 2.7 min for
#   preterm scans at TR 2 s), the fetal-preterm differences were unchanged: 1,248 edges were significant, 1,206 of the
#   1,404 primary edges (85.9%) remained significant, and all 1,404 had the same direction. Effect sizes were
#   preserved (Pearson r = 0.99 and Spearman rho = 0.98 with the primary betas; regression slope 1.00; median beta ratio
#   0.996). Before truncation, fetal scans had 129.9 (SD 47.8) usable volumes and preterm scans 216.7 (SD 53.7). The
#   truncated connectivity correlated 0.93 (fetal) and 0.84 (preterm) with the full-data connectivity within scans."
vol_equal_n    <- 80          # volumes per scan after truncation (the fetal minimum)
vol_equal_mode <- "first"     # "first": the first n usable volumes; "random": a random set of n usable volumes (seeded)
if (RUN_EQUAL_VOLUME_CONTROL) {
  if (length(ts_files) == 0) {                                          # ts_files and read_ts come from section 3e
    message("No time series found in data/rawTS or data/rawRS; skipping 9k.")
  } else {
    stopifnot(!anyDuplicated(names(ts_files)))
    miss <- setdiff(scanInfo$longID, names(ts_files))
    if (length(miss) > 0) stop(length(miss), " scans have no time series file, e.g. ", paste(head(miss, 3), collapse = ", "))

    ts_roi <- gsub("-", ".", roiData$FullLabel[1:111], fixed = TRUE)      # column order of the .1D files = first 111 rows of the label key
    parts  <- strsplit(allConn, "_", fixed = TRUE)
    ei <- match(sapply(parts, `[`, 1), ts_roi); ej <- match(sapply(parts, `[`, 2), ts_roi)
    stopifnot(!anyNA(ei), !anyNA(ej))
    used <- sort(unique(c(ei, ej))); ci <- match(ei, used); cj <- match(ej, used)
    fz <- function(x) { cm <- suppressWarnings(cor(x[, used])); atanh(pmin(pmax(cm[cbind(ci, cj)], -0.999999), 0.999999)) }   # Fisher z for the 3,003 edges
    ts_of <- function(id) { m <- read_ts(ts_files[[id]]); stopifnot(ncol(m) == 111); m[rowSums(abs(m)) > 0, , drop = FALSE] }   # usable volumes only

    ## check 1: full time series must reproduce the stored unharmonized connectivity
    chk_ids <- head(scanInfo$longID[scanInfo$scanType == "Fetal"], 3); chk_ids <- c(chk_ids, head(scanInfo$longID[scanInfo$scanType == "Preterm"], 3))
    chk <- t(sapply(chk_ids, function(id) c(usable_volumes = nrow(ts_of(id)), finalVol = scanInfo$finalVol[scanInfo$longID == id],
                                            r_with_stored = round(cor(fz(ts_of(id)), as.numeric(combFCData[id, allConn])), 4))))
    print(chk)
    if (any(chk[, "r_with_stored"] < 0.99)) stop("Recomputed connectivity does not match the stored values (r < .99): check the ROI column order / preprocessing before using this control.")

    trunc_key <- paste(data_key, "truncFC", vol_equal_n, vol_equal_mode, length(ts_files), sum(file.size(ts_files[scanInfo$longID])))
    tr <- cached("truncFC", key = trunc_key, {
      out <- matrix(NA_real_, nrow(scanInfo), length(allConn), dimnames = list(scanInfo$longID, allConn)); n_use <- integer(nrow(scanInfo)); fid <- numeric(nrow(scanInfo))
      for (k in seq_len(nrow(scanInfo))) {
        x <- ts_of(scanInfo$longID[k]); n_use[k] <- nrow(x)
        if (nrow(x) < vol_equal_n) stop(scanInfo$longID[k], " has only ", nrow(x), " usable volumes")
        if (vol_equal_mode == "random") { set.seed(2026 + k); sel <- sort(sample(nrow(x), vol_equal_n)) } else sel <- seq_len(vol_equal_n)
        out[k, ] <- fz(x[sel, , drop = FALSE]); fid[k] <- cor(out[k, ], fz(x))        # fidelity: truncated vs. full-data connectivity
      }
      list(fc = as.data.frame(out), usable = n_use, fidelity = fid)
    })
    truncFC <- tr$fc
    ## check 2 and fidelity
    print(c(scans_where_usable_ne_finalVol = sum(tr$usable != scanInfo$finalVol)))
    print(do.call(rbind, tapply(scanInfo$finalVol - tr$usable, scanInfo$scanType,                    # finalVol minus the all-zero (censored) rows
                                function(x) c(scans_differing = sum(x != 0), min = min(x), median = median(x), max = max(x)))))
    print(do.call(rbind, tapply(tr$fidelity, scanInfo$scanType, function(x) c(mean = mean(x), min = min(x), sd = sd(x)))))

    ## fetal vs. preterm on the truncated connectivity (scanner is a covariate; values are not harmonized)
    res_eq_fp <- with_fc(truncFC, run_fetal(scanInfo, paste0("equalVol", vol_equal_n, vol_equal_mode)), paste0("truncFC", vol_equal_n))
    print(cmp_fetal(res_eq_fp, n_scans = nrow(scanInfo)))
    if (exists("fetalRes_raw")) print(c(r_beta_with_unharmonized_full = round(cor(fetalRes_raw$scanTypeb, res_eq_fp$scanTypeb, use = "complete.obs"), 3)))

    ## ---- Tables for the supplement ----
    # (1) Descriptives: scans, infants, usable volumes and minutes before and after truncation, and the fidelity of the
    #     truncated connectivity to the full-data connectivity, by group.
    fmt_ms <- function(x, d = 1) sprintf(paste0("%.", d, "f (%.", d, "f)"), mean(x), sd(x))
    table_s_equal_desc <- do.call(rbind, lapply(c("Fetal", "Preterm"), function(g) {
      i  <- scanInfo$scanType == g
      tr_s <- TR_SEC[[g]]
      data.frame(Group = g, Scans = sum(i), Infants = length(unique(scanInfo$shortID[i])), TR_s = tr_s,
                 `Usable volumes before, mean (SD)` = fmt_ms(tr$usable[i]), `Usable volumes before, range` = paste(range(tr$usable[i]), collapse = "-"),
                 `Minutes before, mean (SD)` = fmt_ms(tr$usable[i] * tr_s / 60),
                 `Volumes after` = vol_equal_n, `Minutes after` = round(vol_equal_n * tr_s / 60, 2),
                 `Fidelity r, mean (SD)` = fmt_ms(tr$fidelity[i], 2), `Fidelity r, minimum` = round(min(tr$fidelity[i]), 2),
                 check.names = FALSE)
    }))
    print(table_s_equal_desc, row.names = FALSE)

    # (2) Results: the primary analysis, the unharmonized full-data analysis (9h, if it was run) and the equal-volume
    #     analysis, each against the primary result. beta ratio = median over the primary significant edges of
    #     (beta in this analysis / primary beta); slope = regression slope of this analysis' betas on the primary betas.
    table_s_equal <- rbind(ctrl_row(fetalRes, "Primary: harmonized, all volumes"),
                           if (exists("fetalRes_raw")) ctrl_row(fetalRes_raw, "Unharmonized, all volumes"),
                           ctrl_row(res_eq_fp, paste0("Unharmonized, ", vol_equal_n, " volumes per scan")))
    print(t(table_s_equal), quote = FALSE)                                    # transposed for reading in the console
    si_table(table_s_equal,      "equal_volume", "equalVolume_results")
    si_table(table_s_equal_desc, "equal_volume", "equalVolume_descriptives")
    eq_list <- c(list(`Unharmonized, all volumes` = if (exists("fetalRes_raw")) fetalRes_raw), setNames(list(res_eq_fp), paste0("Unharmonized, ", vol_equal_n, " volumes")))
    eq_list <- eq_list[!vapply(eq_list, is.null, NA)]
    d_fid <- data.frame(scanType = scanInfo$scanType, fidelity = tr$fidelity)
    dm <- facet_dims(length(eq_list))
    si_fig("equal_volume", "betas_and_fidelity", dm[["w"]] + 3, dm[["h"]],
           function() print(ggpubr::ggarrange(beta_scatter(eq_list), dist_plot(d_fid, "fidelity", "Correlation of truncated with full-data connectivity"),
                                              nrow = 1, widths = c(dm[["w"]], 3), labels = c("A", "B"))))
  }
}

###############################################################################
# 10. SAVE TABLES AND RECORD SESSION
###############################################################################

if (SAVE_TABLES) {
  write.csv(cbind(Row = rownames(Table1_ByScan),        Table1_ByScan),        file.path(figureDir, "Table1_ByScan.csv"),        row.names = FALSE)
  write.csv(cbind(Row = rownames(Table1_ByParticipant), Table1_ByParticipant), file.path(figureDir, "Table1_ByParticipant.csv"), row.names = FALSE)
  write.csv(fetalRes,        file.path(figureDir, "Data_FetalVsPreterm_AllEdges.csv"),        row.names = FALSE)
  write.csv(table2,          file.path(figureDir, "Table2_HomotopicFetalVsPreterm.csv"),               row.names = FALSE)
  write.csv(intxRes,         file.path(figureDir, "Data_PmaByPrematurity_AllEdges.csv"),      row.names = FALSE)
  write.csv(tab3_display,    file.path(figureDir, "Table3_PmaByPrematurity.csv"),                         row.names = FALSE)
  write.csv(slopesMain_wide, file.path(figureDir, "Data_PmaByPrematurity_SlopesSigEdges.csv"),    row.names = FALSE)
}

###############################################################################
# 11. SEX DIFFERENCES (exploratory)
#     Repeats the two main comparisons with sex in the model, emphasizing sex effects:
#       11a  descriptives: sex by group, birth group and scanner
#       11b  fetal vs. preterm:  FC ~ PMA + scanner + group * sex   (sex main effect, group x sex), plus the
#            fetal-vs-preterm scan run separately in females and males
#       11c  within preterm:     FC ~ scanner + PMA * birth group + sex + PMA:sex + birth group:sex
#            (sex main effect, PMA x sex, birth group x sex; the PMA x birth group test is carried along),
#            the three-way PMA x birth group x sex test, and the primary PMA x birth group scan adjusted for sex
#            and within each sex (11c); 11d: focused 1-df test on the 16 primary edges of whether the MPT - EPT and
#            MPT - VPT slope contrasts differ between males and females
#     All edge-wise tests are corrected across the 3,003 edges (BH-FDR) separately for each term. These are
#     exploratory: the sample is small for interactions (especially the three-way test), so each summary table
#     also gives the raw-p count against the 5% expected by chance. Sex is coded from goodData$sex (0 = female,
#     1 = male, as in Table 1; character codes beginning F / M are also accepted); check the printed coding table.
#     PMA is centered in the preterm models so the sex coefficient is the male - female difference at the mean PMA.
#     Output numbers: "sex_fetal_preterm" and "sex_preterm" in SI_ID.
###############################################################################

if (RUN_SEX_ANALYSES) {
  sex_cols <- c(Female = "#cc79a7", Male = "#009e73")
  snm <- function(x) gsub("[^A-Za-z0-9]+", "_", x)                       # term / coefficient name -> column name stem

  ## ---- coding and scan tables ------------------------------------------------------------------------------
  sex_raw <- scanInfo$sex
  sexM <- if (is.numeric(sex_raw)) ifelse(sex_raw %in% c(0, 1), sex_raw, NA) else {
    s1 <- toupper(substr(trimws(as.character(sex_raw)), 1, 1)); ifelse(s1 == "M", 1, ifelse(s1 == "F", 0, NA)) }
  print(table(raw = sex_raw, coded_male = sexM, useNA = "ifany"))        # check: 0 = female, 1 = male
  si_sex <- scanInfo; si_sex$sexF <- factor(ifelse(sexM == 1, "Male", "Female"), c("Female", "Male"))
  message("Scans without a sex value (dropped): ", sum(is.na(si_sex$sexF)))
  si_sex <- droplevels(si_sex[!is.na(si_sex$sexF), ])
  stopifnot(all(tapply(as.character(si_sex$sexF), si_sex$shortID, function(x) length(unique(x))) == 1))   # one sex per infant
  si_sex$grp <- factor(ifelse(si_sex$scanType == "Fetal", "Fetal", as.character(si_sex$BirthCat)), c("Fetal", "EPT", "VPT", "MPT"))
  print(table(si_sex$scanType, si_sex$sexF))
  stopifnot(all(table(si_sex$scanType, si_sex$sexF) > 0))
  pt_sex <- droplevels(si_sex[si_sex$scanType == "Preterm", ])
  pt_sex$ScanPMAc <- pt_sex$ScanPMA - mean(pt_sex$ScanPMA)
  pt_f <- droplevels(pt_sex[pt_sex$sexF == "Female", ]); pt_m <- droplevels(pt_sex[pt_sex$sexF == "Male", ])

  ## ---- helpers -----------------------------------------------------------------------------------------------
  # One lmer per edge: Type II Wald chi-square for each of `terms` and the estimate / SE of each of `coefs`
  # (coefficient names as in summary(lmer)); FDR-adjusted p per term. Cached; parallel over edges.
  sex_scan <- function(d0, fml, terms, coefs = character(0), tag) {
    fml_re <- update(fml, . ~ . + (1 | shortID))
    nm <- c(unlist(lapply(terms, function(t) paste0(snm(t), c("_chi", "_df", "_p")))),
            unlist(lapply(coefs, function(k) paste0(snm(k), c("_b", "_se")))))
    chk <- local({ d <- d0; d$FCz <- as.numeric(goodFCData[d$longID, allConn[1]])        # the names must exist, or every edge would silently be NA
      m <- suppressMessages(suppressWarnings(lmer(fml_re, data = d)))
      list(an = rownames(car::Anova(m, type = "II", test = "Chisq")), co = rownames(summary(m)$coefficients)) })
    stopifnot(all(terms %in% chk$an), all(coefs %in% chk$co))
    key <- paste(data_key, "sexScan", tag, paste(deparse(fml), collapse = ""), paste(terms, collapse = "|"), paste(coefs, collapse = "|"), table_key(d0))
    cached(paste0("sexScan_", tag), key = key, {
      out <- par_lapply(allConn, function(cn) {
        d <- d0; d$FCz <- as.numeric(goodFCData[d$longID, cn])
        tryCatch({
          m  <- suppressMessages(suppressWarnings(lmer(fml_re, data = d)))
          an <- car::Anova(m, type = "II", test = "Chisq"); sm <- summary(m)$coefficients
          setNames(c(unlist(lapply(terms, function(t) c(an[t, "Chisq"], an[t, "Df"], an[t, "Pr(>Chisq)"]))),
                     unlist(lapply(coefs, function(k) c(sm[k, "Estimate"], sm[k, "Std. Error"])))), nm)
        }, error = function(e) setNames(rep(NA_real_, length(nm)), nm))
      })
      res <- data.frame(Connection = allConn, do.call(rbind, out), row.names = NULL)
      for (t in terms) res[[paste0(snm(t), "_pAdj")]] <- p.adjust(res[[paste0(snm(t), "_p")]], "fdr")
      res
    })
  }
  sex_term_summary <- function(res, terms, labels = terms) do.call(rbind, lapply(seq_along(terms), function(i) {
    p <- res[[paste0(snm(terms[i]), "_p")]]; pa <- res[[paste0(snm(terms[i]), "_pAdj")]]
    data.frame(Term = labels[i], `Edges tested` = sum(!is.na(p)), `Raw p < .05` = sum(p < .05, na.rm = TRUE),
               `Expected by chance` = round(0.05 * sum(!is.na(p))), `FDR p < .05` = sum(pa < .05, na.rm = TRUE),
               `Smallest FDR p` = signif(suppressWarnings(min(pa, na.rm = TRUE)), 3), check.names = FALSE)
  }))
  sex_top <- function(res, term, coef = NULL, n = 20) {                  # the n edges with the smallest p for one term
    s <- snm(term); o <- order(res[[paste0(s, "_p")]])[seq_len(n)]
    out <- data.frame(Connection = res$Connection[o], chi = round(res[[paste0(s, "_chi")]][o], 2),
                      p = signif(res[[paste0(s, "_p")]][o], 3), pAdj = signif(res[[paste0(s, "_pAdj")]][o], 3))
    if (!is.null(coef)) { k <- snm(coef); out$beta <- round(res[[paste0(k, "_b")]][o], 4); out$se <- round(res[[paste0(k, "_se")]][o], 4) }
    out
  }
  sex_volcano <- function(res, term, coef, xlab) {
    s <- snm(term); k <- snm(coef)
    d <- data.frame(b = res[[paste0(k, "_b")]], y = -log10(res[[paste0(s, "_p")]]),
                    sig = factor(!is.na(res[[paste0(s, "_pAdj")]]) & res[[paste0(s, "_pAdj")]] < .05, c(TRUE, FALSE), c("FDR p < .05", "Not significant")))
    ggplot(d, aes(b, y, color = sig)) + geom_vline(xintercept = 0, linewidth = 0.3) +
      geom_point(size = 0.5, alpha = 0.5, stroke = 0) +
      scale_color_manual(values = c("FDR p < .05" = "#d55e00", "Not significant" = "#8a8a8a"), name = NULL, drop = FALSE) +
      labs(x = xlab, y = "-log10 p") + si_theme()
  }
  sex_qq <- function(res, terms, labels) {                              # observed vs. expected -log10 p (a bend above the line = more small p than chance)
    d <- do.call(rbind, lapply(seq_along(terms), function(i) {
      p <- sort(res[[paste0(snm(terms[i]), "_p")]]); p <- p[!is.na(p)]
      data.frame(Term = labels[i], exp = -log10(ppoints(length(p))), obs = -log10(p)) }))
    d$Term <- factor(d$Term, labels)
    ggplot(d, aes(exp, obs)) + geom_abline(linetype = "dashed", linewidth = 0.3) + geom_point(size = 0.5, stroke = 0) +
      facet_wrap(~Term, ncol = min(length(labels), 2)) + labs(x = "Expected -log10 p", y = "Observed -log10 p") + si_theme()
  }
  sex_edge_plot <- function(d, edges, by_type = FALSE) {                # FC vs. PMA by sex for chosen edges
    long <- do.call(rbind, lapply(edges, function(cn) data.frame(d[c("ScanPMA", "sexF", "scanType")], Edge = cn, FC = as.numeric(goodFCData[d$longID, cn]))))
    long$Edge <- factor(long$Edge, edges)
    ggplot(long, aes(ScanPMA, FC, color = sexF)) + geom_point(alpha = 0.35, size = 0.6, stroke = 0) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.7) +
      scale_color_manual(values = sex_cols, name = NULL) +
      (if (by_type) facet_grid(scanType ~ Edge, scales = "free") else facet_wrap(~Edge, scales = "free_y", ncol = 3)) +
      labs(x = "PMA at scan (weeks)", y = "Connectivity (Fisher z)") + si_theme()
  }

  ## ---- 11a. Descriptives -------------------------------------------------------------------------------------
  sex_desc <- do.call(rbind, lapply(levels(si_sex$grp), function(g) {
    d <- si_sex[si_sex$grp == g, ]; inf <- d[!duplicated(d$shortID), ]
    data.frame(Group = g, Infants = nrow(inf), `Female infants` = sum(inf$sexF == "Female"), `Male infants` = sum(inf$sexF == "Male"),
               `Male, %` = round(100 * mean(inf$sexF == "Male"), 1), Scans = nrow(d), `Female scans` = sum(d$sexF == "Female"), `Male scans` = sum(d$sexF == "Male"),
               `PMA at scan, female, mean (SD)` = sprintf("%.1f (%.1f)", mean(d$ScanPMA[d$sexF == "Female"]), sd(d$ScanPMA[d$sexF == "Female"])),
               `PMA at scan, male, mean (SD)`   = sprintf("%.1f (%.1f)", mean(d$ScanPMA[d$sexF == "Male"]),   sd(d$ScanPMA[d$sexF == "Male"])),
               `3T scans, female, %` = if (g == "Fetal") NA else round(100 * mean(d$scanner[d$sexF == "Female"] == levels(si_sex$scanner)[2]), 1),
               `3T scans, male, %`   = if (g == "Fetal") NA else round(100 * mean(d$scanner[d$sexF == "Male"]   == levels(si_sex$scanner)[2]), 1),
               check.names = FALSE)
  }))
  print(sex_desc, row.names = FALSE)
  inf_sex <- si_sex[!duplicated(si_sex$shortID), ]
  tryCatch(print(chisq.test(table(inf_sex$scanType, inf_sex$sexF))), error = function(e) message("Sex by group test skipped: ", e$message))
  tryCatch(print(chisq.test(table(inf_sex$BirthCat[inf_sex$scanType == "Preterm"], inf_sex$sexF[inf_sex$scanType == "Preterm"]))), error = function(e) message("Sex by birth group test skipped: ", e$message))
  tryCatch(print(t.test(ScanPMA ~ sexF, data = pt_sex)), error = function(e) message("PMA by sex test skipped: ", e$message))
  si_table(sex_desc, "sex_fetal_preterm", "sex_descriptives")

  ## ---- 11b. Fetal vs. preterm ----------------------------------------------------------------------------------
  # Model: FC ~ PMA + scanner + group * sex. Male - female difference within the fetal group = sexFMale; the change of that
  # difference in preterm scans = scanTypePreterm:sexFMale (the interaction tested by scanType:sexF).
  sexfp_terms <- c("sexF", "scanType:sexF"); sexfp_labels <- c("Sex (male vs. female)", "Group x sex")
  res_sex_fp <- sex_scan(si_sex, FCz ~ ScanPMA + scanner + scanType * sexF, sexfp_terms,
                         coefs = c("sexFMale", "scanTypePreterm:sexFMale"), tag = "fetalPreterm")
  sexfp_sum <- sex_term_summary(res_sex_fp, sexfp_terms, sexfp_labels); print(sexfp_sum, row.names = FALSE)
  sexfp_top_main <- sex_top(res_sex_fp, "sexF", "sexFMale");                         print(head(sexfp_top_main, 10))
  sexfp_top_intx <- sex_top(res_sex_fp, "scanType:sexF", "scanTypePreterm:sexFMale"); print(head(sexfp_top_intx, 10))
  # Direction of the FDR-significant edges (sexFMale > 0 = higher FC in males)
  for (k in c("sexF", "scanType:sexF")) {
    sg <- res_sex_fp[[paste0(snm(k), "_pAdj")]] < .05 & !is.na(res_sex_fp[[paste0(snm(k), "_pAdj")]])
    bk <- res_sex_fp[[paste0(snm(if (k == "sexF") "sexFMale" else "scanTypePreterm:sexFMale"), "_b")]]
    message(k, ": ", sum(sg), " FDR-significant edges (", sum(sg & bk > 0), " with a positive coefficient, ", sum(sg & bk < 0), " negative)")
  }
  # Fetal-vs-preterm scan run separately in each sex (the primary scan model); compared with the primary result and with each other
  sex_f <- droplevels(si_sex[si_sex$sexF == "Female", ]); sex_m <- droplevels(si_sex[si_sex$sexF == "Male", ])
  res_sex_f <- run_fetal(sex_f, "sexFemale"); res_sex_m <- run_fetal(sex_m, "sexMale")
  print(cmp_fetal(res_sex_f, n_scans = nrow(sex_f))); print(cmp_fetal(res_sex_m, n_scans = nrow(sex_m)))
  sex_beta_r <- c(r_female_vs_male_betas = round(cor(res_sex_f$scanTypeb, res_sex_m$scanTypeb, use = "complete.obs"), 3),
                  rho_female_vs_male_betas = round(cor(res_sex_f$scanTypeb, res_sex_m$scanTypeb, method = "spearman", use = "complete.obs"), 3),
                  both_FDR_sig = sum(res_sex_f$scanTypep.Adj < .05 & res_sex_m$scanTypep.Adj < .05, na.rm = TRUE),
                  female_only_FDR_sig = sum(res_sex_f$scanTypep.Adj < .05 & !(res_sex_m$scanTypep.Adj < .05), na.rm = TRUE),
                  male_only_FDR_sig = sum(res_sex_m$scanTypep.Adj < .05 & !(res_sex_f$scanTypep.Adj < .05), na.rm = TRUE),
                  same_direction_of_primary_sig_edges = sum(sigFetal & sign(res_sex_f$scanTypeb) == sign(res_sex_m$scanTypeb), na.rm = TRUE))
  print(sex_beta_r)

  si_table(sexfp_sum,      "sex_fetal_preterm", "sexTests_summary")
  si_table(sexfp_top_main, "sex_fetal_preterm", "sexMainEffect_topEdges")
  si_table(sexfp_top_intx, "sex_fetal_preterm", "groupBySex_topEdges")
  si_table(data.frame(Measure = names(sex_beta_r), Value = unname(sex_beta_r)), "sex_fetal_preterm", "bySex_femaleVsMale_agreement")
  si_fetal("sex_fetal_preterm", "bySex", list(`Females only` = res_sex_f, `Males only` = res_sex_m),
           list(`Females only` = nrow(sex_f), `Males only` = nrow(sex_m)))                              # results table + betas vs. primary
  d_fm <- data.frame(f = res_sex_f$scanTypeb, m = res_sex_m$scanTypeb, sig = sigFetal)
  p_fm <- ggplot(d_fm, aes(f, m)) + geom_abline(linetype = "dashed", linewidth = 0.3) +
    geom_point(data = d_fm[!d_fm$sig, ], color = "#8a8a8a", alpha = 0.2, size = 0.4, stroke = 0) +
    geom_point(data = d_fm[d_fm$sig, ], color = col_preterm, alpha = 0.35, size = 0.5, stroke = 0) +
    labs(x = "Females: group difference (beta)", y = "Males: group difference (beta)") + si_theme()
  top_fp <- head(sexfp_top_intx$Connection, 4)
  si_fig("sex_fetal_preterm", "sexEffects", 8.2, 6.6, function()
    print(ggpubr::ggarrange(
      ggpubr::ggarrange(sex_volcano(res_sex_fp, "sexF", "sexFMale", "Male - female (beta, fetal group)"),
                        sex_volcano(res_sex_fp, "scanType:sexF", "scanTypePreterm:sexFMale", "Group x sex (beta)"),
                        p_fm, nrow = 1, labels = c("A", "B", "C"), common.legend = TRUE, legend = "bottom"),
      sex_edge_plot(si_sex, top_fp, by_type = TRUE), ncol = 1, heights = c(1, 1.15), labels = c("", "D"))))
  si_fig("sex_fetal_preterm", "PvalueQQ", 4.6, 2.9, function() print(sex_qq(res_sex_fp, sexfp_terms, sexfp_labels)))

  ## ---- 11c. Within preterm ---------------------------------------------------------------------------------------
  print(group_counts(pt_f)); print(group_counts(pt_m))
  pts_terms  <- c("sexF", "ScanPMAc:sexF", "BirthCat:sexF", "ScanPMAc:BirthCat")
  pts_labels <- c("Sex (male vs. female)", "PMA x sex", "Birth group x sex", "PMA x birth group")
  res_sex_pt <- sex_scan(pt_sex, FCz ~ scanner + ScanPMAc * BirthCat + sexF + ScanPMAc:sexF + BirthCat:sexF, pts_terms,
                         coefs = c("sexFMale", "ScanPMAc:sexFMale"), tag = "preterm")
  pts_sum <- sex_term_summary(res_sex_pt, pts_terms, pts_labels); print(pts_sum, row.names = FALSE)
  pts_top_main <- sex_top(res_sex_pt, "sexF", "sexFMale");              print(head(pts_top_main, 10))
  pts_top_pma  <- sex_top(res_sex_pt, "ScanPMAc:sexF", "ScanPMAc:sexFMale"); print(head(pts_top_pma, 10))
  pts_top_bc   <- sex_top(res_sex_pt, "BirthCat:sexF");                 print(head(pts_top_bc, 10))
  # Three-way PMA x birth group x sex (does the Fig. 3 slope pattern differ by sex?)  Exploratory: small cells.
  res_sex_pt3 <- sex_scan(pt_sex, FCz ~ scanner + ScanPMAc * BirthCat * sexF, "ScanPMAc:BirthCat:sexF", tag = "preterm3way")
  pts3_sum <- sex_term_summary(res_sex_pt3, "ScanPMAc:BirthCat:sexF", "PMA x birth group x sex"); print(pts3_sum, row.names = FALSE)
  pts3_top <- sex_top(res_sex_pt3, "ScanPMAc:BirthCat:sexF"); print(head(pts3_top, 10))
  # The 16 primary PMA x prematurity edges: do any of the sex terms reach nominal significance there?
  prim_sex <- data.frame(Connection = sigEdges,
                         `Sex, p` = signif(res_sex_pt$sexF_p[match(sigEdges, allConn)], 3),
                         `PMA x sex, p` = signif(res_sex_pt$ScanPMAc_sexF_p[match(sigEdges, allConn)], 3),
                         `Birth group x sex, p` = signif(res_sex_pt$BirthCat_sexF_p[match(sigEdges, allConn)], 3),
                         `PMA x birth group x sex, p` = signif(res_sex_pt3$ScanPMAc_BirthCat_sexF_p[match(sigEdges, allConn)], 3), check.names = FALSE)
  print(prim_sex, row.names = FALSE)
  # The primary scan adjusted for sex, and run within each sex
  res_pt_sexadj <- wald_scan(pt_sex, FCz ~ scanner + sexF + ScanPMA * BirthCat, "ScanPMA:BirthCat")
  res_pt_f      <- wald_scan(pt_f, fml_main, "ScanPMA:BirthCat")
  res_pt_m      <- wald_scan(pt_m, fml_main, "ScanPMA:BirthCat")
  cmp_pt_sexadj <- compare_sets(res_pt_sexadj); cmp_pt_f <- compare_sets(res_pt_f); cmp_pt_m <- compare_sets(res_pt_m)
  print(rbind(`Adjusted for sex` = cmp_pt_sexadj$summary, Female = cmp_pt_f$summary, Male = cmp_pt_m$summary))
  slopes_f <- tryCatch(slope_wide(slope_table(pt_f, sigEdges, fml_main)), error = function(e) { message("Female slopes skipped: ", e$message); NULL })
  slopes_m <- tryCatch(slope_wide(slope_table(pt_m, sigEdges, fml_main)), error = function(e) { message("Male slopes skipped: ", e$message); NULL })

  si_table(rbind(pts_sum, pts3_sum), "sex_preterm", "sexTests_summary")
  si_table(pts_top_main, "sex_preterm", "sexMainEffect_topEdges")
  si_table(pts_top_pma,  "sex_preterm", "pmaBySex_topEdges")
  si_table(pts_top_bc,   "sex_preterm", "birthGroupBySex_topEdges")
  si_table(pts3_top,     "sex_preterm", "threeWay_topEdges")
  si_table(prim_sex,     "sex_preterm", "sexTerms_on16PrimaryEdges")
  si_intx("sex_preterm", "pmaByPrematurity_bySex",
          list(`Adjusted for sex` = res_pt_sexadj, `Females only` = res_pt_f, `Males only` = res_pt_m),
          list(`Adjusted for sex` = cmp_pt_sexadj, `Females only` = cmp_pt_f, `Males only` = cmp_pt_m))
  if (!is.null(slopes_f)) si_table(slopes_f, "sex_preterm", "pmaByPrematurity_slopesFemales")
  if (!is.null(slopes_m)) si_table(slopes_m, "sex_preterm", "pmaByPrematurity_slopesMales")
  top_pt <- head(pts_top_pma$Connection, 6)
  si_fig("sex_preterm", "sexEffects", 8.2, 7.4, function()
    print(ggpubr::ggarrange(
      ggpubr::ggarrange(sex_volcano(res_sex_pt, "sexF", "sexFMale", "Male - female (beta, at mean PMA)"),
                        sex_volcano(res_sex_pt, "ScanPMAc:sexF", "ScanPMAc:sexFMale", "PMA x sex (beta per week)"),
                        nrow = 1, labels = c("A", "B"), common.legend = TRUE, legend = "bottom"),
      sex_edge_plot(pt_sex, top_pt), ncol = 1, heights = c(1, 1.7), labels = c("", "C"))))
  si_fig("sex_preterm", "PvalueQQ", 6.2, 5.0, function() print(sex_qq(res_sex_pt, pts_terms, pts_labels)))

  ## ---- 11d. Focused test on the 16 primary edges: does the Fig. 3 slope pattern differ by sex? ---------------------
  # One model per edge: FC ~ scanner + PMA * birth group * sex (random intercept per infant). From it: the PMA slope in each
  # birth group and sex, the male - female difference of each slope (3 contrasts), and the difference-in-differences that
  # tests the Fig. 3 pattern directly: (MPT slope - EPT slope) and (MPT slope - VPT slope) in males minus the same in
  # females (1 df each). FDR is applied across the 16 edges for each contrast family. This is a targeted follow-up on
  # edges chosen from the primary analysis, not a discovery scan, so report it as such.
  sex_focus <- function(cn) {
    d <- pt_sex; d$FCz <- as.numeric(goodFCData[d$longID, cn])
    m <- lme(FCz ~ scanner + ScanPMAc * BirthCat * sexF, random = ~ 1 | shortID, data = d, method = "REML", control = ctrl)
    b <- nlme::fixef(m); V <- vcov(m)
    P <- "ScanPMAc"; PV <- "ScanPMAc:BirthCatVPT"; PM <- "ScanPMAc:BirthCatMPT"; PS <- "ScanPMAc:sexFMale"
    PVS <- "ScanPMAc:BirthCatVPT:sexFMale"; PMS <- "ScanPMAc:BirthCatMPT:sexFMale"
    stopifnot(all(c(P, PV, PM, PS, PVS, PMS) %in% names(b)))
    w <- function(...) { a <- c(...); x <- setNames(numeric(length(b)), names(b)); x[names(a)] <- a; x }
    L <- list(
      slope_EPT_F = w(setNames(1, P)),                         slope_EPT_M = w(setNames(c(1, 1), c(P, PS))),
      slope_VPT_F = w(setNames(c(1, 1), c(P, PV))),            slope_VPT_M = w(setNames(c(1, 1, 1, 1), c(P, PV, PS, PVS))),
      slope_MPT_F = w(setNames(c(1, 1), c(P, PM))),            slope_MPT_M = w(setNames(c(1, 1, 1, 1), c(P, PM, PS, PMS))),
      diff_EPT    = w(setNames(1, PS)),                        # male - female slope, by birth group
      diff_VPT    = w(setNames(c(1, 1), c(PS, PVS))),
      diff_MPT    = w(setNames(c(1, 1), c(PS, PMS))),
      MPTvEPT_F   = w(setNames(1, PM)),                        MPTvEPT_M = w(setNames(c(1, 1), c(PM, PMS))),
      MPTvEPT_diff = w(setNames(1, PMS)),                      # focused test 1: (MPT - EPT) in males minus females
      MPTvVPT_F   = w(setNames(c(1, -1), c(PM, PV))),          MPTvVPT_M = w(setNames(c(1, -1, 1, -1), c(PM, PV, PMS, PVS))),
      MPTvVPT_diff = w(setNames(c(1, -1), c(PMS, PVS))))       # focused test 2: (MPT - VPT) in males minus females
    dfd <- min(m$fixDF$X[c(P, PV, PM, PS, PVS, PMS)]); tq <- qt(.975, dfd)
    do.call(rbind, lapply(names(L), function(k) {
      x <- L[[k]]; est <- sum(x * b); se <- sqrt(drop(t(x) %*% V %*% x))
      data.frame(Connection = cn, term = k, est = est, se = se, lo = est - tq * se, hi = est + tq * se, p = 2 * stats::pt(-abs(est / se), dfd))
    }))
  }
  focus_long <- do.call(rbind, lapply(sigEdges, function(cn) tryCatch(sex_focus(cn), error = function(e) { message(cn, ": ", e$message); NULL })))
  focus_long$pAdj <- ave(focus_long$p, focus_long$term, FUN = function(p) p.adjust(p, "fdr"))      # across the 16 edges, per contrast
  fterm <- function(k) focus_long[focus_long$term == k, ]
  cell  <- function(r) sprintf("%s [%s, %s]", fmt3(r$est), fmt3(r$lo), fmt3(r$hi))
  focus_wide <- data.frame(Connection = fterm("slope_MPT_F")$Connection,
                           `MPT slope, female` = cell(fterm("slope_MPT_F")), `MPT slope, male` = cell(fterm("slope_MPT_M")),
                           `MPT - EPT, female` = cell(fterm("MPTvEPT_F")), `MPT - EPT, male` = cell(fterm("MPTvEPT_M")),
                           `Male - female (MPT - EPT)` = cell(fterm("MPTvEPT_diff")),
                           `p` = signif(fterm("MPTvEPT_diff")$p, 3), `FDR p` = signif(fterm("MPTvEPT_diff")$pAdj, 3),
                           `MPT - VPT, female` = cell(fterm("MPTvVPT_F")), `MPT - VPT, male` = cell(fterm("MPTvVPT_M")),
                           `Male - female (MPT - VPT)` = cell(fterm("MPTvVPT_diff")),
                           `p ` = signif(fterm("MPTvVPT_diff")$p, 3), `FDR p ` = signif(fterm("MPTvVPT_diff")$pAdj, 3), check.names = FALSE)
  focus_sum <- do.call(rbind, lapply(c("MPTvEPT_diff", "MPTvVPT_diff", "diff_EPT", "diff_VPT", "diff_MPT"), function(k) {
    r <- fterm(k)
    data.frame(Contrast = c(MPTvEPT_diff = "(MPT - EPT slope): male - female", MPTvVPT_diff = "(MPT - VPT slope): male - female",
                            diff_EPT = "EPT slope: male - female", diff_VPT = "VPT slope: male - female", diff_MPT = "MPT slope: male - female")[[k]],
               Edges = nrow(r), `p < .05` = sum(r$p < .05), `FDR p < .05` = sum(r$pAdj < .05), `Smallest p` = signif(min(r$p), 3),
               `Smallest FDR p` = signif(min(r$pAdj), 3), `Median estimate (male - female)` = round(median(r$est), 4),
               `Edges with male > female` = sum(r$est > 0), check.names = FALSE)
  }))
  print(focus_sum, row.names = FALSE)
  print(focus_wide[, c(1, 2, 3, 6, 7, 8, 11, 12, 13)], row.names = FALSE)
  si_table(focus_sum,  "sex_preterm", "focused16_summary")
  si_table(focus_wide, "sex_preterm", "focused16_slopesAndContrasts")
  si_table(focus_long, "sex_preterm", "focused16_allEstimates")
  # Figure: A. the MPT - EPT slope contrast in each sex; B. its male - female difference (the focused test), per edge
  fl <- rbind(transform(fterm("MPTvEPT_F"), Sex = "Female"), transform(fterm("MPTvEPT_M"), Sex = "Male"))
  fl$Connection <- factor(fl$Connection, rev(sigEdges))
  fd <- fterm("MPTvEPT_diff"); fd$Connection <- factor(fd$Connection, rev(sigEdges))
  p_fa <- ggplot(fl, aes(est, Connection, color = Sex)) + geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6), size = 0.25, linewidth = 0.4) +
    scale_color_manual(values = sex_cols, name = NULL) + labs(x = "MPT - EPT slope (Fisher z per week)", y = NULL) + si_theme()
  p_fb <- ggplot(fd, aes(est, Connection)) + geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.25, linewidth = 0.4, color = "#444444") +
    labs(x = "Male - female difference in (MPT - EPT) slope", y = NULL) + si_theme() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
  si_fig("sex_preterm", "focused16_MPTvsEPT", 7.0, 4.4, function()
    print(ggpubr::ggarrange(p_fa, p_fb, nrow = 1, widths = c(1.55, 1), labels = c("A", "B"))))
}

###############################################################################
# 12. FOUR-GROUP MODEL: Fetal, EPT, VPT, MPT
#     FC ~ PMA at scan + scanner + group (4 levels, Fetal = reference), random intercept per infant, all 3,003 edges.
#     Reported per edge: the omnibus group effect (3-df Wald chi-square) and six pairwise contrasts (each preterm
#     group vs. fetal, and EPT / VPT / MPT against each other; Wald t on the model coefficients). FDR is applied across
#     the 3,003 edges separately for the omnibus test and for each contrast. A second model adds PMA x group (3 df) to
#     ask whether the PMA slope differs among the four groups.
#     CAUTION: group, PMA range and scanner overlap heavily (fetal scans are all 1.5T and span younger ages than most
#     preterm scans), so the fetal contrasts rest on the PMA and scanner adjustment, as in the two-group model of
#     section 4. The preterm-only contrasts are not affected by this.
#     Output number: "four_group" in SI_ID.
###############################################################################

if (RUN_FOUR_GROUP) {
  si_4 <- scanInfo
  si_4$grp <- factor(ifelse(si_4$scanType == "Fetal", "Fetal", as.character(si_4$BirthCat)), c("Fetal", "EPT", "VPT", "MPT"))
  stopifnot(!anyNA(si_4$grp))
  print(rbind(scans = table(si_4$grp), infants = tapply(as.character(si_4$shortID), si_4$grp, function(x) length(unique(x)))))

  four_group_edge <- function(cn, data) {
    cts <- c("EPT_vs_Fetal", "VPT_vs_Fetal", "MPT_vs_Fetal", "VPT_vs_EPT", "MPT_vs_EPT", "MPT_vs_VPT")
    nm  <- c(unlist(lapply(cts, function(k) paste0(k, c("_b", "_se", "_p")))), "omni_chi", "omni_p", "pmaXgrp_chi", "pmaXgrp_p")
    d <- data; d$FC <- as.numeric(goodFCData[d$longID, cn])
    tryCatch({
      m <- lme(FC ~ ScanPMA + scanner + grp, random = ~ 1 | shortID, data = d, control = ctrl)
      b <- nlme::fixef(m); V <- vcov(m); cf <- c("grpEPT", "grpVPT", "grpMPT"); stopifnot(all(cf %in% names(b)))
      w <- function(a) { x <- setNames(numeric(length(b)), names(b)); x[names(a)] <- a; x }
      Lc <- list(EPT_vs_Fetal = w(c(grpEPT = 1)), VPT_vs_Fetal = w(c(grpVPT = 1)), MPT_vs_Fetal = w(c(grpMPT = 1)),
                 VPT_vs_EPT = w(c(grpVPT = 1, grpEPT = -1)), MPT_vs_EPT = w(c(grpMPT = 1, grpEPT = -1)), MPT_vs_VPT = w(c(grpMPT = 1, grpVPT = -1)))
      dfd <- min(m$fixDF$X[cf])
      ct <- unlist(lapply(cts, function(k) { x <- Lc[[k]]; est <- sum(x * b); se <- sqrt(drop(t(x) %*% V %*% x)); c(est, se, 2 * pt(-abs(est / se), dfd)) }))
      wald3 <- function(b, V, coefs) { L <- matrix(0, 3, length(b), dimnames = list(NULL, names(b))); for (i in 1:3) L[i, coefs[i]] <- 1
        lb <- L %*% b; chi <- drop(t(lb) %*% solve(L %*% V %*% t(L)) %*% lb); c(chi, pchisq(chi, 3, lower.tail = FALSE)) }
      omni <- wald3(b, V, cf)
      m2 <- lme(FC ~ scanner + ScanPMA * grp, random = ~ 1 | shortID, data = d, control = ctrl)
      ix <- wald3(nlme::fixef(m2), vcov(m2), paste0("ScanPMA:grp", c("EPT", "VPT", "MPT")))
      setNames(c(ct, omni, ix), nm)
    }, error = function(e) setNames(rep(NA_real_, length(nm)), nm))
  }
  res4 <- cached("fourGroup", key = paste(data_key, "fourGroup", table_key(si_4)),
                 data.frame(Connection = allConn, do.call(rbind, par_lapply(allConn, four_group_edge, data = si_4)), row.names = NULL))
  cts4 <- c("EPT_vs_Fetal", "VPT_vs_Fetal", "MPT_vs_Fetal", "VPT_vs_EPT", "MPT_vs_EPT", "MPT_vs_VPT")
  for (k in c(cts4, "omni", "pmaXgrp")) res4[[paste0(k, "_pAdj")]] <- p.adjust(res4[[paste0(k, "_p")]], "fdr")

  lab4 <- c(omni = "Omnibus group effect (3 df)", EPT_vs_Fetal = "EPT vs. fetal", VPT_vs_Fetal = "VPT vs. fetal", MPT_vs_Fetal = "MPT vs. fetal",
            VPT_vs_EPT = "VPT vs. EPT", MPT_vs_EPT = "MPT vs. EPT", MPT_vs_VPT = "MPT vs. VPT", pmaXgrp = "PMA x group (3 df)")
  sum4 <- do.call(rbind, lapply(names(lab4), function(k) {
    p <- res4[[paste0(k, "_p")]]; pa <- res4[[paste0(k, "_pAdj")]]; sg <- !is.na(pa) & pa < .05
    bb <- if (k %in% cts4) res4[[paste0(k, "_b")]] else NULL
    data.frame(Test = lab4[[k]], `Raw p < .05` = sum(p < .05, na.rm = TRUE), `FDR p < .05` = sum(sg),
               `First group higher` = if (is.null(bb)) NA else sum(sg & bb > 0), `Second group higher` = if (is.null(bb)) NA else sum(sg & bb < 0),
               `Overlap with the 1,404 fetal-vs-preterm edges` = sum(sg & sigFetal),
               `Smallest FDR p` = signif(suppressWarnings(min(pa, na.rm = TRUE)), 3), check.names = FALSE)
  }))
  print(sum4, row.names = FALSE)   # contrasts are "A vs. B" = A - B: "First group higher" counts FDR-significant edges with A > B

  # Size of the difference from fetal in each preterm group, on the 1,404 primary fetal-vs-preterm edges
  # (graded by birth group? a larger median beta for a group = a larger departure from fetal)
  gr4 <- do.call(rbind, lapply(c("EPT_vs_Fetal", "VPT_vs_Fetal", "MPT_vs_Fetal"), function(k) {
    bb <- res4[[paste0(k, "_b")]][sigFetal]
    data.frame(Contrast = lab4[[k]], `Median beta, primary edges` = round(median(bb, na.rm = TRUE), 4),
               `Edges with same sign as primary` = sum(sign(bb) == sign(fetalRes$scanTypeb[sigFetal]), na.rm = TRUE),
               `Correlation of beta with primary beta` = round(cor(res4[[paste0(k, "_b")]], fetalRes$scanTypeb, use = "complete.obs"), 3), check.names = FALSE)
  }))
  print(gr4, row.names = FALSE)
  # Fig. 3 edges: the preterm group contrasts on the 16 primary PMA x prematurity edges (adjusted for PMA, not group x PMA)
  e16 <- res4[match(sigEdges, res4$Connection), c("Connection", paste0(rep(c("MPT_vs_EPT", "MPT_vs_VPT", "VPT_vs_EPT"), each = 2), c("_b", "_p")), "pmaXgrp_p")]
  e16[-1] <- lapply(e16[-1], signif, 3); print(e16, row.names = FALSE)

  si_table(sum4, "four_group", "fourGroup_summary")
  si_table(gr4,  "four_group", "fourGroup_departureFromFetal")
  si_table(e16,  "four_group", "fourGroup_on16PrimaryEdges")
  si_table(res4, "four_group", "fourGroup_allEdges")

  # Figure: A. FDR-significant edges per test, split by direction; B. each preterm group's beta vs. fetal against the primary beta
  d_cnt <- do.call(rbind, lapply(cts4, function(k) { sg <- !is.na(res4[[paste0(k, "_pAdj")]]) & res4[[paste0(k, "_pAdj")]] < .05; bb <- res4[[paste0(k, "_b")]]
    parts <- strsplit(lab4[[k]], " vs. ", fixed = TRUE)[[1]]
    data.frame(Contrast = lab4[[k]], Higher = c(parts[1], parts[2]), n = c(sum(sg & bb > 0), -sum(sg & bb < 0))) }))
  d_cnt$Contrast <- factor(d_cnt$Contrast, rev(lab4[cts4]))
  p4a <- ggplot(d_cnt, aes(n, Contrast, fill = n > 0)) + geom_col(width = 0.65) + geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_text(aes(label = abs(n), hjust = ifelse(n > 0, -0.2, 1.2)), size = 2.6) +
    scale_fill_manual(values = c(`TRUE` = col_preterm, `FALSE` = col_fetal), guide = "none") +
    scale_x_continuous(expand = expansion(mult = 0.18), labels = abs) +
    labs(x = "FDR-significant edges (left: second group higher; right: first group higher)", y = NULL) + si_theme()
  d_b4 <- do.call(rbind, lapply(c("EPT_vs_Fetal", "VPT_vs_Fetal", "MPT_vs_Fetal"), function(k)
    data.frame(Contrast = lab4[[k]], x = fetalRes$scanTypeb, y = res4[[paste0(k, "_b")]], sig = sigFetal)))
  d_b4$Contrast <- factor(d_b4$Contrast, lab4[c("EPT_vs_Fetal", "VPT_vs_Fetal", "MPT_vs_Fetal")])
  p4b <- ggplot(d_b4, aes(x, y)) + geom_abline(linetype = "dashed", linewidth = 0.3) +
    geom_point(data = d_b4[!d_b4$sig, ], color = "#8a8a8a", alpha = 0.2, size = 0.4, stroke = 0) +
    geom_point(data = d_b4[d_b4$sig, ], color = col_preterm, alpha = 0.35, size = 0.5, stroke = 0) +
    facet_wrap(~Contrast, nrow = 1) + labs(x = "Primary group difference (beta, preterm - fetal)", y = "Beta vs. fetal, by birth group") + si_theme()
  si_fig("four_group", "fourGroup", 8.2, 3.4, function()
    print(ggpubr::ggarrange(p4a, p4b, nrow = 1, widths = c(1, 1.5), labels = c("A", "B"))))
}

###############################################################################
# 13. TABLE 1 GROUP COMPARISONS
#     (A) Preterm (all) vs. fetal and (B) EPT vs. VPT vs. MPT, for the variables in Table 1.
#     Infant-level variables (scans per infant, birth PMA and weight, sex, APGAR, clinical items, ventilation days) are
#     tested with one value per infant. Scan-level variables (PMA at scan, volumes, minutes, mean FD) use a linear mixed
#     model with a random intercept per infant, so infants with several scans are not counted as independent.
#       continuous         Welch t-test (A) / Welch ANOVA with Games-Howell-style pairwise Welch t-tests, Holm (B)
#       skewed or ordinal  Wilcoxon rank-sum (A) / Kruskal-Wallis with pairwise Wilcoxon, Holm (B)   [scans per infant, APGAR, ventilation days]
#       yes/no             Fisher exact test; pairwise Fisher tests, Holm (B)
#       scan level         lme(): Wald t for the group difference (A) / F test, then pairwise contrasts, Holm (B)
#     FDR (Benjamini-Hochberg) is applied across the variables within each table. Clinical items and ventilation days
#     exist only for preterm infants, so they appear in (B) only. PMA at birth defines the preterm groups, so its (B)
#     test is shown for completeness and is not informative. All p-values are exploratory descriptors of the sample.
#     Output number: "table1_tests" in SI_ID.
###############################################################################

if (RUN_TABLE1_TESTS) {
  inf_t$cohort   <- factor(ifelse(inf_t$grp == "Fetal", "Fetal", "Preterm"), c("Fetal", "Preterm"))
  scan_df$cohort <- factor(ifelse(scan_df$grp == "Fetal", "Fetal", "Preterm"), c("Fetal", "Preterm"))
  inf_t$male_yes <- ifelse(is.na(inf_t$male), NA, inf_t$male == 1)
  yn_items <- setdiff(clin_items, c("Birth weight (g)", "APGAR 1 minute", "APGAR 5 minute", "Days on ventilation"))
  for (k in yn_items) inf_t[[paste0(k, "_yes")]] <- as_yes(inf_t[[k]])
  pre_inf  <- droplevels(inf_t[inf_t$grp != "Fetal", ]);   pre_scan <- droplevels(scan_df[scan_df$grp != "Fetal", ])
  pre_inf$grp <- factor(as.character(pre_inf$grp), c("EPT", "VPT", "MPT")); pre_scan$grp <- factor(as.character(pre_scan$grp), c("EPT", "VPT", "MPT"))
  ci_txt <- function(est, lo, hi, d = 2) sprintf(paste0("%.", d, "f [%.", d, "f, %.", d, "f]"), est, lo, hi)
  sg     <- function(p) signif(p, 3)

  # variable specs: label, data column, kind (cont / nonpar / yn / scan), digits for summaries, preterm-only flag
  specs <- list(
    list("Scans per infant",              "n_scans",            "nonpar", 1, FALSE),
    list("PMA at birth, wk (defines the preterm groups)", "birthPMA", "cont", 2, FALSE),
    list("Birth weight, g",               "Birth weight (g)",   "cont",   0, FALSE),
    list("Male sex",                      "male_yes",           "yn",     0, FALSE),
    list("APGAR, 1 minute",               "APGAR 1 minute",     "nonpar", 2, FALSE),
    list("APGAR, 5 minute",               "APGAR 5 minute",     "nonpar", 2, FALSE),
    lapply(yn_items, function(k) list(k, paste0(k, "_yes"), "yn", 0, TRUE)),
    list(list("Days on ventilation",      "Days on ventilation", "nonpar", 2, TRUE)),
    list(list("PMA at scan, wk",          "pma",                "scan",   2, FALSE)),
    list(list("Usable volumes per scan",  "vol",                "scan",   1, FALSE)),
    list(list("Minutes of data per scan", "mins",               "scan",   1, FALSE)),
    list(list("Mean framewise displacement, mm", "fd",          "scan",   3, FALSE)))
  specs <- c(specs[1:6], specs[[7]], specs[8:length(specs)])            # flatten the clinical items
  specs <- lapply(specs, function(x) if (is.list(x[[1]])) x[[1]] else x)

  summ <- function(x, kind, d) if (kind == "yn") npc(x) else msd(x, d)
  adj  <- function(df) { df$`FDR p` <- NA_real_; i <- !is.na(df$p); df$`FDR p`[i] <- sg(p.adjust(df$p[i], "fdr")); df }

  ## ---- (A) preterm (all) vs. fetal ---------------------------------------------------------------------------
  rowA <- function(sp) {
    lab <- sp[[1]]; col <- sp[[2]]; kind <- sp[[3]]; dg <- sp[[4]]
    base <- data.frame(Variable = lab, Level = if (kind == "scan") "Scan (mixed model)" else "Infant", check.names = FALSE, stringsAsFactors = FALSE)
    if (sp[[5]]) return(cbind(base, Test = "preterm only", Fetal = "-", Preterm = "-", `Preterm - fetal [95% CI]` = "-", Statistic = "-", p = NA_real_))
    d <- if (kind == "scan") scan_df else inf_t; x <- d[[col]]; g <- d$cohort; ok <- !is.na(x); x <- x[ok]; g <- droplevels(g[ok]); d <- d[ok, ]
    cells <- c(Fetal = summ(x[g == "Fetal"], kind, dg), Preterm = summ(x[g == "Preterm"], kind, dg))
    out <- tryCatch(switch(kind,
      cont = { tt <- t.test(x ~ g); list("Welch t-test", ci_txt(mean(x[g == "Preterm"]) - mean(x[g == "Fetal"]), -tt$conf.int[2], -tt$conf.int[1], dg),
                                         sprintf("t(%.1f) = %.2f", tt$parameter, -tt$statistic), tt$p.value) },
      nonpar = { wt <- suppressWarnings(wilcox.test(x ~ g, conf.int = TRUE)); list("Wilcoxon rank-sum", ci_txt(-wt$estimate, -wt$conf.int[2], -wt$conf.int[1], dg + 1),
                                         sprintf("W = %.0f", wt$statistic), wt$p.value) },
      yn = { tb <- table(g, factor(x, c(FALSE, TRUE))); ft <- fisher.test(tb); list("Fisher exact", sprintf("OR %.2f [%.2f, %.2f]", ft$estimate, ft$conf.int[1], ft$conf.int[2]), "-", ft$p.value) },
      scan = { m <- lme(x ~ cohort, random = ~ 1 | shortID, data = transform(d, x = x), control = ctrl); tt <- summary(m)$tTable["cohortPreterm", ]
               q <- qt(.975, tt["DF"]); list("Mixed model (random intercept per infant)", ci_txt(tt["Value"], tt["Value"] - q * tt["Std.Error"], tt["Value"] + q * tt["Std.Error"], dg + 1),
                                             sprintf("t(%d) = %.2f", tt["DF"], tt["t-value"]), tt["p-value"]) }),
      error = function(e) list("failed", "-", "-", NA_real_))
    cbind(base, Test = out[[1]], Fetal = cells[["Fetal"]], Preterm = cells[["Preterm"]], `Preterm - fetal [95% CI]` = out[[2]], Statistic = out[[3]], p = sg(unname(out[[4]])))
  }
  tabA <- adj(do.call(rbind, lapply(specs, rowA)))
  print(tabA, row.names = FALSE)

  ## ---- (B) EPT vs. VPT vs. MPT ---------------------------------------------------------------------------------
  pair_fisher <- function(x, g) { pr <- list(c("EPT", "VPT"), c("EPT", "MPT"), c("VPT", "MPT"))
    p <- sapply(pr, function(z) { i <- g %in% z; fisher.test(table(droplevels(g[i]), factor(x[i], c(FALSE, TRUE))))$p.value }); p.adjust(p, "holm") }
  rowB <- function(sp) {
    lab <- sp[[1]]; col <- sp[[2]]; kind <- sp[[3]]; dg <- sp[[4]]
    d <- if (kind == "scan") pre_scan else pre_inf; x <- d[[col]]; g <- d$grp; ok <- !is.na(x); x <- x[ok]; g <- droplevels(g[ok]); d <- d[ok, ]
    base <- data.frame(Variable = lab, Level = if (kind == "scan") "Scan (mixed model)" else "Infant", check.names = FALSE, stringsAsFactors = FALSE)
    cells <- sapply(levels(pre_inf$grp), function(l) summ(x[g == l], kind, dg))
    out <- tryCatch(switch(kind,
      cont = { ao <- oneway.test(x ~ g); pm <- pairwise.t.test(x, g, pool.sd = FALSE, p.adjust.method = "holm")$p.value
               list("Welch ANOVA", sprintf("F(%.0f, %.1f) = %.2f", ao$parameter[1], ao$parameter[2], ao$statistic), ao$p.value, c(pm["VPT", "EPT"], pm["MPT", "EPT"], pm["MPT", "VPT"])) },
      nonpar = { kw <- kruskal.test(x ~ g); pm <- suppressWarnings(pairwise.wilcox.test(x, g, p.adjust.method = "holm", exact = FALSE))$p.value
               list("Kruskal-Wallis", sprintf("H(2) = %.2f", kw$statistic), kw$p.value, c(pm["VPT", "EPT"], pm["MPT", "EPT"], pm["MPT", "VPT"])) },
      yn = { ft <- fisher.test(table(g, factor(x, c(FALSE, TRUE)))); list("Fisher exact", "-", ft$p.value, pair_fisher(x, g)) },
      scan = { m <- lme(x ~ grp, random = ~ 1 | shortID, data = transform(d, x = x), control = ctrl); an <- anova(m)["grp", ]
               b <- nlme::fixef(m); V <- vcov(m); df2 <- m$fixDF$X[c("grpVPT", "grpMPT")]
               cp <- function(w, dfx) { est <- sum(w * b); se <- sqrt(drop(t(w) %*% V %*% w)); 2 * pt(-abs(est / se), dfx) }
               w0 <- setNames(numeric(length(b)), names(b)); wv <- w0; wv["grpVPT"] <- 1; wm <- w0; wm["grpMPT"] <- 1; wvm <- w0; wvm["grpMPT"] <- 1; wvm["grpVPT"] <- -1
               pp <- p.adjust(c(cp(wv, df2[1]), cp(wm, df2[2]), cp(wvm, min(df2))), "holm")
               list("Mixed model (random intercept per infant)", sprintf("F(%d, %d) = %.2f", an[["numDF"]], an[["denDF"]], an[["F-value"]]), an[["p-value"]], pp) }),
      error = function(e) list("failed", "-", NA_real_, rep(NA_real_, 3)))
    data.frame(base, Test = out[[1]], EPT = cells[["EPT"]], VPT = cells[["VPT"]], MPT = cells[["MPT"]], Statistic = out[[2]], p = sg(unname(out[[3]])),
               `EPT vs VPT (Holm p)` = sg(unname(out[[4]][1])), `EPT vs MPT (Holm p)` = sg(unname(out[[4]][2])), `VPT vs MPT (Holm p)` = sg(unname(out[[4]][3])),
               check.names = FALSE, stringsAsFactors = FALSE)
  }
  tabB <- adj(do.call(rbind, lapply(specs, rowB)))
  print(tabB, row.names = FALSE)

  si_table(tabA, "table1_tests", "PretermVsFetal")
  si_table(tabB, "table1_tests", "BirthGroups")
}

sessionInfo()
