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
#   4  Fetal vs. preterm, all edges (Fig. 2: matrices, chords, homotopic bars; Table 2; SI Fig. 2 mean matrices)
#   5  Preterm only: PMA x prematurity, all edges (Fig. 3, Table 3)
#   6  Fig. 3 Panel B (group trajectories)
#   7  Fig. 3 Panel A (BrainNet Viewer edge file; 7b: the same picture in R with rgl)
#   8  Group slopes and pairwise contrasts for the significant edges
#   9  Sensitivity analyses and controls
#  10  Save tables, session info
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
# install.packages(c("ggplot2","nlme","lmerTest","car","emmeans","DescTools",
#                    "circlize","ggpubr","scales","ragg","svglite"))
# devtools::install_github("jcbeer/longCombat")
library(ggplot2)
library(nlme)       # lme(): linear mixed models with a random intercept
library(lmerTest)   # lmer() with p-values (used for Table 3, with car)
library(car)        # Anova(): Type II Wald chi-square tests for lmer models
library(emmeans)    # emtrends() for group slopes (section 8)
library(DescTools)  # FisherZ / FisherZInv
library(longCombat) # longitudinal ComBat
library(circlize)   # chord diagrams (Fig. 2)
library(ggpubr)     # ggarrange() for the multi-panel Fig. 3B

## ---- Switches --------------------------------------------------------------
SAVE_FIGS                 <- FALSE  # TRUE: write figure files to figureDir
FIG_FORMATS               <- c("tiff", "svg", "pdf")   # formats written when SAVE_FIGS is TRUE: any of "tiff","png","svg","pdf"
SAVE_TABLES               <- FALSE  # TRUE: write .csv result tables to figureDir
RUN_QA_PLOTS              <- TRUE   # scan-QA figures (section 1)
RUN_HARMONIZATION_CHECKS  <- TRUE   # scanner effect / PCA / PMA preservation (section 2)
RUN_SENSITIVITY           <- TRUE   # section 9: five lmer scans over all edges (slow serially; fast with N_CORES, cached)
RUN_EXTRA_CONTROLS        <- TRUE   # sections 9f-9h: age-matched fetal vs. preterm, postnatal-age check, unharmonized data
RUN_DISTANCE_CONTROL      <- TRUE   # section 9i: do fetal-vs-preterm differences depend on the distance between ROIs?
RUN_QC_COVARIATE_CONTROL  <- TRUE   # section 9e: finalVol covariate and similar-volume reruns (add tSNR / motion via qc_more)
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
#   "new"      : brighter blue (#2a78d6; passes the color-blind/contrast checks), and the
#                preterm heatmap and the Fig. 3B subgroup ramp are built from that same blue
#   "cyan"     : like "new", but the lightest Fig. 3B group is a cyan (#12b5cf) instead of a pale blue
#                (the pale blue reads gray and has poor contrast on white)
#   "periwinkle": navy -> periwinkle (#4a63d9) -> cyan; the most separated of the blue options for the three groups
#   "original" : earlier palette (navy #08415c, indigo heatmap)  -> set this to revert
COLOR_SCHEME <- "cyan"
pal <- switch(COLOR_SCHEME,
  new = list(fetal = "#cc2936", preterm = "#2a78d6",
             cols_pt      = c(EPT = "#0d366b", VPT = "#2a78d6", MPT = "#86b6ef"),   # dark -> light blue
             heat_preterm = c("#cde2fb", "#2a78d6", "#0d366b"),                     # low, mid, high
             heat_fetal   = c("#f4b1ac", "#cc2936", "#650f16")),
  cyan = list(fetal = "#cc2936", preterm = "#2a78d6",
             cols_pt      = c(EPT = "#0d366b", VPT = "#2a78d6", MPT = "#12b5cf"),
             heat_preterm = c("#cde2fb", "#2a78d6", "#0d366b"),
             heat_fetal   = c("#f4b1ac", "#cc2936", "#650f16")),
  periwinkle = list(fetal = "#cc2936", preterm = "#2a78d6",
             cols_pt      = c(EPT = "#14306b", VPT = "#4a63d9", MPT = "#12b5cf"),
             heat_preterm = c("#cde2fb", "#2a78d6", "#0d366b"),
             heat_fetal   = c("#f4b1ac", "#cc2936", "#650f16")),
  original = list(fetal = "#cc2936", preterm = "#08415c",
                  cols_pt      = c(EPT = "#08415c", VPT = "#2f7fb0", MPT = "#8cc2e0"),
                  heat_preterm = c("#adb3cf", "#373d7e", "#15085c"),
                  heat_fetal   = c("#f4b1ac", "#cc2936", "#650f16")))
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
#   pdf      : vector via cairo_pdf.
# install.packages(c("ragg", "svglite")) once for the best output.
show_fig <- function(file, w, h, draw, res = 600) {
  if (!SAVE_FIGS) { draw(); return(invisible(NULL)) }
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

# p-value formatting for tables: "<.001" or ".034"
fmt_p <- function(p) ifelse(p < .001, "<.001", sub("^0", "", sprintf("%.3f", p)))
# 3-decimal formatting that never prints a negative zero ("-0.000" -> "0.000")
fmt3 <- function(x) sub("^-(0\\.000)$", "\\1", sprintf("%.3f", x))

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
        pct_usable_lt37  = mean(scanQA$finalVol[scanQA$scanPMA <  37] > 119) * 100,
        pct_usable_gt37  = mean(scanQA$finalVol[scanQA$scanPMA >  37] > 119) * 100), 1)
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

  # (b) usable scans only, by scanner
  p_qa_good <- ggplot(subset(scanQA, good4Min == 1), aes(x = shortID, y = scanPMA)) +
    scale_y_continuous(limits = c(24, 46), breaks = seq(20, 60, 2), expand = c(0, 0)) +
    scale_shape_manual(values = c(21, 23)) +
    scale_fill_manual(values = c("#0096FF", "#90EE90")) +
    scale_color_manual(values = c("#00008B", "#095c00")) +
    geom_line(aes(group = shortID), linewidth = 1) +
    geom_point(aes(shape = as.factor(scanner), color = as.factor(scanner), fill = as.factor(scanner)), size = 3, stroke = 2) +
    xlab("Participants") + ylab("PMA at Scan") +
    theme_paper() + theme(axis.text.x = element_blank())

  # (c) salvage histogram: PMA of usable vs. excluded scans
  p_qa_hist <- ggplot(scanQA, aes(x = scanPMA)) +
    scale_x_continuous(limits = c(24, 46), breaks = seq(20, 60, 2), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 30), breaks = seq(0, 40, 2), expand = c(0, 0)) +
    geom_histogram(data = subset(scanQA, good4Min == 1), color = "#095c00", fill = "#7ef472aa", bins = 30) +
    geom_histogram(data = subset(scanQA, good4Min == 0), color = "#911717", fill = "#c75b5baa", bins = 30) +
    xlab("PMA at Scan") + theme_paper() + theme(axis.title.y = element_blank())

  # (d) volumes per scan with the 120-volume cutoff
  p_qa_vol <- ggplot(scanQA, aes(x = finalVol)) +
    geom_vline(xintercept = 120, linetype = "dashed") +
    geom_histogram(color = "#000000", fill = "#000055aa", bins = 30) +
    xlab("Final Volumes") + theme_paper() + theme(axis.title.y = element_blank())

  show_fig("qaPlot_goodOnly.tif", 18, 8, function() print(p_qa_good))
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
# Same style as qaPlot_goodOnly (section 1); both panels share the PMA axis so the groups are comparable.
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
  show_fig("qaPlot_goodOnly_stacked", 18, 12, function() print(p_qa_stack))
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
  show_fig("harmonizationScannerDif.tif", 6, 6, function() print(p_scan_d))

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
  show_fig("harmonizationPreservePMA.tif", 6, 6, function() print(p_bio))
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
                               Class1 = roiData$Class[i1], Class2 = roiData$Class[i2],
                               Module1 = roiData$Module[i1], Module2 = roiData$Module[i2],
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

# Table 1 data-amount rows (scan level): volumes retained and minutes of data per scan, mean (SD).
# Minutes = volumes x TR / 60 (TR = 2 s preterm, 3 s fetal). Columns: fetal, all preterm, then preterm groups.
TR_SEC  <- c(Fetal = 3, Preterm = 2)
scan_min <- scanInfo$finalVol * TR_SEC[as.character(scanInfo$scanType)] / 60   # kept out of scanInfo so cache keys do not change
scan_cols <- list(Fetal = scanInfo$scanType == "Fetal", `All preterm` = scanInfo$scanType == "Preterm",
                  EPT = scanInfo$BirthCat %in% "EPT", VPT = scanInfo$BirthCat %in% "VPT", MPT = scanInfo$BirthCat %in% "MPT")
table1_scan <- rbind(
  `Scans, n`                          = sapply(scan_cols, sum),
  `Volumes per scan, mean (SD)`       = sapply(scan_cols, function(i) ms(scanInfo$finalVol[i])),
  `Minutes of data per scan, mean (SD)` = sapply(scan_cols, function(i) sprintf("%.1f (%.1f)", mean(scan_min[i]), sd(scan_min[i]))))
table1_scan

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
show_fig("Fig2_matrix_fetal",   6, 6, function() print(p_tri_fetal))
show_fig("Fig2_matrix_preterm", 6, 6, function() print(p_tri_preterm))

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
show_fig("Fig2_chord_fetal",   9, 9, function() draw_chord(chord_col_fetal))
show_fig("Fig2_chord_preterm", 9, 9, function() draw_chord(chord_col_preterm))

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
# One plot: rows = bilateral ROIs in lobe order (Precentral at the top), Fetal above Preterm.
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
homoH$Group <- factor(homoH$Group, levels = c("Preterm", "Fetal"))      # dodge order: Fetal ends up on top

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
  guides(color = guide_legend(reverse = TRUE), fill = guide_legend(reverse = TRUE)) +   # legend reads Fetal, Preterm
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
show_fig("Fig2_homotopic_horizontal", 6.2, 10, function() print(p_homo_h))

## ---- 4f. SI Fig. 2: mean connectivity matrices, one per group (Fetal, EPT, VPT, MPT) ----
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
  show_fig("SIFig2_connectivity_matrices", 3 * si2_cell_w_in, 2 * si2_cell_h_in, draw_si2)   # the only file written for SI Fig. 2
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
#    extrapolation). Layout: panel_ncol columns with one shared legend below.
###############################################################################

# Appearance (tune here)
show_ci         <- FALSE   # TRUE: 95% bands
y_span          <- 2.4     # same vertical range (Fisher z) for every panel, centered on each panel's data
pt_size         <- 0.75    # point size
pt_alpha        <- 0.45    # point opacity
line_w          <- 0.9     # colored line width
halo_w          <- 0.2     # extra width of the white halo under each line
panel_ncol      <- 4       # panels per row; 16 edges fit a 4x4 grid exactly (legend goes below the grid)
# Size: design the figure at the size it will be printed, so the text sizes below are the printed sizes.
panelB_w_in      <- 7.5    # export width (inches); a bit wider than tall leaves room for the brain image above
panelB_row_h_in  <- 1.3    # height of each row of panels (inches); the plots are square, so this only needs to fit the square plus its x-axis label (+ 0.55 for the legend when it sits below the grid)
axis_text_size   <- 5.5    # tick labels (points)
axis_title_size  <- 6      # axis titles (points), both axes
legend_text_size <- 7.5    # legend (points)
panelB_labels    <- "short"   # y-axis names: "short" = display names with hemisphere, e.g. "Superior Orbital (L)"; "long" = atlas long names

panel_data <- function(cn) {
  d <- pt; d$FCz <- as.numeric(goodFCData[d$longID, cn])
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

make_panel_b <- function(edges, ncol = panel_ncol, legend = TRUE) {
  pd        <- lapply(edges, panel_data)
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
      scale_x_continuous(breaks = seq(25, 45, 5)) +
      scale_y_continuous(breaks = scales::breaks_width(0.5)) +
      coord_cartesian(xlim = c(24, 49), ylim = ylim) +
      labs(x = "PMA at scan (weeks)", y = ylab_for(cn), color = NULL, fill = NULL) +
      guides(color = guide_legend(override.aes = list(alpha = 1, size = 2, linewidth = 1.2)), fill = "none") +
      theme_classic(base_size = 8) +
      theme(aspect.ratio = 1,                                            # square plotting area
            axis.line = element_line(color = "#000000", linewidth = .5),
            axis.ticks = element_line(color = "#000000", linewidth = .4),
            axis.text = element_text(size = axis_text_size, color = "#000000"),
            axis.title.x = element_text(size = axis_title_size, color = "#000000"),
            axis.title.y = element_text(size = axis_title_size, color = "#000000"),
            plot.margin = margin(4, 8, 2, 2),
            panel.grid.major.y = element_line(color = "#e6e6e6", linewidth = .3),
            legend.position = "bottom", legend.text = element_text(size = legend_text_size))
  }
  plots <- lapply(seq_along(edges), function(j) edge_plot(edges[j], lab_order[j]))
  nrow_fig <- ceiling(length(plots) / ncol)
  # One shared legend below the grid (taken from the first panel); legend = FALSE leaves it out
  list(fig = ggpubr::ggarrange(plotlist = plots, ncol = ncol, nrow = nrow_fig, align = "hv",
                               common.legend = legend, legend = if (legend) "bottom" else "none"),
       nrow = nrow_fig)
}

# Legend-only plot for the empty cells beside Panel C
pt_legend_plot <- function() {
  d <- data.frame(g = factor(rep(names(cols_pt), each = 2), levels = names(cols_pt)), x = rep(1:2, 3), y = 1)
  p <- ggplot(d, aes(x, y, color = g, group = g)) + geom_line(linewidth = line_w) +
    scale_color_manual(values = cols_pt, labels = grp_labels, name = NULL) +
    guides(color = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(linewidth = 1.4))) +
    theme_void() + theme(legend.position = "left", legend.text = element_text(size = legend_text_size, color = "#000000"),
                         legend.key.width = unit(0.45, "in"), legend.key.height = unit(0.24, "in"), legend.margin = margin(0, 0, 0, 0))
  ggpubr::as_ggplot(ggpubr::get_legend(p))
}
# draws the legend in the right-hand part of a Panel C area of the given size (columns beyond the three plots)
draw_pt_legend <- function(w_in, h_in, x0_in = 0, y0_in = 0, first_col = 3, ncols = panel_ncol_split) {
  print(pt_legend_plot(), vp = grid::viewport(x = grid::unit(x0_in + w_in * first_col / ncols, "in"), y = grid::unit(y0_in, "in"),
                                              width = grid::unit(w_in * (ncols - first_col) / ncols, "in"), height = grid::unit(h_in, "in"), just = c("left", "bottom")))
}

# Split into B and C? TRUE: Panel B = the connections where MPT differs from VPT and EPT (which look alike), Panel C = the
# three oddballs (listed below), both panel_ncol_split columns wide; the legend sits in the empty cells to the right of C. FALSE: one Panel B, as before.
panelB_split     <- TRUE
panelC_edges     <- c("PreCG.L_ORBmid.R", "SPG.L_HES.R", "SMG.L_SMG.R")
panel_ncol_split <- 5
if (panelB_split) {
  stopifnot(all(panelC_edges %in% sigEdges))
  edgesB <- setdiff(sigEdges, panelC_edges)                                      # keeps the order set in edge_order
  panelB <- make_panel_b(edgesB,       ncol = panel_ncol_split, legend = FALSE)
  panelC <- make_panel_b(panelC_edges, ncol = panel_ncol_split, legend = FALSE)
  show_fig("Figure3_PanelB", panelB_w_in, panelB_row_h_in * panelB$nrow, function() print(panelB$fig))                 # width x height in inches
  show_fig("Figure3_PanelC", panelB_w_in, panelB_row_h_in * panelC$nrow,
           function() { print(panelC$fig); draw_pt_legend(panelB_w_in, panelB_row_h_in * panelC$nrow) })
} else {
  panelB <- make_panel_b(sigEdges)
  show_fig("Figure3_PanelB", panelB_w_in, panelB_row_h_in * panelB$nrow + 0.55, function() print(panelB$fig))           # width x height in inches
}


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
  edgeMat <- write_brainnet_edges(sigEdges, nodeFile, file.path(figureDir, "bainNetPlotting", "brainDiff_Update.edge"))
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
# The bounding boxes are printed so you can confirm the mesh and nodes overlap.
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
bn_inset_zoom     <- NULL             # NULL = estimated from the box size; set a number (smaller = more magnified) to tune by hand
bn_box_col        <- "#000000"        # color of the box and the magnification frame
bn_layout         <- "row"            # "row": all images side by side (a wide, short figure that sits above the interaction plots); "grid": two rows
bn_panel_h_in     <- 2.3              # height of every image in the combined figure, inches
bn_pad_x_in       <- 0.2              # white space between images, inches (room for ROI labels; raise it if you need more)
bn_pad_y_in       <- 0.1              # white space above and below the images, inches
bn_save_strip     <- FALSE            # TRUE: also save the image strip (Panel A) on its own
bn_gap_ab_in      <- 0.1              # space between Panel A and Panel B in the combined figure, inches
bn_letter_pt      <- 14               # panel letters A and B (points)
bn_panel_left_in  <- 0.22             # Panels B and C are indented by this much (inches) so the letters do not sit on the y-axis labels
bn_letter_in      <- 0.03             # distance of the letters from the top-left corner of each panel, inches
bn_save_views     <- FALSE            # TRUE: also write each view as its own PNG (needs SAVE_FIGS)

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
  print(rbind(mesh = apply(xyz, 2, range)[2, ] - apply(xyz, 2, range)[1, ],                 # extent (mm) in x, y, z
              nodes = apply(bn_nodes[, c("x", "y", "z")], 2, range)[2, ] - apply(bn_nodes[, c("x", "y", "z")], 2, range)[1, ]))
  print(rbind(mesh_min = apply(xyz, 2, min), mesh_max = apply(xyz, 2, max),
              node_min = apply(bn_nodes[, c("x", "y", "z")], 2, min), node_max = apply(bn_nodes[, c("x", "y", "z")], 2, max)))
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
  f_png <- function(nm) if (bn_save_views && SAVE_FIGS) file.path(figureDir, paste0("Fig3_PanelA_", nm, ".png")) else file.path(tempdir(), paste0("bn_", nm, ".png"))
  rgl::open3d(windowRect = c(20, 20, 20 + bn_png_px[1], 20 + bn_png_px[2]))      # the saved PNG has the size of this window
  for (v in bn_panels) snap(v, f_png(v), bn_zoom[[v]], box = bn_inset && v == bn_inset_view)
  if (bn_inset) {
    # Calibrate the magnification from the dorsal picture itself: pixels per mm there (brain width in pixels / width in mm),
    # then choose the zoom at which the rectangle just fits the window. Zoom is inversely proportional to scale.
    im_d  <- png::readPNG(f_png(bn_inset_view))[, , 1:3, drop = FALSE]
    nw_d  <- im_d[, , 1] < .985 | im_d[, , 2] < .985 | im_d[, , 3] < .985
    px_mm <- diff(range(which(colSums(nw_d) > 0))) / diff(range(xyz[, 1]))
    px_mm_t <- min(dim(im_d)[2] / (2 * hx), dim(im_d)[1] / (2 * hy))                 # target pixels per mm in the magnified picture
    inset_zoom <- if (is.null(bn_inset_zoom)) bn_zoom[[bn_inset_view]] * px_mm / px_mm_t else bn_inset_zoom
    message("Magnification zoom: ", round(inset_zoom, 3), " (set bn_inset_zoom to override)")
    snap(bn_inset_view, f_png("inset"), inset_zoom, center_on = to_rgl(matrix(box_ctr, 1))[1, ])
  }

  # ---- combine into one figure: the views in bn_panels, then the magnification (bn_layout = "row"; "grid" gives 2 rows) ----
  # Each image is trimmed to the brain and surrounded by bn_pad_x_in and bn_pad_y_in of white space for labeling the ROIs by hand.
  # The magnification is calibrated from the dorsal image so that the box in the dorsal view and the magnified panel cover the same region.
  if (!requireNamespace("png", quietly = TRUE)) {
    message("Package png not installed (install.packages('png')); the views were saved separately instead of combined.")
    for (nm in c(bn_panels, if (bn_inset) "inset")) file.copy(f_png(nm), file.path(figureDir, paste0("Fig3_PanelA_", nm, ".png")), overwrite = TRUE)
  } else {
    rd <- function(f) png::readPNG(f)[, , 1:3, drop = FALSE]
    trim <- function(im, m = 6) { nw <- im[, , 1] < .99 | im[, , 2] < .99 | im[, , 3] < .99          # trim to the drawing, keep m px of white (the faint outer surface must not be clipped)
                           r <- range(which(rowSums(nw) > 0)); cc <- range(which(colSums(nw) > 0))
                           r <- c(max(1, r[1] - m), min(dim(im)[1], r[2] + m)); cc <- c(max(1, cc[1] - m), min(dim(im)[2], cc[2] + m))
                           im[r[1]:r[2], cc[1]:cc[2], , drop = FALSE] }
    centre_rect <- function(im, w, h) { h <- min(h, dim(im)[1]); w <- min(w, dim(im)[2])                    # w x h pixels from the middle of the window
                                        r0 <- floor((dim(im)[1] - h) / 2); c0 <- floor((dim(im)[2] - w) / 2)
                                        im[r0 + seq_len(h), c0 + seq_len(w), , drop = FALSE] }
    ims   <- setNames(lapply(bn_panels, function(v) trim(rd(f_png(v)))), bn_panels)
    if (bn_inset) ims$inset <- centre_rect(rd(f_png("inset")), round(2 * hx * px_mm_t), round(2 * hy * px_mm_t))
    H     <- bn_panel_h_in; padx <- bn_pad_x_in; pady <- bn_pad_y_in
    hgt   <- setNames(rep(H, length(ims)), names(ims)); if (bn_inset) hgt[["inset"]] <- H * bn_inset_scale
    wid   <- sapply(names(ims), function(nm) hgt[[nm]] * dim(ims[[nm]])[2] / dim(ims[[nm]])[1])
    rows  <- if (bn_layout == "row") list(c(bn_panels, if (bn_inset) "inset"))
             else list(bn_panels[-length(bn_panels)], c(bn_panels[length(bn_panels)], if (bn_inset) "inset"))
    row_w <- function(r) sum(wid[r]) + (length(r) + 1) * padx
    W     <- max(sapply(rows, row_w)); Hfig <- length(rows) * H + (length(rows) + 1) * pady
    # Draws the image strip (Panel A) in the current viewport; sc scales every size (1 = natural size in inches).
    draw_combined <- function(sc = 1) {
      grid::grid.newpage()
      for (k in seq_along(rows)) {
        r <- rows[[k]]; x <- (W - row_w(r)) / 2 + padx; y0 <- Hfig - k * (H + pady)               # row k from the top
        for (nm in r) {
          yy <- y0 + (H - hgt[[nm]]) / 2                                                          # centred vertically in the row
          grid::grid.raster(ims[[nm]], x = grid::unit(x * sc, "in"), y = grid::unit(yy * sc, "in"), width = grid::unit(wid[[nm]] * sc, "in"),
                            height = grid::unit(hgt[[nm]] * sc, "in"), just = c("left", "bottom"), interpolate = TRUE)
          if (nm == "inset") grid::grid.rect(x = grid::unit(x * sc, "in"), y = grid::unit(yy * sc, "in"), width = grid::unit(wid[[nm]] * sc, "in"),
                                             height = grid::unit(hgt[[nm]] * sc, "in"), just = c("left", "bottom"), gp = grid::gpar(col = bn_box_col, fill = NA, lwd = 1.5))
          x <- x + wid[[nm]] + padx
        }
      }
    }
    if (bn_save_strip) show_fig("Fig3_PanelA_brain", W, Hfig, draw_combined)

    # Figure 3, Panels A + B: the image strip on top (scaled to the width of Panel B), the interaction plots below.
    # Panel B is drawn at its own design size (section 6), so its text sizes are unchanged.
    sc_a    <- panelB_w_in / W
    a_h_in  <- Hfig * sc_a
    b_h_in  <- if (panelB_split) panelB_row_h_in * panelB$nrow else panelB_row_h_in * panelB$nrow + 0.55   # Panel B (without legend when split)
    c_h_in  <- if (panelB_split) panelB_row_h_in * panelC$nrow else 0                                        # Panel C (its legend sits beside the plots)
    fig3_h  <- a_h_in + bn_gap_ab_in + b_h_in + c_h_in
    draw_fig3 <- function() {
      grid::grid.newpage()
      grid::pushViewport(grid::viewport(x = 0, y = grid::unit(fig3_h, "in"), width = grid::unit(panelB_w_in, "in"), height = grid::unit(a_h_in, "in"), just = c("left", "top")))
      draw_strip <- function() {                                                                  # Panel A (same drawing, scaled)
        for (k in seq_along(rows)) {
          r <- rows[[k]]; x <- (W - row_w(r)) / 2 + padx; y0 <- Hfig - k * (H + pady)
          for (nm in r) {
            yy <- y0 + (H - hgt[[nm]]) / 2
            grid::grid.raster(ims[[nm]], x = grid::unit(x * sc_a, "in"), y = grid::unit(yy * sc_a, "in"), width = grid::unit(wid[[nm]] * sc_a, "in"),
                              height = grid::unit(hgt[[nm]] * sc_a, "in"), just = c("left", "bottom"), interpolate = TRUE)
            if (nm == "inset") grid::grid.rect(x = grid::unit(x * sc_a, "in"), y = grid::unit(yy * sc_a, "in"), width = grid::unit(wid[[nm]] * sc_a, "in"),
                                               height = grid::unit(hgt[[nm]] * sc_a, "in"), just = c("left", "bottom"), gp = grid::gpar(col = bn_box_col, fill = NA, lwd = 1))
            x <- x + wid[[nm]] + padx
          }
        }
      }
      draw_strip()
      grid::popViewport()
      letter <- function(txt, y_top_in) grid::grid.text(txt, x = grid::unit(bn_letter_in, "in"), y = grid::unit(y_top_in, "in") - grid::unit(0.02, "in"),
                                                        just = c("left", "top"), gp = grid::gpar(fontsize = bn_letter_pt, fontface = "bold"))
      bx <- bn_panel_left_in; bw <- panelB_w_in - bx                                              # Panels B and C start bx in from the left, leaving room for the letters
      print(panelB$fig, vp = grid::viewport(x = grid::unit(bx, "in"), y = grid::unit(c_h_in, "in"), width = grid::unit(bw, "in"), height = grid::unit(b_h_in, "in"), just = c("left", "bottom")))
      if (panelB_split) {
        print(panelC$fig, vp = grid::viewport(x = grid::unit(bx, "in"), y = 0, width = grid::unit(bw, "in"), height = grid::unit(c_h_in, "in"), just = c("left", "bottom")))
        draw_pt_legend(bw, c_h_in, x0_in = bx)
      }
      letter("A", fig3_h); letter("B", c_h_in + b_h_in); if (panelB_split) letter("C", c_h_in)
    }
    show_fig("Fig3_PanelsAB", panelB_w_in, fig3_h, draw_fig3)
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


###############################################################################
# 9. SENSITIVITY ANALYSES AND CONTROLS
#    Each repeats the PMA x prematurity scan under a different specification and
#    compares the FDR-significant set with the primary result (`sigEdges`).
#
#    OVERVIEW (what each check tests and what it showed; details above each block)
#      9a  clinical cutoffs ........ 9 of the edges FDR-significant; 8 of the 16 hold (cutoff-dependent set)
#      9b  continuous birth GA ..... 0 significant (the effect is a step, not a linear gradient)
#      9c  scans after 32 wk PMA ... 6 of 16 hold; same slope pattern (less power)
#      9d  within/between infant ... omnibus not significant (underpowered); within-infant slopes match
#      9e  data amount (finalVol) .. robust (direction kept for >= 97% of edges); tSNR / motion not yet run
#      9f  age-matched fetal/preterm  robust (1,174-1,354 of 1,404 edges stay significant; all same direction)
#      9g  postnatal age ........... robust (14 of 16 edges hold when adjusted); (C) PMA vs. postnatal clock: see 9g
#      9h  unharmonized data ....... fetal-vs-preterm fully robust; PMA x prematurity 9 of 16 hold
#      9i  distance between ROIs ... group differences unrelated to distance
###############################################################################

compare_sets <- function(res_alt) {
  s <- res_alt$Connection[res_alt$pAdj < .05 & !is.na(res_alt$pAdj)]
  summary_vec <- c(raw_p05 = sum(res_alt$p < .05, na.rm = TRUE), fdr05 = length(s),
                   of_primary = sum(sigEdges %in% s), n_primary = length(sigEdges), newly_sig = sum(!s %in% sigEdges))
  list(summary = summary_vec, dropped = setdiff(sigEdges, s), entered = setdiff(s, sigEdges),
       primary_edges = res_alt[match(sigEdges, res_alt$Connection), c("Connection", "chi", "p", "pAdj")])
}

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

## ---- 9e. Fetal vs. preterm with data-quality covariates and similar data amounts -----
# WHAT IT TESTS: whether the fetal-vs-preterm group difference is explained by differences in the amount
#   or quality of data between fetal and preterm scans. finalVol = number of volumes retained per scan.
#   (A) finalVol as a covariate (raw); (B) finalVol centered within group (so the covariate is not
#   collinear with group); (C) similar volumes: restrict to the finalVol range both groups share, and
#   1:1 match fetal to preterm scans on finalVol; (D) all QC covariates together once tSNR / motion are
#   added (add columns to scanInfo and list them in qc_more).
# STATUS: finalVol RUN; tSNR / motion not yet available. Result summaries use the same columns as 9f:
#   fdr05, still_sig_of_primary (of the 1,404 primary edges), same_direction_of_primary, r_beta_with_primary.
# NOTE: this is a proxy for "equal amounts of data". The strict version (truncate every scan to the same number
#   of volumes and recompute FC) needs the ROI time series.
# DATA: finalVol differs by group: fetal 80-281 (median 117, 192 scans), preterm 120-296 (median 220, 346
#   scans); correlation with group r = 0.63. (An earlier minimum of 32 was a data-entry error, corrected to 132.)
#   102 of the 192 fetal scans (53%) have < 120 volumes while every preterm scan has >= 120, because the
#   criterion is 4 min of data in both groups: >= 120 volumes at TR = 2 s (preterm) and >= 80 volumes at
#   TR = 3 s (fetal). Fetal scans therefore have fewer volumes by design, and both groups have the same
#   duration of data; the volume controls below address the difference in volume count.
# RESULT (each vs. the primary 1,404 significant edges; 'same direction' = of those 1,404):
#   A raw finalVol covariate ........ 1,113 significant; 1,029 still significant; 1,403 same direction; r = 0.962
#   B finalVol centered in group .... 1,426 significant; 1,387 still significant; 1,404 same direction; r = 0.999
#   C1 shared range (120-281) ....... 375 scans (90 fetal); 1,164 significant; 1,083 still; 1,403 same; r = 0.962
#   C2 1:1 matched on finalVol ...... 92 pairs (165.0 vs 167.7 volumes, SD 48.5 vs 46.5); 603 significant;
#                                      581 still; 1,373 same direction; r = 0.910
# INTERPRETATION: the group difference is not explained by the amount of data: direction is preserved for at
#   least 97.8% of edges and the betas correlate 0.91-0.999. Fewer significant edges in A (collinearity inflates
#   standard errors by about 30%) and C (90-92 fetal scans) reflect precision, not a changed effect. B equals the
#   primary result by construction (centering leaves the group coefficient unchanged); it only shows that
#   volumes matter little within groups. This covers data AMOUNT, not data QUALITY: tSNR and motion are still open.
# REPORT: "Results were unchanged when including finalVol as a covariate (1,113 significant edges; 1,403 of
#   1,404 primary edges in the same direction) or restricting to scans with >= 120 volumes in both groups
#   (1,164 significant edges; 1,403 of 1,404 in the same direction)."
if (RUN_QC_COVARIATE_CONTROL) {
  stopifnot("finalVol" %in% names(scanInfo))
  print(tapply(scanInfo$finalVol, scanInfo$scanType, function(x)
    c(n = sum(!is.na(x)), missing = sum(is.na(x)), min = min(x, na.rm = TRUE), median = median(x, na.rm = TRUE), max = max(x, na.rm = TRUE))))
  stopifnot(!anyNA(scanInfo$finalVol))                                 # every scan needs a value
  print(c(r_finalVol_with_group = cor(scanInfo$finalVol, as.numeric(scanInfo$scanType == "Preterm"))))   # collinearity with group
  scanInfo$finalVol_c <- scanInfo$finalVol - ave(scanInfo$finalVol, scanInfo$scanType)                   # group-mean-centered

  ## A. raw finalVol as a covariate
  res_qc_raw <- run_fetal(scanInfo, "finalVol", extra = "finalVol")
  print(cmp_fetal(res_qc_raw, n_scans = nrow(scanInfo)))
  ## B. finalVol centered within group
  res_qc_cen <- run_fetal(scanInfo, "finalVolC", extra = "finalVol_c")
  print(cmp_fetal(res_qc_cen, n_scans = nrow(scanInfo)))

  ## C. similar volumes between fetal and preterm scans
  vol_rng <- lapply(split(scanInfo$finalVol, scanInfo$scanType), range)
  vol_lo <- max(sapply(vol_rng, `[`, 1)); vol_hi <- min(sapply(vol_rng, `[`, 2))
  si_vol_range <- droplevels(scanInfo[scanInfo$finalVol >= vol_lo & scanInfo$finalVol <= vol_hi, ])
  print(c(vol_lo = vol_lo, vol_hi = vol_hi)); print(table(si_vol_range$scanType))
  res_vol_range <- run_fetal(si_vol_range, "volRange")
  print(cmp_fetal(res_vol_range, n_scans = nrow(si_vol_range)))
  si_vol_match <- droplevels(match_scans(scanInfo, "finalVol", caliper = 0.1 * sd(scanInfo$finalVol)))
  print(c(matched_pairs = sum(si_vol_match$scanType == "Fetal")))
  print(tapply(si_vol_match$finalVol, si_vol_match$scanType, function(x) c(mean = mean(x), sd = sd(x))))   # balance on volumes
  res_vol_match <- run_fetal(si_vol_match, "volMatched")
  print(cmp_fetal(res_vol_match, n_scans = nrow(si_vol_match)))

  ## D. every QC covariate together (add e.g. scanInfo$tSNR_c <- ... and set qc_more <- c("tSNR_c"))
  qc_more <- character(0)
  if (length(qc_more) > 0) {
    stopifnot(all(qc_more %in% names(scanInfo)), !anyNA(scanInfo[qc_more]))
    res_qc_all <- run_fetal(scanInfo, "allQC", extra = c("finalVol_c", qc_more))
    print(cmp_fetal(res_qc_all, n_scans = nrow(scanInfo)))
  }
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
  # Temporarily swap the connectivity matrix (and the cache key) used by the scan functions
  with_fc <- function(fc, code) {
    old <- list(goodFCData = goodFCData, data_key = data_key)
    on.exit({ goodFCData <<- old$goodFCData; data_key <<- old$data_key }, add = TRUE)
    goodFCData <<- fc
    data_key   <<- paste(old$data_key, "unharmonized", sum(as.matrix(fc)))
    force(code)
  }
  fetalRes_raw <- with_fc(rawFC, run_fetal(scanInfo, "unharmonized"))
  print(cmp_fetal(fetalRes_raw, n_scans = nrow(scanInfo)))
  res_raw <- with_fc(rawFC, wald_scan(pt, fml_main, "ScanPMA:BirthCat"))
  cmp_raw <- compare_sets(res_raw); print(cmp_raw)
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

  dir_sig <- ifelse(fetalRes$scanTypeb[sigFetal] > 0, "Preterm greater", "Fetal greater")      # distance of the significant edges, by direction
  print(tapply(edge_dist[sigFetal], dir_sig, function(x) round(c(n = length(x), median = median(x), q25 = quantile(x, .25, names = FALSE), q75 = quantile(x, .75, names = FALSE)), 1)))
  print(wilcox.test(edge_dist[sigFetal] ~ dir_sig))

  ## ---- (2) sensitivity: distances from every template age ----
  rho_by_ga <- data.frame(GA = ga_avail,
                          rho_beta = sapply(as.character(ga_avail), function(g) cor(Dedge[, g], fetalRes$scanTypeb, method = "spearman")))
  rho_by_ga$rho_beta <- round(rho_by_ga$rho_beta, 3); rownames(rho_by_ga) <- NULL
  print(rho_by_ga)

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
  show_fig("SuppFig_distance", 3.4, 3.0, function() print(p_dist))

  dc <- scanInfo; dc$coupling <- coupling(subsets[["All edges"]])
  p_coup <- ggplot(dc, aes(ScanPMA, coupling, color = scanType)) +
    geom_point(alpha = 0.5, size = 0.9, stroke = 0) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.8) +
    scale_color_manual(values = c(Fetal = col_fetal, Preterm = col_preterm)) +
    labs(x = "PMA at scan (weeks)", y = "Distance dependence of FC (Spearman rho)", color = NULL) +
    theme_classic(base_size = 9) + theme(axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4), legend.position = "bottom")
  show_fig("SuppFig_distance_coupling", 3.4, 3.2, function() print(p_coup))
}

###############################################################################
# 10. SAVE TABLES AND RECORD SESSION
###############################################################################

if (SAVE_TABLES) {
  write.csv(fetalRes,        file.path(figureDir, "fetalVsPreterm_allEdges.csv"),        row.names = FALSE)
  write.csv(table2,          file.path(figureDir, "Table2_homotopic.csv"),               row.names = FALSE)
  write.csv(intxRes,         file.path(figureDir, "pmaByPrematurity_allEdges.csv"),      row.names = FALSE)
  write.csv(tab3_display,    file.path(figureDir, "Table3.csv"),                         row.names = FALSE)
  write.csv(slopesMain_wide, file.path(figureDir, "slopesAndContrasts_sigEdges.csv"),    row.names = FALSE)
}
sessionInfo()
