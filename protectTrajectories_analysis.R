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
#   4  Fetal vs. preterm, all edges (Fig. 2A-B, Table 2)
#   5  Preterm only: PMA x prematurity, all edges (Fig. 3, Table 3)
#   6  Fig. 3 Panel B (group trajectories)
#   7  Fig. 3 Panel A (BrainNet Viewer edge file)
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
#                    "circlize","ggpubr","scales","reshape2"))
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
SAVE_FIGS                 <- FALSE  # TRUE: write .tif files to figureDir
SAVE_TABLES               <- FALSE  # TRUE: write .csv result tables to figureDir
RUN_QA_PLOTS              <- TRUE   # scan-QA figures (section 1)
RUN_HARMONIZATION_CHECKS  <- TRUE   # scanner effect / PCA / PMA preservation (section 2)
RUN_SENSITIVITY           <- TRUE   # section 9: five lmer scans over all edges (slow serially; fast with N_CORES, cached)
RUN_QC_COVARIATE_CONTROL  <- FALSE  # needs finalVol / tSNR for fetal scans (section 9e)
N_CORES                   <- max(1, parallel::detectCores() - 1, na.rm = TRUE)  # cores for the edge-wise loops; 1 = serial
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

# Colors: fetal = red, preterm = navy; preterm subgroups = navy ramp
col_fetal   <- "#cc2936"
col_preterm <- "#08415c"
cols_pt     <- c(EPT = "#08415c", VPT = "#2f7fb0", MPT = "#8cc2e0")
grp_labels  <- c(EPT = "Extremely Premature", VPT = "Very Premature", MPT = "Moderately Premature")

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

# Draw a figure on screen and, if SAVE_FIGS, also write it to a .tif file.
# `draw` is a function that prints a ggplot / draws base graphics.
show_fig <- function(file, w, h, draw, res = 600) {
  if (SAVE_FIGS) {
    tiff(file.path(figureDir, file), units = "in", width = w, height = h, res = res, compression = "lzw")
    on.exit(dev.off())
  }
  draw()
}

# p-value formatting for tables: "<.001" or ".034"
fmt_p <- function(p) ifelse(p < .001, "<.001", sub("^0", "", sprintf("%.3f", p)))

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

# A scan is usable if it has > 119 volumes (>= 4 min of data at TR = 2 s)
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

## ---- 4b. Fig. 2A: connectivity matrices (edge beta, significant edges only) -
cp <- data.frame(Conn = fetalRes$Connection,
                 ROI1 = sapply(strsplit(fetalRes$Connection, "_", fixed = TRUE), `[`, 1),
                 ROI2 = sapply(strsplit(fetalRes$Connection, "_", fixed = TRUE), `[`, 2),
                 Region1 = columnMembership$Lobe1[match(fetalRes$Connection, columnMembership$Connection)],
                 B = fetalRes$scanTypeb, pAdj = fetalRes$scanTypep.Adj, stringsAsFactors = FALSE)
roiOrder <- unlist(lapply(c("Frontal", "Temporal", "Parietal", "Occipital"),
                          function(l) unique(cp$ROI1[cp$Region1 == l])))
heat_lim <- ceiling(max(abs(fetalRes$scanTypeb[sigFetal])) * 100) / 100   # symmetric color limit

# Mirror each edge to fill both triangles, keeping only edges that satisfy `keep`
heat_data <- function(keep) {
  keep[is.na(keep)] <- FALSE
  d <- cp; d$B[!keep] <- NA
  d2 <- d; names(d2)[names(d2) %in% c("ROI1", "ROI2")] <- c("ROI2", "ROI1")
  rbind(d, d2)
}
heat_plot <- function(d, low, mid, high, midpoint, limits) {
  ggplot(d, aes(x = factor(ROI1, levels = roiOrder, ordered = TRUE),
                y = factor(ROI2, levels = roiOrder, ordered = TRUE), fill = B)) +
    scale_y_discrete(limits = rev) +
    scale_fill_gradient2(low = low, mid = mid, high = high, na.value = "#ffffff00", midpoint = midpoint, limits = limits) +
    geom_tile() + xlab("") + ylab("") +
    theme_paper(grid_y = FALSE, legend = "bottom") + theme(axis.text.x = element_text(angle = 90))
}
p_heat_preterm <- heat_plot(heat_data(cp$B > 0 & cp$pAdj < .05), "#adb3cf", "#373d7e", "#15085c", heat_lim / 2, c(0, heat_lim))
p_heat_fetal   <- heat_plot(heat_data(cp$B < 0 & cp$pAdj < .05), "#650f16", "#cc2936", "#f4b1ac", -heat_lim / 2, c(-heat_lim, 0))
show_fig("preterm_connectivityHeatMap.tiff", 12, 12, function() print(p_heat_preterm))
show_fig("term_connectivityHeatMap.tiff",    12, 12, function() print(p_heat_fetal))

## ---- 4c. Fig. 2B: chord diagrams -------------------------------------------
# ROI order and gaps: right hemisphere then left (mirrored), grouped by lobe
lobeOrder <- c("Frontal", "Temporal", "Parietal", "Occipital")
ord       <- order(factor(columnMembershipCortOnly$Lobe1, levels = lobeOrder))
bilatROI  <- unique(sub("-[LR]$", "", columnMembershipCortOnly$ROI1[ord]))   # bilateral ROI names in lobe order
roiList   <- c(paste0(bilatROI, "-R"), rev(paste0(bilatROI, "-L")))
roi_hemi  <- substr(roiList, nchar(roiList), nchar(roiList))
roi_lobe  <- roiData$Lobe[match(roiList, roiData$FullLabel)]
gapSize   <- sapply(seq_along(roiList), function(i) {
  if (i == length(roiList) || roi_hemi[i] != roi_hemi[i + 1]) return(10)      # hemisphere break
  if (roi_lobe[i] == roi_lobe[i + 1]) 2 else 5                                # within lobe vs. between lobes
})

# Sector colors by lobe (frontal, temporal, parietal/cingulate, occipital)
lobe_cols <- list("#ef476f" = c("PreCG","SFGdor","ORBsupb","MFG","ORBmid","IFGoperc","IFGtriang","ORBinf","ROL","SMA","OLF","SFGmed","ORBmed","REC","ACG"),
                  "#ffd166" = c("INS","PHG","HES","STG","TPOsup","MTG","TPOmid","ITG"),
                  "#06d6a0" = c("MCG","PCG","PoCG","SPG","IPL","SMG","ANG","PCUN","PCL"),
                  "#118ab2" = c("CAL","CUN","LING","SOG","MOG","IOG","FFG"))
roiColor <- unlist(lapply(names(lobe_cols), function(cl)
  setNames(rep(cl, 2 * length(lobe_cols[[cl]])), c(paste0(lobe_cols[[cl]], "-L"), paste0(lobe_cols[[cl]], "-R")))))

chord_df <- data.frame(From = columnMembershipCortOnly$ROI1, To = columnMembershipCortOnly$ROI2,
                       Value = as.numeric(sigFetal[match(allConn, fetalRes$Connection)]), stringsAsFactors = FALSE)
b_dir    <- fetalRes$scanTypeb[match(allConn, fetalRes$Connection)] > 0
chord_col_all     <- ifelse(b_dir, col_preterm, col_fetal)
chord_col_fetal   <- ifelse(b_dir, "#00000000", col_fetal)
chord_col_preterm <- ifelse(b_dir, col_preterm, "#00000000")

draw_chord <- function(cols) {
  circos.clear()
  circos.par(start.degree = 86, gap.after = gapSize)
  chordDiagram(chord_df, order = roiList, grid.col = roiColor, col = cols, symmetric = TRUE, scale = TRUE,
               annotationTrack = "grid", annotationTrackHeight = c(0.03, 0.03))
  circos.track(track.index = 1, bg.border = NA, panel.fun = function(x, y)
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(-.25, .5)))
  circos.clear()
}
show_fig("chordPlot2Groups_UPDATE.tiff", 16, 8, function() {
  par(mfrow = c(1, 2)); draw_chord(chord_col_fetal); draw_chord(chord_col_preterm); par(mfrow = c(1, 1)) })

## ---- 4d. Table 2 + homotopic bar plot --------------------------------------
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
p_homo <- ggplot(homoPlot, aes(x = Connection, y = Mean, group = Group)) +
  scale_color_manual(values = c(col_fetal, col_preterm)) +
  scale_fill_manual(values = paste0(c(col_fetal, col_preterm), "aa")) +
  geom_bar(width = .8, stat = "identity", position = position_dodge2(), aes(color = Group, fill = Group), linewidth = 1) +
  geom_errorbar(width = .8, position = position_dodge2(width = .2), aes(ymin = Mean - SE, ymax = Mean + SE)) +
  geom_hline(yintercept = 0, linewidth = 1.5) +
  xlab("") + ylab("Estimated marginal mean connectivity (Fisher z)") +
  theme_paper(grid_x = TRUE, legend = "bottom") + theme(axis.text.x = element_text(angle = 60, hjust = 1))
show_fig("Homotopic2Groups_UPDATE.tiff", 19.5, 6, function() print(p_homo))


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
  `PMA B` = sprintf("%.3f", tab3$PMA_B), `PMA p` = fmt_p(tab3$PMA_p),
  `Scanner B` = sprintf("%.3f", tab3$Scanner_B), `Scanner p` = fmt_p(tab3$Scanner_p),
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
pt_size         <- 0.9     # point size
pt_alpha        <- 0.45    # point opacity
line_w          <- 1.1     # colored line width
halo_w          <- 0.2     # extra width of the white halo under each line
axis_title_size <- 6.5     # used for BOTH axis titles
panel_ncol      <- 4       # panels per row; 16 edges fit a 4x4 grid exactly (legend goes below the grid)

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

ylab_for <- function(cn) {                                            # two lines: ROI 1 -, ROI 2
  p <- strsplit(cn, "_", fixed = TRUE)[[1]]
  paste0(full_name(p[1]), " \u2013\n", full_name(p[2]))
}

make_panel_b <- function(edges) {
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
      theme(axis.line = element_line(color = "#000000", linewidth = .5),
            axis.ticks = element_line(color = "#000000", linewidth = .4),
            axis.text = element_text(size = 6.5, color = "#000000"),
            axis.title.x = element_text(size = axis_title_size, color = "#000000"),
            axis.title.y = element_text(size = axis_title_size, color = "#000000"),
            plot.margin = margin(4, 8, 2, 2),
            panel.grid.major.y = element_line(color = "#e6e6e6", linewidth = .3),
            legend.position = "bottom", legend.text = element_text(size = 7))
  }
  plots <- lapply(seq_along(edges), function(j) edge_plot(edges[j], lab_order[j]))
  nrow_fig <- ceiling(length(plots) / panel_ncol)
  # One shared legend below the grid (taken from the first panel)
  list(fig = ggpubr::ggarrange(plotlist = plots, ncol = panel_ncol, nrow = nrow_fig, align = "hv",
                               common.legend = TRUE, legend = "bottom"),
       nrow = nrow_fig)
}

panelB <- make_panel_b(sigEdges)
show_fig("Figure3_PanelB.tiff", 2.6 * panel_ncol, 2.5 * panelB$nrow + 0.5, function() print(panelB$fig))   # width x height in inches


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
  fmt <- function(r) sprintf("%.3f [%.3f, %.3f]", r$est, r$lo, r$hi)
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
  # Changes group membership of the infants born at exactly 28.0 and 32.0 weeks.
  # Harmonization is NOT rerun (it used the paper's cutoffs).
  pt_clin <- pt; pt_clin$BirthCat <- birth_cat(pt_clin$BirthPMA, rule = "clinical")
  print(group_counts(pt_clin))
  res_clin <- wald_scan(pt_clin, fml_main, "ScanPMA:BirthCat")
  cmp_clin <- compare_sets(res_clin); print(cmp_clin)

  ## ---- 9b. Birth GA as a continuous moderator (no cutoffs; 1 df) ------------
  pt_cont <- pt; pt_cont$GA <- pt_cont$BirthPMA - mean(pt_cont$BirthPMA)
  res_cont <- wald_scan(pt_cont, FCz ~ scanner + ScanPMA * GA, "ScanPMA:GA")
  cmp_cont <- compare_sets(res_cont); print(cmp_cont)

  ## ---- 9c. Scans after 32 weeks PMA only ------------------------------------
  # Drops the earliest EPT/VPT scans (all MPT scans are already later).
  pt32 <- droplevels(pt[pt$ScanPMA > 32, ])
  print(c(scans = nrow(pt32), infants = nlevels(pt32$shortID), dropped = nrow(pt) - nrow(pt32)))
  print(table(pt32$BirthCat, pt32$scanner))                            # every group should still have both scanners
  res_gt32 <- wald_scan(pt32, fml_main, "ScanPMA:BirthCat")
  cmp_gt32 <- compare_sets(res_gt32); print(cmp_gt32)
  slopes_gt32_wide <- slope_wide(slope_table(pt32, sigEdges, fml_main)); print(slopes_gt32_wide)

  ## ---- 9d. Within- vs. between-infant PMA effects ---------------------------
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

## ---- 9e. Fetal vs. preterm with data-quality covariates (pending data) ------
# Adds covariates to the fetal-vs-preterm model and compares the group effect.
# Requires the covariate for BOTH fetal and preterm scans. finalVol is already in
# scanInfo (fetal values must be present); add tSNR / motion as new columns, e.g.:
#   scanInfo$tSNR <- tsnr_table$tSNR[match(scanInfo$longID, tsnr_table$longID)]
# Re-centering within group (group-mean-centered covariate) is also advisable,
# because raw covariates can be collinear with group.
if (RUN_QC_COVARIATE_CONTROL) {
  qc_covars <- c("finalVol")                                          # add "tSNR", ... when available
  stopifnot(!anyNA(scanInfo[qc_covars]))
  fetalRes_qc <- do.call(rbind, lapply(allConn, fit_fetal_edge, extra = qc_covars))
  fetalRes_qc$scanTypep.Adj <- p.adjust(fetalRes_qc$scanTypep, "fdr")
  print(c(fdr05_primary = sum(sigFetal),
    fdr05_with_covariates = sum(fetalRes_qc$scanTypep.Adj < .05, na.rm = TRUE),
    still_sig = sum(sigFetal & fetalRes_qc$scanTypep.Adj < .05, na.rm = TRUE),
    same_direction = sum(sigFetal & sign(fetalRes_qc$scanTypeb) == sign(fetalRes$scanTypeb), na.rm = TRUE),
    r_beta = cor(fetalRes$scanTypeb, fetalRes_qc$scanTypeb, use = "complete.obs")))
}
# Still to do (need data not yet in this script): equal-data-per-scan reanalysis,
# age-matched fetal-vs-preterm comparison, postnatal-age check, and an
# unharmonized-data sensitivity analysis (scanner as a covariate).


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
