# EarlyAutismMarkers_updated.R
#
# Clean restart of the regional anatomical analysis from EarlyAutismMarkers.qmd.
# Question: do high-risk (M-CHAT Fail) and low-risk (M-CHAT Pass) fetuses differ in
# regional brain growth trajectories (group x scan age)?
#
# Differences from the .qmd:
#   * scan age is CENTERED, so the "group" term is the difference at the average scan age
#     (not at week 0)
#   * Pass is the reference level, so positive coefficients mean HIGH RISK > LOW RISK
#   * one ID helper instead of the repeated if/else loops; no setwd()
#   * results are read from the model's named coefficients, not tTable row numbers
#   * adds a global check (how many ROIs show group x age?) and a permutation test
#   * figures show fitted trajectories over scan age, not bars of group means
#
# NOT included yet: behavioral QA tables and the functional (FC) analyses.


# ---- 1. Setup ----------------------------------------------------------------

library(here)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(MatchIt)
library(nlme)
library(emmeans)

baseDir <- here("data")
figDir  <- here("paper", "figures")

tissueDir   <- file.path(baseDir, "anatomical", "tissueVolumes_Individual")
regionalDir <- file.path(baseDir, "anatomical", "regionalVolumes_Individual")

# Scans excluded from the regional models (documented here so tables and figures agree).
# The old figure code also dropped 'fetus_02178_01'; the old table code did not.
# Decide once and keep it in one place:
excludeScans <- c("rescue_0175_01", "rescue_0184_01")   # add "fetus_02178_01" if intended

lmeCtl <- lmeControl(msMaxIter = 200, tolerance = 1e-4, msTol = 1e-5, opt = "optim")


# ---- 2. Helpers --------------------------------------------------------------

# Study-specific raw study_id -> short image ID (e.g. 'fetus_00556', 'pl_10011', 'rescue_0142')
makeShortID <- function(study_id, study) {
  mapply(function(id, proj) {
    id <- as.character(id)
    if (proj == "ABC") {
      paste0("fetus_00", id)
    } else if (proj == "CHD R01") {
      paste0("fetus_0", id)
    } else if (proj == "Placenta") {
      paste0("pl_1", sprintf("%04d", as.numeric(substr(id, 4, nchar(id)))))
    } else {
      n <- as.numeric(substr(id, 5, nchar(id)))
      paste0("rescue_", sprintf(ifelse(n < 100, "%03d", "%04d"), n))
    }
  }, study_id, study, USE.NAMES = FALSE)
}

# Short ID + scan timepoint ('ScanAge1' -> '_01')
makeLongID <- function(shortID, scanTimepoint) {
  tp <- as.character(scanTimepoint)
  paste0(shortID, "_0", substr(tp, nchar(tp), nchar(tp)))
}


# ---- 3. Build the matched sample --------------------------------------------

testData <- read.csv(file.path(baseDir, "anatomical", "combinedASDData.csv"), stringsAsFactors = FALSE)
testFetal <- subset(testData, PretermStatus == "Term" & CHDStatus != "CHD" & AnyFetal == TRUE)

# One row per scan
fetalMchatMelt <- melt(testFetal, measure.vars = c("ScanAge1", "ScanAge2"),
                       id.vars = c("study_id", "mchat_result", "study"))
fetalMchatMelt <- subset(fetalMchatMelt, mchat_result != "" & value != "")
rownames(fetalMchatMelt) <- seq_len(nrow(fetalMchatMelt))
colnames(fetalMchatMelt) <- c("study_id", "mchat_result", "study", "scanTimepoint", "scanAge")

# 1:1 match of M-CHAT Fail to Pass on scan age and study
set.seed(1000)
fetalMchatMatch <- matchit(factor(mchat_result, levels = c("Pass", "Fail")) ~ scanAge + study,
                           data = fetalMchatMelt)
mm <- fetalMchatMatch$match.matrix
matchedRows <- c(as.vector(mm), rownames(mm)[!is.na(mm)])
fetalMatched <- fetalMchatMelt[matchedRows[!is.na(matchedRows)], ]

fetalMatched$shortID <- makeShortID(fetalMatched$study_id, fetalMatched$study)
fetalMatched$scanTimepoint <- as.character(fetalMatched$scanTimepoint)
fetalMatched$longID <- makeLongID(fetalMatched$shortID, fetalMatched$scanTimepoint)
rownames(fetalMatched) <- fetalMatched$longID

# Hand swaps: replace Placenta-cohort Pass scans (high-risk pregnancies) with matched controls
HANDMATCH <- TRUE
if (HANDMATCH) {
  removeList <- c("pl_10011_01", "pl_10027_01", "pl_10029_01", "pl_10029_02", "pl_10051_01",
                  "pl_10091_01", "pl_10106_01", "pl_10124_01", "pl_10125_01", "fetus_02153_01")
  addListShort <- c("2120", "2121", "2020", "RES-142", "556", "2172", "2107", "RES-193", "2104", "RES-042")
  addListFull  <- c("fetus_02120_02", "fetus_02121_02", "fetus_02020_01", "rescue_0142_02",
                    "fetus_00556_01", "fetus_02172_01", "fetus_02107_02", "rescue_0193_01",
                    "fetus_02104_02", "rescue_0042_01")

  addFrame <- fetalMchatMelt[fetalMchatMelt$study_id %in% addListShort, ]
  addFrame$shortID <- makeShortID(addFrame$study_id, addFrame$study)
  addFrame$scanTimepoint <- as.character(addFrame$scanTimepoint)
  addFrame$longID <- makeLongID(addFrame$shortID, addFrame$scanTimepoint)
  addFrame <- addFrame[addFrame$longID %in% addListFull, ]
  rownames(addFrame) <- addFrame$longID

  fetalMatched <- fetalMatched[!(fetalMatched$longID %in% removeList), ]
  fetalMatched <- rbind(fetalMatched, addFrame)
}

# Drop poor-quality placenta-cohort scans and their matched controls
dropList <- c("pl_10028_02", "pl_10149_01", "pl_10149_02", "pl_10151_01",
              "rescue_071_02", "fetus_02075_01", "rescue_017_02", "fetus_02070_02")
fetalMatched <- fetalMatched[!(fetalMatched$longID %in% dropList), ]

# Quick balance checks
with(fetalMatched, {
  cat("Children  (Pass / Fail):", length(unique(shortID[mchat_result == "Pass"])),
      length(unique(shortID[mchat_result == "Fail"])), "\n")
  cat("Scans     (Pass / Fail):", sum(mchat_result == "Pass"), sum(mchat_result == "Fail"), "\n")
})
print(t.test(scanAge ~ mchat_result, data = fetalMatched))


# ---- 4. Read regional volumes -----------------------------------------------

regionalFiles <- list.files(regionalDir, full.names = TRUE)
regionalIDs   <- substr(basename(regionalFiles), 1, nchar(basename(regionalFiles)) - 18)

regionalLabels <- read.table(file.path(baseDir, "anatomical", "regionalROILabels.txt"))
regionalAll <- as.data.frame(matrix(NA_real_, nrow = nrow(fetalMatched), ncol = 87))
rownames(regionalAll) <- fetalMatched$longID
colnames(regionalAll) <- t(regionalLabels)

for (i in seq_len(nrow(regionalAll))) {
  f <- regionalFiles[regionalIDs %in% rownames(regionalAll)[i]]
  if (length(f) == 0) next
  tmp <- read.table(f)
  regionalAll[i, 1:(ncol(tmp) - 2)] <- tmp[2, 3:ncol(tmp)]
}
regionalAll[] <- lapply(regionalAll, function(x) as.numeric(as.character(x)))

# Total brain volume (excludes CSF / white-matter-adjacent labels as in the original), then ROI set
tbv <- rowSums(regionalAll[, c(1:48, 51:82)], na.rm = TRUE)
regionalData <- regionalAll[, 1:82]
for (drop in c("WM", "LatVen", "CRB", "CorpusCallosum")) {
  regionalData <- regionalData[, grep(drop, colnames(regionalData), invert = TRUE)]
}

# Missing scans (no volume file) show up here
cat("Scans with no regional volume file:", sum(rowSums(!is.na(regionalData)) == 0), "\n")


# ---- 5. Per-ROI model: risk group x scan age --------------------------------

# Long data for one ROI; scan age centered at the ROI sample's mean
getROIData <- function(roi) {
  d <- data.frame(longID = rownames(regionalData),
                  Volume = regionalData[[roi]],
                  TBV = tbv,
                  stringsAsFactors = FALSE)
  d <- merge(fetalMatched[, c("shortID", "longID", "scanAge", "mchat_result")], d, by = "longID")
  d <- d[!(d$longID %in% excludeScans), ]
  d$VolNormed <- d$Volume / d$TBV
  d <- d[is.finite(d$VolNormed), ]
  d$scanAge <- as.numeric(d$scanAge)
  d$scanAgeC <- d$scanAge - mean(d$scanAge)
  d$Risk <- factor(ifelse(d$mchat_result == "Fail", "High Risk", "Low Risk"),
                   levels = c("Low Risk", "High Risk"))      # Low Risk = reference
  d
}

fitROI <- function(d) {
  lme(VolNormed ~ Risk * scanAgeC, random = ~ 1 | factor(shortID), data = d, control = lmeCtl)
}

roiNames <- colnames(regionalData)
resultsList <- lapply(roiNames, function(roi) {
  d <- getROIData(roi)
  m <- tryCatch(fitROI(d), error = function(e) NULL)
  if (is.null(m)) return(data.frame(ROI = roi))
  tt <- summary(m)$tTable
  data.frame(ROI = roi,
             n.scans = nrow(d),
             LowRisk.Mean  = mean(d$VolNormed[d$Risk == "Low Risk"]),
             HighRisk.Mean = mean(d$VolNormed[d$Risk == "High Risk"]),
             # difference at the average scan age (High - Low)
             Grp.B = tt["RiskHigh Risk", "Value"],
             Grp.p = tt["RiskHigh Risk", "p-value"],
             # growth slope in Low Risk, and how much steeper it is in High Risk
             Age.B = tt["scanAgeC", "Value"],
             Age.p = tt["scanAgeC", "p-value"],
             Itx.B = tt["RiskHigh Risk:scanAgeC", "Value"],
             Itx.p = tt["RiskHigh Risk:scanAgeC", "p-value"])
})
results <- do.call(rbind, resultsList)
results <- results[!is.na(results$Itx.p), ]

results$Grp.p.Adj <- p.adjust(results$Grp.p, method = "fdr")
results$Age.p.Adj <- p.adjust(results$Age.p, method = "fdr")
results$Itx.p.Adj <- p.adjust(results$Itx.p, method = "fdr")

cat("\nROIs with group difference at average age, FDR < .05:\n")
print(subset(results, Grp.p.Adj < .05))
cat("\nROIs with group x age interaction, FDR < .05:\n")
print(subset(results, Itx.p.Adj < .05))
cat("\nROIs with group x age interaction, uncorrected p < .05:\n")
print(subset(results, Itx.p < .05)[, c("ROI", "Itx.B", "Itx.p", "Itx.p.Adj")])


# ---- 6. Global check: more interactions than chance? -------------------------

nObs  <- sum(results$Itx.p < .05)
nExp  <- nrow(results) * .05
cat(sprintf("\nROIs with Itx.p < .05: %d observed vs %.1f expected by chance (of %d)\n",
            nObs, nExp, nrow(results)))
cat("Sign of Itx.B among those:", table(sign(subset(results, Itx.p < .05)$Itx.B)), "\n")

# Permutation test (accounts for correlated ROIs): shuffle risk labels across CHILDREN,
# keeping each child's scans together, and recount ROIs with Itx.p < .05.
# This refits every ROI for every permutation -- slow. Turn on when ready.
RUN_PERMUTATION <- FALSE
nPerm <- 500

if (RUN_PERMUTATION) {
  roiData <- lapply(setNames(results$ROI, results$ROI), getROIData)
  childRisk <- unique(roiData[[1]][, c("shortID", "Risk")])

  countHits <- function(riskByChild) {
    sum(sapply(roiData, function(d) {
      d$Risk <- riskByChild[as.character(d$shortID)]
      p <- tryCatch(summary(fitROI(d))$tTable["RiskHigh Risk:scanAgeC", "p-value"],
                    error = function(e) NA_real_)
      !is.na(p) && p < .05
    }))
  }

  set.seed(2026)
  permHits <- replicate(nPerm, {
    shuffled <- setNames(sample(childRisk$Risk), childRisk$shortID)
    countHits(factor(shuffled, levels = levels(childRisk$Risk)))
  })
  cat(sprintf("Permutation p (>= %d ROIs with Itx.p < .05): %.4f\n", nObs,
              (sum(permHits >= nObs) + 1) / (nPerm + 1)))
}


# ---- 7. Figures: fitted trajectories by risk group --------------------------

figROIs <- c(HPC_L = "Hippocampus (Left)", HPC_R = "Hippocampus (Right)",
             Amyg_L = "Amygdala (Left)", FusiformAnt_GM_R = "Fusiform (Right)",
             Caudate_R = "Caudate (Right)", Subthalamic_L = "Subthalamic (Left)",
             Lentiform_L = "Pallidum (Left)")

plotTrajectory <- function(roi, label) {
  d <- getROIData(roi)
  m <- fitROI(d)
  newd <- expand.grid(scanAgeC = seq(min(d$scanAgeC), max(d$scanAgeC), length.out = 50),
                      Risk = levels(d$Risk))
  newd$pred <- predict(m, newd, level = 0)
  newd$scanAge <- newd$scanAgeC + mean(d$scanAge)
  q <- results$Itx.p.Adj[results$ROI == roi]
  ggplot() +
    geom_point(data = d, aes(scanAge, VolNormed, color = Risk), alpha = .5, size = 2) +
    geom_line(data = newd, aes(scanAge, pred, color = Risk), linewidth = 1.5) +
    scale_color_manual(values = c("#4a7fc1", "#c24a4a"), name = NULL) +
    labs(x = "Scan age (weeks)", y = paste(label, "volume / TBV"),
         subtitle = sprintf("group x age: p = %.3f, FDR q = %.2f",
                            results$Itx.p[results$ROI == roi], q)) +
    theme_classic(base_size = 13) +
    theme(legend.position = "bottom")
}

trajPlots <- Map(plotTrajectory, names(figROIs), figROIs)

#tiff(file.path(figDir, "regionalTrajectories.tiff"), height = 8, width = 16, units = "in", res = 600, compression = "lzw")
ggarrange(plotlist = trajPlots, nrow = 2, ncol = 4, common.legend = TRUE, legend = "bottom")
#dev.off()
