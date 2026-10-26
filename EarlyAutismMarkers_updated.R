# EarlyAutismMarkers_updated.R
#
# Clean restart of EarlyAutismMarkers.qmd. Same analyses, rebuilt to be faster and to
# answer the right question for the group comparisons:
#   do high-risk (M-CHAT Fail) and low-risk (M-CHAT Pass) fetuses differ in brain
#   structure / connectivity, and in how it changes with scan age (group x age)?
#
# What changed from the .qmd
#   * Scan age is CENTERED (CENTER_AGE). The "group" coefficient is then the difference at
#     the average scan age instead of at week 0 (an extrapolation that inflated the old
#     group effects). Set CENTER_AGE <- FALSE to reproduce the old behaviour.
#   * Low risk (Pass / TD) is always the reference level, so Grp.B = High - Low.
#   * One model function for every analysis (volumes, tissue, thickness, connectivity).
#     It fits the additive and interaction models once per ROI and returns everything the
#     four old tables (raw/normed x additive/interaction) contained, using named coefficients.
#   * Per-ROI / per-connection fits run in parallel (parallel::mclapply; single core on Windows).
#   * Data are merged once, not once per ROI. The 19,900 whole-brain connection models are
#     fit once, and every network table is a subset of that fit (only the FDR correction is
#     redone within each network) -- the old code refit the same models for each network.
#   * The 40,000-iteration connectogram loop is vectorized.
#   * No setwd(), no repeated if/else ID loops, no one-block-per-ROI copy/paste figures.
#   * Partial eta-squared (the source of hundreds of warnings) is dropped.
#
# What was dropped / simplified (tell me if you want any back)
#   * The older bar-chart regional figures (raw means; and the "Alt" versions): replaced by
#     fitted-trajectory figures.
#   * Network INTERACTION tables: the old code used model-adjusted means/SEs (emmeans) for CRB,
#     DMN, SomMot, SubCort, Auditory and Visual, but raw means for Limbic. All networks now use
#     the model-adjusted means/SEs (TD.Mean etc.); raw ones are kept as TD.Mean.Raw etc.
#     The trend-connection figure therefore plots adjusted means, as the old figure did for
#     all networks except Limbic.
#   * Cerebellar power-analysis scraps (hard-coded eta2_to_f2 values) and unused objects.
#   * Cortical thickness and tissue models now use the same exclusion list as the regional
#     models (excludeScans) instead of having none / a different one.
#
# NOTE: this script has not been run against the data. Expect to fix small typos on first run.


# ---- 1. Setup ----------------------------------------------------------------

library(here)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(MatchIt)
library(nlme)
library(circlize)
library(parallel)

baseDir <- here("data")
figDir  <- here("paper", "figures")

tissueDir   <- file.path(baseDir, "anatomical", "tissueVolumes_Individual")
regionalDir <- file.path(baseDir, "anatomical", "regionalVolumes_Individual")

# Options can be pre-set before source()-ing this file (the verification script does this);
# otherwise these defaults are used.
opt <- function(name, default) {
  env <- parent.frame()
  if (exists(name, envir = env, inherits = FALSE)) get(name, envir = env, inherits = FALSE) else default
}

CENTER_AGE      <- opt("CENTER_AGE", TRUE)       # center scan age so group effect = difference at the average age
HANDMATCH       <- opt("HANDMATCH", TRUE)        # hand swaps of Placenta-cohort Pass scans in the anatomical sample
WRITE_OUTPUT    <- opt("WRITE_OUTPUT", TRUE)     # write the functional CSV tables (as the .qmd did)
RUN_PERMUTATION <- opt("RUN_PERMUTATION", FALSE) # global permutation test for regional group x age (slow)
nPerm           <- opt("nPerm", 500)
INCLUDE_CEREBELLAR_PARCELLATION <- opt("INCLUDE_CEREBELLAR_PARCELLATION", FALSE)
outDir          <- opt("outDir", file.path(baseDir, "functional", "output"))   # where functional CSVs go

# Scans excluded from the anatomical models. (The old figure code also dropped 'fetus_02178_01';
# the old table code did not; the old tissue code dropped all three; old cortical thickness none.)
# Decide once, here. Tissue and thickness can be overridden separately.
excludeScans  <- opt("excludeScans", c("rescue_0175_01", "rescue_0184_01"))
excludeTissue <- opt("excludeTissue", excludeScans)
excludeCT     <- opt("excludeCT", excludeScans)

lmeCtl <- lmeControl(msMaxIter = 200, tolerance = 1e-4, msTol = 1e-5, opt = "optim")

nCores <- if (.Platform$OS.type == "windows") 1 else max(1, detectCores() - 1)
pLapply <- function(X, FUN) {
  if (nCores > 1) mclapply(X, FUN, mc.cores = nCores) else lapply(X, FUN)
}

if (WRITE_OUTPUT) dir.create(outDir, showWarnings = FALSE, recursive = TRUE)


# ---- 2. Helpers --------------------------------------------------------------

# Study-specific study_id -> short image ID ('fetus_00556', 'pl_10011', 'rescue_0142')
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

# Group comparison printout used for the demographics tables
compareGroups <- function(df, var, group = "mchat_result", levs = c("Pass", "Fail")) {
  a <- df[[var]][df[[group]] == levs[1]]
  b <- df[[var]][df[[group]] == levs[2]]
  print(t.test(a, b))
  cat(sprintf("%s: %s %.3f (%.3f) | %s %.3f (%.3f)\n", var,
              levs[1], mean(a, na.rm = TRUE), sd(a, na.rm = TRUE),
              levs[2], mean(b, na.rm = TRUE), sd(b, na.rm = TRUE)))
}

# --- the one model function used for every ROI / connection ---
# d needs: Risk (factor, 'Low Risk' reference), scanAgeC, shortID (factor), and the outcome column.
statNames <- c("LowRisk.Mean", "HighRisk.Mean", "LowRisk.SE", "HighRisk.SE",
               "Add.Grp.B", "Add.Grp.p", "Add.Age.B", "Add.Age.p",
               "Itx.Grp.B", "Itx.Grp.p", "Itx.Age.B", "Itx.Age.p", "Itx.B", "Itx.p", "n",
               "AdjLow.Mean", "AdjHigh.Mean", "AdjLow.SE", "AdjHigh.SE")
emptyStats <- setNames(rep(NA_real_, length(statNames)), statNames)

fitModels <- function(d, y) {
  d$y <- d[[y]]
  tryCatch({
    add <- summary(lme(y ~ Risk + scanAgeC, random = ~ 1 | shortID, data = d, control = lmeCtl))$tTable
    mi  <- lme(y ~ Risk * scanAgeC, random = ~ 1 | shortID, data = d, control = lmeCtl)
    itx <- summary(mi)$tTable
    # Model-adjusted group means (what emmeans(model, "Grp") returned) at the sample's average age,
    # from the interaction model's fixed effects: order is (Intercept), High Risk, age, High Risk:age
    xbar <- mean(d$scanAgeC)
    beta <- fixef(mi)
    V <- mi$varFix
    adj <- function(L) c(sum(L * beta), sqrt(drop(t(L) %*% V %*% L)))
    adjLow  <- adj(c(1, 0, xbar, 0))
    adjHigh <- adj(c(1, 1, xbar, xbar))
    lo <- d$y[d$Risk == "Low Risk"]
    hi <- d$y[d$Risk == "High Risk"]
    se <- function(x) sd(x) / sqrt(length(x))
    c(LowRisk.Mean = mean(lo), HighRisk.Mean = mean(hi), LowRisk.SE = se(lo), HighRisk.SE = se(hi),
      Add.Grp.B = add["RiskHigh Risk", "Value"], Add.Grp.p = add["RiskHigh Risk", "p-value"],
      Add.Age.B = add["scanAgeC", "Value"],      Add.Age.p = add["scanAgeC", "p-value"],
      Itx.Grp.B = itx["RiskHigh Risk", "Value"], Itx.Grp.p = itx["RiskHigh Risk", "p-value"],
      Itx.Age.B = itx["scanAgeC", "Value"],      Itx.Age.p = itx["scanAgeC", "p-value"],
      Itx.B = itx["RiskHigh Risk:scanAgeC", "Value"], Itx.p = itx["RiskHigh Risk:scanAgeC", "p-value"],
      n = nrow(d),
      AdjLow.Mean = adjLow[1], AdjHigh.Mean = adjHigh[1], AdjLow.SE = adjLow[2], AdjHigh.SE = adjHigh[2])
  }, error = function(e) emptyStats)
}

# Interaction p only (cheaper; used by the permutation test)
fitItxP <- function(d, y) {
  d$y <- d[[y]]
  tryCatch(summary(lme(y ~ Risk * scanAgeC, random = ~ 1 | shortID, data = d,
                       control = lmeCtl))$tTable["RiskHigh Risk:scanAgeC", "p-value"],
           error = function(e) NA_real_)
}

adjustP <- function(df) {
  for (p in c("Add.Grp.p", "Add.Age.p", "Itx.Grp.p", "Itx.Age.p", "Itx.p")) {
    df[[paste0(p, ".Adj")]] <- p.adjust(df[[p]], method = "fdr")
  }
  df
}

stackStats <- function(rows, ids) {
  rows <- lapply(rows, function(r) if (is.numeric(r) && length(r) == length(statNames)) r else emptyStats)
  df <- data.frame(ROI = ids, do.call(rbind, rows), row.names = NULL)
  adjustP(df)
}

reportSig <- function(res, label) {
  cat("\n====", label, "====\n")
  cat("Group difference at average age (additive), FDR < .05:\n")
  print(subset(res, Add.Grp.p.Adj < .05)[, c("ROI", "Add.Grp.B", "Add.Grp.p", "Add.Grp.p.Adj")])
  cat("Group difference at average age (interaction model), FDR < .05:\n")
  print(subset(res, Itx.Grp.p.Adj < .05)[, c("ROI", "Itx.Grp.B", "Itx.Grp.p", "Itx.Grp.p.Adj")])
  cat("Group x age interaction, FDR < .05:\n")
  print(subset(res, Itx.p.Adj < .05)[, c("ROI", "Itx.B", "Itx.p", "Itx.p.Adj")])
  nObs <- sum(res$Itx.p < .05, na.rm = TRUE)
  nTot <- sum(!is.na(res$Itx.p))
  cat(sprintf("Group x age p < .05 (uncorrected): %d observed vs %.1f expected by chance (of %d)\n",
              nObs, nTot * .05, nTot))
}

# --- anatomical (one row per scan) data builders ---
getAnatData <- function(base, valMat, roi, denom = NULL) {
  d <- base
  d$Value <- valMat[match(d$longID, rownames(valMat)), roi]
  ok <- is.finite(d$Value)
  if (!is.null(denom)) {
    d$ValueNormed <- d$Value / denom[match(d$longID, names(denom))]
    ok <- ok & is.finite(d$ValueNormed)
  }
  d <- d[ok, ]
  d$scanAgeC <- if (CENTER_AGE) d$scanAge - mean(d$scanAge) else d$scanAge
  d$shortID <- factor(d$shortID)
  d
}

# Fits every ROI in valMat; returns list(raw=, normed=) (normed only when denom is given)
analyzeAnatomy <- function(base, valMat, denom = NULL) {
  rois <- colnames(valMat)
  out <- pLapply(rois, function(roi) {
    d <- getAnatData(base, valMat, roi, denom)
    r <- list(raw = fitModels(d, "Value"))
    if (!is.null(denom)) r$normed <- fitModels(d, "ValueNormed")
    r
  })
  res <- list(raw = stackStats(lapply(out, function(o) o$raw), rois))
  if (!is.null(denom)) res$normed <- stackStats(lapply(out, function(o) o$normed), rois)
  res
}


# ---- 3. Behavioral QA: which assessments exist for which child ---------------

assessmentNames <- c("Bayley", "ITSEA", "MCHAT", "Vineland", "SRS", "CBCL", "BRIEF")
behavData <- setNames(lapply(assessmentNames, function(a)
  read.csv(file.path(baseDir, "behavioral", paste0(a, ".csv")))), assessmentNames)

allID <- unique(unlist(lapply(behavData, function(d) d$study_id)))
assessmentSummary <- data.frame(ID = allID, Study = NA_character_, TimePoint = NA)
for (a in assessmentNames) {                     # later assessments overwrite Study/TimePoint, as before
  d <- behavData[[a]]
  idx <- match(allID, d$study_id)
  has <- !is.na(idx)
  assessmentSummary$Study[has] <- as.character(d$study[idx[has]])
  if (a != "BRIEF") assessmentSummary$TimePoint[has] <- d$timepoint[idx[has]]   # BRIEF timepoint intentionally skipped
  assessmentSummary[[a]] <- has
}
# write.table(assessmentSummary, file.path(baseDir, "behavioral", "AssessmentSummary.csv"), col.names = T, row.names = F, sep = ",")

# Super-wide merge of all behavioral data (with a presence flag for each assessment)
for (a in assessmentNames) behavData[[a]][[if (a == "Vineland") "VINELAND" else a]] <- TRUE
assessmentDataSuperWide <- Reduce(function(x, y) merge(x, y, all = TRUE), behavData)
# write.table(assessmentDataSuperWide, file.path(baseDir, "behavioral", "AssessmentSummarySuperWide.csv"), col.names = T, row.names = F, sep = ",")


# ---- 4. Anatomical sample (matched M-CHAT Pass / Fail) ------------------------

testData <- read.csv(file.path(baseDir, "anatomical", "combinedASDData.csv"), stringsAsFactors = FALSE)
testFetal <- subset(testData, PretermStatus == "Term" & CHDStatus != "CHD" & AnyFetal == TRUE)

# One row per scan
fetalMchatMelt <- melt(testFetal, measure.vars = c("ScanAge1", "ScanAge2"),
                       id.vars = c("study_id", "mchat_result", "study"))
fetalMchatMelt <- subset(fetalMchatMelt, mchat_result != "" & value != "")
rownames(fetalMchatMelt) <- seq_len(nrow(fetalMchatMelt))
colnames(fetalMchatMelt) <- c("study_id", "mchat_result", "study", "scanTimepoint", "scanAge")

# 1:1 match of Fail to Pass on scan age and study
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

with(fetalMatched, {
  cat("Children (Pass / Fail):", length(unique(shortID[mchat_result == "Pass"])),
      length(unique(shortID[mchat_result == "Fail"])), "\n")
  cat("Scans    (Pass / Fail):", sum(mchat_result == "Pass"), sum(mchat_result == "Fail"), "\n")
})
print(t.test(scanAge ~ mchat_result, data = fetalMatched))

# --- Image file availability ---
tissueFiles   <- list.files(tissueDir)
regionalFiles <- list.files(regionalDir)
tissueScanIDs   <- substr(tissueFiles,   1, nchar(tissueFiles)   - 16)    # e.g. fetus_00556_01
regionalScanIDs <- substr(regionalFiles, 1, nchar(regionalFiles) - 18)
cat("Matched scans missing tissue file:  ", setdiff(fetalMatched$longID, tissueScanIDs), "\n")
cat("Matched scans missing regional file:", setdiff(fetalMatched$longID, regionalScanIDs), "\n")

# Timepoints available per child (all term, non-CHD fetal M-CHAT scans)
testFetalIDs <- makeShortID(testFetal$study_id, testFetal$study)
anatTimePoints <- data.frame(ID = testFetalIDs, mchat_result = testFetal$mchat_result,
  tissueFiles   = vapply(testFetalIDs, function(i) sum(substr(tissueScanIDs,   1, nchar(tissueScanIDs)   - 3) == i), numeric(1)),
  regionalFiles = vapply(testFetalIDs, function(i) sum(substr(regionalScanIDs, 1, nchar(regionalScanIDs) - 3) == i), numeric(1)),
  row.names = NULL)

# --- Demographics for the anatomical sample ---
demoDataPaper <- merge(fetalMatched, testData[, c("study_id", "Sex", "BirthGA", "X18MLTFAge", "birthWeight", "Race")],
                       all.x = TRUE, all.y = FALSE, by = "study_id")
compareGroups(demoDataPaper, "BirthGA")
compareGroups(demoDataPaper, "birthWeight")
compareGroups(demoDataPaper, "X18MLTFAge")
print(chisq.test(table(demoDataPaper$mchat_result, demoDataPaper$Sex)[, c("m", "f")]))

# Race table: cell corrections below were hard-coded in the original -- verify they still apply
try({
  contingency_table <- table(demoDataPaper$mchat_result, demoDataPaper$Race)
  contingency_table[1, 2] <- 8
  contingency_table[2, 1] <- 0
  contingency_table[2, 8] <- 3
  contingency_table <- contingency_table[1:2, 2:9]
  print(chisq.test(contingency_table))
})

# Analysis base: one row per scan, read once and reused by every anatomical analysis
makeBase <- function(exclude) {
  b <- data.frame(longID = fetalMatched$longID, shortID = fetalMatched$shortID,
                  scanAge = as.numeric(fetalMatched$scanAge), mchat_result = fetalMatched$mchat_result,
                  stringsAsFactors = FALSE)
  b <- b[!(b$longID %in% exclude), ]
  b$Risk <- factor(ifelse(b$mchat_result == "Fail", "High Risk", "Low Risk"), levels = c("Low Risk", "High Risk"))
  b
}
base <- makeBase(excludeScans)
tissueBase <- makeBase(excludeTissue)
ctBase <- makeBase(excludeCT)


# ---- 5. Tissue volumes -------------------------------------------------------

tissueRaw <- read.csv(file.path(baseDir, "anatomical", "tissueVolumes_12122025_compiled.csv"))
tissueCols <- c("CorticalGrayMatter", "CorticalWhiteMatter", "Cerebellum", "SubcorticalGray",
                "Brainstem", "AmygdalaHippocampus")
tissueTBV <- setNames(rowSums(tissueRaw[, tissueCols], na.rm = TRUE), tissueRaw$longID)
tissueMat <- as.matrix(tissueRaw[, tissueCols])
rownames(tissueMat) <- tissueRaw$longID
storage.mode(tissueMat) <- "numeric"

tissueRes <- analyzeAnatomy(tissueBase, tissueMat, tissueTBV)
reportSig(tissueRes$raw,    "Tissue volumes (raw)")
reportSig(tissueRes$normed, "Tissue volumes (normalized to TBV)")


# ---- 6. Regional volumes -----------------------------------------------------

regionalPaths <- file.path(regionalDir, regionalFiles)
regionalLabels <- read.table(file.path(baseDir, "anatomical", "regionalROILabels.txt"))

regionalAll <- as.data.frame(matrix(NA_real_, nrow = nrow(fetalMatched), ncol = 87))
rownames(regionalAll) <- fetalMatched$longID
colnames(regionalAll) <- t(regionalLabels)

readRegionalRow <- function(id) {
  f <- regionalPaths[regionalScanIDs == id]
  if (length(f) == 0) return(NULL)
  tmp <- read.table(f[1], nrows = 2)
  as.numeric(as.character(unlist(tmp[2, 3:ncol(tmp)])))
}
regionalRows <- pLapply(rownames(regionalAll), readRegionalRow)
for (i in seq_along(regionalRows)) {
  v <- regionalRows[[i]]
  if (!is.null(v)) regionalAll[i, seq_along(v)] <- v
}
regionalAll[] <- lapply(regionalAll, as.numeric)

# Total brain volume uses the original label ranges, then drop CSF / WM / cerebellar / CC labels
regionalTBV <- setNames(rowSums(regionalAll[, c(1:48, 51:82)], na.rm = TRUE), rownames(regionalAll))
regionalData <- regionalAll[, 1:82]
for (drop in c("WM", "LatVen", "CRB", "CorpusCallosum")) {
  regionalData <- regionalData[, grep(drop, colnames(regionalData), invert = TRUE)]
}

if (INCLUDE_CEREBELLAR_PARCELLATION) {
  vermisVol <- read.csv(here("_Earlier_Iteration_NeurodevelopmentalRisk", "data", "anatomicalData",
                             "vermisParcVolumes_071026.csv"), stringsAsFactors = FALSE)
  vermisVol$totalVermis <- vermisVol$supVerm + vermisVol$infVerm + vermisVol$midVerm
  regionalData$fullID <- rownames(regionalData)
  regionalData <- merge(regionalData, vermisVol, all = TRUE, by = "fullID")
  rownames(regionalData) <- regionalData$fullID
  regionalData$fullID <- NULL
  warning("Vermis TBV adjustment is not carried over; regionalTBV is unchanged.")
}
regionalMat <- as.matrix(regionalData)

cat("Scans with no regional volume file:", sum(rowSums(!is.na(regionalMat)) == 0), "\n")

regionalRes <- analyzeAnatomy(base, regionalMat, regionalTBV)
reportSig(regionalRes$raw,    "Regional volumes (raw)")
reportSig(regionalRes$normed, "Regional volumes (normalized to TBV)")

# Global check: more group x age interactions than chance? Permutation test accounts for
# correlated ROIs: shuffle risk labels across CHILDREN (scans stay together), recount hits.
if (RUN_PERMUTATION) {
  observed <- sum(regionalRes$normed$Itx.p < .05, na.rm = TRUE)
  permData <- lapply(setNames(colnames(regionalMat), colnames(regionalMat)),
                     function(r) getAnatData(base, regionalMat, r, regionalTBV))
  childRisk <- unique(base[, c("shortID", "Risk")])

  set.seed(2026)
  perms <- replicate(nPerm, setNames(sample(as.character(childRisk$Risk)), childRisk$shortID), simplify = FALSE)
  permHits <- unlist(pLapply(perms, function(rk) {
    sum(vapply(permData, function(d) {
      d$Risk <- factor(rk[as.character(d$shortID)], levels = levels(base$Risk))
      p <- fitItxP(d, "ValueNormed")
      !is.na(p) && p < .05
    }, logical(1)))
  }))
  cat(sprintf("Permutation p (>= %d ROIs with Itx.p < .05): %.4f\n", observed,
              (sum(permHits >= observed) + 1) / (nPerm + 1)))
}

# --- Figures: fitted trajectories by risk group ---
figROIs <- c(HPC_L = "Hippocampus (Left)", HPC_R = "Hippocampus (Right)",
             Amyg_L = "Amygdala (Left)", FusiformAnt_GM_R = "Fusiform (Right)",
             Caudate_R = "Caudate (Right)", Subthalamic_L = "Subthalamic (Left)",
             Lentiform_L = "Pallidum (Left)")

plotTrajectory <- function(roi, label) {
  d <- getAnatData(base, regionalMat, roi, regionalTBV)
  m <- lme(ValueNormed ~ Risk * scanAgeC, random = ~ 1 | shortID, data = d, control = lmeCtl)
  shift <- if (CENTER_AGE) mean(d$scanAge) else 0
  newd <- expand.grid(scanAgeC = seq(min(d$scanAgeC), max(d$scanAgeC), length.out = 50),
                      Risk = levels(d$Risk))
  newd$pred <- predict(m, newd, level = 0)
  newd$scanAge <- newd$scanAgeC + shift
  row <- regionalRes$normed[regionalRes$normed$ROI == roi, ]
  ggplot() +
    geom_point(data = d, aes(scanAge, ValueNormed, color = Risk), alpha = .5, size = 2) +
    geom_line(data = newd, aes(scanAge, pred, color = Risk), linewidth = 1.5) +
    scale_color_manual(values = c("#4a7fc1", "#c24a4a"), name = NULL) +
    labs(x = "Scan age (weeks)", y = paste(label, "volume / TBV"),
         subtitle = sprintf("group x age: p = %.3f, FDR q = %.2f", row$Itx.p, row$Itx.p.Adj)) +
    theme_classic(base_size = 13) + theme(legend.position = "bottom")
}
trajPlots <- Map(plotTrajectory, names(figROIs), figROIs)

#tiff(file.path(figDir, "regionalTrajectories.tiff"), height = 8, width = 16, units = "in", res = 600, compression = "lzw")
ggarrange(plotlist = trajPlots, nrow = 2, ncol = 4, common.legend = TRUE, legend = "bottom")
#dev.off()


# ---- 7. Cortical thickness ---------------------------------------------------

ctDataRaw <- read.csv(file.path(baseDir, "anatomical", "CT_Ave_Thickness.csv"), header = FALSE)
ctIDList  <- read.csv(file.path(baseDir, "anatomical", "CT_SubList.1D"), header = FALSE)
rownames(ctDataRaw) <- substr(ctIDList$V1, 1, nchar(ctIDList$V1) - 14)
ctROILabel <- read.csv(file.path(baseDir, "anatomical", "CTROILabels.txt"), header = FALSE)
colnames(ctDataRaw) <- ctROILabel$V1
ctMat <- as.matrix(ctDataRaw[, c(1:70, 79:90)])
globalCT <- setNames(as.numeric(t(read.csv(file.path(baseDir, "anatomical", "CT_Global_Thickness.csv"), header = FALSE))),
                     rownames(ctDataRaw))

ctRes <- analyzeAnatomy(ctBase, ctMat, globalCT)
reportSig(ctRes$raw,    "Cortical thickness (raw)")
reportSig(ctRes$normed, "Cortical thickness (normalized to global CT)")


# ---- 8. Functional sample ----------------------------------------------------

asdDemo <- read.csv(file.path(baseDir, "functional", "asdDemo_expand.csv"), stringsAsFactors = FALSE)
tdDemo  <- read.csv(file.path(baseDir, "functional", "tdDemo.csv"), stringsAsFactors = FALSE)

asdDemo <- asdDemo[asdDemo$ID != "rescue_0187_01", ]   # (old note: why is this excluded?)
rownames(asdDemo) <- seq_len(nrow(asdDemo))

asdDemo$Grp <- "ASD"
tdDemo$Grp  <- "TD"
combDemo <- rbind(asdDemo, tdDemo)

# Match TD to ASD on scan GA. (Kept exactly as in the .qmd: control rows are looked up in tdDemo by
# match.matrix row names. Check that this pulls the intended TD rows; NA rows are dropped below.)
matchFrame <- MatchIt::matchit(as.factor(Grp) ~ ScanGA, data = combDemo)
tdDemoMatch <- tdDemo[matchFrame$match.matrix, ]
tdDemoMatch <- tdDemoMatch[complete.cases(tdDemoMatch), ]
tdDemoMatch <- tdDemoMatch[order(as.numeric(row.names(tdDemoMatch))), ]

dropPtList <- c("none")
keepRows <- as.numeric(rownames(asdDemo[!(asdDemo$ID %in% dropPtList), ]))
tdDemoFinal  <- tdDemoMatch[keepRows, ]
asdDemoFinal <- asdDemo[keepRows, ]

cat("Functional sample scan GA: TD", mean(tdDemoFinal$ScanGA), sd(tdDemoFinal$ScanGA),
    "| ASD", mean(asdDemoFinal$ScanGA), sd(asdDemoFinal$ScanGA), "\n")
print(t.test(tdDemoFinal$ScanGA, asdDemoFinal$ScanGA))
print(chisq.test(rbind(table(factor(tdDemoFinal$Sex,  levels = c("m", "f"))),
                       table(factor(asdDemoFinal$Sex, levels = c("m", "f"))))))

tdBirthAge  <- c(40.00, 39.00, 40.29, 39.86, 39.71, 38.71, 40.43, 39.86, 38.88, 40.57, 38.43, 39.71, 39.71, 38.43, 39.43, 39.00, 40.14, 38.14, 38.71)
asdBirthAge <- c(37.71, 39.57, 38.86, 38.86, 39.29, 38.29, 38.14, 39.71, 38.00, 39.00, 39.00, 41.00, 40.86, 38.86, 38.86, 38.00, 38.00, 39.43, 37.29)
print(t.test(tdBirthAge, asdBirthAge))
cat("Birth GA: TD", mean(tdBirthAge), sd(tdBirthAge), "| ASD", mean(asdBirthAge), sd(asdBirthAge), "\n")

# Final exclusions after matching
tdDemoFinal  <- tdDemoFinal[!(tdDemoFinal$ID %in% c("fetus_00748_02", "fetus_00556_01")), ]
asdDemoFinal <- asdDemoFinal[asdDemoFinal$ID != "pl_0028_02", ]

# Demographics keyed to the main data file
fcDemoTemp <- rbind(asdDemoFinal, tdDemoFinal)
idParts <- strsplit(fcDemoTemp$ID, "_", fixed = TRUE)
idNum <- vapply(idParts, function(p) as.numeric(p[2]), numeric(1))
fcDemoTemp$study_id <- ifelse(substr(fcDemoTemp$ID, 1, 1) == "f", as.character(idNum),
                       ifelse(substr(fcDemoTemp$ID, 1, 1) == "r", paste0("RES-", sprintf("%03d", idNum)),
                                                                  paste0("PL-",  sprintf("%03d", idNum))))
testDataFC <- read.csv(file.path(baseDir, "functional", "combinedASDData.csv"), stringsAsFactors = FALSE)
demoDataPaperFC <- merge(fcDemoTemp, testDataFC[, c("study_id", "Sex", "BirthGA", "X18MLTFAge", "birthWeight", "Race", "mchat_result")],
                         all.x = TRUE, all.y = FALSE, by = "study_id")
print(chisq.test(table(demoDataPaperFC$mchat_result, demoDataPaperFC$Sex.x)))
compareGroups(demoDataPaperFC, "birthWeight")
compareGroups(demoDataPaperFC, "BirthGA")
print(chisq.test(table(demoDataPaperFC$mchat_result, demoDataPaperFC$Race)))

# Functional data
allTDFC  <- read.csv(file.path(baseDir, "functional", "UseableData", "tdCorrSummary.csv"), row.names = 1)
allASDFC <- read.csv(file.path(baseDir, "functional", "UseableData", "asdCorrSummary_Expand.csv"), row.names = 1)
rownames(allTDFC) <- substr(rownames(allTDFC), 1, nchar(rownames(allTDFC)) - 28)
tdFC  <- as.matrix(allTDFC[rownames(allTDFC) %in% tdDemoFinal$ID, ])
asdFC <- as.matrix(allASDFC[rownames(allASDFC) %in% asdDemoFinal$ID, ])
rm(allTDFC, allASDFC)
stopifnot(identical(colnames(tdFC), colnames(asdFC)))

# Overlap with the anatomical sample
cat("Functional ASD scans not in anatomical sample:", setdiff(asdDemoFinal$ID, fetalMatched$longID), "\n")
cat("Anatomical Fail scans not in functional sample:",
    setdiff(fetalMatched$longID[fetalMatched$mchat_result == "Fail"], asdDemoFinal$ID), "\n")


# ---- 9. Whole-brain connectivity (every connection fit once) ------------------

tdD  <- tdDemoFinal[match(rownames(tdFC),  tdDemoFinal$ID), ]
asdD <- asdDemoFinal[match(rownames(asdFC), asdDemoFinal$ID), ]
fcDemo <- rbind(data.frame(ID = tdD$ID,  ScanGA = tdD$ScanGA,  Risk = "Low Risk",  stringsAsFactors = FALSE),
                data.frame(ID = asdD$ID, ScanGA = asdD$ScanGA, Risk = "High Risk", stringsAsFactors = FALSE))
fcDemo$Risk <- factor(fcDemo$Risk, levels = c("Low Risk", "High Risk"))
fcDemo$shortID <- factor(substr(fcDemo$ID, 1, nchar(fcDemo$ID) - 3))
fcDemo$scanAgeC <- if (CENTER_AGE) fcDemo$ScanGA - mean(fcDemo$ScanGA) else fcDemo$ScanGA
fcMat <- rbind(tdFC, asdFC)

fitConnection <- function(i) {
  y <- fcMat[, i]
  if (anyNA(y)) return(emptyStats)             # connections with a missing value are skipped, as before
  d <- fcDemo
  d$FC <- y
  fitModels(d, "FC")
}
fcRows <- pLapply(seq_len(ncol(fcMat)), fitConnection)
fcAll <- stackStats(fcRows, colnames(fcMat))
fcAll$TD.Mean  <- fcAll$LowRisk.Mean
fcAll$ASD.Mean <- fcAll$HighRisk.Mean

cat("\nWhole-brain: connections with group difference (interaction model) FDR < .05:",
    sum(fcAll$Itx.Grp.p.Adj < .05, na.rm = TRUE), "\n")
cat("Whole-brain: connections with group difference (additive model) FDR < .05:",
    sum(fcAll$Add.Grp.p.Adj < .05, na.rm = TRUE), "\n")
cat("Whole-brain: connections with group x age FDR < .05:",
    sum(fcAll$Itx.p.Adj < .05, na.rm = TRUE), "\n")

if (WRITE_OUTPUT) {
  write.table(fcAll, file.path(outDir, "WholeBrain_WholeTable.csv"), col.names = TRUE, row.names = FALSE, sep = ",")
  write.table(subset(fcAll, Itx.Grp.p.Adj < .05), file.path(outDir, "WholeBrain_P05FDR.csv"),
              col.names = TRUE, row.names = FALSE, sep = ",")
}
sigConn <- fcAll$ROI[which(fcAll$Itx.Grp.p.Adj < .05)]


# ---- 10. Connectogram and lobe-by-lobe summary --------------------------------

fullFCList <- read.csv(file.path(baseDir, "functional", "100PerHemiLabels_Decimal.1D"), header = FALSE)[[1]]
NetList <- read.csv(file.path(baseDir, "functional", "r200_Net.csv"))
ROIOrder <- NetList[order(NetList$OrderNew), ]$Name_Abrv

# Sector colors by lobe, in the order of the 200 labels
roiColorSpec <- data.frame(
  name = c("ACG.L","IFGtriang.L","MFG.L","ORBinf.L","ORBmid.L","ORBsupb.L","PreCG.L","REC.L","SFGdor.L","SFGmed.L",
           "INS.L","IPL.L","MCG.L","PCG.L","PCUN.L","PoCG.L","SMG.L","SPG.L",
           "FFG.L","ITG.L","MTG.L","STG.L","TPOmid.L",
           "CAL.L","CUN.L","LING.L","IOG.L","MOG.L","SOG.L",
           "AMY.L","CAU.L","PAL.L","PHG.L","PUT.L","THA.L",
           "CRB.L","BS","CRB.R",
           "THA.R","PUT.R","PHG.R","PAL.R","HIP.R","AMY.R",
           "SOG.R","MOG.R","IOG.R","LING.R","CUN.R","CAL.R",
           "TPOsup.R","TPOmid.R","STG.R","MTG.R","ITG.R","FFG.R",
           "SPG.R","SMG.R","PoCG.R","PCUN.R","PCL.R","MCG.R","IPL.R","INS.R","HES.R",
           "SMA.R","SFGmed.R","SFGdor.R","PreCG.R","ORBsupb.R","ORBmid.R","ORBmed.R","ORBinf.R","MFG.R","IFGtriang.R","IFGoperc.R","ACG.R"),
  n    = c(3,3,4,2,2,2,5,1,6,2,
           4,3,4,1,3,3,2,4,
           1,3,5,4,1,
           2,3,2,2,3,1,
           1,3,1,2,1,2,
           5,8,5,
           3,3,1,1,1,1,
           2,3,2,2,2,2,
           1,1,3,4,3,2,
           4,1,6,3,2,3,4,2,1,
           2,2,4,3,4,1,1,2,4,2,1,2),
  col  = c(rep("#efa973", 10), rep("#d9c97e", 8), rep("#a9dd8c", 5), rep("#79f09f", 6), rep("#89d4c1", 6),
           "#96a5d7", "#a67de8", "#96a5d7",
           rep("#89d4c1", 6), rep("#79f09f", 6), rep("#a9dd8c", 6), rep("#d9c97e", 9), rep("#efa973", 12)),
  stringsAsFactors = FALSE)
ROIColor <- setNames(rep(roiColorSpec$col, roiColorSpec$n), rep(roiColorSpec$name, roiColorSpec$n))
stopifnot(length(ROIColor) == 200)

# Every ordered ROI pair (i = 'To', j = 'From'), vectorized (was a 40,000-iteration loop with a subset() each)
nLab <- length(fullFCList)
chord <- expand.grid(j = seq_len(nLab), i = seq_len(nLab))      # j varies fastest, as in the old loop
chord$ToAbbrev   <- fullFCList[chord$i]
chord$FromAbbrev <- fullFCList[chord$j]
chord$ConnFull   <- paste0(fullFCList[chord$i], ".", chord$i, "_", fullFCList[chord$j], ".", chord$j)
meanDiff <- (fcAll$TD.Mean - fcAll$ASD.Mean)[match(chord$ConnFull, fcAll$ROI)]
meanDiff[is.na(meanDiff)] <- 0
isSig <- chord$ConnFull %in% sigConn
chord$ASDConn <- as.numeric(meanDiff < 0)
chord$TDConn  <- as.numeric(meanDiff > 0)
chord$ASDCol <- ifelse(chord$ASDConn == 1 & isSig, "#a60000aa", "#00000000")
chord$TDCol  <- ifelse(chord$TDConn  == 1 & isSig, "#0050bfaa", "#00000000")

gapSizes <- c(rep(2, 9), 5, rep(2, 7), 5, rep(2, 4), 5, rep(2, 5), 5, rep(2, 5), 5, 5, 5, 5,
              rep(2, 5), 5, rep(2, 5), 5, rep(2, 5), 5, rep(2, 8), 5, rep(2, 11), 13)

drawChord <- function(valueCol, colCol) {
  circos.clear()
  circos.par(start.degree = 81, gap.after = gapSizes)
  chordDiagram(chord[, c("ToAbbrev", "FromAbbrev", valueCol)],
               order = ROIOrder, grid.col = ROIColor, grid.border = "#000000",
               col = chord[[colCol]], symmetric = TRUE, scale = TRUE,
               link.border = chord[[colCol]], link.lwd = 3,
               annotationTrack = c("grid"), annotationTrackHeight = c(0.03, 0.03))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(-.25, .5))
  }, bg.border = NA)
}

# tiff(file.path(figDir, "chordPlot_ASD_newColors.tiff"), units = "in", height = 8, width = 8, res = 600, compression = "lzw")
drawChord("ASDConn", "ASDCol")
# dev.off()
# tiff(file.path(figDir, "chordPlot_TD_newColors.tiff"), units = "in", height = 8, width = 8, res = 600, compression = "lzw")
drawChord("TDConn", "TDCol")
# dev.off()

print(sort(table(c(chord$ToAbbrev[chord$ASDCol == "#a60000aa"], chord$FromAbbrev[chord$ASDCol == "#a60000aa"]))))
print(sort(table(c(chord$ToAbbrev[chord$TDCol  == "#0050bfaa"], chord$FromAbbrev[chord$TDCol  == "#0050bfaa"]))))

# Connection counts between lobes
lobeNames <- c("Frontal", "Parietal", "Temporal", "Occipital", "Subcortical", "Cerebellum", "Brainstem")
connRoiName <- function(x) {                       # 'CAU.L.39' -> 'CAU.L'; brainstem labels -> 'BS'
  nm <- vapply(strsplit(x, ".", fixed = TRUE), function(p) paste(p[1:2], collapse = "."), character(1))
  ifelse(substr(nm, 1, 2) == "BS", "BS", nm)
}
connParts <- strsplit(fcAll$ROI, "_", fixed = TRUE)
roi1 <- vapply(connParts, `[`, character(1), 1)
roi2 <- vapply(connParts, `[`, character(1), 2)
lobe1 <- match(NetList$Lobe[match(connRoiName(roi1), NetList$Name_Abrv)], lobeNames)
lobe2 <- match(NetList$Lobe[match(connRoiName(roi2), NetList$Name_Abrv)], lobeNames)

lobeCounts <- function(keep) {                     # symmetric lobe x lobe counts (diagonal counted twice, as before)
  tb <- table(factor(lobe1[keep], levels = seq_along(lobeNames)), factor(lobe2[keep], levels = seq_along(lobeNames)))
  m <- as.matrix(tb) + t(as.matrix(tb))
  dimnames(m) <- list(lobeNames, lobeNames)
  as.data.frame(m)
}
isFDR  <- fcAll$ROI %in% sigConn
tdHigh <- fcAll$TD.Mean > fcAll$ASD.Mean
allConnDesc        <- lobeCounts(rep(TRUE, nrow(fcAll)) & !is.na(lobe1) & !is.na(lobe2))
allConnDescSig     <- lobeCounts(isFDR & !is.na(lobe1) & !is.na(lobe2))
allConnDescSig.TD  <- lobeCounts(isFDR & tdHigh  & !is.na(lobe1) & !is.na(lobe2))
allConnDescSig.ASD <- lobeCounts(isFDR & !tdHigh & !is.na(lobe1) & !is.na(lobe2))

plotLobeTable <- function(counts, high, limits = NULL) {
  m <- counts
  m$ROWNAMES <- rownames(counts)
  long <- melt(m, id.vars = "ROWNAMES")
  colnames(long) <- c("Row", "Col", "Num")
  long$Row <- factor(long$Row, levels = lobeNames)
  long$Col <- factor(long$Col, levels = rev(lobeNames))
  ggplot(long, aes(x = Row, y = Col)) +
    scale_fill_gradient(low = "#ffffff", high = high, limits = limits) +
    geom_tile(aes(fill = Num), color = "#000000") +
    geom_text(aes(label = Num), color = "#000000") +
    theme(axis.title.x = element_text(color = "black", size = 16),
          axis.text.x = element_text(color = "black", size = 12), axis.ticks.x = element_blank(),
          axis.title.y = element_text(color = "black", size = 16),
          axis.text.y = element_text(color = "black", size = 12), axis.ticks.y = element_blank(),
          axis.line.x = element_blank(), panel.grid.major = element_blank(),
          panel.background = element_rect(fill = "white"), legend.position = "right",
          legend.text = element_text(color = "black", size = 12), legend.title = element_blank())
}
# tiff(file.path(figDir, "tdWholeBrainTableFig.tiff"), width = 8.95, height = 8, units = "in", res = 600, compression = "lzw")
plotLobeTable(allConnDescSig.TD, "#0050bf")
# dev.off()
# tiff(file.path(figDir, "asdWholeBrainTableFig.tiff"), width = 8.95, height = 8, units = "in", res = 600, compression = "lzw")
plotLobeTable(allConnDescSig.ASD, "#a60000", limits = c(0, 25))
# dev.off()


# ---- 11. Network analyses (subsets of the whole-brain fit, re-corrected within network) ----

roiNum <- function(x) vapply(strsplit(x, ".", fixed = TRUE), function(p) suppressWarnings(as.numeric(p[3])), numeric(1))
num1 <- roiNum(roi1)
num2 <- roiNum(roi2)
inNums <- function(nums) num1 %in% nums & num2 %in% nums

networks <- list(
  CRB      = grepl("CRB", roi1) & grepl("CRB", roi2),
  Limbic   = inNums(c(5, 9, 20, 29, 56, 61, 67, 68, 81, 86, 94, 108, 120, 133, 141, 143, 155, 160, 188, 196, 198, 200,
                      43, 86, 125)),      # includes parahippocampal gyrus
  DMN      = inNums(c(4, 11, 18, 22, 23, 34, 40, 42, 48, 53, 70, 85, 88, 96, 100, 102, 104, 106, 123, 124, 128, 137,
                      150, 152, 156, 163, 168, 171, 176, 178)),       # hippocampus / parahippocampal removed
  SomMot   = grepl("PreCG|PoCG|SMA", roi1) & grepl("PreCG|PoCG|SMA", roi2),
  SubCort  = inNums(c(31, 32, 39, 51, 80, 92, 97, 107, 111, 134, 139, 145, 151)),
  Auditory = inNums(c(88, 94, 2, 12, 48, 171, 107, 139, 152, 185)),
  Visual   = inNums(c(181, 174, 153, 114, 149, 195, 167, 105, 116, 194, 197, 162, 131, 186, 135,
                      72, 60, 19, 52, 35, 75, 25, 15, 42, 30, 91, 87, 98, 57)))

networkFits <- lapply(names(networks), function(nm) {
  net <- fcAll[networks[[nm]], ]
  net <- adjustP(net)                               # FDR within this network only

  # additive-model table (old fcGLM) and interaction-model table (old fcIntxGLM)
  add <- data.frame(ROI = net$ROI, TD.Mean = net$TD.Mean, TD.SE = net$LowRisk.SE,
                    ASD.Mean = net$ASD.Mean, ASD.SE = net$HighRisk.SE,
                    Grp.B = net$Add.Grp.B, Grp.p = net$Add.Grp.p, GA.B = net$Add.Age.B, GA.p = net$Add.Age.p,
                    Grp.p.adjust = net$Add.Grp.p.Adj, GA.p.adjust = net$Add.Age.p.Adj)
  itx <- data.frame(ROI = net$ROI, TD.Mean = net$AdjLow.Mean, TD.SE = net$AdjLow.SE,
                    ASD.Mean = net$AdjHigh.Mean, ASD.SE = net$AdjHigh.SE,
                    TD.Mean.Raw = net$TD.Mean, ASD.Mean.Raw = net$ASD.Mean,
                    TD.SE.Raw = net$LowRisk.SE, ASD.SE.Raw = net$HighRisk.SE,
                    Grp.B = net$Itx.Grp.B, Grp.p = net$Itx.Grp.p, GA.B = net$Itx.Age.B, GA.p = net$Itx.Age.p,
                    Intx.B = net$Itx.B, Intx.p = net$Itx.p,
                    Grp.p.adjust = net$Itx.Grp.p.Adj, GA.p.adjust = net$Itx.Age.p.Adj, Intx.p.adjust = net$Itx.p.Adj)

  cat(sprintf("\n== %s network (%d connections) ==\n", nm, nrow(net)))
  cat("Group difference FDR < .05 (additive):", sum(add$Grp.p.adjust < .05, na.rm = TRUE),
      "| (interaction model):", sum(itx$Grp.p.adjust < .05, na.rm = TRUE),
      "| group x age:", sum(itx$Intx.p.adjust < .05, na.rm = TRUE), "\n")

  if (WRITE_OUTPUT) {
    wt <- function(df, suffix) write.table(df, file.path(outDir, paste0(nm, "_", suffix, ".csv")),
                                           col.names = TRUE, row.names = FALSE, sep = ",")
    wt(add, "WholeTable")
    wt(subset(add, Grp.p < .05), "P05")
    wt(subset(add, Grp.p < .01), "P01")
    wt(subset(add, Grp.p < .001), "P001")
    wt(subset(add, Grp.p.adjust < .05), "P05FDR")
    wt(subset(itx, Grp.p.adjust < .05 | Intx.p.adjust < .05), "Intx_P05FDR")
    wt(subset(itx, Grp.p.adjust < .1 | Intx.p.adjust < .1), "Intx_P10FDR")
    wt(itx, "Intx_Wholebrain")
  }
  list(add = add, itx = itx)
})
names(networkFits) <- names(networks)

# Connection-by-connection bar plots, ordered by group p-value
barTheme <- theme(axis.title.x = element_text(color = "black", size = 16), axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(), axis.title.y = element_text(color = "black", size = 16),
                  axis.text.y = element_text(color = "black", size = 12), axis.ticks.y = element_blank(),
                  axis.line.x = element_line(color = "black", linewidth = 1),
                  axis.line.y = element_line(color = "black", linewidth = 1),
                  panel.grid.major = element_line(color = "dark gray"),
                  panel.background = element_rect(fill = "white"), legend.position = "none")

plotNetworkBars <- function(add, label, ylim = NULL) {
  long <- melt(add, id.vars = "ROI", measure.vars = c("TD.Mean", "ASD.Mean"), variable.name = "Group", value.name = "FC")
  pOrder <- add$ROI[order(add$Grp.p, decreasing = TRUE)]
  p <- ggplot(long, aes(x = factor(ROI, levels = pOrder), y = FC, group = factor(Group))) +
    scale_fill_manual(values = c("#d88080cc", "#8cb0e3cc")) +
    scale_color_manual(values = c("#bb0000", "#0000bb")) +
    geom_hline(yintercept = 0, linewidth = 1.5) +
    geom_bar(stat = "identity", position = "dodge", aes(color = Group, fill = Group)) +
    ylab("Connectivity Strength") + xlab(paste(label, "Connections")) + barTheme
  if (!is.null(ylim)) p <- p + scale_y_continuous(limits = ylim, breaks = seq(-2, 2, .2), expand = c(0, 0))
  p
}
plotNetworkBars(networkFits$CRB$add,     "Cerebellar",  c(0, 1.2))
plotNetworkBars(networkFits$Limbic$add,  "Limbic",      c(-.4, 1.4))
plotNetworkBars(networkFits$DMN$add,     "Default Mode", c(-.6, 1))
plotNetworkBars(networkFits$SomMot$add,  "Sensorimotor", c(-.3, 1.1))
plotNetworkBars(networkFits$SubCort$add, "Subcortical", c(-.2, 1.2))
plotNetworkBars(networkFits$Auditory$add, "Auditory")
plotNetworkBars(networkFits$Visual$add,  "Visual")


# ---- 12. Trend-level connections across the four main networks -----------------

trendAll <- do.call(rbind, lapply(c("DMN", "Limbic", "SubCort", "CRB"), function(nm) {
  itx <- networkFits[[nm]]$itx
  out <- subset(itx, Grp.p.adjust < .1 | Intx.p.adjust < .1)
  if (nrow(out)) out$Net <- nm
  out
}))

plotOrder <- c("CAU.L.39_THA.R.139", "PAL.L.97_PUT.R.111", "CAU.L.39_THA.L.80", "THA.L.51_PUT.L.92",
               "PAL.L.97_PAL.R.107", "PAL.L.97_THA.R.139", "PAL.R.107_PUT.R.145", "TPOmid.L.9_ITG.R.188",
               "REC.L.20_TPOmid.R.108", "REC.L.20_ORBsupb.R.155", "ITG.L.56_ITG.R.188", "SFGdor.L.67_AMY.R.133",
               "SFGdor.L.67_ORBsupb.R.196", "ORBsupb.R.143_ITG.R.188", "ITG.R.188_ORBsupb.R.196",
               "ITG.R.188_ORBsupb.R.198", "SFGmed.L.18_ACG.R.123", "MTG.L.22_ACG.L.34", "MTG.L.22_ACG.R.156",
               "PCUN.L.23_STG.L.88", "MTG.L.96_STG.R.102", "SFGmed.R.128_MTG.R.150", "PHG.L.43_PCUN.R.178",
               "CRB.L.13_CRB.R.142", "CRB.L.13_CRB.R.192", "CRB.L.17_CRB.L.99", "CRB.L.17_CRB.R.192",
               "CRB.L.54_CRB.R.192", "CRB.L.99_CRB.R.142", "CRB.L.99_CRB.R.192", "CRB.R.129_CRB.R.192")

trendFrame <- rbind(
  data.frame(Connection = trendAll$ROI, Mean = trendAll$TD.Mean,  SE = trendAll$TD.SE,  Grp = "Low Risk"),
  data.frame(Connection = trendAll$ROI, Mean = trendAll$ASD.Mean, SE = trendAll$ASD.SE, Grp = "High Risk"))
trendFrame$Connection <- factor(trendFrame$Connection, levels = plotOrder)
trendFrame$Grp <- factor(trendFrame$Grp, levels = c("Low Risk", "High Risk"))

# Note: the old plot silently dropped (NA) any trend connection that is not in plotOrder.
cat("Trend connections not in plotOrder:", setdiff(trendAll$ROI, plotOrder), "\n")

ggplot(data = trendFrame, aes(x = Connection, y = Mean, group = Grp)) +
  scale_fill_manual(values = c("#d88080cc", "#8cb0e3cc")) +
  scale_color_manual(values = c("#bb0000", "#0000bb")) +
  geom_hline(yintercept = 0, linewidth = 1) +
  geom_bar(stat = "identity", position = position_dodge2(0.2), aes(color = Grp, fill = Grp)) +
  ylab("Connectivity Strength") + xlab("Connections") +
  theme(axis.title.x = element_text(color = "black", size = 16),
        axis.text.x = element_text(color = "black", size = 4, angle = 90), axis.ticks.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 16),
        axis.text.y = element_text(color = "black", size = 12), axis.ticks.y = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 1),
        axis.line.y = element_line(color = "black", linewidth = 1),
        panel.grid.major = element_line(color = "#aaaaaa"), panel.background = element_rect(fill = "white"),
        legend.position = "right", legend.title = element_blank())
