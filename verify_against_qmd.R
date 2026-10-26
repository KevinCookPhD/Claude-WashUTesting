# verify_against_qmd.R
#
# Checks that EarlyAutismMarkers_updated.R reproduces the results of the original
# EarlyAutismMarkers.qmd BEFORE you trust any of the new numbers.
#
# The new script differs from the old on purpose in two ways, and this script undoes both:
#   * scan age centering  -> run with CENTER_AGE <- FALSE
#   * reference group     -> old: Fail (ASD) is the reference, so Grp.B = Pass - Fail
#                            new: Low Risk (Pass) is the reference, so Grp.B = High - Low
#     The mapping below accounts for this (sign flips; the "scan age" slope in an
#     interaction model is the reference group's slope, so old = new age slope + interaction).
#
# STEP 1 (old code). BACK UP the old functional outputs first, then run the old .qmd
# through the end of the Regional, Tissue, Cortical Thickness, and Whole Brain / network
# chunks, and run this in that same session:
#
#   file.copy(file.path(baseDir, "functional", "output"), baseDir, recursive = TRUE)   # then rename to output_old
#   oldResults <- list(
#     fetalMatched  = fetalMatched,  tdDemoFinal = tdDemoFinal, asdDemoFinal = asdDemoFinal,
#     tissueRaw     = list(add = glmTissueVolTable,       itx = glmTissueVolTable.Itx),
#     tissueNormed  = list(add = glmTissueNormedVolTable, itx = glmTissueNormedVolTable.Itx),
#     regRaw        = list(add = glmRegionalVolTable,       itx = glmRegionalVolTable.Itx),
#     regNormed     = list(add = glmRegionalNormedVolTable, itx = glmRegionalNormedVolTable.Itx),
#     ctRaw         = list(add = glmCTTable,       itx = glmCTTable.Itx),
#     ctNormed      = list(add = glmCTNormedTable, itx = glmCTNormedTable.Itx))
#   saveRDS(oldResults, file.path(baseDir, "old_results.rds"))
#
# STEP 2 (new code). Edit the two paths below if needed, then source this file.
# It runs the new script in "match the old code" mode (no centering, old exclusions,
# output written to a separate folder) and prints a comparison.


# ---- comparison helpers ------------------------------------------------------

# Compare one quantity between an old and a new table (rows matched on ROI).
compareOne <- function(label, oldVals, newVals, relTol = 1e-4) {
  both <- !is.na(oldVals) & !is.na(newVals)
  naMismatch <- sum(is.na(oldVals) != is.na(newVals))
  d <- abs(oldVals[both] - newVals[both])
  scale <- pmax(abs(oldVals[both]), abs(newVals[both]), 1e-12)
  rel <- d / scale
  data.frame(quantity = label, n = sum(both), naMismatch = naMismatch,
             maxAbsDiff = if (length(d)) max(d) else NA,
             maxRelDiff = if (length(rel)) max(rel) else NA,
             ok = naMismatch == 0 && (length(rel) == 0 || max(rel) < relTol))
}

# old, new: data frames with an ROI column. map: list(label = list(oldFun, newFun)).
compareTables <- function(old, new, map, title, relTol = 1e-4) {
  common <- intersect(old$ROI, new$ROI)
  cat(sprintf("\n--- %s: %d ROIs in old, %d in new, %d in common ---\n",
              title, length(unique(old$ROI)), length(unique(new$ROI)), length(common)))
  if (length(setdiff(old$ROI, new$ROI)) || length(setdiff(new$ROI, old$ROI))) {
    cat("  only in old:", head(setdiff(old$ROI, new$ROI), 8), "\n")
    cat("  only in new:", head(setdiff(new$ROI, old$ROI), 8), "\n")
  }
  o <- old[match(common, old$ROI), ]
  n <- new[match(common, new$ROI), ]
  res <- do.call(rbind, lapply(names(map), function(k)
    compareOne(k, map[[k]][[1]](o), map[[k]][[2]](n), relTol)))
  res[, c("maxAbsDiff", "maxRelDiff")] <- signif(res[, c("maxAbsDiff", "maxRelDiff")], 3)
  print(res, row.names = FALSE)
  invisible(res)
}

# Old interaction-model table <-> new stats table (new = fitModels() output with Low Risk reference)
mapItx <- function(hasMeans = TRUE) {
  m <- list(
    Grp.B   = list(function(o) o$Grp.B,   function(n) -n$Itx.Grp.B),
    Grp.p   = list(function(o) o$Grp.p,   function(n) n$Itx.Grp.p),
    Age.B   = list(function(o) o$scanGA.B, function(n) n$Itx.Age.B + n$Itx.B),   # old reference = Fail slope
    Itx.B   = list(function(o) o$Itx.B,   function(n) -n$Itx.B),
    Itx.p   = list(function(o) o$Itx.p,   function(n) n$Itx.p),
    Grp.p.Adj = list(function(o) o$Grp.p.Adj, function(n) n$Itx.Grp.p.Adj),
    Itx.p.Adj = list(function(o) o$Itx.p.Adj, function(n) n$Itx.p.Adj))
  if (hasMeans) m <- c(list(
    TD.Mean  = list(function(o) o$TD.Mean,  function(n) n$LowRisk.Mean),
    ASD.Mean = list(function(o) o$ASD.Mean, function(n) n$HighRisk.Mean)), m)
  m
}
# Old additive-model table <-> new stats table
mapAdd <- function(hasGrpB = TRUE) {
  m <- list(
    Grp.p      = list(function(o) o$Grp.p,     function(n) n$Add.Grp.p),
    Age.B      = list(function(o) o$scanGA.B,  function(n) n$Add.Age.B),
    Age.p      = list(function(o) o$scanGA.p,  function(n) n$Add.Age.p),
    Grp.p.Adj  = list(function(o) o$Grp.p.Adj, function(n) n$Add.Grp.p.Adj),
    Age.p.Adj  = list(function(o) o$scanGA.p.Adj, function(n) n$Add.Age.p.Adj))
  if (hasGrpB) m <- c(list(Grp.B = list(function(o) o$Grp.B, function(n) -n$Add.Grp.B)), m)
  m
}

# Old whole-brain CSV (fcIntxGLM) <-> new fcAll
mapWholeBrain <- list(
  TD.Mean  = list(function(o) o$TD.Mean,  function(n) n$LowRisk.Mean),
  ASD.Mean = list(function(o) o$ASD.Mean, function(n) n$HighRisk.Mean),
  TD.SE    = list(function(o) o$TD.SE,    function(n) n$LowRisk.SE),
  ASD.SE   = list(function(o) o$ASD.SE,   function(n) n$HighRisk.SE),
  Grp.B    = list(function(o) o$Grp.B,    function(n) -n$Itx.Grp.B),
  Grp.p    = list(function(o) o$Grp.p,    function(n) n$Itx.Grp.p),
  GA.B     = list(function(o) o$GA.B,     function(n) n$Itx.Age.B + n$Itx.B),
  Intx.B   = list(function(o) o$Intx.B,   function(n) -n$Itx.B),
  Intx.p   = list(function(o) o$Intx.p,   function(n) n$Itx.p),
  nIncluded = list(function(o) o$nIncluded, function(n) n$n),
  Grp.p.adjust  = list(function(o) o$Grp.p.adjust,  function(n) n$Itx.Grp.p.Adj),
  Intx.p.adjust = list(function(o) o$Intx.p.adjust, function(n) n$Itx.p.Adj))

compareSets <- function(title, old, new) {
  a <- sort(unique(old)); b <- sort(unique(new))
  cat(sprintf("%-34s old %4d | new %4d | only old %3d | only new %3d  %s\n", title, length(a), length(b),
              length(setdiff(a, b)), length(setdiff(b, a)),
              if (identical(a, b)) "IDENTICAL" else "DIFFERENT"))
  if (!identical(a, b)) {
    cat("   only in old:", head(setdiff(a, b), 10), "\n")
    cat("   only in new:", head(setdiff(b, a), 10), "\n")
  }
}


# ---- main (skipped when this file is sourced just for its helper functions) ----

if (!exists("VERIFY_HELPERS_ONLY")) {

  library(here)
  baseDir     <- here("data")
  oldResults  <- readRDS(file.path(baseDir, "old_results.rds"))
  oldOutDir   <- file.path(baseDir, "functional", "output_old")     # backup of the OLD outputs
  newOutDir   <- file.path(baseDir, "functional", "output_new")     # new script writes here

  # Run the new script in "match the old code" mode
  newEnv <- new.env()
  newEnv$CENTER_AGE    <- FALSE                                          # old code did not center
  newEnv$WRITE_OUTPUT  <- TRUE
  newEnv$outDir        <- newOutDir
  newEnv$excludeScans  <- c("rescue_0175_01", "rescue_0184_01")          # old regional table code
  newEnv$excludeTissue <- c("rescue_0175_01", "fetus_02178_01", "rescue_0184_01")   # old tissue code
  newEnv$excludeCT     <- character(0)                                   # old CT code excluded nothing
  # Where the new script lives. Set newScriptPath before sourcing this file to point anywhere, e.g.
  #   newScriptPath <- "C:/Users/Kevin/Downloads/EarlyAutismMarkers_updated.R"
  newScript <- if (exists("newScriptPath")) newScriptPath else
    c(here("code", "EarlyAutismMarkers_updated.R"), here("EarlyAutismMarkers_updated.R"))[
      file.exists(c(here("code", "EarlyAutismMarkers_updated.R"), here("EarlyAutismMarkers_updated.R")))][1]
  if (is.na(newScript)) newScript <- here("code", "EarlyAutismMarkers_updated.R")
  if (!file.exists(newScript)) stop("Can't find the new script at: ", newScript,
                                    "\nDownload it from the GitHub branch (or set newScriptPath).")
  source(newScript, local = newEnv)

  cat("\n\n=================== VERIFICATION ===================\n")
  cat("(ok = same to within 0.01% after accounting for the intentional sign/reference changes)\n")

  # 1. Samples
  cat("\n--- Samples ---\n")
  compareSets("Anatomical scans (longID)", oldResults$fetalMatched$longID, newEnv$fetalMatched$longID)
  compareSets("  Fail scans", oldResults$fetalMatched$longID[oldResults$fetalMatched$mchat_result == "Fail"],
              newEnv$fetalMatched$longID[newEnv$fetalMatched$mchat_result == "Fail"])
  compareSets("Functional ASD IDs", oldResults$asdDemoFinal$ID, newEnv$asdDemoFinal$ID)
  compareSets("Functional TD IDs",  oldResults$tdDemoFinal$ID,  newEnv$tdDemoFinal$ID)

  # 2. Anatomical models (old: raw / normed x additive / interaction)
  for (nm in c("regNormed", "regRaw", "tissueNormed", "tissueRaw", "ctNormed", "ctRaw")) {
    newRes <- switch(nm,
      regNormed = newEnv$regionalRes$normed, regRaw = newEnv$regionalRes$raw,
      tissueNormed = newEnv$tissueRes$normed, tissueRaw = newEnv$tissueRes$raw,
      ctNormed = newEnv$ctRes$normed, ctRaw = newEnv$ctRes$raw)
    compareTables(oldResults[[nm]]$itx, newRes, mapItx(), paste(nm, "- interaction model"))
    # old regional raw additive table assigned Grp.p to Grp.B by mistake (stray '<-'); skip Grp.B there
    compareTables(oldResults[[nm]]$add, newRes, mapAdd(hasGrpB = nm != "regRaw"), paste(nm, "- additive model"))
  }

  # 3. Whole-brain connectivity
  oldWB <- read.csv(file.path(oldOutDir, "WholeBrain_WholeTable.csv"), stringsAsFactors = FALSE)
  compareTables(oldWB, newEnv$fcAll, mapWholeBrain, "Whole-brain connectivity")

  # 4. Networks (old CSVs written by the old code)
  for (net in names(newEnv$networkFits)) {
    f1 <- file.path(oldOutDir, paste0(net, "_WholeTable.csv"))
    f2 <- file.path(oldOutDir, paste0(net, "_Intx_Wholebrain.csv"))
    if (file.exists(f1)) {
      oldAdd <- read.csv(f1, stringsAsFactors = FALSE)
      newAdd <- newEnv$networkFits[[net]]$add
      compareTables(oldAdd, newAdd, list(
        Grp.B = list(function(o) o$Grp.B, function(n) -n$Grp.B),
        Grp.p = list(function(o) o$Grp.p, function(n) n$Grp.p),
        GA.B  = list(function(o) o$GA.B,  function(n) n$GA.B),
        GA.p  = list(function(o) o$GA.p,  function(n) n$GA.p),
        TD.Mean = list(function(o) o$TD.Mean, function(n) n$TD.Mean),
        Grp.p.adjust = list(function(o) o$Grp.p.adjust, function(n) n$Grp.p.adjust)),
        paste(net, "network - additive"))
    }
    if (file.exists(f2)) {
      oldItx <- read.csv(f2, stringsAsFactors = FALSE)
      newItx <- newEnv$networkFits[[net]]$itx
      m <- list(
        Grp.B = list(function(o) o$Grp.B, function(n) -n$Grp.B),
        Grp.p = list(function(o) o$Grp.p, function(n) n$Grp.p),
        GA.B  = list(function(o) o$GA.B,  function(n) n$GA.B + n$Intx.B),
        Intx.B = list(function(o) o$Intx.B, function(n) -n$Intx.B),
        Intx.p = list(function(o) o$Intx.p, function(n) n$Intx.p),
        Grp.p.adjust  = list(function(o) o$Grp.p.adjust,  function(n) n$Grp.p.adjust),
        Intx.p.adjust = list(function(o) o$Intx.p.adjust, function(n) n$Intx.p.adjust))
      # old Limbic interaction table used RAW means; every other network used model-adjusted (emmeans) means
      suf <- if (net == "Limbic") ".Raw" else ""
      meanMap <- list(
        TD.Mean  = list(function(o) o$TD.Mean,  function(n) n[[paste0("TD.Mean", suf)]]),
        ASD.Mean = list(function(o) o$ASD.Mean, function(n) n[[paste0("ASD.Mean", suf)]]),
        TD.SE    = list(function(o) o$TD.SE,    function(n) n[[paste0("TD.SE", suf)]]),
        ASD.SE   = list(function(o) o$ASD.SE,   function(n) n[[paste0("ASD.SE", suf)]]))
      m <- c(meanMap, m)
      compareTables(oldItx, newItx, m, paste(net, "network - interaction"))
    }
  }

  cat("\nAnything marked ok = FALSE needs a look before you trust the new numbers.\n",
      "Common legitimate causes: a different sample (section 1), or old code that differs by design\n",
      "(see the header of EarlyAutismMarkers_updated.R).\n")
}
