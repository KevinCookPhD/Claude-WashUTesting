# PMA CENTERING CHECK ----
#
# How much do the results depend on where postmenstrual age (PMA = scan age) is centered?
#
# Whenever a model has a group x age interaction, the group term is the group difference AT THE AGE WHERE THE AGE VARIABLE IS ZERO.
# With no centering that is age 0 (a fetus nobody scanned); with mean centering it is the average scan age. This script refits the
# volume and functional models under five centerings and compares them:
#     none  = scan age as recorded (the current analysis)
#     mean  = centered at the mean scan age of the analysis sample
#     wk20  = centered at 20 weeks,  wk28 = centered at 28 weeks (start of the third trimester),  wk40 = centered at 40 weeks
# (the weeks assume scan age is recorded in weeks)
#
# What changes with centering: only the GROUP term of the interaction models (its estimate, p-value and so the FDR result).
# The interaction term, the age term and the main-effects models do not change (the script checks the interaction term).
#
# HOW TO RUN: first run EarlyAutismMarkers_edit.R (the cached fits make that quick) so that the data and helper functions exist in
# the R session, then run this whole file. Slow fits are cached in cacheDir, so a second run is quick.
#
# OUTPUT: a comparison table printed in the console, and the figure figDir/PMACheck/PMACheck_network_interactions
# (one panel per significant or trending connection, for every network).

# Settings ----

# Which centerings to run (names below); the first one is the reference the others are compared with
centerings <- list(none = 0, mean = NA, wk20 = 20, wk28 = 28, wk40 = 40)      # NA = the sample mean
centeringLabels <- c(none = 'No centering', mean = 'Mean centered', wk20 = '20 weeks', wk28 = '28 weeks', wk40 = '40 weeks', main = 'No interaction term')   # 'main' = the model without an interaction (same under any centering)

# The whole-brain models (19,900 connections) take long: about as long as the original run for each centering. TRUE adds them.
includeWholeBrain <- FALSE

# Which centering decides what counts as significant or trending in the figure
figureSelectionCentering <- 'none'

# File types for the PMACheck figure (it is large, so a tiff at 600 dpi would be huge: pdf is the default). 'tiff', 'pdf', 'svg', 'eps'
pmaFigureTypes <- c('pdf')

# Region of significance (see the end of the script): at which ages do the groups differ? The groups are compared at every half week
# between the youngest age present in BOTH groups and the oldest age present in BOTH groups.
# 'fdr' = p-values are FDR-corrected across the connections (or ROIs) of a network at each age; 'none' = uncorrected
regionOfSignificanceMethod <- 'fdr'
regionOfSignificanceAlpha <- .05
regionMaxPanels <- 24                     # most panels drawn per network in the region-of-significance figure

# Networks to analyze: the a priori networks set in EarlyAutismMarkers_edit.R (networksOfInterest); all of networkSpecs if that is missing
pmaNetworks <- if (exists('networksOfInterest')) networksOfInterest else names(networkSpecs)

# Make sure the data and functions from EarlyAutismMarkers_edit.R are available ----

pmaNeeded <- c('fetalMatched', 'tissueData', 'tissueDataRaw', 'tissueExclude', 'regionalData', 'regionalDataRaw', 'regionalExclude',
               'prepVolumeAnalysis', 'fitRoiModels', 'runRoiModels', 'addFDR', 'addFDRfc',
               'tdFC', 'asdFC', 'tdDemoFinal', 'asdDemoFinal', 'connROIs', 'networkSpecs', 'fitNetworkModels',
               'cached', 'showFigure', 'connectionLabel', 'groupColors', 'figDir')
pmaMissing <- pmaNeeded[!vapply(pmaNeeded, exists, logical(1))]
if (length(pmaMissing) > 0) {
  stop('Run EarlyAutismMarkers_edit.R first. These objects are missing: ', paste(pmaMissing, collapse = ', '))
}
if (includeWholeBrain && !all(vapply(c('runWholeBrainModels', 'nCores'), exists, logical(1)))) {
  stop('The whole-brain models need runWholeBrainModels and nCores from EarlyAutismMarkers_edit.R')
}
if (!all(names(centerings) %in% names(centeringLabels))) { stop('Every centering needs a label in centeringLabels') }

# Fit the models under each centering ----

# x minus the center (NA = the mean of x)
centerBy <- function(x, center) x - ifelse(is.na(center), mean(x, na.rm = TRUE), center)

# Volume models (tissue and regional): scan age is centered in each ROI's analysis data (mean = that ROI's sample mean)
fitVolumeCentering <- function(label) {
  center <- centerings[[label]]
  centered <- function(prepFun) function(i) { d <- prepFun(i); d$scanAge <- centerBy(d$scanAge, center); d }
  cached(paste0('pmacheck_volume_', label),
         inputs = list(tissueData, tissueDataRaw, regionalData, regionalDataRaw, fetalMatched, tissueExclude, regionalExclude, center,
                       deparse(prepVolumeAnalysis), deparse(fitRoiModels), deparse(runRoiModels)),
         compute = function() {
           list(tissue = runRoiModels(colnames(tissueData),
                                      centered(function(i) prepVolumeAnalysis(tissueData[,i, drop=FALSE], fetalMatched, tissueDataRaw, excludeIDs=tissueExclude)),
                                      rawVar='Volume', normVar='VolNormed'),
                regional = runRoiModels(colnames(regionalData),
                                        centered(function(i) prepVolumeAnalysis(regionalData[,i, drop=FALSE], fetalMatched, regionalDataRaw,
                                                                                excludeIDs=regionalExclude, dropInf=TRUE)),
                                        rawVar='Volume', normVar='VolNormed'))
         })
}

# Functional demographics with scan age centered (mean = the mean of both groups together)
centeredDemos <- function(label) {
  center <- centerings[[label]]
  if (is.na(center)) { center <- mean(c(tdDemoFinal$ScanGA, asdDemoFinal$ScanGA), na.rm = TRUE) }
  td <- tdDemoFinal;  td$ScanGA <- td$ScanGA - center
  asd <- asdDemoFinal;  asd$ScanGA <- asd$ScanGA - center
  list(td = td, asd = asd)
}

# Network models for every network
fitNetworkCentering <- function(label) {
  demos <- centeredDemos(label)
  fits <- lapply(pmaNetworks, function(net) {
    both <- networkSpecs[[net]]$match(connROIs$ROI1) & networkSpecs[[net]]$match(connROIs$ROI2)
    if (sum(both) == 0) { return(NULL) }
    tdNetwork <- tdFC[, both, drop = FALSE]
    asdNetwork <- asdFC[, both, drop = FALSE]
    models <- cached(paste0('pmacheck_network_', net, '_', label),
                     inputs = list(tdNetwork, asdNetwork, demos$td, demos$asd, deparse(fitNetworkModels)),
                     compute = function() fitNetworkModels(tdNetwork, asdNetwork, demos$td, demos$asd))
    list(fcGLM = addFDRfc(models$fcGLM), fcIntxGLM = addFDRfc(models$fcIntxGLM))
  })
  names(fits) <- pmaNetworks
  fits[!vapply(fits, is.null, logical(1))]
}

# Whole-brain models (optional)
fitWholeBrainCentering <- function(label) {
  demos <- centeredDemos(label)
  models <- cached(paste0('pmacheck_wholebrain_', label),
                   inputs = list(tdFC, asdFC, demos$td, demos$asd, deparse(fitWholeBrainChunk)),
                   compute = function() runWholeBrainModels(tdFC, asdFC, demos$td, demos$asd, nCores))
  list(fcGLM = addFDRfc(models$fcGLM), fcIntxGLM = addFDRfc(models$fcIntxGLM))
}

message('Fitting volume models...')
volumeFits <- lapply(setNames(names(centerings), names(centerings)), function(label) {
  message('  ', centeringLabels[[label]]); fitVolumeCentering(label)
})
message('Fitting network models...')
networkFits <- lapply(setNames(names(centerings), names(centerings)), function(label) {
  message('  ', centeringLabels[[label]]); fitNetworkCentering(label)
})
if (includeWholeBrain) {
  message('Fitting whole-brain models...')
  wholeBrainFits <- lapply(setNames(names(centerings), names(centerings)), fitWholeBrainCentering)
}

# Gather the interaction-model tables: one list per analysis, one table per centering ----
# (each table: group term in Grp.p.adj / Grp.B, interaction term in Itx.p.adj; TD.Mean and ASD.Mean for the direction)

# itxAdj = NULL for a model without an interaction term (its interaction columns are NA)
standardTable <- function(tbl, grpAdj, itxAdj = NULL) {
  data.frame(ROI = tbl$ROI, TD.Mean = tbl$TD.Mean, ASD.Mean = tbl$ASD.Mean, Grp.B = tbl$Grp.B,
             Grp.p.adj = tbl[[grpAdj]],
             Itx.B = if (is.null(itxAdj)) NA_real_ else tbl[[intersect(c('Itx.B', 'Intx.B'), colnames(tbl))[1]]],
             Itx.p.adj = if (is.null(itxAdj)) NA_real_ else tbl[[itxAdj]],
             stringsAsFactors = FALSE)
}

analyses <- list()
volumeFamilies <- list('Tissue, raw volume' = c('tissue', 'RawItx'), 'Tissue, normalized' = c('tissue', 'NormedItx'),
                       'Regional, raw volume' = c('regional', 'RawItx'), 'Regional, normalized' = c('regional', 'NormedItx'))
for (fam in names(volumeFamilies)) {
  analyses[[fam]] <- lapply(volumeFits, function(f) standardTable(addFDR(f[[volumeFamilies[[fam]][1]]][[volumeFamilies[[fam]][2]]]), 'Grp.p.Adj', 'Itx.p.Adj'))
}
for (net in names(networkFits[[1]])) {
  analyses[[paste('Network:', net)]] <- lapply(networkFits, function(f) standardTable(f[[net]]$fcIntxGLM, 'Grp.p.adjust', 'Intx.p.adjust'))
}
if (includeWholeBrain) {
  analyses[['Whole brain']] <- lapply(wholeBrainFits, function(f) standardTable(f$fcIntxGLM, 'Grp.p.adjust', 'Intx.p.adjust'))
}

# The model WITHOUT an interaction term (FC or volume ~ Group + age). Centering does not change it, so one copy is added to every analysis
# as the extra column 'main'. (Volume group means in it are the raw means; network group means are raw means too.)
for (fam in names(volumeFamilies)) {
  mainName <- sub('Itx', '', volumeFamilies[[fam]][2])
  analyses[[fam]]$main <- standardTable(addFDR(volumeFits[[1]][[volumeFamilies[[fam]][1]]][[mainName]]), 'Grp.p.Adj')
}
for (net in names(networkFits[[1]])) {
  analyses[[paste('Network:', net)]]$main <- standardTable(networkFits[[1]][[net]]$fcGLM, 'Grp.p.adjust')
}
if (includeWholeBrain) { analyses[['Whole brain']]$main <- standardTable(wholeBrainFits[[1]]$fcGLM, 'Grp.p.adjust') }
allModels <- c(names(centerings), 'main')

# Comparison table ----

countsAt <- function(p, cutoff) sum(p < cutoff, na.rm = TRUE)
summaryRows <- do.call(rbind, lapply(names(analyses), function(an) {
  do.call(rbind, lapply(allModels, function(label) {
    tbl <- analyses[[an]][[label]]
    trending <- which(tbl$Grp.p.adj < .10)
    # a group effect whose beta points the opposite way from the group means
    opposes <- sum((tbl$Grp.B[trending] > 0) != (tbl$TD.Mean[trending] > tbl$ASD.Mean[trending]), na.rm = TRUE)
    data.frame(analysis = an, centering = label, n = nrow(tbl),
               group05 = countsAt(tbl$Grp.p.adj, .05), group10 = countsAt(tbl$Grp.p.adj, .10),
               interaction05 = countsAt(tbl$Itx.p.adj, .05), interaction10 = countsAt(tbl$Itx.p.adj, .10),
               betaOpposesMeans = opposes, stringsAsFactors = FALSE)
  }))
}))

# One cell per analysis x centering, written "count at p.adj < .05 / count at p.adj < .10"
widen <- function(first, second = NULL, valueLabel) {
  out <- data.frame(analysis = unique(summaryRows$analysis), stringsAsFactors = FALSE)
  out$tested <- summaryRows$n[match(out$analysis, summaryRows$analysis)]
  for (label in allModels) {
    rows <- summaryRows[summaryRows$centering == label, ]
    rows <- rows[match(out$analysis, rows$analysis), ]
    cells <- if (is.null(second)) as.character(rows[[first]]) else paste(rows[[first]], '/', rows[[second]])
    cells[is.na(rows[[first]])] <- '-'                 # the no-interaction model has no interaction term
    out[[centeringLabels[[label]]]] <- cells
  }
  out
}

oldWidth <- options(width = 220)       # keep each table on one line; restored below
cat('\n=== PMA centering check ===\n')
cat('\nGROUP TERM: number of ROIs/connections with FDR p < .05 / < .10 (this is the term that changes with centering)\n')
groupSummary <- widen('group05', 'group10')
print(groupSummary, row.names = FALSE)

cat('\nINTERACTION TERM: number with FDR p < .05 / < .10 (should be identical in every column)\n')
interactionSummary <- widen('interaction05', 'interaction10')
print(interactionSummary, row.names = FALSE)

cat('\nDIRECTION: among ROIs/connections with group FDR p < .10, how many have a group beta that points the OPPOSITE way from the\n',
    'group means (TD.Mean vs ASD.Mean). Mean centering should bring this to (near) 0; large numbers mean the beta is describing an age far from the data.\n', sep = '')
directionSummary <- widen('betaOpposesMeans')
print(directionSummary, row.names = FALSE)

centeredRows <- summaryRows[summaryRows$centering != 'main', ]
interactionChanged <- any(centeredRows$interaction05 != ave(centeredRows$interaction05, centeredRows$analysis, FUN = function(x) x[1]) |
                          centeredRows$interaction10 != ave(centeredRows$interaction10, centeredRows$analysis, FUN = function(x) x[1]))
cat('\nCheck: the interaction term', ifelse(interactionChanged, 'DID change with centering (unexpected; look into it).', 'is unchanged by centering, as expected.'), '\n')

options(oldWidth)

# Which ROIs/connections change group-effect status (significant < .05, trending < .10, not significant) with centering ----

statusOf <- function(p) ifelse(!is.na(p) & p < .05, 'sig', ifelse(!is.na(p) & p < .10, 'trend', 'ns'))
statusChanges <- do.call(rbind, lapply(names(analyses), function(an) {
  tables <- analyses[[an]]
  status <- sapply(allModels, function(label) statusOf(tables[[label]]$Grp.p.adj))
  if (is.null(dim(status))) { status <- matrix(status, nrow = 1, dimnames = list(NULL, allModels)) }
  changed <- apply(status[, names(centerings), drop = FALSE], 1, function(x) length(unique(x)) > 1)      # change across the centerings
  if (!any(changed)) { return(NULL) }
  out <- data.frame(analysis = an, ROI = tables[[1]]$ROI[changed], stringsAsFactors = FALSE)
  for (label in allModels) { out[[paste0('status_', label)]] <- status[changed, label] }
  for (label in allModels) { out[[paste0('betaSign_', label)]] <- ifelse(tables[[label]]$Grp.B[changed] > 0, '+', '-') }
  out
}))
cat('\nROIs/connections whose group-effect status (sig = p.adj < .05, trend = < .10, ns) differs between centerings:',
    ifelse(is.null(statusChanges), 0, nrow(statusChanges)), '\n')
if (!is.null(statusChanges)) {
  print(table(statusChanges$analysis))
  cat('(the full list is in the object statusChanges; betaSign + means Low Risk higher at that centering)\n')
}

if (exists('saveTables') && isTRUE(saveTables)) {
  pmaTableDir <- file.path(tableDir, 'PMACheck')
  dir.create(pmaTableDir, showWarnings = FALSE, recursive = TRUE)
  write.csv(groupSummary, file.path(pmaTableDir, 'PMACheck_group_term_counts.csv'), row.names = FALSE)
  write.csv(interactionSummary, file.path(pmaTableDir, 'PMACheck_interaction_term_counts.csv'), row.names = FALSE)
  write.csv(directionSummary, file.path(pmaTableDir, 'PMACheck_direction_conflicts.csv'), row.names = FALSE)
  if (!is.null(statusChanges)) { write.csv(statusChanges, file.path(pmaTableDir, 'PMACheck_status_changes.csv'), row.names = FALSE) }
}

# Figure: PMA x Group interaction for the significant and trending connections of every network ----

# A connection is shown if its group or interaction term has FDR p < .10 (the same rule as the P10FDR tables), at the centering
# chosen in figureSelectionCentering. Lines are the fixed-effect fit of FC ~ Group * PMA (random intercept per participant);
# points are the individual scans. Raw scan age is on the x axis, so the figure does not depend on the centering.

pmaConnectionData <- function(roi) {
  one <- function(fc, demo) {
    d <- data.frame(ID = rownames(fc), FC = as.numeric(fc[[roi]]), stringsAsFactors = FALSE)
    merge(d, demo[, c('ID', 'ScanGA', 'Grp')], by = 'ID')
  }
  d <- rbind(one(tdFC, tdDemoFinal), one(asdFC, asdDemoFinal))
  d <- d[!is.na(d$FC), ]
  d$shortID <- substr(d$ID, 1, nchar(d$ID) - 3)
  d$Grp <- factor(d$Grp, levels = c('ASD', 'TD'))
  d
}

pValueText <- function(p) ifelse(p < .001, '<.001', sub('^0', '', sprintf('%.3f', p)))

pmaNetworkPanel <- function(net, tbl) {
  keep <- tbl[!is.na(tbl$Grp.p.adj) & (tbl$Grp.p.adj < .10 | tbl$Itx.p.adj < .10), ]
  if (nrow(keep) == 0) { return(NULL) }
  keep <- keep[order(pmin(keep$Grp.p.adj, keep$Itx.p.adj)), ]
  names <- connectionLabel(keep$ROI)
  keep$panel <- paste0(gsub('\n', ' ', names), '\nGroup p.adj = ', pValueText(keep$Grp.p.adj),
                              ', Interaction p.adj = ', pValueText(keep$Itx.p.adj),
                              ifelse(pmin(keep$Grp.p.adj, keep$Itx.p.adj) < .05, '  (significant)', '  (trending)'))
  keep$panel <- factor(keep$panel, levels = unique(keep$panel))

  lmeCtrl <- nlme::lmeControl(msMaxIter = 200, tolerance = 1e-4, msTol = 1e-5, opt = 'optim')
  points <- list(); lines <- list()
  for (i in seq_len(nrow(keep))) {
    d <- pmaConnectionData(keep$ROI[i])
    fit <- nlme::lme(FC ~ Grp * ScanGA, random = ~1 | factor(shortID), data = d, control = lmeCtrl)
    grid <- do.call(rbind, lapply(levels(d$Grp), function(g) {
      r <- range(d$ScanGA[d$Grp == g])
      data.frame(Grp = factor(g, levels = levels(d$Grp)), ScanGA = seq(r[1], r[2], length.out = 40))
    }))
    grid$FC <- predict(fit, grid, level = 0)
    d$panel <- keep$panel[i]; grid$panel <- keep$panel[i]
    points[[i]] <- d[, c('panel', 'Grp', 'ScanGA', 'FC')]
    lines[[i]] <- grid[, c('panel', 'Grp', 'ScanGA', 'FC')]
  }
  points <- do.call(rbind, points); lines <- do.call(rbind, lines)

  title <- if (exists('networkNodeLabels') && net %in% names(networkNodeLabels)) paste(networkNodeLabels[[net]], 'Network') else paste(net, 'Network')
  plotObj <- ggplot2::ggplot(mapping = ggplot2::aes(x = ScanGA, y = FC, color = Grp)) +
    ggplot2::geom_point(data = points, alpha = .45, size = 1.1) +
    ggplot2::geom_line(data = lines, linewidth = 1.1) +
    ggplot2::facet_wrap(~panel, ncol = 4, scales = 'free_y') +
    ggplot2::scale_color_manual(values = c(TD = groupColors[['low']], ASD = groupColors[['high']]),
                                labels = c(TD = 'Low Risk', ASD = 'High Risk'), name = NULL) +
    ggplot2::labs(title = paste0(title, ' (', sum(pmin(keep$Grp.p.adj, keep$Itx.p.adj) < .05), ' significant, ',
                                 sum(pmin(keep$Grp.p.adj, keep$Itx.p.adj) >= .05), ' trending)'),
                  x = 'Postmenstrual age at scan (weeks)', y = 'Connectivity strength') +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 7.5, hjust = 0), strip.background = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 16, face = 'bold'), legend.position = 'bottom',
                   panel.grid.major = ggplot2::element_line(color = '#e6e6e6'))
  list(plot = plotObj, rows = ceiling(nrow(keep) / 4))
}

selectionTables <- lapply(names(networkFits[[figureSelectionCentering]]), function(net) {
  analyses[[paste('Network:', net)]][[figureSelectionCentering]]
})
names(selectionTables) <- names(networkFits[[figureSelectionCentering]])

panels <- lapply(names(selectionTables), function(net) pmaNetworkPanel(net, selectionTables[[net]]))
names(panels) <- names(selectionTables)
panels <- panels[!vapply(panels, is.null, logical(1))]
if (length(panels) == 0) {
  message('No network has a significant or trending connection, so no PMA figure was made.')
} else {
  rowsPerNetwork <- vapply(panels, `[[`, numeric(1), 'rows')
  pmaFigure <- ggpubr::ggarrange(plotlist = lapply(panels, `[[`, 'plot'), ncol = 1, nrow = length(panels),
                                 heights = rowsPerNetwork + .6, common.legend = TRUE, legend = 'bottom')
  pmaWidth <- 15
  pmaHeight <- 2.9 * sum(rowsPerNetwork) + 0.9 * length(panels)

  oldSave <- saveFigures; oldTypes <- figureFileTypes
  saveFigures <- TRUE; figureFileTypes <- pmaFigureTypes
  showFigure(pmaFigure, 'PMACheck_network_interactions', width = pmaWidth, height = pmaHeight, subfolder = 'PMACheck')
  saveFigures <- oldSave; figureFileTypes <- oldTypes
  message('Figure saved in ', file.path(figDir, 'PMACheck'))
}

# Region of significance: at which ages do the groups differ? ----
#
# The model is  outcome ~ Low-Risk indicator * age  with a random intercept per participant, so the Low minus High Risk difference at a
# given age is  (group beta) + (interaction beta) * age.  Instead of testing it at one chosen age (what centering does), this tests it
# at every half week of age from the youngest to the oldest age present in BOTH groups, and reports the ages where it is significant.
# Nothing outside the observed ages is tested, so there is no extrapolation.
#
# Per connection/ROI this gives: the crossover age (where the difference is zero, if it falls inside the tested range) and the ages
# where Low Risk is higher and where High Risk is higher, both with the p-values uncorrected and FDR-corrected across the
# connections/ROIs of that network at each age. regionOfSignificanceMethod picks which one the figure uses.

# Low minus High Risk difference at each age for one connection/ROI.
#   d: data frame with columns y (outcome), age, low (1 = Low Risk, 0 = High Risk) and id (participant)
differenceAcrossAge <- function(d) {
  d <- d[complete.cases(d[, c('y', 'age', 'low', 'id')]), ]
  lo <- max(tapply(d$age, d$low, min)); hi <- min(tapply(d$age, d$low, max))
  if (length(unique(d$low)) < 2 || hi < lo) { return(NULL) }
  ages <- seq(ceiling(lo * 2) / 2, floor(hi * 2) / 2, by = .5)
  if (length(ages) < 2) { return(NULL) }
  d$id <- factor(d$id)
  fit <- tryCatch(nlme::lme(y ~ low * age, random = ~1 | id, data = d,
                            control = nlme::lmeControl(msMaxIter = 200, tolerance = 1e-4, msTol = 1e-5, opt = 'optim')),
                  error = function(e) NULL)
  if (is.null(fit)) { return(NULL) }
  b <- nlme::fixef(fit); v <- fit$varFix; tt <- summary(fit)$tTable
  estimate <- b[['low']] + b[['low:age']] * ages
  se <- sqrt(v['low', 'low'] + 2 * ages * v['low', 'low:age'] + ages^2 * v['low:age', 'low:age'])
  df <- tt['low', 'DF']
  crossover <- -b[['low']] / b[['low:age']]
  list(curve = data.frame(age = ages, estimate = estimate, se = se, df = df, p = 2 * pt(-abs(estimate / se), df)),
       crossover = ifelse(crossover >= min(ages) & crossover <= max(ages), crossover, NA_real_),
       range = c(min(ages), max(ages)), interactionP = tt['low:age', 'p-value'])
}

# Region of significance for a set of connections/ROIs. dataFun(roi) returns the data frame for differenceAcrossAge.
regionOfSignificance <- function(rois, dataFun) {
  results <- lapply(rois, function(roi) differenceAcrossAge(dataFun(roi)))
  names(results) <- rois
  results <- results[!vapply(results, is.null, logical(1))]
  if (length(results) == 0) { return(NULL) }
  # FDR across the connections/ROIs at each age (ages are matched by value; a connection with no data at an age does not enter it)
  allAges <- sort(unique(unlist(lapply(results, function(r) r$curve$age))))
  pMatrix <- sapply(results, function(r) r$curve$p[match(allAges, r$curve$age)])
  if (is.null(dim(pMatrix))) { pMatrix <- matrix(pMatrix, nrow = length(allAges)) }
  adjusted <- t(apply(pMatrix, 1, p.adjust, method = 'fdr'))
  if (ncol(pMatrix) == 1) { adjusted <- matrix(adjusted, ncol = 1) }
  for (j in seq_along(results)) {
    results[[j]]$curve$pAdjusted <- adjusted[match(results[[j]]$curve$age, allAges), j]
  }
  results
}

# Contiguous runs of ages as text, e.g. 'wk 31-37' or 'wk 20-22, 30-33'; '' if none
ageRunsText <- function(ages, flag) {
  if (!any(flag)) { return('') }
  hit <- ages[flag]; runStart <- c(TRUE, diff(hit) > .5 + 1e-9); runId <- cumsum(runStart)
  paste0('wk ', paste(vapply(split(hit, runId), function(a) if (length(a) == 1) format(a) else paste0(format(min(a)), '-', format(max(a))), character(1)), collapse = ', '))
}

rosTable <- function(family, results) {
  do.call(rbind, lapply(names(results), function(roi) {
    r <- results[[roi]]; cv <- r$curve
    sigNone <- cv$p < regionOfSignificanceAlpha; sigFdr <- cv$pAdjusted < regionOfSignificanceAlpha
    data.frame(analysis = family, ROI = roi, testedWeeks = paste0(r$range[1], '-', r$range[2]),
               interactionP = r$interactionP, crossoverWeek = round(r$crossover, 1),
               lowHigher_uncorrected = ageRunsText(cv$age, sigNone & cv$estimate > 0),
               highHigher_uncorrected = ageRunsText(cv$age, sigNone & cv$estimate < 0),
               lowHigher_fdr = ageRunsText(cv$age, sigFdr & cv$estimate > 0),
               highHigher_fdr = ageRunsText(cv$age, sigFdr & cv$estimate < 0),
               nAgesSig_fdr = sum(sigFdr), nAgesTested = nrow(cv), stringsAsFactors = FALSE)
  }))
}

# The data for each analysis family
rosVolumeData <- function(column, data, raw, exclude, dropInf) {
  function(roi) {
    d <- prepVolumeAnalysis(data[, roi, drop = FALSE], fetalMatched, raw, excludeIDs = exclude, dropInf = dropInf)
    data.frame(y = d$VolNormed, age = d$scanAge, low = as.numeric(d$mchat_result == 'Pass'), id = d$shortID)
  }
}
rosNetworkData <- function(roi) {
  one <- function(fc, demo) { x <- merge(data.frame(ID = rownames(fc), y = as.numeric(fc[[roi]]), stringsAsFactors = FALSE), demo[, c('ID', 'ScanGA')], by = 'ID'); x }
  td <- one(tdFC, tdDemoFinal); td$low <- 1
  asd <- one(asdFC, asdDemoFinal); asd$low <- 0
  d <- rbind(td, asd)
  data.frame(y = d$y, age = d$ScanGA, low = d$low, id = substr(d$ID, 1, nchar(d$ID) - 3))
}

message('Region of significance...')
rosResults <- list()
rosResults[['Regional, normalized']] <- cached('pmacheck_ros_regional_normalized',
  inputs = list(regionalData, regionalDataRaw, fetalMatched, regionalExclude, deparse(differenceAcrossAge), deparse(regionOfSignificance)),
  compute = function() regionOfSignificance(colnames(regionalData), rosVolumeData(NULL, regionalData, regionalDataRaw, regionalExclude, TRUE)))
for (net in names(networkFits[[1]])) {
  both <- networkSpecs[[net]]$match(connROIs$ROI1) & networkSpecs[[net]]$match(connROIs$ROI2)
  rosResults[[paste('Network:', net)]] <- cached(paste0('pmacheck_ros_network_', net),
    inputs = list(tdFC[, both, drop = FALSE], asdFC[, both, drop = FALSE], tdDemoFinal, asdDemoFinal,
                  deparse(differenceAcrossAge), deparse(regionOfSignificance), deparse(rosNetworkData)),
    compute = function() regionOfSignificance(colnames(tdFC)[both], rosNetworkData))
}
rosResults <- rosResults[!vapply(rosResults, is.null, logical(1))]

regionTable <- do.call(rbind, lapply(names(rosResults), function(an) rosTable(an, rosResults[[an]])))
rownames(regionTable) <- NULL

sigColumnLow <- ifelse(regionOfSignificanceMethod == 'fdr', 'lowHigher_fdr', 'lowHigher_uncorrected')
sigColumnHigh <- ifelse(regionOfSignificanceMethod == 'fdr', 'highHigher_fdr', 'highHigher_uncorrected')
regionSummary <- do.call(rbind, lapply(split(regionTable, regionTable$analysis), function(t) {
  lowAny <- t[[sigColumnLow]] != ''; highAny <- t[[sigColumnHigh]] != ''
  data.frame(analysis = t$analysis[1], tested = nrow(t),
             anyAgeSignificant = sum(lowAny | highAny), lowHigherOnly = sum(lowAny & !highAny), highHigherOnly = sum(highAny & !lowAny),
             bothDirections = sum(lowAny & highAny),
             medianCrossoverWeek = ifelse(all(is.na(t$crossoverWeek)), NA, median(t$crossoverWeek, na.rm = TRUE)),
             stringsAsFactors = FALSE)
}))
rownames(regionSummary) <- NULL
cat('\nREGION OF SIGNIFICANCE (', ifelse(regionOfSignificanceMethod == 'fdr', 'FDR-corrected across each network at every age', 'uncorrected'),
    ', alpha = ', regionOfSignificanceAlpha, ')\n', sep = '')
cat('anyAgeSignificant = connections/ROIs where the groups differ at at least one age; bothDirections = the difference changes sign\n',
    'with both signs significant; medianCrossoverWeek = median age where the difference is zero (inside the tested ages only).\n', sep = '')
oldWidth <- options(width = 220)
print(regionSummary, row.names = FALSE)
options(oldWidth)
cat('(per-connection ages are in the object regionTable)\n')

if (exists('saveTables') && isTRUE(saveTables)) {
  write.csv(regionTable, file.path(pmaTableDir, 'PMACheck_region_of_significance.csv'), row.names = FALSE)
  write.csv(regionSummary, file.path(pmaTableDir, 'PMACheck_region_of_significance_summary.csv'), row.names = FALSE)
}

# Figure: the Low minus High Risk difference across age, for the connections/ROIs that differ at some age
pmaRegionPanel <- function(family, results) {
  sigColumn <- ifelse(regionOfSignificanceMethod == 'fdr', 'pAdjusted', 'p')
  nSig <- vapply(results, function(r) sum(r$curve[[sigColumn]] < regionOfSignificanceAlpha), numeric(1))
  keep <- names(sort(nSig[nSig > 0], decreasing = TRUE))
  if (length(keep) == 0) { return(NULL) }
  keep <- head(keep, regionMaxPanels)
  labels <- if (grepl('^Network', family)) gsub('\n', ' ', connectionLabel(keep)) else keep
  curves <- do.call(rbind, lapply(seq_along(keep), function(i) {
    cv <- results[[keep[i]]]$curve
    cv$panel <- paste0(labels[i], '\n', ifelse(is.na(results[[keep[i]]]$crossover), 'no crossover in range',
                                              paste0('crossover at wk ', round(results[[keep[i]]]$crossover, 1))))
    cv$lower <- cv$estimate - qt(.975, cv$df) * cv$se; cv$upper <- cv$estimate + qt(.975, cv$df) * cv$se
    cv$sig <- ifelse(cv[[sigColumn]] < regionOfSignificanceAlpha, ifelse(cv$estimate > 0, 'Low Risk higher', 'High Risk higher'), NA)
    cv
  }))
  curves$panel <- factor(curves$panel, levels = unique(curves$panel))
  yLabel <- if (grepl('^Network', family)) 'Low minus High Risk connectivity' else 'Low minus High Risk normalized volume'
  plotObj <- ggplot2::ggplot(curves, ggplot2::aes(x = age, y = estimate)) +
    ggplot2::geom_hline(yintercept = 0, color = '#444444', linewidth = .5) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = '#bbbbbb', alpha = .5) +
    ggplot2::geom_line(linewidth = .8, color = '#333333') +
    ggplot2::geom_point(data = curves[!is.na(curves$sig), ], ggplot2::aes(color = sig), size = 1.8) +
    ggplot2::facet_wrap(~panel, ncol = 4, scales = 'free_y') +
    ggplot2::scale_color_manual(values = c('Low Risk higher' = groupColors[['low']], 'High Risk higher' = groupColors[['high']]),
                                limits = c('Low Risk higher', 'High Risk higher'), drop = FALSE, name = NULL) +
    ggplot2::labs(title = paste0(family, ' (', length(keep), ifelse(sum(nSig > 0) > length(keep), paste0(' of ', sum(nSig > 0)), ''),
                                 ' connections/ROIs with a significant age range)'),
                  x = 'Postmenstrual age at scan (weeks)', y = yLabel) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 7.5, hjust = 0), strip.background = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 16, face = 'bold'), legend.position = 'bottom',
                   panel.grid.major = ggplot2::element_line(color = '#e6e6e6'))
  list(plot = plotObj, rows = ceiling(length(keep) / 4))
}

regionPanels <- lapply(names(rosResults), function(an) pmaRegionPanel(an, rosResults[[an]]))
names(regionPanels) <- names(rosResults)
regionPanels <- regionPanels[!vapply(regionPanels, is.null, logical(1))]
if (length(regionPanels) == 0) {
  message('No connection or ROI differs between the groups at any age, so no region-of-significance figure was made.')
} else {
  regionRows <- vapply(regionPanels, `[[`, numeric(1), 'rows')
  regionFigure <- ggpubr::ggarrange(plotlist = lapply(regionPanels, `[[`, 'plot'), ncol = 1, nrow = length(regionPanels),
                                    heights = regionRows + .6, common.legend = TRUE, legend = 'bottom')
  oldSave <- saveFigures; oldTypes <- figureFileTypes
  saveFigures <- TRUE; figureFileTypes <- pmaFigureTypes
  showFigure(regionFigure, 'PMACheck_region_of_significance', width = 15, height = 2.9 * sum(regionRows) + 0.9 * length(regionPanels), subfolder = 'PMACheck')
  saveFigures <- oldSave; figureFileTypes <- oldTypes
  message('Region-of-significance figure saved in ', file.path(figDir, 'PMACheck'))
}
