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
centeringLabels <- c(none = 'No centering', mean = 'Mean centered', wk20 = '20 weeks', wk28 = '28 weeks', wk40 = '40 weeks')

# The whole-brain models (19,900 connections) take long: about as long as the original run for each centering. TRUE adds them.
includeWholeBrain <- FALSE

# Which centering decides what counts as significant or trending in the figure
figureSelectionCentering <- 'none'

# File types for the PMACheck figure (it is large, so a tiff at 600 dpi would be huge: pdf is the default). 'tiff', 'pdf', 'svg', 'eps'
pmaFigureTypes <- c('pdf')

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
  fits <- lapply(names(networkSpecs), function(net) {
    both <- networkSpecs[[net]]$match(connROIs$ROI1) & networkSpecs[[net]]$match(connROIs$ROI2)
    if (sum(both) == 0) { return(NULL) }
    tdNetwork <- tdFC[, both, drop = FALSE]
    asdNetwork <- asdFC[, both, drop = FALSE]
    models <- cached(paste0('pmacheck_network_', net, '_', label),
                     inputs = list(tdNetwork, asdNetwork, demos$td, demos$asd, deparse(fitNetworkModels)),
                     compute = function() fitNetworkModels(tdNetwork, asdNetwork, demos$td, demos$asd))
    list(fcGLM = addFDRfc(models$fcGLM), fcIntxGLM = addFDRfc(models$fcIntxGLM))
  })
  names(fits) <- names(networkSpecs)
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

standardTable <- function(tbl, grpAdj, itxAdj) {
  data.frame(ROI = tbl$ROI, TD.Mean = tbl$TD.Mean, ASD.Mean = tbl$ASD.Mean, Grp.B = tbl$Grp.B,
             Grp.p.adj = tbl[[grpAdj]], Itx.B = tbl[[intersect(c('Itx.B', 'Intx.B'), colnames(tbl))[1]]], Itx.p.adj = tbl[[itxAdj]],
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

# Comparison table ----

countsAt <- function(p, cutoff) sum(p < cutoff, na.rm = TRUE)
summaryRows <- do.call(rbind, lapply(names(analyses), function(an) {
  do.call(rbind, lapply(names(centerings), function(label) {
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
  for (label in names(centerings)) {
    rows <- summaryRows[summaryRows$centering == label, ]
    rows <- rows[match(out$analysis, rows$analysis), ]
    out[[centeringLabels[[label]]]] <- if (is.null(second)) as.character(rows[[first]]) else paste(rows[[first]], '/', rows[[second]])
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

interactionChanged <- any(summaryRows$interaction05 != ave(summaryRows$interaction05, summaryRows$analysis, FUN = function(x) x[1]) |
                          summaryRows$interaction10 != ave(summaryRows$interaction10, summaryRows$analysis, FUN = function(x) x[1]))
cat('\nCheck: the interaction term', ifelse(interactionChanged, 'DID change with centering (unexpected; look into it).', 'is unchanged by centering, as expected.'), '\n')

options(oldWidth)

# Which ROIs/connections change group-effect status (significant < .05, trending < .10, not significant) with centering ----

statusOf <- function(p) ifelse(!is.na(p) & p < .05, 'sig', ifelse(!is.na(p) & p < .10, 'trend', 'ns'))
statusChanges <- do.call(rbind, lapply(names(analyses), function(an) {
  tables <- analyses[[an]]
  status <- sapply(names(centerings), function(label) statusOf(tables[[label]]$Grp.p.adj))
  if (is.null(dim(status))) { status <- matrix(status, nrow = 1, dimnames = list(NULL, names(centerings))) }
  changed <- apply(status, 1, function(x) length(unique(x)) > 1)
  if (!any(changed)) { return(NULL) }
  out <- data.frame(analysis = an, ROI = tables[[1]]$ROI[changed], stringsAsFactors = FALSE)
  for (label in names(centerings)) { out[[paste0('status_', label)]] <- status[changed, label] }
  for (label in names(centerings)) { out[[paste0('betaSign_', label)]] <- ifelse(tables[[label]]$Grp.B[changed] > 0, '+', '-') }
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
