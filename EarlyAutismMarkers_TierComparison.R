# TIER 2 COMPARISON: ORIGINAL vs REVISED NETWORK MEMBERSHIP ----
#
# Lists the ROI-to-ROI connections that gained or lost significance when the network membership changed from the original lists to the revised ones
# (networkSpecsRevised in EarlyAutismMarkers_edit.R). Both versions are fitted here (and cached separately, so neither overwrites the other), connection by connection,
# with the same model as the Tier 2 analysis in EarlyAutismMarkers_edit.R:  connectivity ~ Group + scan age, random intercept per participant.
#
# A connection can change status for two reasons, and the table says which:
#   membership  the connection is inside the network in only one of the two versions (a parcel was added or removed)
#   FDR         the connection is in both versions, but the FDR correction covers a different set of connections, so its adjusted p moved
#
# HOW TO RUN: first run EarlyAutismMarkers_edit.R (so the data and the network definitions are in the R session), then run this whole file.
# OUTPUT: the lists printed in the console, and tableDir/connectivityTiers/tier2_membership_comparison.csv (only if saveTables is TRUE).

# Settings ----

# Networks to compare (the ones that pass the Tier 1 gate in both versions, plus DMN for the exploratory supplement)
comparisonNetworks <- c('CRB', 'Limbic', 'SubCort', 'DMN')
cmpCores <- if (exists('nCores')) nCores else 1

cmpNeeded <- c('tdFC', 'asdFC', 'tdDemoFinal', 'asdDemoFinal', 'connROIs', 'networkSpecs', 'cached', 'connectionLabel', 'tableDir')
cmpMissing <- cmpNeeded[!vapply(cmpNeeded, exists, logical(1))]
if (length(cmpMissing) > 0) { stop('Run EarlyAutismMarkers_edit.R first. These objects are missing: ', paste(cmpMissing, collapse = ', ')) }

# The two sets of network definitions ----

if (!exists('networkSpecsOriginal') || !exists('networkSpecsRevised')) {
  stop('networkSpecsOriginal and networkSpecsRevised are not defined; run EarlyAutismMarkers_edit.R through the network definitions (section 3b) first')
}
originalSpecs <- networkSpecsOriginal
revisedSpecs <- networkSpecsRevised
specsOf <- list(original = originalSpecs, revised = revisedSpecs)

# The same connection model as the Tiers script (copied so this file stands alone) ----

cmpFitChunk <- function(tdChunk, asdChunk, tdBase, asdBase) {
  ctrl <- nlme::lmeControl(msMaxIter = 200, tolerance = 1e-4, msTol = 1e-5, opt = 'optim')
  out <- matrix(NA_real_, ncol(tdChunk), 5)       # TD.Mean, ASD.Mean, Grp.B, Grp.p, nIncluded
  for (j in seq_len(ncol(tdChunk))) {
    tempTD <- tdBase; tempTD$FC <- tdChunk[tdBase$rowIdx, j]
    tempASD <- asdBase; tempASD$FC <- asdChunk[asdBase$rowIdx, j]
    d <- rbind(tempTD, tempASD)
    d$shortID <- substr(d$ID, 1, nchar(d$ID) - 3)
    if (is.na(sum(d$FC))) { next }
    fit <- nlme::lme(FC ~ Grp + ScanGA, random = ~1 | factor(shortID), data = d, control = ctrl)
    tt <- summary(fit)$tTable
    out[j, ] <- c(mean(tempTD$FC), mean(tempASD$FC), tt[2, 1], tt[2, 5], nrow(d))
  }
  out
}
cmpFit <- function(tdMat, asdMat) {
  tdBase <- merge(data.frame(ID = rownames(tdMat), rowIdx = seq_len(nrow(tdMat)), Grp = 'TD', stringsAsFactors = FALSE), tdDemoFinal)
  asdBase <- merge(data.frame(ID = rownames(asdMat), rowIdx = seq_len(nrow(asdMat)), Grp = 'ASD', stringsAsFactors = FALSE), asdDemoFinal)
  tdMat <- as.matrix(tdMat); asdMat <- as.matrix(asdMat)
  chunks <- parallel::splitIndices(ncol(tdMat), max(1, min(ncol(tdMat), cmpCores * 4)))
  chunkData <- lapply(chunks, function(cols) list(td = tdMat[, cols, drop = FALSE], asd = asdMat[, cols, drop = FALSE]))
  if (cmpCores > 1 && length(chunkData) > 1) {
    cl <- parallel::makeCluster(cmpCores); on.exit(parallel::stopCluster(cl))
    parallel::clusterExport(cl, c('tdBase', 'asdBase'), envir = environment())
    parallel::clusterExport(cl, 'cmpFitChunk', envir = globalenv())
    worker <- function(d) cmpFitChunk(d$td, d$asd, tdBase, asdBase); environment(worker) <- globalenv()
    pieces <- parallel::parLapply(cl, chunkData, worker)
  } else {
    pieces <- lapply(chunkData, function(d) cmpFitChunk(d$td, d$asd, tdBase, asdBase))
  }
  res <- data.frame(ROI = colnames(tdMat), do.call(rbind, pieces), stringsAsFactors = FALSE)
  colnames(res) <- c('ROI', 'TD.Mean', 'ASD.Mean', 'Grp.B', 'Grp.p', 'nIncluded')
  res
}

# Fit every network under both versions ----

withinFits <- function(version, net) {
  specs <- specsOf[[version]]
  mask <- specs[[net]]$match(connROIs$ROI1) & specs[[net]]$match(connROIs$ROI2)
  if (sum(mask) == 0) {
    message('No connections in ', net, ' (', version, ')')
    return(data.frame(ROI = character(0), TD.Mean = numeric(0), ASD.Mean = numeric(0), Grp.B = numeric(0), Grp.p = numeric(0), nIncluded = numeric(0), Grp.p.adj = numeric(0)))
  }
  tdN <- tdFC[, mask, drop = FALSE]; asdN <- asdFC[, mask, drop = FALSE]
  fits <- cached(paste0('tiercmp_', version, '_', net),
                 inputs = list(tdN, asdN, tdDemoFinal, asdDemoFinal, deparse(cmpFitChunk), deparse(cmpFit)),
                 compute = function() cmpFit(tdN, asdN))
  fits$Grp.p.adj <- p.adjust(fits$Grp.p, method = 'fdr')       # FDR across the connections of this network in this version
  fits
}
statusOf <- function(p) ifelse(is.na(p), 'not in network', ifelse(p < .05, 'significant', ifelse(p < .10, 'trending', 'ns')))

comparison <- do.call(rbind, lapply(comparisonNetworks, function(net) {
  message('Fitting ', net, ' (original and revised)...')
  o <- withinFits('original', net); r <- withinFits('revised', net)
  both <- merge(o[, c('ROI', 'Grp.B', 'Grp.p', 'Grp.p.adj')], r[, c('ROI', 'Grp.B', 'Grp.p', 'Grp.p.adj')], by = 'ROI', all = TRUE, suffixes = c('.original', '.revised'))
  data.frame(network = net, both, stringsAsFactors = FALSE)
}))
comparison$status.original <- statusOf(comparison$Grp.p.adj.original)
comparison$status.revised <- statusOf(comparison$Grp.p.adj.revised)
comparison$direction <- ifelse(is.na(comparison$Grp.B.revised), ifelse(comparison$Grp.B.original > 0, 'Low-risk higher', 'High-risk higher'),
                               ifelse(comparison$Grp.B.revised > 0, 'Low-risk higher', 'High-risk higher'))
comparison$label <- gsub('\n', ' ', connectionLabel(comparison$ROI))
comparison$reason <- ifelse(comparison$status.original == 'not in network', 'membership: added in the revised network',
                     ifelse(comparison$status.revised == 'not in network', 'membership: removed in the revised network', 'FDR: same connection, different correction set'))
rank <- c('not in network' = 0, ns = 1, trending = 2, significant = 3)
comparison$change <- ifelse(comparison$status.original == comparison$status.revised, 'same',
                     ifelse(rank[comparison$status.revised] > rank[comparison$status.original], 'gained', 'lost'))

# Report ----

changed <- comparison[comparison$change != 'same', ]
changed <- changed[order(changed$network, changed$change, changed$Grp.p.original), ]
showCols <- c('network', 'ROI', 'label', 'direction', 'status.original', 'status.revised', 'Grp.p.original', 'Grp.p.revised', 'Grp.p.adj.original', 'Grp.p.adj.revised', 'reason')
cmpWidth <- options(width = 250)
cat('\n=== Tier 2 connections whose significance changed: original -> revised membership ===\n')
cat('Status: significant = FDR p < .05, trending = < .10, ns = not significant; "not in network" = the connection is not part of that network in that version.\n')
cat('Only ', 'significant/trending vs ns or absent counts as a change; ns -> not in network (and the reverse) is not shown.\n', sep = '')
meaningful <- changed[!(changed$status.original %in% c('ns', 'not in network') & changed$status.revised %in% c('ns', 'not in network')), ]
cat('\nGained significance or trend (', sum(meaningful$change == 'gained'), '):\n', sep = '')
if (any(meaningful$change == 'gained')) print(format(meaningful[meaningful$change == 'gained', showCols], digits = 3), row.names = FALSE) else cat('  none\n')
cat('\nLost significance or trend (', sum(meaningful$change == 'lost'), '):\n', sep = '')
if (any(meaningful$change == 'lost')) print(format(meaningful[meaningful$change == 'lost', showCols], digits = 3), row.names = FALSE) else cat('  none\n')
cat('\nSummary by network:\n')
print(do.call(rbind, lapply(comparisonNetworks, function(net) {
  d <- comparison[comparison$network == net, ]
  data.frame(network = net, significantOriginal = sum(d$status.original == 'significant'), significantRevised = sum(d$status.revised == 'significant'),
             trendingOriginal = sum(d$status.original == 'trending'), trendingRevised = sum(d$status.revised == 'trending'),
             unchangedSigOrTrend = sum(d$status.original == d$status.revised & d$status.original %in% c('significant', 'trending')),
             stringsAsFactors = FALSE)
})), row.names = FALSE)
options(cmpWidth)

if (exists('saveTables') && isTRUE(saveTables)) {
  dir.create(file.path(tableDir, 'connectivityTiers'), showWarnings = FALSE, recursive = TRUE)
  write.csv(comparison[order(comparison$network, comparison$Grp.p.revised), ], file.path(tableDir, 'connectivityTiers', 'tier2_membership_comparison.csv'), row.names = FALSE)
}
