# CONNECTIVITY IN THREE TIERS ----
#
# A prototype of the new connectivity analysis, kept separate from EarlyAutismMarkers_edit.R until it is approved.
# Every model here is  connectivity ~ Group + scan age  (no interaction) with a random intercept per participant.
# The group estimate is Low minus High Risk, adjusted for age (positive = Low Risk higher).
#
#   TIER 1  Network averages: the average connectivity within each a priori network and between each pair of them
#           (one score per scan). FDR across the networks (within) and across the network pairs (between), separately.
#   TIER 2  ROI-to-ROI connections, ONLY inside the networks (and pairs) that are significant in Tier 1.
#           FDR across the connections of each network (or across all Tier 2 connections; see tier2FDR).
#   TIER 3  Exploratory: every ROI to every other ROI (whole brain), FDR across all connections.
#
# HOW TO RUN: first run EarlyAutismMarkers_edit.R (so the data and functions are in the R session), then run this whole file.
# Fits are cached in cacheDir, so a second run is quick.
#
# OUTPUT: figures in figDir/connectivityTiers, tables in tableDir/connectivityTiers (tables only if saveTables is TRUE).

# Settings ----

# A Tier 1 result opens Tier 2 for that network (or pair) when its FDR-adjusted p is below this
tier1Gate <- .05

# Also look at the ROI-to-ROI connections BETWEEN two networks when that pair is significant in Tier 1
tier2IncludeBetween <- TRUE

# 'withinNetwork' = FDR across the connections of each gated network separately (as in the earlier network analyses);
# 'allTier2' = one FDR across every Tier 2 connection together
tier2FDR <- 'withinNetwork'

# Networks that did NOT pass the Tier 1 gate but are still taken through Tier 2 for exploratory reporting (for example in a supplement).
# They are marked 'exploratory' (gated = FALSE) in the tables and figure names, and they never change the FDR of the gated networks.
# Use character(0) for none.
tier2ExploratoryNetworks <- c('DMN')

# The whole-brain exploratory tier fits ~20,000 models; FALSE skips it
runTier3 <- TRUE

# How the connections of a network are combined in Tier 1: 'mean' = plain average; 'absmean' = average of absolute values
networkScore <- 'mean'

# Use the revised network membership (the parcel-by-parcel decisions in networkSpecs_revised.R) instead of the original lists in networkSpecs.
# FALSE runs the original lists (needed if you want to compare the two). The file is looked for in baseDir/functional/atlasParcellattion,
# then in the working folder.
useRevisedNetworks <- TRUE
revisedNetworksFile <- if (exists('baseDir') && file.exists(file.path(baseDir, 'functional', 'atlasParcellattion', 'networkSpecs_revised.R')))
  file.path(baseDir, 'functional', 'atlasParcellattion', 'networkSpecs_revised.R') else 'networkSpecs_revised.R'

# Subfolder (inside the figure and table folders) for everything this script writes
tiersFolder <- 'connectivityTiers'

# Cores for the model fits (uses nCores from EarlyAutismMarkers_edit.R when it exists)
tierCores <- if (exists('nCores')) nCores else 1

# Make sure the data and functions from EarlyAutismMarkers_edit.R are available ----

tierNeeded <- c('tdFC', 'asdFC', 'tdDemoFinal', 'asdDemoFinal', 'connROIs', 'networkSpecs', 'networksOfInterest', 'groupColors', 'networkBarFill',
                'showFigure', 'plotNetworkBars', 'addFDRfc', 'cached', 'connectionLabel', 'figDir', 'tableDir')
tierMissing <- tierNeeded[!vapply(tierNeeded, exists, logical(1))]
if (length(tierMissing) > 0) { stop('Run EarlyAutismMarkers_edit.R first (current version). These objects are missing: ', paste(tierMissing, collapse = ', ')) }
if (!networkScore %in% c('mean', 'absmean')) { stop("networkScore must be 'mean' or 'absmean'") }
if (!tier2FDR %in% c('withinNetwork', 'allTier2')) { stop("tier2FDR must be 'withinNetwork' or 'allTier2'") }
if (useRevisedNetworks) {
  if (!file.exists(revisedNetworksFile)) { stop('Revised network file not found: ', revisedNetworksFile, '. Put networkSpecs_revised.R there or set revisedNetworksFile.') }
  source(revisedNetworksFile)            # replaces networkSpecs for the four networks and keeps the old ones in networkSpecsOriginal
} else if (exists('networkSpecsOriginal')) {
  networkSpecs <- networkSpecsOriginal   # go back to the original lists if the revised ones were loaded earlier in this session
}
message('Tier networks: ', paste(networksOfInterest, collapse = ', '), if (useRevisedNetworks) ' (revised membership)' else ' (original membership)')
if (!exists('resultTables')) { resultTables <- list() }

# Helpers ----

networkLabel <- function(net) {
  if (exists('networkNodeLabels') && net %in% names(networkNodeLabels)) networkNodeLabels[[net]] else net
}
starFor <- function(p) ifelse(p < .05, '*', ifelse(p < .10, '+', ''))
fmtP <- function(p) ifelse(p < .001, '<.001', sub('^0', '', sprintf('%.3f', p)))

saveTierTable <- function(tbl, name) {
  resultTables[[name]] <<- tbl
  if (exists('saveTables') && isTRUE(saveTables)) {
    dir.create(file.path(tableDir, tiersFolder), showWarnings = FALSE, recursive = TRUE)
    write.csv(tbl, file.path(tableDir, tiersFolder, paste0(name, '.csv')), row.names = FALSE)
  }
  invisible(tbl)
}

# The ROI-to-ROI model, for a set of connections (columns). Self-contained so it can run on a worker.
# tdBase / asdBase: each group's demographics with rowIdx = the participant's row in the FC matrix.
# Connections with any missing value are left as NA (as in the earlier whole-brain analysis).
fitMainChunk <- function(tdChunk, asdChunk, tdBase, asdBase) {
  ctrl <- nlme::lmeControl(msMaxIter = 200, tolerance = 1e-4, msTol = 1e-5, opt = 'optim')
  out <- matrix(NA_real_, ncol(tdChunk), 9)       # TD.Mean, ASD.Mean, TD.SE, ASD.SE, Grp.B, Grp.p, GA.B, GA.p, nIncluded
  for (j in seq_len(ncol(tdChunk))) {
    tempTD <- tdBase; tempTD$FC <- tdChunk[tdBase$rowIdx, j]
    tempASD <- asdBase; tempASD$FC <- asdChunk[asdBase$rowIdx, j]
    d <- rbind(tempTD, tempASD)
    d$shortID <- substr(d$ID, 1, nchar(d$ID) - 3)
    if (is.na(sum(d$FC))) { next }
    fit <- nlme::lme(FC ~ Grp + ScanGA, random = ~1 | factor(shortID), data = d, control = ctrl)
    tt <- summary(fit)$tTable
    out[j, ] <- c(mean(tempTD$FC), mean(tempASD$FC), sd(tempTD$FC) / sqrt(nrow(tempTD)), sd(tempASD$FC) / sqrt(nrow(tempASD)),
                  tt[2, 1], tt[2, 5], tt[3, 1], tt[3, 5], nrow(d))
  }
  out
}

# Run fitMainChunk over connections, split across cores. Returns one row per connection.
fitConnections <- function(tdMat, asdMat, tdDemo, asdDemo, cores = 1) {
  tdBase <- merge(data.frame(ID = rownames(tdMat), rowIdx = seq_len(nrow(tdMat)), Grp = 'TD', stringsAsFactors = FALSE), tdDemo)
  asdBase <- merge(data.frame(ID = rownames(asdMat), rowIdx = seq_len(nrow(asdMat)), Grp = 'ASD', stringsAsFactors = FALSE), asdDemo)
  tdMat <- as.matrix(tdMat); asdMat <- as.matrix(asdMat)
  chunks <- parallel::splitIndices(ncol(tdMat), max(1, min(ncol(tdMat), cores * 4)))
  chunkData <- lapply(chunks, function(cols) list(td = tdMat[, cols, drop = FALSE], asd = asdMat[, cols, drop = FALSE]))
  if (cores > 1 && length(chunkData) > 1) {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterExport(cl, c('tdBase', 'asdBase'), envir = environment())
    parallel::clusterExport(cl, 'fitMainChunk', envir = globalenv())
    worker <- function(d) fitMainChunk(d$td, d$asd, tdBase, asdBase)
    environment(worker) <- globalenv()
    pieces <- parallel::parLapply(cl, chunkData, worker)
  } else {
    pieces <- lapply(chunkData, function(d) fitMainChunk(d$td, d$asd, tdBase, asdBase))
  }
  res <- data.frame(ROI = colnames(tdMat), do.call(rbind, pieces), stringsAsFactors = FALSE)
  colnames(res) <- c('ROI', 'TD.Mean', 'ASD.Mean', 'TD.SE', 'ASD.SE', 'Grp.B', 'Grp.p', 'GA.B', 'GA.p', 'nIncluded')
  res
}

# TIER 1: network averages ----

inNetwork <- lapply(setNames(networksOfInterest, networksOfInterest), function(net) {
  list(roi1 = networkSpecs[[net]]$match(connROIs$ROI1), roi2 = networkSpecs[[net]]$match(connROIs$ROI2))
})
withinMask <- function(net) inNetwork[[net]]$roi1 & inNetwork[[net]]$roi2
# one region in each network, and not also a within-network connection (some regions belong to two networks)
betweenMask <- function(a, b) {
  ((inNetwork[[a]]$roi1 & inNetwork[[b]]$roi2) | (inNetwork[[b]]$roi1 & inNetwork[[a]]$roi2)) & !withinMask(a) & !withinMask(b)
}

scoreData <- function(mask) {
  one <- function(fc, demo, grp) {
    values <- as.matrix(fc[, mask, drop = FALSE])
    if (networkScore == 'absmean') { values <- abs(values) }
    d <- merge(data.frame(ID = rownames(fc), FC = rowMeans(values, na.rm = TRUE), stringsAsFactors = FALSE), demo[, c('ID', 'ScanGA')], by = 'ID')
    d$Grp <- grp
    d
  }
  d <- rbind(one(tdFC, tdDemoFinal, 'TD'), one(asdFC, asdDemoFinal, 'ASD'))
  d <- d[is.finite(d$FC), ]
  d$shortID <- substr(d$ID, 1, nchar(d$ID) - 3)
  d$Grp <- factor(d$Grp, levels = c('ASD', 'TD'))        # ASD (High Risk) is the reference, so the estimate is Low minus High Risk
  attr(d, 'nConnections') <- sum(mask)
  d
}

fitScore <- function(d) {
  fit <- nlme::lme(FC ~ Grp + ScanGA, random = ~1 | factor(shortID), data = d,
                   control = nlme::lmeControl(msMaxIter = 200, tolerance = 1e-4, msTol = 1e-5, opt = 'optim'))
  tt <- summary(fit)$tTable
  adj <- as.data.frame(emmeans::emmeans(fit, 'Grp', data = d))
  data.frame(nConnections = attr(d, 'nConnections'), nScans = nrow(d), nParticipants = length(unique(d$shortID)),
             LowRisk.Mean = adj$emmean[adj$Grp == 'TD'], LowRisk.SE = adj$SE[adj$Grp == 'TD'],
             HighRisk.Mean = adj$emmean[adj$Grp == 'ASD'], HighRisk.SE = adj$SE[adj$Grp == 'ASD'],
             Diff.LowMinusHigh = tt['GrpTD', 'Value'], Diff.SE = tt['GrpTD', 'Std.Error'], Diff.df = tt['GrpTD', 'DF'],
             Diff.inSD = tt['GrpTD', 'Value'] / sd(d$FC), Grp.p = tt['GrpTD', 'p-value'], GA.B = tt['ScanGA', 'Value'], GA.p = tt['ScanGA', 'p-value'],
             stringsAsFactors = FALSE)
}

runScores <- function(masks, ids) {
  masks <- masks[vapply(masks, sum, numeric(1)) > 0]
  frames <- lapply(masks, scoreData)
  fits <- do.call(rbind, lapply(names(frames), function(nm) cbind(ids[[nm]], fitScore(frames[[nm]]), stringsAsFactors = FALSE)))
  fits$Grp.p.adj <- p.adjust(fits$Grp.p, method = 'fdr')
  fits$GA.p.adj <- p.adjust(fits$GA.p, method = 'fdr')
  rownames(fits) <- NULL
  list(fits = fits, frames = frames)
}

withinMasks <- setNames(lapply(networksOfInterest, withinMask), networksOfInterest)
withinIds <- setNames(lapply(networksOfInterest, function(n) data.frame(network = n, label = networkLabel(n), stringsAsFactors = FALSE)), networksOfInterest)
tier1Within <- runScores(withinMasks, withinIds)

pairs <- combn(networksOfInterest, 2, simplify = FALSE)
pairNames <- vapply(pairs, paste, character(1), collapse = ' x ')
betweenMasks <- setNames(lapply(pairs, function(p) betweenMask(p[1], p[2])), pairNames)
betweenIds <- setNames(lapply(pairs, function(p) data.frame(network = paste(p, collapse = ' x '),
                                                             label = paste(networkLabel(p[1]), 'x', networkLabel(p[2])), stringsAsFactors = FALSE)), pairNames)
tier1Between <- runScores(betweenMasks, betweenIds)

showColumns <- c('label', 'nConnections', 'nScans', 'LowRisk.Mean', 'HighRisk.Mean', 'Diff.LowMinusHigh', 'Diff.inSD', 'Grp.p', 'Grp.p.adj', 'GA.p')
tierWidth <- options(width = 220)
cat('\n=== TIER 1: network averages (score = ', networkScore, ') ===\n', sep = '')
cat('Diff.LowMinusHigh: Low minus High Risk, age-adjusted (positive = Low Risk higher); Diff.inSD: in SD of the scores.\n',
    'p.adj: FDR within the table (', nrow(tier1Within$fits), ' networks; ', nrow(tier1Between$fits), ' network pairs). Gate for Tier 2: p.adj < ', tier1Gate, '\n', sep = '')
cat('\nWITHIN networks\n'); print(format(tier1Within$fits[, showColumns], digits = 3), row.names = FALSE)
cat('\nBETWEEN networks\n'); print(format(tier1Between$fits[, showColumns], digits = 3), row.names = FALSE)
options(tierWidth)
saveTierTable(tier1Within$fits, 'tier1_within')
saveTierTable(tier1Between$fits, 'tier1_between')

# Tier 1 figures
barsFor <- function(fits, xAngle = 0) {
  labels <- fits$label
  barData <- rbind(data.frame(label = labels, Group = 'Low Risk', Mean = fits$LowRisk.Mean, SE = fits$LowRisk.SE),
                   data.frame(label = labels, Group = 'High Risk', Mean = fits$HighRisk.Mean, SE = fits$HighRisk.SE))
  barData$label <- factor(barData$label, levels = labels); barData$Group <- factor(barData$Group, levels = c('Low Risk', 'High Risk'))
  starData <- data.frame(label = factor(labels, levels = labels), star = starFor(fits$Grp.p.adj),
                         y = pmax(fits$LowRisk.Mean + fits$LowRisk.SE, fits$HighRisk.Mean + fits$HighRisk.SE))
  starData$y <- starData$y + .04 * diff(range(c(0, barData$Mean + barData$SE)))
  ggplot2::ggplot(barData, ggplot2::aes(x = label, y = Mean, fill = Group, color = Group)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(.8), width = .72, linewidth = .8) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Mean - SE, ymax = Mean + SE), position = ggplot2::position_dodge(.8), width = .25, color = 'black', linewidth = .6) +
    ggplot2::geom_text(data = starData, ggplot2::aes(x = label, y = y, label = star), inherit.aes = FALSE, size = 7) +
    ggplot2::scale_fill_manual(values = c('Low Risk' = networkBarFill[['low']], 'High Risk' = networkBarFill[['high']]), name = NULL) +
    ggplot2::scale_color_manual(values = c('Low Risk' = groupColors[['low']], 'High Risk' = groupColors[['high']]), name = NULL) +
    ggplot2::labs(x = NULL, y = paste0('Average connectivity', ifelse(networkScore == 'absmean', ' (absolute values)', '')),
                  caption = '* p.adj < .05, + p.adj < .10 (FDR within this figure); age-adjusted means +/- SE') +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(legend.position = 'bottom', axis.text.y = ggplot2::element_text(color = 'black'),
                   axis.text.x = ggplot2::element_text(color = 'black', angle = xAngle, hjust = ifelse(xAngle == 0, .5, 1), vjust = ifelse(xAngle == 0, 1, .5)))
}
showFigure(barsFor(tier1Within$fits), 'tier1_within_bars', width = 9, height = 6, subfolder = tiersFolder)
showFigure(barsFor(tier1Between$fits, xAngle = 45), 'tier1_between_bars', width = 11, height = 7, subfolder = tiersFolder)

cells <- rbind(
  data.frame(row = tier1Within$fits$network, col = tier1Within$fits$network, diff = tier1Within$fits$Diff.LowMinusHigh, star = starFor(tier1Within$fits$Grp.p.adj), stringsAsFactors = FALSE),
  do.call(rbind, lapply(seq_along(pairs), function(i) {
    r <- tier1Between$fits[tier1Between$fits$network == pairNames[i], ]
    if (nrow(r) == 0) { return(NULL) }
    data.frame(row = pairs[[i]], col = rev(pairs[[i]]), diff = r$Diff.LowMinusHigh, star = starFor(r$Grp.p.adj), stringsAsFactors = FALSE)
  })))
cells$row <- factor(cells$row, levels = networksOfInterest, labels = vapply(networksOfInterest, networkLabel, character(1)))
cells$col <- factor(cells$col, levels = rev(networksOfInterest), labels = rev(vapply(networksOfInterest, networkLabel, character(1))))
cells$text <- paste0(sprintf('%.3f', cells$diff), '\n', cells$star)
limit <- max(abs(cells$diff))
tier1_matrix <- ggplot2::ggplot(cells, ggplot2::aes(x = row, y = col, fill = diff)) +
  ggplot2::geom_tile(color = 'white', linewidth = 1) + ggplot2::geom_text(ggplot2::aes(label = text), size = 4, lineheight = .9) +
  ggplot2::scale_fill_gradient2(low = groupColors[['high']], mid = 'white', high = groupColors[['low']], midpoint = 0, limits = c(-limit, limit), name = 'Low minus High\nRisk difference') +
  ggplot2::labs(x = NULL, y = NULL, caption = 'Diagonal = within the network; off-diagonal = between two networks. Blue = Low Risk higher, red = High Risk higher.\n* p.adj < .05, + p.adj < .10 (FDR within the within and between families separately); age-adjusted') +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(color = 'black', angle = 45, hjust = 1), axis.text.y = ggplot2::element_text(color = 'black'))
showFigure(tier1_matrix, 'tier1_matrix', width = 8, height = 7, subfolder = tiersFolder)

# TIER 2: ROI-to-ROI connections, only in the networks (pairs) that are significant in Tier 1 ----

gatedWithin <- tier1Within$fits$network[tier1Within$fits$Grp.p.adj < tier1Gate]
gatedBetween <- if (tier2IncludeBetween) tier1Between$fits$network[tier1Between$fits$Grp.p.adj < tier1Gate] else character(0)
cat('\n=== TIER 2: ROI-to-ROI connections ===\n')
cat('Networks significant in Tier 1 (p.adj <', tier1Gate, '):', ifelse(length(gatedWithin) == 0, 'none', paste(gatedWithin, collapse = ', ')), '\n')
if (tier2IncludeBetween) { cat('Network pairs significant in Tier 1:', ifelse(length(gatedBetween) == 0, 'none', paste(gatedBetween, collapse = ', ')), '\n') }

exploratoryWithin <- setdiff(intersect(tier2ExploratoryNetworks, networksOfInterest), gatedWithin)
if (length(exploratoryWithin) > 0) { cat('Exploratory (not gated by Tier 1, reported separately):', paste(exploratoryWithin, collapse = ', '), '\n') }
tier2Units <- c(setNames(lapply(gatedWithin, withinMask), gatedWithin),
                setNames(lapply(gatedBetween, function(nm) { p <- strsplit(nm, ' x ', fixed = TRUE)[[1]]; betweenMask(p[1], p[2]) }), gatedBetween),
                setNames(lapply(exploratoryWithin, withinMask), exploratoryWithin))
tier2 <- NULL
if (length(tier2Units) == 0) {
  message('No network or pair passes the Tier 1 gate, so Tier 2 is empty.')
} else {
  tier2Parts <- lapply(names(tier2Units), function(unit) {
    mask <- tier2Units[[unit]]
    if (sum(mask) == 0) { message('No connections in ', unit, '; skipped'); return(NULL) }
    tdN <- tdFC[, mask, drop = FALSE]; asdN <- asdFC[, mask, drop = FALSE]
    fits <- cached(paste0('tiers_t2_', gsub('[^A-Za-z0-9]+', '_', unit)),
                   inputs = list(tdN, asdN, tdDemoFinal, asdDemoFinal, deparse(fitMainChunk), deparse(fitConnections)),
                   compute = function() fitConnections(tdN, asdN, tdDemoFinal, asdDemoFinal, tierCores))
    cbind(unit = unit, fits, stringsAsFactors = FALSE)
  })
  tier2 <- do.call(rbind, tier2Parts); rownames(tier2) <- NULL
  tier2$gated <- !tier2$unit %in% exploratoryWithin
  # FDR: within each network, or (allTier2) across all the gated connections together; exploratory networks are always corrected on their own
  adjustTier2 <- function(p) {
    out <- rep(NA_real_, length(p))
    pooled <- if (tier2FDR == 'allTier2') tier2$gated else rep(FALSE, length(p))
    if (any(pooled)) { out[pooled] <- p.adjust(p[pooled], method = 'fdr') }
    for (u in unique(tier2$unit[!pooled])) { i <- which(tier2$unit == u & !pooled); out[i] <- p.adjust(p[i], method = 'fdr') }
    out
  }
  tier2$Grp.p.adj <- adjustTier2(tier2$Grp.p)
  tier2$GA.p.adj <- adjustTier2(tier2$GA.p)
  tier2$label <- gsub('\n', ' ', connectionLabel(tier2$ROI))
  tier2$direction <- ifelse(tier2$Grp.B > 0, 'Low Risk higher', 'High Risk higher')
  tier2$status <- ifelse(tier2$Grp.p.adj < .05, 'significant', ifelse(tier2$Grp.p.adj < .10, 'trending', 'ns'))

  summaryTier2 <- do.call(rbind, lapply(split(tier2, tier2$unit), function(d) data.frame(unit = d$unit[1], gated = d$gated[1], tested = nrow(d), sig05 = sum(d$Grp.p.adj < .05, na.rm = TRUE),
                                                                                       trend10 = sum(d$Grp.p.adj >= .05 & d$Grp.p.adj < .10, na.rm = TRUE),
                                                                                       lowHigher = sum(d$Grp.p.adj < .10 & d$Grp.B > 0, na.rm = TRUE),
                                                                                       highHigher = sum(d$Grp.p.adj < .10 & d$Grp.B < 0, na.rm = TRUE), stringsAsFactors = FALSE)))
  rownames(summaryTier2) <- NULL
  cat('\nConnections tested and significant/trending (FDR ', ifelse(tier2FDR == 'allTier2', 'across all Tier 2 connections', 'within each network'), ')\n', sep = '')
  print(summaryTier2, row.names = FALSE)
  saveTierTable(tier2, 'tier2_connections')
  saveTierTable(tier2[tier2$status != 'ns', ], 'tier2_significant_or_trending')

  # bar chart per network (the same chart as the earlier network analyses); exploratory networks get '_exploratory' in the file name
  for (unit in c(gatedWithin, exploratoryWithin)) {
    d <- tier2[tier2$unit == unit, ]
    showFigure(plotNetworkBars(d[, c('ROI', 'TD.Mean', 'ASD.Mean', 'Grp.p')], networkSpecs[[unit]]),
               paste0('tier2_', unit, '_bars', ifelse(unit %in% exploratoryWithin, '_exploratory', '')), width = 12, height = 6, subfolder = tiersFolder)
  }
}

# TIER 3: exploratory whole brain, every ROI to every other ROI ----

tier3 <- NULL
if (runTier3) {
  cat('\n=== TIER 3: exploratory whole-brain connectivity ===\n')
  tier3Fits <- cached('tiers_t3_wholebrain',
                      inputs = list(tdFC, asdFC, tdDemoFinal, asdDemoFinal, deparse(fitMainChunk), deparse(fitConnections)),
                      compute = function() {
                        message('Fitting ', ncol(tdFC), ' connections on ', tierCores, ' core(s)...')
                        fitConnections(tdFC, asdFC, tdDemoFinal, asdDemoFinal, tierCores)
                      })
  tier3 <- tier3Fits
  tier3$Grp.p.adj <- p.adjust(tier3$Grp.p, method = 'fdr')
  tier3$GA.p.adj <- p.adjust(tier3$GA.p, method = 'fdr')
  tier3$direction <- ifelse(tier3$Grp.B > 0, 'Low Risk higher', 'High Risk higher')
  cat('Connections tested:', sum(!is.na(tier3$Grp.p)), 'of', nrow(tier3), '\n')
  cat('Group effect, FDR p < .05:', sum(tier3$Grp.p.adj < .05, na.rm = TRUE), ' (Low Risk higher:', sum(tier3$Grp.p.adj < .05 & tier3$Grp.B > 0, na.rm = TRUE),
      ', High Risk higher:', sum(tier3$Grp.p.adj < .05 & tier3$Grp.B < 0, na.rm = TRUE), ')\n')
  cat('Group effect, FDR p < .10:', sum(tier3$Grp.p.adj < .10, na.rm = TRUE), '\n')
  top <- head(tier3[order(tier3$Grp.p), c('ROI', 'TD.Mean', 'ASD.Mean', 'Grp.B', 'Grp.p', 'Grp.p.adj')], 15)
  tierWidth <- options(width = 200); cat('\nSmallest p-values:\n'); print(format(top, digits = 3), row.names = FALSE); options(tierWidth)
  saveTierTable(tier3, 'tier3_wholebrain_all')
  saveTierTable(tier3[which(tier3$Grp.p.adj < .10), ], 'tier3_wholebrain_significant_or_trending')
}
