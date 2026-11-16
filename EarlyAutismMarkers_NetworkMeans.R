# NETWORK-LEVEL CONNECTIVITY (exploratory) ----
#
# Instead of testing every ROI-to-ROI connection separately, this averages connectivity across all the connections of
#   - a network ("within": both regions in the same network, the connections used by the network analyses), and
#   - a pair of networks ("between": one region in each network)
# into one score per scan, then compares the groups on that score. One test per network or network pair, so far fewer tests to correct for.
#
# A connection counts as "between" networks A and B when one region is in A and the other in B, and it is NOT also a within-A or
# within-B connection (some regions belong to two networks, e.g. thalamus regions in both the subcortical and auditory networks).
#
# For each score it fits  (1) score ~ Group + age   (primary: the average group difference)
#                         (2) score ~ Group * age   (age centered at the sample mean: does the difference change with age?)
# with a random intercept per participant (the same models as the connection-level analyses). FDR is applied across the networks
# (within family) and across the network pairs (between family) separately.
#
# HOW TO RUN: first run EarlyAutismMarkers_edit.R (so the data and functions are in the R session), then run this whole file.
#
# OUTPUT (figures in figDir/networkAverages, tables in tableDir/networkAverages; tables only if saveTables is TRUE):
#   networkAverages_within_bars / _within_by_age     adjusted group means and score-by-age plots for the networks
#   networkAverages_between_bars / _between_by_age   the same for the network pairs
#   networkAverages_matrix                           Low minus High Risk difference for every network and pair in one matrix
#   tables: networkAverages_within, networkAverages_between

# Settings ----

# How the connections are combined: 'mean' = plain average of the connection values;
#   'absmean' = average of their absolute values (strength regardless of sign)
networkScore <- 'mean'

# Networks to analyze (names in networkSpecs); NULL = all of them
networkMeanNetworks <- NULL

# Subfolder (inside the figure and table folders) for everything this script writes
networkAveragesFolder <- 'networkAverages'

# Make sure the data and functions from EarlyAutismMarkers_edit.R are available ----

nmNeeded <- c('tdFC', 'asdFC', 'tdDemoFinal', 'asdDemoFinal', 'connROIs', 'networkSpecs', 'groupColors', 'showFigure', 'figDir', 'tableDir')
nmMissing <- nmNeeded[!vapply(nmNeeded, exists, logical(1))]
if (length(nmMissing) > 0) { stop('Run EarlyAutismMarkers_edit.R first. These objects are missing: ', paste(nmMissing, collapse = ', ')) }
if (!networkScore %in% c('mean', 'absmean')) { stop("networkScore must be 'mean' or 'absmean'") }
if (is.null(networkMeanNetworks)) { networkMeanNetworks <- names(networkSpecs) }

# Which connections belong to a network or a network pair ----

inNetwork <- lapply(setNames(networkMeanNetworks, networkMeanNetworks), function(net) {
  list(roi1 = networkSpecs[[net]]$match(connROIs$ROI1), roi2 = networkSpecs[[net]]$match(connROIs$ROI2))
})
withinMask <- function(net) inNetwork[[net]]$roi1 & inNetwork[[net]]$roi2
betweenMask <- function(a, b) {
  mask <- (inNetwork[[a]]$roi1 & inNetwork[[b]]$roi2) | (inNetwork[[b]]$roi1 & inNetwork[[a]]$roi2)
  mask & !withinMask(a) & !withinMask(b)
}

networkLabel <- function(net) {
  if (exists('networkNodeLabels') && net %in% names(networkNodeLabels)) networkNodeLabels[[net]] else net
}

# One score per scan ----

scoreData <- function(mask) {
  if (sum(mask) == 0) { return(NULL) }
  one <- function(fc, demo, grp) {
    values <- as.matrix(fc[, mask, drop = FALSE])
    if (networkScore == 'absmean') { values <- abs(values) }
    d <- data.frame(ID = rownames(fc), FC = rowMeans(values, na.rm = TRUE), stringsAsFactors = FALSE)
    d <- merge(d, demo[, c('ID', 'ScanGA')], by = 'ID')
    d$Grp <- grp
    d
  }
  d <- rbind(one(tdFC, tdDemoFinal, 'TD'), one(asdFC, asdDemoFinal, 'ASD'))
  d <- d[is.finite(d$FC), ]
  d$shortID <- substr(d$ID, 1, nchar(d$ID) - 3)
  d$Grp <- factor(d$Grp, levels = c('ASD', 'TD'))        # ASD (High Risk) is the reference, so the group estimate is Low minus High Risk
  attr(d, 'nConnections') <- sum(mask)
  d
}

# Fit both models for one score ----

fitScore <- function(d) {
  ctrl <- nlme::lmeControl(msMaxIter = 200, tolerance = 1e-4, msTol = 1e-5, opt = 'optim')
  main <- nlme::lme(FC ~ Grp + ScanGA, random = ~1 | factor(shortID), data = d, control = ctrl)
  d$ScanGAc <- d$ScanGA - mean(d$ScanGA)
  intx <- nlme::lme(FC ~ Grp * ScanGAc, random = ~1 | factor(shortID), data = d, control = ctrl)
  tm <- summary(main)$tTable; ti <- summary(intx)$tTable
  adj <- as.data.frame(emmeans::emmeans(main, 'Grp', data = d))
  data.frame(nConnections = attr(d, 'nConnections'), nScans = nrow(d), nParticipants = length(unique(d$shortID)),
             nLow = sum(d$Grp == 'TD'), nHigh = sum(d$Grp == 'ASD'),
             LowRisk.Mean = adj$emmean[adj$Grp == 'TD'], LowRisk.SE = adj$SE[adj$Grp == 'TD'],
             HighRisk.Mean = adj$emmean[adj$Grp == 'ASD'], HighRisk.SE = adj$SE[adj$Grp == 'ASD'],
             Diff.LowMinusHigh = tm['GrpTD', 'Value'], Diff.SE = tm['GrpTD', 'Std.Error'], Diff.df = tm['GrpTD', 'DF'],
             Diff.inSD = tm['GrpTD', 'Value'] / sd(d$FC),
             Grp.p = tm['GrpTD', 'p-value'], GA.B = tm['ScanGA', 'Value'], GA.p = tm['ScanGA', 'p-value'],
             Intx.B = ti['GrpTD:ScanGAc', 'Value'], Intx.p = ti['GrpTD:ScanGAc', 'p-value'],
             stringsAsFactors = FALSE)
}

# A family of scores: masks is a named list (name -> logical vector over connections). Returns the fitted table and the data frames.
runFamily <- function(masks, idColumns) {
  frames <- lapply(masks, scoreData)
  frames <- frames[!vapply(frames, is.null, logical(1))]
  fits <- do.call(rbind, lapply(names(frames), function(nm) cbind(idColumns[[nm]], fitScore(frames[[nm]]), stringsAsFactors = FALSE)))
  fits$Grp.p.adj <- p.adjust(fits$Grp.p, method = 'fdr')
  fits$Intx.p.adj <- p.adjust(fits$Intx.p, method = 'fdr')
  fits$GA.p.adj <- p.adjust(fits$GA.p, method = 'fdr')
  rownames(fits) <- NULL
  list(fits = fits, frames = frames)
}

withinMasks <- lapply(setNames(networkMeanNetworks, networkMeanNetworks), withinMask)
withinId <- lapply(names(withinMasks), function(net) data.frame(network = net, label = networkLabel(net), stringsAsFactors = FALSE))
names(withinId) <- names(withinMasks)
within <- runFamily(withinMasks, withinId)

pairs <- combn(networkMeanNetworks, 2, simplify = FALSE)
pairNames <- vapply(pairs, function(p) paste(p, collapse = ' x '), character(1))
betweenMasks <- setNames(lapply(pairs, function(p) betweenMask(p[1], p[2])), pairNames)
betweenId <- setNames(lapply(pairs, function(p) data.frame(network = paste(p, collapse = ' x '),
                                                            label = paste(networkLabel(p[1]), 'x', networkLabel(p[2])), stringsAsFactors = FALSE)), pairNames)
between <- runFamily(betweenMasks, betweenId)

# Print and save the tables ----

showColumns <- c('label', 'nConnections', 'nScans', 'LowRisk.Mean', 'HighRisk.Mean', 'Diff.LowMinusHigh', 'Diff.inSD', 'Grp.p', 'Grp.p.adj', 'GA.p', 'Intx.p', 'Intx.p.adj')
nmWidth <- options(width = 220)
cat('\n=== Network-level connectivity (average of connections; score = ', networkScore, ') ===\n', sep = '')
cat('Diff.LowMinusHigh: Low minus High Risk, adjusted for age (positive = Low Risk higher). Diff.inSD: that difference in SD of the scores.\n',
    'p.adj: FDR within the table (', nrow(within$fits), ' networks; ', nrow(between$fits), ' network pairs). Intx.p is from the model with age centered at the mean.\n', sep = '')
cat('\nWITHIN networks\n')
print(format(within$fits[, showColumns], digits = 3), row.names = FALSE)
cat('\nBETWEEN networks\n')
print(format(between$fits[, showColumns], digits = 3), row.names = FALSE)
options(nmWidth)

saveNetworkAveragesTable <- function(tbl, name) {
  resultTables[[name]] <<- tbl
  if (exists('saveTables') && isTRUE(saveTables)) {
    dir.create(file.path(tableDir, networkAveragesFolder), showWarnings = FALSE, recursive = TRUE)
    write.csv(tbl, file.path(tableDir, networkAveragesFolder, paste0(name, '.csv')), row.names = FALSE)
  }
  invisible(tbl)
}
if (!exists('resultTables')) { resultTables <- list() }
saveNetworkAveragesTable(within$fits, 'networkAverages_within')
saveNetworkAveragesTable(between$fits, 'networkAverages_between')

# Figures ----

starFor <- function(p) ifelse(p < .05, '*', ifelse(p < .10, '+', ''))
scoreAxisLabel <- paste0('Average connectivity', ifelse(networkScore == 'absmean', ' (absolute values)', ''))

# Bars: adjusted group means for every score in a family
familyBars <- function(fits, xAngle = 0) {
  labels <- fits$label
  barData <- rbind(data.frame(label = labels, Group = 'Low Risk', Mean = fits$LowRisk.Mean, SE = fits$LowRisk.SE),
                   data.frame(label = labels, Group = 'High Risk', Mean = fits$HighRisk.Mean, SE = fits$HighRisk.SE))
  barData$label <- factor(barData$label, levels = labels)
  barData$Group <- factor(barData$Group, levels = c('Low Risk', 'High Risk'))
  starData <- data.frame(label = factor(labels, levels = labels),
                         y = pmax(fits$LowRisk.Mean + fits$LowRisk.SE, fits$HighRisk.Mean + fits$HighRisk.SE), star = starFor(fits$Grp.p.adj))
  starData$y <- starData$y + .04 * diff(range(c(0, barData$Mean + barData$SE)))
  ggplot2::ggplot(barData, ggplot2::aes(x = label, y = Mean, fill = Group, color = Group)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(.8), width = .72, linewidth = .8) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Mean - SE, ymax = Mean + SE), position = ggplot2::position_dodge(.8), width = .25, color = 'black', linewidth = .6) +
    ggplot2::geom_text(data = starData, ggplot2::aes(x = label, y = y, label = star), inherit.aes = FALSE, size = 7) +
    ggplot2::scale_fill_manual(values = c('Low Risk' = networkBarFill[['low']], 'High Risk' = networkBarFill[['high']]), name = NULL) +
    ggplot2::scale_color_manual(values = c('Low Risk' = groupColors[['low']], 'High Risk' = groupColors[['high']]), name = NULL) +
    ggplot2::labs(x = NULL, y = scoreAxisLabel, caption = '* p.adj < .05, + p.adj < .10 (FDR within this figure); age-adjusted means +/- SE') +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(legend.position = 'bottom', axis.text.y = ggplot2::element_text(color = 'black'),
                   axis.text.x = ggplot2::element_text(color = 'black', angle = xAngle, hjust = ifelse(xAngle == 0, .5, 1), vjust = ifelse(xAngle == 0, 1, .5)))
}

# Score against age: one panel per score, points per scan and the fitted line per group (interaction model)
familyByAge <- function(family) {
  fits <- family$fits; stripNames <- character(0); ageLines <- list(); agePoints <- list()
  for (i in seq_len(nrow(fits))) {
    d <- family$frames[[fits$network[i]]]
    fmt <- function(p) ifelse(p < .001, '<.001', sub('^0', '', sprintf('%.3f', p)))
    strip <- paste0(fits$label[i], ' (', fits$nConnections[i], ' connections)\nGroup p.adj = ', fmt(fits$Grp.p.adj[i]), ', Interaction p.adj = ', fmt(fits$Intx.p.adj[i]))
    stripNames <- c(stripNames, strip)
    d$ScanGAc <- d$ScanGA - mean(d$ScanGA)
    fit <- nlme::lme(FC ~ Grp * ScanGAc, random = ~1 | factor(shortID), data = d,
                     control = nlme::lmeControl(msMaxIter = 200, tolerance = 1e-4, msTol = 1e-5, opt = 'optim'))
    grid <- do.call(rbind, lapply(levels(d$Grp), function(g) {
      r <- range(d$ScanGA[d$Grp == g])
      data.frame(Grp = factor(g, levels = levels(d$Grp)), ScanGA = seq(r[1], r[2], length.out = 40))
    }))
    grid$ScanGAc <- grid$ScanGA - mean(d$ScanGA)
    grid$FC <- predict(fit, grid, level = 0)
    grid$panel <- strip; d$panel <- strip
    ageLines[[i]] <- grid[, c('panel', 'Grp', 'ScanGA', 'FC')]; agePoints[[i]] <- d[, c('panel', 'Grp', 'ScanGA', 'FC')]
  }
  ageLines <- do.call(rbind, ageLines); agePoints <- do.call(rbind, agePoints)
  ageLines$panel <- factor(ageLines$panel, levels = stripNames); agePoints$panel <- factor(agePoints$panel, levels = stripNames)
  ggplot2::ggplot(mapping = ggplot2::aes(x = ScanGA, y = FC, color = Grp)) +
    ggplot2::geom_point(data = agePoints, alpha = .5, size = 1.1) +
    ggplot2::geom_line(data = ageLines, linewidth = 1.1) +
    ggplot2::facet_wrap(~panel, ncol = 4, scales = 'free_y') +
    ggplot2::scale_color_manual(values = c(TD = groupColors[['low']], ASD = groupColors[['high']]), labels = c(TD = 'Low Risk', ASD = 'High Risk'), name = NULL) +
    ggplot2::labs(x = 'Postmenstrual age at scan (weeks)', y = scoreAxisLabel) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 8, hjust = 0), strip.background = ggplot2::element_blank(), legend.position = 'bottom',
                   panel.grid.major = ggplot2::element_line(color = '#e6e6e6'))
}

showFigure(familyBars(within$fits), 'networkAverages_within_bars', width = 9, height = 6, subfolder = networkAveragesFolder)
showFigure(familyByAge(within), 'networkAverages_within_by_age', width = 14, height = 2.2 + 2.6 * ceiling(nrow(within$fits) / 4), subfolder = networkAveragesFolder)
showFigure(familyBars(between$fits, xAngle = 45), 'networkAverages_between_bars', width = 14, height = 8, subfolder = networkAveragesFolder)
showFigure(familyByAge(between), 'networkAverages_between_by_age', width = 14, height = 2.2 + 2.6 * ceiling(nrow(between$fits) / 4), subfolder = networkAveragesFolder)

# One matrix: Low minus High Risk difference for every network (diagonal = within) and network pair (off the diagonal)
cells <- rbind(
  data.frame(row = within$fits$network, col = within$fits$network, diff = within$fits$Diff.LowMinusHigh, star = starFor(within$fits$Grp.p.adj), stringsAsFactors = FALSE),
  do.call(rbind, lapply(seq_along(pairs), function(i) {
    r <- between$fits[between$fits$network == pairNames[i], ]
    if (nrow(r) == 0) { return(NULL) }
    data.frame(row = c(pairs[[i]][1], pairs[[i]][2]), col = c(pairs[[i]][2], pairs[[i]][1]), diff = r$Diff.LowMinusHigh, star = starFor(r$Grp.p.adj), stringsAsFactors = FALSE)
  })))
cells$row <- factor(cells$row, levels = networkMeanNetworks, labels = vapply(networkMeanNetworks, networkLabel, character(1)))
cells$col <- factor(cells$col, levels = rev(networkMeanNetworks), labels = rev(vapply(networkMeanNetworks, networkLabel, character(1))))
cells$text <- paste0(sprintf('%.3f', cells$diff), '\n', cells$star)
limit <- max(abs(cells$diff))
networkAverages_matrix <- ggplot2::ggplot(cells, ggplot2::aes(x = row, y = col, fill = diff)) +
  ggplot2::geom_tile(color = 'white', linewidth = 1) +
  ggplot2::geom_text(ggplot2::aes(label = text), size = 3.6, lineheight = .9) +
  ggplot2::scale_fill_gradient2(low = groupColors[['high']], mid = 'white', high = groupColors[['low']], midpoint = 0, limits = c(-limit, limit),
                                name = 'Low minus High\nRisk difference') +
  ggplot2::labs(x = NULL, y = NULL, caption = 'Diagonal = within the network; off-diagonal = between two networks. Blue = Low Risk higher, red = High Risk higher.\n* p.adj < .05, + p.adj < .10 (FDR within the within-network and between-network families separately); age-adjusted') +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(color = 'black', angle = 45, hjust = 1), axis.text.y = ggplot2::element_text(color = 'black'))
showFigure(networkAverages_matrix, 'networkAverages_matrix', width = 9, height = 7.5, subfolder = networkAveragesFolder)
