# NETWORK-LEVEL CONNECTIVITY (exploratory) ----
#
# Instead of testing every ROI-to-ROI connection separately, this averages the connectivity of all connections WITHIN a network
# (both regions in the network, the same connections used by the network analyses) into one score per scan, then compares the groups on
# that score. One test per network, so far fewer tests to correct for.
#
# For each network it fits  (1) network score ~ Group + age   (primary: the average group difference)
#                           (2) network score ~ Group * age   (age centered at the sample mean: does the difference change with age?)
# with a random intercept per participant (the same models as the connection-level analyses), FDR-corrected across networks.
#
# HOW TO RUN: first run EarlyAutismMarkers_edit.R (so the data and functions are in the R session), then run this whole file.
# OUTPUT: a table printed in the console (and saved as network_means_models if saveTables is TRUE) and two figures:
#   network_mean_bars    adjusted group means (+/- SE) for every network
#   network_mean_by_age  the network score against scan age, per network, with the fitted line for each group

# Settings ----

# How the connections are combined: 'mean' = plain average of the connection values;
#   'absmean' = average of their absolute values (strength regardless of sign)
networkScore <- 'mean'

# Networks to analyze (names in networkSpecs); NULL = all of them
networkMeanNetworks <- NULL

# Make sure the data and functions from EarlyAutismMarkers_edit.R are available ----

nmNeeded <- c('tdFC', 'asdFC', 'tdDemoFinal', 'asdDemoFinal', 'connROIs', 'networkSpecs', 'groupColors', 'showFigure', 'saveTable', 'figDir')
nmMissing <- nmNeeded[!vapply(nmNeeded, exists, logical(1))]
if (length(nmMissing) > 0) { stop('Run EarlyAutismMarkers_edit.R first. These objects are missing: ', paste(nmMissing, collapse = ', ')) }
if (!networkScore %in% c('mean', 'absmean')) { stop("networkScore must be 'mean' or 'absmean'") }
if (is.null(networkMeanNetworks)) { networkMeanNetworks <- names(networkSpecs) }

# One score per scan for a network ----

networkMeanData <- function(net) {
  both <- networkSpecs[[net]]$match(connROIs$ROI1) & networkSpecs[[net]]$match(connROIs$ROI2)
  if (sum(both) == 0) { return(NULL) }
  one <- function(fc, demo, grp) {
    values <- as.matrix(fc[, both, drop = FALSE])
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
  attr(d, 'nConnections') <- sum(both)
  d
}

# Fit both models for one network ----

fitNetworkMean <- function(d) {
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

networkMeanFrames <- lapply(setNames(networkMeanNetworks, networkMeanNetworks), networkMeanData)
networkMeanFrames <- networkMeanFrames[!vapply(networkMeanFrames, is.null, logical(1))]
networkMeanFits <- do.call(rbind, lapply(names(networkMeanFrames), function(net) {
  cbind(network = net, fitNetworkMean(networkMeanFrames[[net]]), stringsAsFactors = FALSE)
}))
networkMeanFits$Grp.p.adj <- p.adjust(networkMeanFits$Grp.p, method = 'fdr')
networkMeanFits$Intx.p.adj <- p.adjust(networkMeanFits$Intx.p, method = 'fdr')
networkMeanFits$GA.p.adj <- p.adjust(networkMeanFits$GA.p, method = 'fdr')
rownames(networkMeanFits) <- NULL

cat('\n=== Network-level connectivity (average of the within-network connections; score = ', networkScore, ') ===\n', sep = '')
cat('Diff.LowMinusHigh: Low minus High Risk, adjusted for age (positive = Low Risk higher). Diff.inSD: that difference in SD of the scores.\n',
    'Grp.p.adj / Intx.p.adj / GA.p.adj: FDR across the ', nrow(networkMeanFits), ' networks. Intx.p is from the model with age centered at the mean.\n', sep = '')
nmWidth <- options(width = 220)
print(format(networkMeanFits[, c('network', 'nConnections', 'nScans', 'LowRisk.Mean', 'HighRisk.Mean', 'Diff.LowMinusHigh', 'Diff.inSD',
                                 'Grp.p', 'Grp.p.adj', 'GA.p', 'Intx.p', 'Intx.p.adj')], digits = 3), row.names = FALSE)
options(nmWidth)
saveTable(networkMeanFits, 'network_means_models')

# Figure 1: adjusted group means per network ----

networkLabel <- function(net) {
  if (exists('networkNodeLabels') && net %in% names(networkNodeLabels)) networkNodeLabels[[net]] else net
}
starFor <- function(p) ifelse(p < .05, '*', ifelse(p < .10, '+', ''))

barData <- rbind(data.frame(network = networkMeanFits$network, Group = 'Low Risk', Mean = networkMeanFits$LowRisk.Mean, SE = networkMeanFits$LowRisk.SE),
                 data.frame(network = networkMeanFits$network, Group = 'High Risk', Mean = networkMeanFits$HighRisk.Mean, SE = networkMeanFits$HighRisk.SE))
barData$network <- factor(barData$network, levels = networkMeanFits$network, labels = vapply(networkMeanFits$network, networkLabel, character(1)))
barData$Group <- factor(barData$Group, levels = c('Low Risk', 'High Risk'))
starData <- data.frame(network = factor(networkMeanFits$network, levels = networkMeanFits$network, labels = vapply(networkMeanFits$network, networkLabel, character(1))),
                       y = pmax(networkMeanFits$LowRisk.Mean + networkMeanFits$LowRisk.SE, networkMeanFits$HighRisk.Mean + networkMeanFits$HighRisk.SE),
                       star = starFor(networkMeanFits$Grp.p.adj))
starData$y <- starData$y + .04 * diff(range(c(0, barData$Mean + barData$SE)))

network_mean_bars <- ggplot2::ggplot(barData, ggplot2::aes(x = network, y = Mean, fill = Group, color = Group)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(.8), width = .72, linewidth = .8) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = Mean - SE, ymax = Mean + SE), position = ggplot2::position_dodge(.8), width = .25, color = 'black', linewidth = .6) +
  ggplot2::geom_text(data = starData, ggplot2::aes(x = network, y = y, label = star), inherit.aes = FALSE, size = 7) +
  ggplot2::scale_fill_manual(values = c('Low Risk' = networkBarFill[['low']], 'High Risk' = networkBarFill[['high']]), name = NULL) +
  ggplot2::scale_color_manual(values = c('Low Risk' = groupColors[['low']], 'High Risk' = groupColors[['high']]), name = NULL) +
  ggplot2::labs(x = NULL, y = paste0('Average connectivity within network', ifelse(networkScore == 'absmean', ' (absolute values)', '')),
                caption = '* p.adj < .05, + p.adj < .10 (FDR across networks); age-adjusted means +/- SE') +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::theme(legend.position = 'bottom', axis.text.x = ggplot2::element_text(color = 'black'), axis.text.y = ggplot2::element_text(color = 'black'))
showFigure(network_mean_bars, 'network_mean_bars', width = 9, height = 6)

# Figure 2: network score against scan age ----

ageLines <- list(); agePoints <- list(); stripNames <- character(0)
for (net in networkMeanFits$network) {
  d <- networkMeanFrames[[net]]
  fitRow <- networkMeanFits[networkMeanFits$network == net, ]
  strip <- paste0(networkLabel(net), ' (', fitRow$nConnections, ' connections)\nGroup p.adj = ', ifelse(fitRow$Grp.p.adj < .001, '<.001', sub('^0', '', sprintf('%.3f', fitRow$Grp.p.adj))),
                  ', Interaction p.adj = ', ifelse(fitRow$Intx.p.adj < .001, '<.001', sub('^0', '', sprintf('%.3f', fitRow$Intx.p.adj))))
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
  ageLines[[net]] <- grid[, c('panel', 'Grp', 'ScanGA', 'FC')]
  agePoints[[net]] <- d[, c('panel', 'Grp', 'ScanGA', 'FC')]
}
ageLines <- do.call(rbind, ageLines); agePoints <- do.call(rbind, agePoints)
ageLines$panel <- factor(ageLines$panel, levels = stripNames); agePoints$panel <- factor(agePoints$panel, levels = stripNames)

network_mean_by_age <- ggplot2::ggplot(mapping = ggplot2::aes(x = ScanGA, y = FC, color = Grp)) +
  ggplot2::geom_point(data = agePoints, alpha = .5, size = 1.3) +
  ggplot2::geom_line(data = ageLines, linewidth = 1.1) +
  ggplot2::facet_wrap(~panel, ncol = 3, scales = 'free_y') +
  ggplot2::scale_color_manual(values = c(TD = groupColors[['low']], ASD = groupColors[['high']]), labels = c(TD = 'Low Risk', ASD = 'High Risk'), name = NULL) +
  ggplot2::labs(x = 'Postmenstrual age at scan (weeks)', y = paste0('Average connectivity within network', ifelse(networkScore == 'absmean', ' (absolute values)', ''))) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(strip.text = ggplot2::element_text(size = 9, hjust = 0), strip.background = ggplot2::element_blank(), legend.position = 'bottom',
                 panel.grid.major = ggplot2::element_line(color = '#e6e6e6'))
showFigure(network_mean_by_age, 'network_mean_by_age', width = 12, height = 3 + 2.6 * ceiling(nrow(networkMeanFits) / 3))
