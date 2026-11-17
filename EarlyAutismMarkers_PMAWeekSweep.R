# PMA CENTERING SWEEP: VOLUMES ----
#
# Refits the tissue and regional volume models (raw and normalized volume; with and without the group x age interaction) with scan age
# (PMA) centered at every week from 20 to 40, plus at the sample mean, and shows how the group-effect results change as the centering
# moves. In the interaction models the group term is the Low vs High Risk difference at the centering age, so this traces the
# difference across ages. (The main-effects model and the interaction term do not depend on the centering; they are shown for reference.)
#
# HOW TO RUN: first run EarlyAutismMarkers_edit.R so the volume data and functions exist in the R session, then run this whole file.
# Fits are cached in cacheDir (one per centering), so a second run is quick.
#
# OUTPUT: a count table in the console, tables in tableDir/PMACheck (only if saveTables is TRUE), and two figures in figDir/PMACheck:
#   PMACheck_volume_week_sweep_counts    how many ROIs are significant (p.adj < .05) or trending (< .10) at each centering week
#   PMACheck_volume_week_sweep_heatmap   every ROI that is significant or trending at some week: its status and direction at every week

# Settings ----

sweepWeeks <- 20:40                # centering ages in weeks (assumes scan age is recorded in weeks)
pmaFigureTypes <- c('pdf')         # file types for the figures ('tiff', 'pdf', 'svg', 'eps'); pdf because the heatmap can be large
sweepFolder <- 'PMACheck'          # subfolder of the figure and table folders

# Make sure the data and functions from EarlyAutismMarkers_edit.R are available ----

sweepNeeded <- c('fetalMatched', 'tissueData', 'tissueDataRaw', 'tissueExclude', 'regionalData', 'regionalDataRaw', 'regionalExclude',
                 'prepVolumeAnalysis', 'fitRoiModels', 'runRoiModels', 'addFDR', 'cached', 'showFigure', 'groupColors', 'figDir', 'tableDir')
sweepMissing <- sweepNeeded[!vapply(sweepNeeded, exists, logical(1))]
if (length(sweepMissing) > 0) { stop('Run EarlyAutismMarkers_edit.R first. These objects are missing: ', paste(sweepMissing, collapse = ', ')) }

# Fit the volume models at each centering ----

# center = a week, or NA for the mean scan age of each ROI's analysis sample
fitVolumesAt <- function(center) {
  label <- ifelse(is.na(center), 'mean', paste0('wk', center))
  centered <- function(prepFun) function(i) {
    d <- prepFun(i)
    d$scanAge <- d$scanAge - ifelse(is.na(center), mean(d$scanAge, na.rm = TRUE), center)
    d
  }
  cached(paste0('pmasweep_volume_', label),
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

centersToRun <- c(sweepWeeks, NA)
sweepLabels <- ifelse(is.na(centersToRun), 'mean', paste0('wk', centersToRun))
message('Fitting volume models at ', length(centersToRun), ' centerings...')
sweepFits <- lapply(seq_along(centersToRun), function(i) { message('  ', sweepLabels[i]); fitVolumesAt(centersToRun[i]) })
names(sweepFits) <- sweepLabels

# One long table: family x ROI x centering ----

families <- list('Tissue, raw volume' = c('tissue', 'Raw'), 'Tissue, normalized' = c('tissue', 'Normed'),
                 'Regional, raw volume' = c('regional', 'Raw'), 'Regional, normalized' = c('regional', 'Normed'))

sweepRows <- do.call(rbind, lapply(names(families), function(fam) {
  part <- families[[fam]][1]; model <- families[[fam]][2]
  do.call(rbind, lapply(seq_along(sweepLabels), function(i) {
    itx <- addFDR(sweepFits[[i]][[part]][[paste0(model, 'Itx')]])
    data.frame(family = fam, ROI = itx$ROI, centering = sweepLabels[i], week = centersToRun[i],
               Grp.B = itx$Grp.B, Grp.p.adj = itx$Grp.p.Adj, Itx.p.adj = itx$Itx.p.Adj,
               TD.Mean = itx$TD.Mean, ASD.Mean = itx$ASD.Mean, stringsAsFactors = FALSE)
  }))
}))
# the model without an interaction (the same at every centering), for reference
mainRows <- do.call(rbind, lapply(names(families), function(fam) {
  part <- families[[fam]][1]; model <- families[[fam]][2]
  main <- addFDR(sweepFits[[1]][[part]][[model]])
  data.frame(family = fam, ROI = main$ROI, Grp.B = main$Grp.B, Grp.p.adj = main$Grp.p.Adj, stringsAsFactors = FALSE)
}))

statusOf <- function(p) ifelse(!is.na(p) & p < .05, 'sig', ifelse(!is.na(p) & p < .10, 'trend', 'ns'))
sweepRows$status <- statusOf(sweepRows$Grp.p.adj)
sweepRows$direction <- ifelse(sweepRows$Grp.B > 0, 'Low', 'High')        # Grp.B > 0 means Low Risk (Pass) higher at that age

# Counts table ----

counts <- do.call(rbind, lapply(split(sweepRows, list(sweepRows$family, sweepRows$centering), drop = TRUE), function(d) {
  data.frame(family = d$family[1], centering = d$centering[1], week = d$week[1], n = nrow(d),
             group05 = sum(d$Grp.p.adj < .05, na.rm = TRUE), group10 = sum(d$Grp.p.adj < .10, na.rm = TRUE),
             interaction05 = sum(d$Itx.p.adj < .05, na.rm = TRUE), interaction10 = sum(d$Itx.p.adj < .10, na.rm = TRUE),
             lowHigher10 = sum(d$Grp.p.adj < .10 & d$direction == 'Low', na.rm = TRUE), highHigher10 = sum(d$Grp.p.adj < .10 & d$direction == 'High', na.rm = TRUE),
             stringsAsFactors = FALSE)
}))
rownames(counts) <- NULL
counts$centering <- factor(counts$centering, levels = sweepLabels)
counts <- counts[order(counts$family, counts$centering), ]

wide <- data.frame(centering = factor(sweepLabels, levels = sweepLabels), stringsAsFactors = FALSE)
for (fam in names(families)) {
  rows <- counts[counts$family == fam, ]
  wide[[fam]] <- paste(rows$group05[match(wide$centering, rows$centering)], '/', rows$group10[match(wide$centering, rows$centering)])
}
mainCounts <- vapply(names(families), function(fam) { d <- mainRows[mainRows$family == fam, ]; paste(sum(d$Grp.p.adj < .05, na.rm = TRUE), '/', sum(d$Grp.p.adj < .10, na.rm = TRUE)) }, character(1))
wide$centering <- as.character(wide$centering)
wide[nrow(wide) + 1, ] <- c('no interaction', mainCounts)
sweepWidth <- options(width = 200)
cat('\n=== Volume group effect by centering week: number of ROIs with FDR p < .05 / < .10 ===\n')
cat('(group term of the interaction models; the last row is the model without an interaction, which does not depend on centering)\n')
print(wide, row.names = FALSE)
cat('\nInteraction term (same at every centering): ',
    paste(names(families), vapply(names(families), function(fam) { d <- counts[counts$family == fam, ]; paste(d$interaction05[1], '/', d$interaction10[1]) }, character(1)), sep = ' = ', collapse = '; '), '\n')
options(sweepWidth)

if (exists('saveTables') && isTRUE(saveTables)) {
  sweepTableDir <- file.path(tableDir, sweepFolder)
  dir.create(sweepTableDir, showWarnings = FALSE, recursive = TRUE)
  write.csv(counts, file.path(sweepTableDir, 'PMACheck_volume_week_sweep_counts.csv'), row.names = FALSE)
  write.csv(sweepRows, file.path(sweepTableDir, 'PMACheck_volume_week_sweep_rois.csv'), row.names = FALSE)
}

# Figures ----

# Where the data are: the mean and the range of scan age in the matched sample (used to mark the plots)
meanAge <- mean(as.numeric(fetalMatched$scanAge), na.rm = TRUE)
ageRange <- range(as.numeric(fetalMatched$scanAge), na.rm = TRUE)

saveSweepFigure <- function(plotObj, name, width, height) {
  oldSave <- saveFigures; oldTypes <- figureFileTypes
  saveFigures <<- TRUE; figureFileTypes <<- pmaFigureTypes
  on.exit({ saveFigures <<- oldSave; figureFileTypes <<- oldTypes })
  showFigure(plotObj, name, width = width, height = height, subfolder = sweepFolder)
}

# 1. Counts across centering weeks
weekCounts <- counts[!is.na(counts$week), ]
countLong <- rbind(data.frame(weekCounts[, c('family', 'week')], count = weekCounts$group05, threshold = 'p.adj < .05'),
                   data.frame(weekCounts[, c('family', 'week')], count = weekCounts$group10, threshold = 'p.adj < .10'))
countLong$family <- factor(countLong$family, levels = names(families))
sweep_counts <- ggplot2::ggplot(countLong, ggplot2::aes(x = week, y = count, color = threshold, linetype = threshold)) +
  ggplot2::annotate('rect', xmin = max(min(sweepWeeks), ageRange[1]), xmax = min(max(sweepWeeks), ageRange[2]), ymin = -Inf, ymax = Inf, fill = '#eeeeee') +
  ggplot2::geom_vline(xintercept = meanAge, linetype = 'dotted', color = '#444444') +
  ggplot2::geom_line(linewidth = 1.1) + ggplot2::geom_point(size = 2) +
  ggplot2::facet_wrap(~family, ncol = 2, scales = 'free_y') +
  ggplot2::scale_x_continuous(breaks = seq(min(sweepWeeks), max(sweepWeeks), by = 2)) +
  ggplot2::scale_color_manual(values = c('p.adj < .05' = '#222222', 'p.adj < .10' = '#888888'), name = NULL) +
  ggplot2::scale_linetype_manual(values = c('p.adj < .05' = 'solid', 'p.adj < .10' = 'dashed'), name = NULL) +
  ggplot2::labs(x = 'Age at which scan age is centered (weeks)', y = 'Number of ROIs with a group difference',
                title = 'Volume group effect as the centering age moves',
                caption = paste0('Group term of the interaction models. Gray band = observed scan ages (', round(ageRange[1], 1), ' to ', round(ageRange[2], 1),
                                 ' weeks); dotted line = mean scan age (', round(meanAge, 1), ').')) +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::theme(legend.position = 'bottom', strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = 'bold', hjust = 0),
                 panel.grid.major = ggplot2::element_line(color = '#eeeeee'))
saveSweepFigure(sweep_counts, 'PMACheck_volume_week_sweep_counts', width = 12, height = 8)

# 2. Heatmap: every ROI that is significant or trending at some week
heat <- sweepRows[!is.na(sweepRows$week), ]
heat$cell <- ifelse(heat$status == 'ns', 'Not significant',
                    paste0(ifelse(heat$direction == 'Low', 'Low Risk higher', 'High Risk higher'), ifelse(heat$status == 'sig', ', p.adj < .05', ', p.adj < .10')))
cellLevels <- c('Low Risk higher, p.adj < .05', 'Low Risk higher, p.adj < .10', 'Not significant', 'High Risk higher, p.adj < .10', 'High Risk higher, p.adj < .05')
heat$cell <- factor(heat$cell, levels = cellLevels)
hits <- unique(heat[heat$status != 'ns', c('family', 'ROI')])
if (nrow(hits) == 0) {
  message('No volume ROI is significant or trending at any centering week, so no heatmap was made.')
} else {
  heat <- merge(heat, hits, by = c('family', 'ROI'))
  nSig <- tapply(heat$status != 'ns', paste(heat$family, heat$ROI), sum)
  heat$rowKey <- paste(heat$family, heat$ROI)
  heat$ROI <- factor(heat$ROI, levels = unique(heat$ROI[order(-nSig[heat$rowKey], heat$ROI)]))
  heat$family <- factor(heat$family, levels = names(families))
  rowsPerFamily <- tapply(hits$ROI, factor(hits$family, levels = names(families)), length)
  sweep_heatmap <- ggplot2::ggplot(heat, ggplot2::aes(x = week, y = ROI, fill = cell)) +
    ggplot2::geom_tile(color = 'white', linewidth = .4) +
    ggplot2::geom_vline(xintercept = meanAge, linetype = 'dotted', color = '#222222') +
    ggplot2::facet_grid(family ~ ., scales = 'free_y', space = 'free_y') +
    ggplot2::scale_x_continuous(breaks = seq(min(sweepWeeks), max(sweepWeeks), by = 2), expand = c(0, 0)) +
    ggplot2::scale_fill_manual(values = c('Low Risk higher, p.adj < .05' = groupColors[['low']], 'Low Risk higher, p.adj < .10' = '#9DBFE0', 'Not significant' = '#f0f0f0',
                                          'High Risk higher, p.adj < .10' = '#E5A3A8', 'High Risk higher, p.adj < .05' = groupColors[['high']]),
                               drop = FALSE, name = NULL) +
    ggplot2::labs(x = 'Age at which scan age is centered (weeks)', y = NULL, title = 'Volume ROIs with a group difference at some centering age',
                  caption = paste0('Dotted line = mean scan age (', round(meanAge, 1), ' weeks); observed scan ages ', round(ageRange[1], 1), ' to ', round(ageRange[2], 1), '.')) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = 'bottom', strip.text.y = ggplot2::element_text(angle = 0, hjust = 0, face = 'bold'),
                   axis.text.y = ggplot2::element_text(color = 'black', size = 8))
  saveSweepFigure(sweep_heatmap, 'PMACheck_volume_week_sweep_heatmap', width = 11, height = 2.5 + .22 * nrow(hits) + .5 * length(unique(hits$family)))
}
