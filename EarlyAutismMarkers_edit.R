# 0. SETUP ----


## Paths and Libraries ----

# HOME LAPTOP
baseDir <- 'C:/Users/Kevin/OneDrive - Childrens National Hospital/ComputerSyncing/04_Projects/01_Active/FASD_FetalASDRisk/data'

figDir <- "C:/Users/Kevin/OneDrive - Childrens National Hospital/ComputerSyncing/04_Projects/01_Active/FASD_FetalASDRisk/paper/figures"

library(ggplot2)
library(ggpubr)
# library(lmerTest)
library(reshape2)
library(MatchIt)
library(nlme)
library(circlize)
library(emmeans)

# effectsize is called as effectsize::... (not loaded with library), so make sure it is installed
if (!requireNamespace('effectsize', quietly = TRUE)) {
  stop("Please install the 'effectsize' package first: install.packages('effectsize')")
}


## Figure Output Settings ----

# Set saveFigures to TRUE to write every figure to figDir (existing files with the same name are overwritten).
saveFigures <- FALSE

# File types to save, in this order. Choose any of: 'tiff', 'pdf', 'svg', 'eps'
figureFileTypes <- c('tiff', 'pdf', 'svg', 'eps')

# The individual pieces of a composite figure (the chord diagrams and heatmaps inside wholebrain_composite) are saved in this
# subfolder of figDir; only the composite itself goes in figDir.
componentFolder <- 'components'

# Every ggplot figure is also kept here by name (e.g. figures[['regional_group_p05']]) so it can be re-shown or re-saved
figures <- list()


## Group Colors ----

# The Low Risk (TD) and High Risk (ASD) colors used in every figure. Lighter tints for bars, points and lines are derived from these.
groupColors <- c(low = '#2166AC', high = '#B2182B')

# Chord diagram links: opacity as a hex alpha ('aa' = 67%, '55' = 33%) and the width of each link's outline.
# Lower the opacity or the outline width for lighter-looking chords.
chordLinkAlpha <- 'aa'
chordLinkBorderWidth <- 3

## Glass Brain Settings ----

# Files for the glass brain figures (the region coordinates and the brain surface are in the same fetal atlas space)
glassBrainDir <- file.path(baseDir, 'functional', 'output', 'brainNet')
roiCoordinateFile <- file.path(glassBrainDir, 'nodeTable.node')                    # one row per region; the first three columns are x y z
brainMeshFile <- file.path(glassBrainDir, 'Segmentation_inner_cortex.mz3')         # brain surface mesh

# The coordinates and mesh are in AFNI orientation (x: positive = left, y: positive = posterior). TRUE converts them to the
# standard orientation (x: right, y: anterior, z: superior) before drawing.
coordinatesAreAFNI <- TRUE

# The top (axial) view: TRUE puts the back of the head at the top of the picture (as in the BrainNet Viewer figures this was
# matched to); FALSE puts the front at the top.
glassAxialPosteriorAtTop <- TRUE

# The four networks shown on the network glass brain: node color and legend name
networkNodeColors <- c(CRB = '#F2B134', DMN = '#E8305E', Limbic = '#14C28D', SubCort = '#1A8BB5')
networkNodeLabels <- c(CRB = 'Cerebellar', DMN = 'Default Mode', Limbic = 'Limbic', SubCort = 'Subcortical')

## Table Output Settings ----

# Set saveTables to TRUE to write the result tables (CSV) to tableDir (existing files with the same name are overwritten).
saveTables <- FALSE

# The tables folder sits next to the figures folder (paper/tables)
tableDir <- file.path(dirname(figDir), 'tables')

# Every exported table is also kept here by name (e.g. resultTables[['tissue_raw_main']]), whether or not it is saved
resultTables <- list()

## Speed Settings ----

# Number of CPU cores for the whole-brain models (1 = no parallel processing).
nCores <- max(1, parallel::detectCores() - 1)

# Slow steps save their results in cacheDir and reload them on later runs. A cache is reused only if its inputs
# and the model code are unchanged; set rerun to TRUE to ignore every cache and recompute everything.
cacheDir <- file.path(baseDir, 'cache')
rerun <- FALSE


## Helper Functions ----

# Save a figure to figDir (or a subfolder of it, e.g. subfolder = componentFolder) in each of figureFileTypes (only if saveFigures is TRUE).
# drawFun is a function with no arguments that draws the figure on the current device.
# eps files are written without transparency (the eps format cannot hold it); svg and pdf keep it.
saveFigure <- function(name, width, height, drawFun, subfolder = NULL) {
  if (!saveFigures) { return(invisible(NULL)) }

  badTypes <- setdiff(figureFileTypes, c('tiff', 'pdf', 'svg', 'eps'))
  if (length(badTypes) > 0) { stop('Unknown figure file type(s): ', paste(badTypes, collapse=', '), ". Use 'tiff', 'pdf', 'svg', 'eps'.") }

  saveDir <- if (is.null(subfolder)) figDir else file.path(figDir, subfolder)
  dir.create(saveDir, showWarnings = FALSE, recursive = TRUE)
  for (type in figureFileTypes) {
    file <- file.path(saveDir, paste0(name, '.', type))
    switch(type,
           tiff = tiff(file, width = width, height = height, units = 'in', res = 600, compression = 'lzw'),
           pdf  = pdf(file, width = width, height = height),
           svg  = svg(file, width = width, height = height),
           eps  = postscript(file, width = width, height = height, horizontal = FALSE, onefile = FALSE, paper = 'special'))
    # the eps device warns once per figure that semi-transparency is unsupported; that is expected, so hide it
    tryCatch(withCallingHandlers(drawFun(), warning = function(w) {
               if (type == 'eps' && grepl('semi-transparency', conditionMessage(w))) invokeRestart('muffleWarning')
             }), finally = dev.off())
  }
  invisible(NULL)
}

# Show a ggplot/ggarrange figure, keep it in figures[[name]], and save it (in `subfolder` of figDir if given). Does nothing for a NULL figure.
showFigure <- function(plotObj, name, width, height, subfolder = NULL) {
  if (is.null(plotObj)) { return(invisible(NULL)) }
  figures[[name]] <<- plotObj
  print(plotObj)
  saveFigure(name, width, height, function() print(plotObj), subfolder)
  invisible(plotObj)
}

# Same for a base-graphics figure (e.g. circlize chord diagram): drawFun draws it.
showBaseFigure <- function(name, width, height, drawFun, subfolder = NULL) {
  drawFun()
  saveFigure(name, width, height, drawFun, subfolder)
  invisible(NULL)
}


# Build the participant ID used in the imaging file names from a study_id and its study.
# Vectorized: pass the study_id and study columns of a data frame.
formatID <- function(study_id, study) {
  study_id <- as.character(study_id)
  study <- as.character(study)
  out <- character(length(study_id))

  abc <- study == 'ABC'
  chd <- study == 'CHD R01'
  pl <- study == 'Placenta'
  rescue <- !(abc | chd | pl)

  out[abc] <- paste0('fetus_00', study_id[abc])
  out[chd] <- paste0('fetus_0', study_id[chd])
  out[pl] <- paste0('pl_1', sprintf("%04d", as.numeric(substr(study_id[pl], 4, nchar(study_id[pl])))))

  rescueNum <- as.numeric(substr(study_id[rescue], 5, nchar(study_id[rescue])))
  out[rescue] <- paste0('rescue_', ifelse(rescueNum < 100, sprintf("%03d", rescueNum), sprintf("%04d", rescueNum)))

  out
}


# Convert a functional-data ID (e.g. 'fetus_02151_02', 'rescue_0175_01', 'pl_10149_01')
# into the study_id used in combinedASDData.csv (2151, 'RES-175', 'PL-149').
# Vectorized: pass the ID column. Returns numeric if every ID is a fetus ID, otherwise character.
formatFCStudyID <- function(ID) {
  ID <- as.character(ID)
  abbrev <- substr(ID, 0, 1)
  num <- as.numeric(vapply(strsplit(ID, '_', fixed = TRUE), function(x) x[2], character(1)))

  out <- ifelse(abbrev == 'f', as.character(num),
         ifelse(abbrev == 'r', paste0('RES-', sprintf("%03d", num)),
                               paste0('PL-', sprintf("%03d", num))))

  if (all(abbrev == 'f')) as.numeric(out) else out
}

# Prepare the analysis data frame for one ROI of a volume measure (tissue or regional).
# roiColumn: one-column data frame of the ROI's volumes (rownames = longID)
# matched:   the matched subject frame (fetalMatched)
# rawData:   data frame with longID and TBV, used for the normalized volume
# excludeIDs: longIDs to drop;  dropInf: drop rows where the normalized volume is Inf
prepVolumeAnalysis <- function(roiColumn, matched, rawData, excludeIDs, dropInf = FALSE) {
  colnames(roiColumn) <- 'Volume'
  roiColumn$Volume <- as.numeric(roiColumn$Volume)

  tempAnalyses <- merge(matched[,c('shortID','longID','scanAge', 'mchat_result')], roiColumn, by="row.names")
  tempAnalyses <- tempAnalyses[complete.cases(tempAnalyses),]
  tempAnalyses$mchat_result <- factor(tempAnalyses$mchat_result)
  tempAnalyses$scanAge <- as.numeric(tempAnalyses$scanAge)

  tempAnalyses <- tempAnalyses[!(tempAnalyses$longID %in% excludeIDs),]
  tempAnalyses <- merge(tempAnalyses, rawData, by='longID', all.x=T, all.y=F)
  tempAnalyses$VolNormed <- tempAnalyses$Volume/tempAnalyses$TBV
  if (dropInf) { tempAnalyses <- tempAnalyses[!tempAnalyses$VolNormed==Inf,] }

  tempAnalyses
}

# Fit the four standard models for one ROI: raw and normalized outcome, each with and
# without the group x scanAge interaction (random intercept per subject).
# Returns a list of four one-row data frames: Raw, RawItx, Normed, NormedItx.
# pEta2 is the group term for the main-effects models and the interaction term for the Itx models.
fitRoiModels <- function(tempAnalyses, roi, rawVar, normVar) {

  lmeCtrl <- lmeControl(msMaxIter=200, tolerance=1e-4, msTol=1e-5, opt='optim')

  fitRow <- function(yVar, interaction) {
    op <- if (interaction) '*' else '+'
    model <- lme(fixed = as.formula(paste(yVar, '~ mchat_result', op, 'scanAge')), data = tempAnalyses,
                 random = ~1|factor(shortID), control = lmeCtrl)
    tTable <- summary(model)$tTable
    eta <- effectsize::eta_squared(model)$Eta2_partial

    passMean <- mean(tempAnalyses[[yVar]][which(tempAnalyses$mchat_result == 'Pass')])
    failMean <- mean(tempAnalyses[[yVar]][which(tempAnalyses$mchat_result == 'Fail')])

    out <- data.frame(ROI = roi, TD.Mean = passMean, ASD.Mean = failMean, TD.ASD.Diff = passMean - failMean,
                      Grp.B = tTable[2,1], Grp.p = tTable[2,5], scanGA.B = tTable[3,1], scanGA.p = tTable[3,5],
                      stringsAsFactors = FALSE)
    if (interaction) {
      out$Itx.B <- tTable[4,1]
      out$Itx.p <- tTable[4,5]
      out$pEta2 <- eta[3]
    } else {
      out$pEta2 <- eta[1]
    }
    out
  }

  list(Raw = fitRow(rawVar, FALSE), RawItx = fitRow(rawVar, TRUE),
       Normed = fitRow(normVar, FALSE), NormedItx = fitRow(normVar, TRUE))
}

# Run fitRoiModels over every ROI. prepFun(i) must return the analysis data frame for ROI number i.
# Returns the four result tables (Raw, RawItx, Normed, NormedItx), one row per ROI.
runRoiModels <- function(roiNames, prepFun, rawVar, normVar) {
  results <- lapply(seq_along(roiNames), function(i) fitRoiModels(prepFun(i), roiNames[i], rawVar, normVar))
  combine <- function(nm) {
    tbl <- do.call(rbind, lapply(results, `[[`, nm))
    rownames(tbl) <- NULL
    tbl
  }
  list(Raw = combine('Raw'), RawItx = combine('RawItx'), Normed = combine('Normed'), NormedItx = combine('NormedItx'))
}


# Read a wide CSV whose first column holds the row names (used for the functional connectivity matrices).
# read.csv is very slow on files with ~20,000 columns, so this uses data.table::fread when it is installed
# (about 30x faster) and falls back to read.csv otherwise. The parsed result is cached, so an unchanged file
# loads in well under a second on later runs. Returns the same data frame as read.csv(file, row.names = 1), except that
# fread can differ from read.csv in the last binary digit (~1e-16) of a few values. Use useFread = FALSE for exact read.csv parsing.
readWideCSV <- function(file, useFread = requireNamespace('data.table', quietly = TRUE)) {
  cached(paste0('csv_', gsub('[^A-Za-z0-9]+', '_', basename(file))),
         inputs = list(unname(tools::md5sum(file)), 'readWideCSV v1'),
         compute = function() {
           if (useFread) {
             raw <- data.table::fread(file, data.table = FALSE, colClasses = list(character = 1))
             out <- raw[, -1, drop = FALSE]
             rownames(out) <- raw[[1]]
             colnames(out) <- make.names(colnames(out), unique = TRUE)
             out
           } else {
             read.csv(file, row.names = 1)
           }
         })
}

# FDR-correct the p-values in an ROI results table (adds Grp.p.Adj, scanGA.p.Adj, and Itx.p.Adj for interaction tables).
addFDR <- function(tbl) {
  tbl$Grp.p.Adj <- p.adjust(tbl$Grp.p, method='fdr')
  tbl$scanGA.p.Adj <- p.adjust(tbl$scanGA.p, method='fdr')
  if ('Itx.p' %in% colnames(tbl)) { tbl$Itx.p.Adj <- p.adjust(tbl$Itx.p, method='fdr') }
  tbl
}

# Print the ROIs that survive FDR (group effect, plus interaction effect for interaction tables).
reportSignificant <- function(tbl) {
  print(subset(tbl, Grp.p.Adj < .05))
  if ('Itx.p.Adj' %in% colnames(tbl)) { print(subset(tbl, Itx.p.Adj < .05)) }
}


# A lighter version of a color: mix it with white (amount 0 = the color, 1 = white), with an optional hex alpha such as '55'
tintColor <- function(color, amount, alpha = '') {
  rgbValues <- col2rgb(color)[, 1]
  mixed <- rgbValues * (1 - amount) + 255 * amount
  paste0(rgb(mixed[1], mixed[2], mixed[3], maxColorValue = 255), alpha)
}

# Tints of the group colors: filled bars, translucent points, fitted lines, and network bars
groupBarFill <- c(low = tintColor(groupColors[['low']], .64), high = tintColor(groupColors[['high']], .64))
groupPointFill <- c(low = tintColor(groupColors[['low']], .64, '55'), high = tintColor(groupColors[['high']], .64, '55'))
groupLineColor <- c(low = tintColor(groupColors[['low']], .2), high = tintColor(groupColors[['high']], .2))
networkBarFill <- c(low = tintColor(groupColors[['low']], .5, 'cc'), high = tintColor(groupColors[['high']], .5, 'cc'))

# Axis numbers with a fixed number of decimals (never scientific notation): the fewest decimals that show every tick exactly
fixedDecimalLabels <- function(breaks) {
  digits <- 0
  while (digits < 8 && any(abs(breaks * 10^digits - round(breaks * 10^digits)) > 1e-8, na.rm = TRUE)) { digits <- digits + 1 }
  formatC(breaks, format = 'f', digits = digits)
}

# A boxed legend (Low Risk / High Risk) drawn as its own panel, like the key in the older hand-formatted figure.
# type = 'group' (filled bars), 'network' (bars with colored outlines) or 'intx' (diamond / circle with a line);
# textSize = size of the 'Low Risk' / 'High Risk' text
riskLegendPanel <- function(type, textSize = 7.5) {
  legendPlot <- ggplot() + xlim(0, 1) + ylim(0, 1) + coord_fixed() + theme_void() +
    annotate('rect', xmin = .08, xmax = .92, ymin = .2, ymax = .8, fill = 'white', colour = 'black', linewidth = 1) +
    annotate('text', x = .40, y = .61, label = 'Low Risk', hjust = 0, size = textSize) +
    annotate('text', x = .40, y = .39, label = 'High Risk', hjust = 0, size = textSize)

  if (type == 'group') {
    legendPlot +
      annotate('rect', xmin = .16, xmax = .32, ymin = .53, ymax = .69, fill = groupBarFill[['low']], colour = 'black', linewidth = 1) +
      annotate('rect', xmin = .16, xmax = .32, ymin = .31, ymax = .47, fill = groupBarFill[['high']], colour = 'black', linewidth = 1)
  } else if (type == 'network') {
    legendPlot +
      annotate('rect', xmin = .16, xmax = .32, ymin = .53, ymax = .69, fill = networkBarFill[['low']], colour = groupColors[['low']], linewidth = 1) +
      annotate('rect', xmin = .16, xmax = .32, ymin = .31, ymax = .47, fill = networkBarFill[['high']], colour = groupColors[['high']], linewidth = 1)
  } else {
    legendPlot +
      annotate('segment', x = .13, xend = .35, y = .61, yend = .61, colour = groupLineColor[['low']], linewidth = 2) +
      annotate('point', x = .24, y = .61, shape = 23, size = 7, fill = groupPointFill[['low']], colour = 'black', stroke = .75) +
      annotate('segment', x = .13, xend = .35, y = .39, yend = .39, colour = groupLineColor[['high']], linewidth = 2) +
      annotate('point', x = .24, y = .39, shape = 21, size = 7, fill = groupPointFill[['high']], colour = 'black', stroke = .75)
  }
}

# Put the boxed Low Risk / High Risk legend to the right of a single plot
addRiskLegend <- function(plotObj, legendType) {
  ggpubr::ggarrange(plotObj, riskLegendPanel(legendType, textSize = 5), ncol = 2, widths = c(5, 1.3), legend = 'none')
}

# Full region names (used on the chord diagram and the connection charts), keyed by the abbreviation without its hemisphere
# suffix (e.g. 'ACG' for 'ACG.L').
# Regions not listed here keep their abbreviation.
roiFullNames <- c(ACG = 'Anterior Cingulate', IFGoperc = 'Inferior Operculum', IFGtriang = 'Inferior Triangle', MFG = 'Middle Frontal',
                  ORBinf = 'Inferior Orbital', ORBmed = 'Medial Orbital', ORBmid = 'Middle Orbital', ORBsupb = 'Superior Orbital',
                  PreCG = 'Precentral', REC = 'Rectus', SFGdor = 'Superior Frontal (Dorsal)', SFGmed = 'Superior Frontal (Medial)',
                  SMA = 'Supplemental Motor',
                  INS = 'Insula', IPL = 'Inferior Parietal Lobule', MCG = 'Middle Cingulate', PCG = 'Posterior Cingulate',
                  PCL = 'Paracentral Lobule', PCUN = 'Precuneus', PoCG = 'Postcentral', SMG = 'Supramarginal',
                  SPG = 'Superior Parietal', HES = "Heschl's",
                  FFG = 'Fusiform', ITG = 'Inferior Temporal', MTG = 'Middle Temporal', STG = 'Superior Temporal',
                  TPOmid = 'Temporal Pole', TPOsup = 'Superior Temporal Pole',
                  CAL = 'Calcarine', CUN = 'Cuneus', LING = 'Lingual', IOG = 'Inferior Occipital', MOG = 'Middle Occipital',
                  SOG = 'Superior Occipital',
                  AMY = 'Amygdala', CAU = 'Caudate', HIP = 'Hippocampus', PAL = 'Pallidum', PHG = 'Parahippocampal',
                  PUT = 'Putamen', THA = 'Thalamus', CRB = 'Cerebellum', BS = 'Brainstem')

roiFullName <- function(sector) {
  abbrev <- sub('\\.[LR]$', '', sector)
  ifelse(abbrev %in% names(roiFullNames), roiFullNames[abbrev], sector)
}

# Label for a connection ID such as 'CAU.L.39_THA.R.139': both regions with hemisphere, on two lines,
# e.g. "Left Caudate -\nRight Thalamus". Connections that would share a label within `ids` also get their region numbers,
# e.g. "Left Cerebellum -\nRight Cerebellum (#13, #142)". Returns a vector the same length as ids.
connectionLabel <- function(ids) {
  regionLabel <- function(region) {
    fields <- strsplit(region, '.', fixed = TRUE)[[1]]            # abbreviation, hemisphere (none for the brainstem), number
    hasHemisphere <- length(fields) >= 3 && fields[2] %in% c('L', 'R')
    name <- roiFullName(if (hasHemisphere) paste(fields[1], fields[2], sep = '.') else fields[1])
    if (hasHemisphere) paste(c(L = 'Left', R = 'Right')[[fields[2]]], name) else name
  }
  regionNumber <- function(region) { fields <- strsplit(region, '.', fixed = TRUE)[[1]]; fields[length(fields)] }

  regions <- strsplit(ids, '_', fixed = TRUE)
  first <- vapply(regions, function(x) regionLabel(x[1]), character(1))
  second <- vapply(regions, function(x) regionLabel(x[2]), character(1))
  base <- paste0(first, ' -\n', second)

  ambiguous <- duplicated(base) | duplicated(base, fromLast = TRUE)
  numbers <- vapply(regions, function(x) paste0(' (#', regionNumber(x[1]), ', #', regionNumber(x[2]), ')'), character(1))
  ifelse(ambiguous, paste0(base, numbers), base)
}

# ---- Glass brain helpers ----

# Read a brain surface mesh in MZ3 format (the format used by Surf Ice / MRIcroGL), in base R. Returns the vertex coordinates
# (one row per vertex) and the triangles (three vertex numbers per row).
readMz3 <- function(file) {
  con <- gzfile(file, 'rb')
  on.exit(close(con))
  magic <- readBin(con, 'integer', 1, size = 2, signed = FALSE)
  attrs <- readBin(con, 'integer', 1, size = 2, signed = FALSE)
  nFace <- readBin(con, 'integer', 1, size = 4)
  nVert <- readBin(con, 'integer', 1, size = 4)
  nSkip <- readBin(con, 'integer', 1, size = 4)
  if (magic != 23117) { stop('Not an MZ3 mesh file: ', file) }
  if (nSkip > 0) { readBin(con, 'raw', nSkip) }
  faces <- matrix(readBin(con, 'integer', nFace * 3, size = 4), ncol = 3, byrow = TRUE) + 1L
  vertices <- matrix(readBin(con, 'numeric', nVert * 3, size = 4), ncol = 3, byrow = TRUE)
  list(vertices = vertices, faces = faces)
}

# Coordinates (columns x, y, z) as a data frame with columns X (right), Y (anterior), Z (superior).
# AFNI orientation has x positive to the left and y positive to the posterior, so those two are flipped.
standardCoordinates <- function(xyz, fromAFNI = coordinatesAreAFNI) {
  xyz <- as.matrix(xyz[, 1:3])        # first three columns only (a BrainNet .node file also has color, size and label columns)
  data.frame(X = xyz[, 1] * ifelse(fromAFNI, -1, 1), Y = xyz[, 2] * ifelse(fromAFNI, -1, 1), Z = xyz[, 3])
}

# How each glass brain view maps the coordinates onto the page. Lateral = left side of the head (front on the left),
# posterior = seen from behind, axial = seen from above or below (see glassAxialPosteriorAtTop).
glassViews <- list(lateral   = list(h = list(axis = 'Y', sign = -1), v = list(axis = 'Z', sign = 1)),
                   posterior = list(h = list(axis = 'X', sign = 1),  v = list(axis = 'Z', sign = 1)),
                   axial     = list(h = list(axis = 'X', sign = 1),  v = list(axis = 'Y', sign = ifelse(glassAxialPosteriorAtTop, -1, 1))))

# The translucent brain behind the nodes: the mesh surface area behind each pixel (smoothed), for every view.
# Folds show through where the surface overlaps itself, which gives the glass look.
glassBackgrounds <- function(mesh, pixels = 300, blur = 1.4, pad = 4) {
  vertices <- standardCoordinates(mesh$vertices)
  tri <- mesh$faces
  a <- as.matrix(vertices[tri[, 2], ] - vertices[tri[, 1], ])
  b <- as.matrix(vertices[tri[, 3], ] - vertices[tri[, 1], ])
  triArea <- .5 * sqrt((a[, 2] * b[, 3] - a[, 3] * b[, 2])^2 + (a[, 3] * b[, 1] - a[, 1] * b[, 3])^2 + (a[, 1] * b[, 2] - a[, 2] * b[, 1])^2)
  vertexArea <- as.vector(tapply(rep(triArea / 3, 3), c(tri), sum))      # a third of each adjacent triangle

  kernel <- dnorm(seq(-4, 4), 0, blur); kernel <- kernel / sum(kernel)
  smooth <- function(m) { m <- t(apply(m, 1, function(r) stats::filter(r, kernel, sides = 2))); m[is.na(m)] <- 0; m }

  lapply(glassViews, function(view) {
    h <- vertices[[view$h$axis]] * view$h$sign
    v <- vertices[[view$v$axis]] * view$v$sign
    hBreaks <- seq(min(h) - pad, max(h) + pad, length.out = pixels)
    vBreaks <- seq(min(v) - pad, max(v) + pad, length.out = pixels)
    pixel <- (findInterval(h, hBreaks, all.inside = TRUE) - 1L) * pixels + findInterval(v, vBreaks, all.inside = TRUE)
    density <- matrix(0, pixels, pixels)                                   # rows = vertical, columns = horizontal
    summed <- tapply(vertexArea, pixel, sum)
    density[as.integer(names(summed))] <- summed
    density <- smooth(t(smooth(t(density))))
    tiles <- expand.grid(h = hBreaks, v = vBreaks)
    tiles$alpha <- as.vector(t(density))
    tiles$alpha <- .7 * pmin(1, (tiles$alpha / quantile(tiles$alpha[tiles$alpha > 0], .97))^.6)   # .7 = how dark the glass is
    tiles[tiles$alpha > 0.005, ]
  })
}

# One glass brain view: the translucent brain, edges, and nodes.
#   view = 'lateral', 'posterior' or 'axial';  backgrounds = the result of glassBackgrounds();  coords = standardCoordinates() of the regions
#   nodes = data frame with the region number (id) and the node color (fill);  edges = data frame with region numbers a and b and dir ('low' or 'high')
glassBrainView <- function(view, backgrounds, coords, nodes, edges = NULL, nodeSize = 2, edgeWidth = 1) {
  vw <- glassViews[[view]]
  hOf <- function(ids) coords[[vw$h$axis]][ids] * vw$h$sign
  vOf <- function(ids) coords[[vw$v$axis]][ids] * vw$v$sign

  p <- ggplot() + geom_raster(data = backgrounds[[view]], aes(h, v, alpha = alpha), fill = '#555555') + scale_alpha_identity()
  if (!is.null(edges) && nrow(edges) > 0) {
    edgeData <- data.frame(x = hOf(edges$a), y = vOf(edges$a), xend = hOf(edges$b), yend = vOf(edges$b), dir = edges$dir)
    p <- p + geom_segment(data = edgeData, aes(x, y, xend = xend, yend = yend, color = dir), linewidth = edgeWidth, alpha = .8, lineend = 'round') +
      scale_color_manual(values = c(low = groupLineColor[['low']], high = groupLineColor[['high']]), guide = 'none')
  }
  nodeData <- data.frame(h = hOf(nodes$id), v = vOf(nodes$id), fill = nodes$fill)
  p + geom_point(data = nodeData, aes(h, v), size = nodeSize + .7, color = '#202020', alpha = .3) +
    geom_point(data = nodeData, aes(h, v, fill = fill), size = nodeSize, shape = 21, color = '#303030', stroke = .4) +
    scale_fill_identity() + coord_fixed() + theme_void()
}

# Load the region coordinates and draw the translucent brain (once; later calls do nothing). Creates roiCoordinates and glassBrainBackgrounds.
# The brain surface is cached, since it does not change.
loadGlassBrain <- function() {
  for (file in c(roiCoordinateFile, brainMeshFile)) {
    if (!file.exists(file)) { stop('Glass brain file not found: ', file, '\nSet roiCoordinateFile and brainMeshFile in the Glass Brain Settings at the top of the script.') }
  }
  if (!exists('roiCoordinates', envir = globalenv())) {
    roiCoordinates <<- standardCoordinates(read.table(roiCoordinateFile))
  }
  if (!exists('glassBrainBackgrounds', envir = globalenv())) {
    glassBrainBackgrounds <<- cached('glass_brain_backgrounds',
                                     inputs = list(unname(tools::md5sum(brainMeshFile)), deparse(glassBackgrounds), deparse(glassViews), coordinatesAreAFNI),
                                     compute = function() glassBackgrounds(readMz3(brainMeshFile)))
  }
  invisible(NULL)
}

# Region numbers at the two ends of each connection ID such as 'CAU.L.39_THA.R.139'
connectionRegions <- function(ids) {
  regionNumbers <- lapply(strsplit(ids, '_', fixed = TRUE), function(x) as.integer(sub('.*\\.', '', x)))
  data.frame(a = vapply(regionNumbers, `[`, integer(1), 1), b = vapply(regionNumbers, `[`, integer(1), 2))
}

# Keep a table in resultTables[[name]] and, if saveTables is TRUE, write it to tableDir/<name>.csv
saveTable <- function(tbl, name) {
  resultTables[[name]] <<- tbl
  if (saveTables) {
    dir.create(tableDir, showWarnings = FALSE, recursive = TRUE)
    write.csv(tbl, file.path(tableDir, paste0(name, '.csv')), row.names = FALSE)
  }
  invisible(tbl)
}

# Save a result table and, next to it, <name>_sigFDR: the rows where the group (or interaction) effect passes FDR at .05.
# Works for the ROI tables (Grp.p.Adj / Itx.p.Adj) and the connectivity tables (Grp.p.adjust / Intx.p.adjust).
saveResultTable <- function(tbl, name) {
  saveTable(tbl, name)
  adjCols <- intersect(c('Grp.p.Adj', 'Itx.p.Adj', 'Grp.p.adjust', 'Intx.p.adjust'), colnames(tbl))
  if (length(adjCols) > 0) {
    passes <- Reduce(`|`, lapply(adjCols, function(col) tbl[[col]] < .05 & !is.na(tbl[[col]])))
    saveTable(tbl[passes,], paste0(name, '_sigFDR'))
  }
  invisible(tbl)
}

# Fingerprint of any R objects (used to tell whether a cached result is still valid)
objectHash <- function(...) {
  tmp <- tempfile()
  saveRDS(list(...), tmp, compress = FALSE)
  hash <- unname(tools::md5sum(tmp))
  unlink(tmp)
  hash
}

# Return compute()'s result, saved in cacheDir under `name`. The saved copy is reused only when `inputs`
# (a list of everything the result depends on) is unchanged and rerun is FALSE.
cached <- function(name, inputs, compute) {
  cacheFile <- file.path(cacheDir, paste0(name, '.rds'))
  key <- objectHash(inputs)

  if (!rerun && file.exists(cacheFile)) {
    saved <- readRDS(cacheFile)
    if (identical(saved$key, key)) {
      message('Loaded cached result: ', name)
      return(saved$value)
    }
    message('Cache out of date, recomputing: ', name)
  }

  value <- compute()
  dir.create(cacheDir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(list(key = key, value = value), cacheFile)
  value
}

# FDR-correct the p-values of a functional connectivity table (Grp.p.adjust, GA.p.adjust, and Intx.p.adjust for the interaction table)
addFDRfc <- function(tbl) {
  tbl$Grp.p.adjust <- p.adjust(tbl$Grp.p, method='fdr')
  tbl$GA.p.adjust <- p.adjust(tbl$GA.p, method='fdr')
  if ('Intx.p' %in% colnames(tbl)) { tbl$Intx.p.adjust <- p.adjust(tbl$Intx.p, method='fdr') }
  tbl
}


# 1. SHARED DATA ----

# (Nothing is shared between the structural and functional analyses: each reads its own data below.)


# 2. STRUCTURAL ANALYSES (run after Setup + Shared Data) ----


## 2a. Subject Selection and Matching ----

### Anatomical Data ----


testData <- read.csv(paste0(baseDir,'/anatomical/combinedASDData.csv'),stringsAsFactors=F)
testFetal <- subset(testData, PretermStatus=='Term' & CHDStatus !='CHD' & AnyFetal==T)

testFetalFail <- subset(testFetal, mchat_result=='Fail')


fetalFailID <- formatID(testFetalFail$study_id, testFetalFail$study)


tissueDir <- file.path(baseDir, "anatomical", "tissueVolumes_Individual")
regionalDir <- file.path(baseDir, "anatomical", "regionalVolumes_Individual")
vermisVol <- read.csv(file.path(dirname(baseDir), "_Earlier_Iteration_NeurodevelopmentalRisk", "data", "anatomicalData", "vermisParcVolumes_071026.csv"), stringsAsFactors=F)
vermisVol$totalVermis <- vermisVol$supVerm + vermisVol$infVerm + vermisVol$midVerm

tissueFiles <- list.files(tissueDir)
fetalFailID[!(fetalFailID %in% substr(tissueFiles, 0, nchar(tissueFiles)-19))]

regionalFiles <- list.files(regionalDir)
fetalFailID[!(fetalFailID %in% substr(regionalFiles, 0, nchar(regionalFiles)-21))]


anatTimePoints <- as.data.frame(matrix(ncol=3, nrow=length(fetalFailID)))
colnames(anatTimePoints) <- c('ID','tissueFiles','regionalFiles')
anatTimePoints$ID <- fetalFailID

for(i in 1:nrow(anatTimePoints)){
  targetID <- anatTimePoints$ID[i]
  anatTimePoints$tissueFiles[i] <- sum(substr(tissueFiles, 0, nchar(tissueFiles)-19) %in% targetID)
  anatTimePoints$regionalFiles[i] <- sum(substr(regionalFiles, 0, nchar(regionalFiles)-21) %in% targetID)
}


### Anatomical Matching ----

# IMPORTANT NOTE - ORIGINALLY WE MATCHED WITHIN COHORT WHEN POSSIBLE.

# TO AVOID ISSUES WITH HIGH RISK PREGNANCIES IN THE PL COHORT THE CODE HAS BEEN WRITTEN TO HARD CODE SOME PARTICIPANT LEVEL SWAPS TO REPLACE THE MCHAT PASSES FROM THE PLACENTA COHORT


fetalMchatMelt <- melt(testFetal, measure.vars=c('ScanAge1', 'ScanAge2'), id.vars=c('study_id', 'mchat_result','study'))
fetalMchatMelt <- subset(fetalMchatMelt, mchat_result!='' & value!='')
rownames(fetalMchatMelt) <- 1:nrow(fetalMchatMelt)
colnames(fetalMchatMelt) <- c('study_id','mchat_result','study','scanTimepoint','scanAge')

# fetalMchatMelt <- fetalMchatMelt[!c(fetalMchatMelt$study == 'Placenta' & fetalMchatMelt$mchat_result == 'Pass'),]

# fetalMchatMelt <- subset(fetalMchatMelt, study_id != 'RES-175' | scanTimepoint != 'ScanAge1')
# fetalMchatMelt <- subset(fetalMchatMelt, study_id != 'RES-184' | scanTimepoint != 'ScanAge1')
# fetalMchatMelt <- subset(fetalMchatMelt, study_id != 'RES-2175' | scanTimepoint != 'ScanAge1')

set.seed(1000)
fetalMchatMatch <- matchit(data=fetalMchatMelt, factor(mchat_result, levels=c('Pass','Fail')) ~ scanAge + study)

fetalMatched <- fetalMchatMelt[c(fetalMchatMatch$match.matrix, rownames(!is.na(fetalMchatMatch$match.matrix))[!is.na(fetalMchatMatch$match.matrix)]),]
fetalMatched <- subset(fetalMatched, !is.na(study_id ))

t.test(subset(fetalMatched, mchat_result=='Fail')$scanAge, subset(fetalMatched, mchat_result=='Pass')$scanAge)

fetalMatched$shortID <- formatID(fetalMatched$study_id, fetalMatched$study)

fetalMatched$scanTimepoint <- as.character(fetalMatched$scanTimepoint)
fetalMatched$longID <- paste0(fetalMatched$shortID,'_0', substr(fetalMatched$scanTimepoint,nchar(fetalMatched$scanTimepoint),1000))
rownames(fetalMatched) <- fetalMatched$longID


fetalMatched[!(fetalMatched$longID %in% substr(tissueFiles, 0, nchar(tissueFiles)-16)),]
fetalMatched[!(fetalMatched$longID %in% substr(regionalFiles, 0, nchar(regionalFiles)-18)),]

# TRUE or FALSE to turn on the hand matching
HANDMATCH <- T

if(HANDMATCH) {
  
  ## HERE IS THE HAND MATCHING CODE
  #View(subset(fetalMchatMelt, mchat_result=='Pass')[!c(subset(fetalMchatMelt, mchat_result=='Pass')$study_id %in% fetalMatched$study_id),])
  
  # WHEN IN DOUBT CAN GO A LITTLE YOUNGER AND STILL BE IN PARITY
  #pl_10011_01 - 38.571 --- fetus_02151_02 - 38.571 - (F)
  #pl_10027_01 - 38.429 --- fetus_02201_02 - 38.429 - (F)
  #pl_10029_01 - 28.429 --- fetus_02020_01 - 28.286 - (M) | others nearby...
  #pl_10029_02 - 31.571 --- rescue_0142_02 - 31.571 - (F) | others nearby (see pl_10125_01 below)
  #pl_10051_01 - 23.571 --- fetus_00556_01 - 23.710 - (M) | others nearby...
  #pl_10091_01 - 28.286 --- fetus_02172_01 - 28.286 - (M) | others nearby...
  #pl_10106_01 - 35.143 --- fetus_02107_02 - 35.143 - (F) | RES-178_02
  #pl_10124_01 - 26.429 --- fetus_02023_01 - 26.429 - (F) | RES-193_01
  #pl_10125_01 - 35.571 --- fetus_02104_02 - 35.571 - (F) | others nearby...
  
  #fetus_02153_01 - 28.286 (F) --  
  
  #Identify who is to be removed and who is to be added
  removeList <- c('pl_10011_01','pl_10027_01','pl_10029_01','pl_10029_02','pl_10051_01','pl_10091_01','pl_10106_01','pl_10124_01','pl_10125_01','fetus_02153_01')
 
  addListShort <- c('2120', #pl_10011_01 match 1
                    #'2151', #pl_10011_01 match 2
                    
                    '2121', #pl_10027_01 match 1
                    #'2201', #pl_10027_01 match 2
                    
                    '2020',
                    'RES-142',
                    '556',
                    '2172',
                    '2107',
                    
                    'RES-193', #pl_10124_01 match 1
                    #'2023', #pl_10124_01 match 2
                    
                    '2104',
                    'RES-042')
  
  addListFull <- c('fetus_02120_02', #pl_10011_01 match 1
                   #'fetus_02151_02', #pl_10011_01 match 2
                   
                   'fetus_02121_02', #pl_10027_01 match 1
                   #'fetus_02201_02', #pl_10027_01 match 2
                   
                   'fetus_02020_01',
                   'rescue_0142_02',
                   
                   'fetus_00556_01', #pl_10051_01 match 1
                   #'RES-217', #pl_10051_01 match 2
                   #'RES-79', #pl_10051_01 match 3
                   
                   'fetus_02172_01',
                   'fetus_02107_02',
                   
                   'rescue_0193_01', #pl_10124_01 match 1
                   #'fetus_02023_01', #pl_10124_01 match 2
                   
                   'fetus_02104_02',
                   'rescue_0042_01')
  
  #Pull those to be added and format dataframe for merging
  shortAddFrame <- fetalMchatMelt[fetalMchatMelt$study_id %in% addListShort,]
  
  shortAddFrame$shortID <- formatID(shortAddFrame$study_id, shortAddFrame$study)
  
  shortAddFrame$scanTimepoint <- as.character(shortAddFrame$scanTimepoint)
  shortAddFrame$longID <- paste0(shortAddFrame$shortID,'_0', substr(shortAddFrame$scanTimepoint,nchar(shortAddFrame$scanTimepoint),1000))
  rownames(shortAddFrame) <- shortAddFrame$longID
  
  shortAddFrame <- shortAddFrame[shortAddFrame$longID %in% addListFull,]

  #Remove the excluded pl participants
  fetalMatched <- fetalMatched[!c(fetalMatched$longID %in% removeList),]
  #Add in the replacements
  fetalMatched <- rbind(fetalMatched, shortAddFrame)
}


#Now drap the bad pl fail scans and their associated controls

dropList <- c('pl_10028_02', 'pl_10149_01', 'pl_10149_02', 'pl_10151_01',
              'rescue_071_02', 'fetus_02075_01', 'rescue_017_02', 'fetus_02070_02')

fetalMatched <- fetalMatched[!c(fetalMatched$longID %in% dropList),]

saveTable(fetalMatched, 'structural_matchedSample')


length(unique(subset(fetalMatched, mchat_result=='Pass')$shortID))
length(unique(subset(fetalMatched, mchat_result=='Fail')$shortID))

length(subset(fetalMatched, mchat_result=='Pass')$shortID)
length(subset(fetalMatched, mchat_result=='Fail')$shortID)

t.test(subset(fetalMatched, mchat_result=='Pass')$scanAge,subset(fetalMatched, mchat_result=='Fail')$scanAge)
mean(subset(fetalMatched, mchat_result=='Pass')$scanAge)
sd(subset(fetalMatched, mchat_result=='Pass')$scanAge)
mean(subset(fetalMatched, mchat_result=='Fail')$scanAge)
sd(subset(fetalMatched, mchat_result=='Fail')$scanAge)


demoDataPaper <- merge(fetalMatched, testData[,c('study_id','Sex','BirthGA','X18MLTFAge','birthWeight', 'Race')], all.x=T, all.y=F, by='study_id')


t.test(subset(demoDataPaper, mchat_result=='Pass')$BirthGA,subset(demoDataPaper, mchat_result=='Fail')$BirthGA)
mean(subset(demoDataPaper, mchat_result=='Pass')$BirthGA, na.rm=T)
sd(subset(demoDataPaper, mchat_result=='Pass')$BirthGA, na.rm=T)
mean(subset(demoDataPaper, mchat_result=='Fail')$BirthGA, na.rm=T)
sd(subset(demoDataPaper, mchat_result=='Fail')$BirthGA, na.rm=T)


xiSqTable <- matrix(ncol=2, nrow=2)
xiSqTable[1,1] <- nrow(subset(demoDataPaper, Sex=='m' & mchat_result=='Pass'))
xiSqTable[1,2] <- nrow(subset(demoDataPaper, Sex=='f' & mchat_result=='Pass'))
xiSqTable[2,1] <- nrow(subset(demoDataPaper, Sex=='m' & mchat_result=='Fail'))
xiSqTable[2,2] <- nrow(subset(demoDataPaper, Sex=='f' & mchat_result=='Fail'))
chisq.test(xiSqTable)


t.test(subset(demoDataPaper, mchat_result=='Pass')$birthWeight,subset(demoDataPaper, mchat_result=='Fail')$birthWeight, paired=F)
mean(subset(demoDataPaper, mchat_result=='Pass')$birthWeight, na.rm=T)
sd(subset(demoDataPaper, mchat_result=='Pass')$birthWeight, na.rm=T)
mean(subset(demoDataPaper, mchat_result=='Fail')$birthWeight, na.rm=T)
sd(subset(demoDataPaper, mchat_result=='Fail')$birthWeight, na.rm=T)


t.test(subset(demoDataPaper, mchat_result=='Pass')$X18MLTFAge,subset(demoDataPaper, mchat_result=='Fail')$X18MLTFAge, paired=F)
mean(subset(demoDataPaper, mchat_result=='Pass')$X18MLTFAge, na.rm=T)
sd(subset(demoDataPaper, mchat_result=='Pass')$X18MLTFAge, na.rm=T)
mean(subset(demoDataPaper, mchat_result=='Fail')$X18MLTFAge, na.rm=T)
sd(subset(demoDataPaper, mchat_result=='Fail')$X18MLTFAge, na.rm=T)


contingency_table <- table(demoDataPaper$mchat_result, demoDataPaper$Race)
contingency_table[1,2] <- 8
contingency_table[2,1] <- 0
contingency_table[2,8] <- 3
contingency_table <- contingency_table[1:2, 2:9]
chisq.test(contingency_table)


## 2b. Tissue Volumes ----


tissueDataRaw <- read.csv(paste0(baseDir,'/anatomical/tissueVolumes_12122025_compiled.csv'))
tissueDataRaw$TBV <- rowSums(tissueDataRaw[,c('CorticalGrayMatter','CorticalWhiteMatter','Cerebellum','SubcorticalGray','Brainstem','AmygdalaHippocampus')], na.rm=T)
tissueData <- tissueDataRaw[,c('longID','CorticalGrayMatter','CorticalWhiteMatter','Cerebellum','SubcorticalGray','Brainstem','AmygdalaHippocampus')]
rownames(tissueData) <- tissueData$longID
tissueData <- subset(tissueData, select = -c(longID))


# Analyze trajectories across tissue ROI
tissueExclude <- c('rescue_0175_01', 'fetus_02178_01','rescue_0184_01')

# The fits are cached: recomputed only if the data, the exclusions or the model code change (or rerun is TRUE)
tissueModels <- cached('structural_tissue_models',
                       inputs = list(tissueData, tissueDataRaw, fetalMatched, tissueExclude,
                                     deparse(prepVolumeAnalysis), deparse(fitRoiModels), deparse(runRoiModels)),
                       compute = function() runRoiModels(colnames(tissueData),
                                                         function(i) prepVolumeAnalysis(tissueData[,i, drop=FALSE], fetalMatched, tissueDataRaw,
                                                                                        excludeIDs=tissueExclude),
                                                         rawVar='Volume', normVar='VolNormed'))

glmTissueVolTable <- tissueModels$Raw
glmTissueVolTable.Itx <- tissueModels$RawItx
glmTissueNormedVolTable <- tissueModels$Normed
glmTissueNormedVolTable.Itx <- tissueModels$NormedItx


# FDR correction and significant ROIs

glmTissueVolTable <- addFDR(glmTissueVolTable)
reportSignificant(glmTissueVolTable)
saveResultTable(glmTissueVolTable, 'tissue_raw_main')

glmTissueVolTable.Itx <- addFDR(glmTissueVolTable.Itx)
reportSignificant(glmTissueVolTable.Itx)
saveResultTable(glmTissueVolTable.Itx, 'tissue_raw_intx')

glmTissueNormedVolTable <- addFDR(glmTissueNormedVolTable)
reportSignificant(glmTissueNormedVolTable)
saveResultTable(glmTissueNormedVolTable, 'tissue_normed_main')

glmTissueNormedVolTable.Itx <- addFDR(glmTissueNormedVolTable.Itx)
reportSignificant(glmTissueNormedVolTable.Itx)
saveResultTable(glmTissueNormedVolTable.Itx, 'tissue_normed_intx')


## 2c. Regional Volumes ----


# Combine the regional volumes into a single data frame
regionalData <- as.data.frame(matrix(ncol = 87, nrow= nrow(fetalMatched)))
rownames(regionalData) <- fetalMatched$longID
regionalLabels<-read.table(paste0(baseDir,'/anatomical/regionalROILabels.txt'))
colnames(regionalData) <- t(regionalLabels)

for(i in 1:nrow(regionalData)){
  
  targetID <- rownames(regionalData)[i]
  targetFile <- regionalFiles[substr(regionalFiles, 0, nchar(regionalFiles)-18) %in% targetID]
  if(length(targetFile) == 0) {next}
  tempData <- read.table(file.path(regionalDir, targetFile))
  ncolData <- ncol(tempData)
  regionalData[i,1:(ncolData-2)] <- tempData[2,3:ncolData]
  
}

# Remove ROI to drop CSF and white matter
regionalDataRaw <- as.data.frame(lapply(regionalData, function(x) as.numeric(as.character(x))))
rownames(regionalDataRaw) <- rownames(regionalData)
regionalDataRaw$TBV <- rowSums(regionalDataRaw[,c(1:48,51:82)], na.rm=T)
regionalData <- regionalData[,1:82]
regionalData <- regionalData[,grep('WM',colnames(regionalData), invert = T)]
regionalData <- regionalData[,grep('LatVen',colnames(regionalData), invert = T)]
regionalData <- regionalData[,grep('CRB',colnames(regionalData), invert = T)]
regionalData <- regionalData[,grep('CorpusCallosum',colnames(regionalData), invert = T)]
regionalDataRaw$longID <- rownames(regionalDataRaw)


# FLAG WHETHER TO INCLUDE THE PARCELLATED CEREBELLUM WITHIN THE APPROACH

IncludeCerebellarParcellation <- F
if(IncludeCerebellarParcellation){
  
  
  regionalData$fullID <- rownames(regionalData)
  regionalData <- merge(regionalData,vermisVol, all=T, by='fullID')
  rownames(regionalData) <- regionalData$fullID
  regionalData$fullID <- NULL
  
  regionalDataRaw <- merge(regionalDataRaw,vermisVol, all=T, by.x='longID', by.y='fullID')
  rownames(regionalDataRaw) <- regionalDataRaw$longID
  
}


# Analyze trajectories across regional ROI
# NOTE: this analysis excludes two participants; the tissue analysis excludes three (fetus_02178_01 is kept here)
regionalExclude <- c('rescue_0175_01', 'rescue_0184_01')

regionalModels <- cached('structural_regional_models',
                         inputs = list(regionalData, regionalDataRaw, fetalMatched, regionalExclude,
                                       deparse(prepVolumeAnalysis), deparse(fitRoiModels), deparse(runRoiModels)),
                         compute = function() runRoiModels(colnames(regionalData),
                                                           function(i) prepVolumeAnalysis(regionalData[,i, drop=FALSE], fetalMatched, regionalDataRaw,
                                                                                          excludeIDs=regionalExclude, dropInf=TRUE),
                                                           rawVar='Volume', normVar='VolNormed'))

glmRegionalVolTable <- regionalModels$Raw
glmRegionalVolTable.Itx <- regionalModels$RawItx
glmRegionalNormedVolTable <- regionalModels$Normed
glmRegionalNormedVolTable.Itx <- regionalModels$NormedItx


# FDR correction and significant ROIs

glmRegionalVolTable <- addFDR(glmRegionalVolTable)
reportSignificant(glmRegionalVolTable)
saveResultTable(glmRegionalVolTable, 'regional_raw_main')

glmRegionalVolTable.Itx <- addFDR(glmRegionalVolTable.Itx)
reportSignificant(glmRegionalVolTable.Itx)
saveResultTable(glmRegionalVolTable.Itx, 'regional_raw_intx')

glmRegionalNormedVolTable <- addFDR(glmRegionalNormedVolTable)
reportSignificant(glmRegionalNormedVolTable)
saveResultTable(glmRegionalNormedVolTable, 'regional_normed_main')

glmRegionalNormedVolTable.Itx <- addFDR(glmRegionalNormedVolTable.Itx)
reportSignificant(glmRegionalNormedVolTable.Itx)
saveResultTable(glmRegionalNormedVolTable.Itx, 'regional_normed_intx')


#### Regional Figure ----

# Figures are built for whichever regions are significant (FDR-adjusted p, from the normalized interaction model table),
# at two thresholds, for two effects:
#   Group main effect: estimated marginal means (bars) for Low vs High Risk
#   Interaction:       normalized volume across gestational age, with a line per risk group
# Change the thresholds here:
figureThresholds <- c(.05, .10)


# Optional hand-tuned axes. ROIs listed here keep these y-axis limits and labels; any other ROI that
# turns up as significant gets a label made from its column name and automatic y-axis limits.
# (yMax = upper y limit; breakMax and breakStep define the tick marks, seq(0, breakMax, breakStep))
regionalFigureROIs <- data.frame(
  roi            = c('HPC_L', 'HPC_R', 'Amyg_L', 'FusiformAnt_GM_R', 'Caudate_R', 'Subthalamic_L', 'Lentiform_L'),
  label          = c('Hippocampus (Left) Volume', 'Hippocampus (Right) Volume', 'Amygdala (Left) Volume',
                     'Fusiform (Right) Volume', 'Caudate (Right) Volume', 'Subthalamic (Left) Volume', 'Pallidum (Left) Volume'),
  groupYMax      = c(.003,  .003,  .002,  .002,  .01,  .0007,  .0105),
  groupBreakMax  = c(.003,  .003,  .003,  .003,  .03,  .03,    .03),
  groupBreakStep = c(.0005, .0005, .0005, .0005, .002, .0002,  .002),
  intxYMax       = c(.0047, .0047, .0037, .003,  .02,  .00125, .02),
  intxBreakMax   = c(.005,  .005,  .005,  .005,  .05,  .005,   .03),
  intxBreakStep  = c(.001,  .001,  .001,  .0005, .0025, .00025, .005),
  stringsAsFactors = FALSE)

# Participants dropped from the figures (the regional stats above drop only the first and third)
regionalFigureExclude <- c('rescue_0175_01', 'fetus_02178_01','rescue_0184_01')

# ROIs whose adjusted p-value (pAdjCol) is below the threshold
significantROIs <- function(tbl, pAdjCol, threshold) {
  tbl$ROI[which(tbl[[pAdjCol]] < threshold)]
}

# Break a long axis label onto two lines so it fits beside the panel
wrapAxisLabel <- function(label, width = 28) paste(strwrap(label, width = width), collapse = '\n')

# Axis label and y-axis settings for one ROI: the hand-tuned values if it is in regionalFigureROIs, otherwise automatic.
# type is 'group' or 'intx'; maxValue is the largest value that must fit on the y axis (used only for automatic limits).
roiAxis <- function(roi, type, maxValue) {
  spec <- regionalFigureROIs[regionalFigureROIs$roi == roi,]
  yMaxCol <- paste0(type, 'YMax'); breakMaxCol <- paste0(type, 'BreakMax'); stepCol <- paste0(type, 'BreakStep')

  if (nrow(spec) == 1) {
    list(label = wrapAxisLabel(spec$label), limits = c(0, spec[[yMaxCol]]), breaks = seq(0, spec[[breakMaxCol]], spec[[stepCol]]))
  } else {
    breaks <- pretty(c(0, maxValue * 1.1), n = 6)
    list(label = wrapAxisLabel(paste(gsub('_', ' ', roi), 'Volume')), limits = c(0, max(breaks)), breaks = breaks)
  }
}

# Shared figure theme. The group (bar) figures have no x-axis title.
themeRegionalFigure <- function(showXTitle = TRUE) {
  theme(axis.title.x = if (showXTitle) element_text(color = 'black', size = 16) else element_blank(),
        axis.text.x = element_text(color = 'black', size = 12),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(color = 'black', size = 16),
        axis.text.y = element_text(color = 'black', size = 12),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(color='black', linewidth = 1),
        axis.line.y = element_line(color='black', linewidth = 1),
        plot.title = element_text(color = 'black', size = 24, hjust=0.5),
        panel.grid.major.x =  element_blank(),
        panel.grid.major.y = element_line(color='#cccccc', linewidth = 1),
        panel.background = element_rect(fill='white'),
        legend.position = 'none',
        legend.text = element_text(color = 'black', size = 12),
        legend.title=element_blank())
}

# Data for one regional ROI's figure
prepRegionalFigureData <- function(roi) {
  prepVolumeAnalysis(regionalData[,roi, drop=FALSE], fetalMatched, regionalDataRaw,
                     excludeIDs = regionalFigureExclude, dropInf = TRUE)
}

# Group main effect figure: model-estimated means (+/- SE) for Low vs High Risk
plotRoiGroup <- function(roi) {
  tempAnalyses <- prepRegionalFigureData(roi)

  model <- lme(data=tempAnalyses, VolNormed ~ mchat_result * scanAge, random= ~1|factor(shortID), control=lmeControl(msMaxIter=200, tolerance=1e-4, msTol=1e-5, opt='optim'))
  meansTable <- emmeans::emmeans(model, ~ mchat_result, data=tempAnalyses)

  meansTablePlot <- summary(meansTable)
  colnames(meansTablePlot)[1] <- 'Autism Risk'
  meansTablePlot$`Autism Risk` <- factor(ifelse(meansTablePlot$`Autism Risk`=='Fail','High Risk','Low Risk'), levels = c('Low Risk','High Risk'))

  axisInfo <- roiAxis(roi, 'group', max(meansTablePlot$emmean + meansTablePlot$SE))

  ggplot(data=meansTablePlot, aes(x=`Autism Risk`, y=emmean)) +
    scale_y_continuous(limits=axisInfo$limits, breaks=axisInfo$breaks, labels=fixedDecimalLabels, expand=c(0,0)) +
    scale_fill_manual(values=c(groupBarFill[['low']], groupBarFill[['high']])) +
    geom_bar(stat='identity',aes(fill=`Autism Risk`), color='#000000', linewidth=1) +
    geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.125, linewidth=2) +
    ylab(axisInfo$label) +
    themeRegionalFigure(showXTitle = FALSE)
}

# Interaction figure: normalized volume by gestational age, with a fitted line per risk group
plotRoiInteraction <- function(roi) {
  tempAnalyses <- prepRegionalFigureData(roi)
  # Low Risk first so it gets the first (blue) color, matching the group figure
  tempAnalyses$`Autism Risk` <- factor(ifelse(tempAnalyses$mchat_result=='Fail','High Risk','Low Risk'), levels = c('Low Risk','High Risk'))

  axisInfo <- roiAxis(roi, 'intx', max(tempAnalyses$VolNormed, na.rm=TRUE))

  ggplot(data=tempAnalyses, aes(x=scanAge, y=VolNormed)) +
    scale_y_continuous(limits=axisInfo$limits, breaks=axisInfo$breaks, labels=fixedDecimalLabels, expand=c(0,0)) +
    scale_x_continuous(limits=c(19.5,40.5), breaks=seq(20,40,5), expand=c(0,0)) +
    scale_shape_manual(values=c(23,21)) +   # Low Risk = diamond, High Risk = circle (as before)
    scale_fill_manual(values=c(groupPointFill[['low']], groupPointFill[['high']])) +
    geom_point(aes(fill=`Autism Risk`, shape=`Autism Risk`), color='#000000', size=4.5, stroke=.75) +
    scale_color_manual(values=c(groupLineColor[['low']], groupLineColor[['high']])) +
    geom_smooth(method='lm', formula='y~x', aes(color=`Autism Risk`), se=F, linewidth=2) +
    ylab(axisInfo$label) + xlab('Gestational Age (Weeks)') +
    themeRegionalFigure(showXTitle = TRUE)
}

# Size (inches) of a figure with one panel per ROI plus the legend: 4 panels per row, 3.5 in per panel
regionalFigureSize <- function(nROIs) {
  nPanels <- nROIs + 1
  list(width = 3.5 * min(4, nPanels), height = 3.5 * ceiling(nPanels / 4))
}

# Arrange one plot per ROI in a grid (4 columns wide) with the legend in the bottom-right cell. Prints a message and draws nothing if no ROI qualifies.
regionalFigure <- function(rois, plotFun, label, legendType) {
  if (length(rois) == 0) {
    message('No regions to plot for: ', label)
    return(invisible(NULL))
  }
  message(label, ': ', length(rois), ' region(s): ', paste(rois, collapse=', '))
  plots <- lapply(rois, plotFun)

  # the legend always sits in the bottom-right cell; any unused cells before it are left blank
  nCols <- min(4, length(plots) + 1)
  nRows <- ceiling((length(plots) + 1) / nCols)
  blankPanel <- ggplot() + theme_void()
  slots <- c(plots, rep(list(blankPanel), nCols * nRows - length(plots) - 1), list(riskLegendPanel(legendType)))

  ggpubr::ggarrange(plotlist = slots, nrow = nRows, ncol = nCols, legend = 'none')
}


## Significant regions

# Group effect and interaction p-values both come from the normalized interaction model table
groupROIs <- lapply(figureThresholds, function(thr) significantROIs(glmRegionalNormedVolTable.Itx, 'Grp.p.Adj', thr))
intxROIs <- lapply(figureThresholds, function(thr) significantROIs(glmRegionalNormedVolTable.Itx, 'Itx.p.Adj', thr))
names(groupROIs) <- names(intxROIs) <- paste0('p', figureThresholds)


## Group main effect figures

for (k in seq_along(figureThresholds)) {
  thrTag <- paste0('p', sprintf('%02d', round(figureThresholds[k]*100)))
  figSize <- regionalFigureSize(length(groupROIs[[k]]))
  showFigure(regionalFigure(groupROIs[[k]], plotRoiGroup, paste('Group effect, p.adj <', figureThresholds[k]), 'group'),
             paste0('regional_group_', thrTag), width = figSize$width, height = figSize$height)
}


## Interaction figures

for (k in seq_along(figureThresholds)) {
  thrTag <- paste0('p', sprintf('%02d', round(figureThresholds[k]*100)))
  figSize <- regionalFigureSize(length(intxROIs[[k]]))
  showFigure(regionalFigure(intxROIs[[k]], plotRoiInteraction, paste('Interaction, p.adj <', figureThresholds[k]), 'intx'),
             paste0('regional_intx_', thrTag), width = figSize$width, height = figSize$height)
}


## 2d. Cortical Thickness ----


ctDataRaw <- read.csv(paste0(baseDir,'/anatomical/CT_Ave_Thickness.csv'), header=F)
ctIDList <- read.csv(paste0(baseDir,'/anatomical/CT_SubList.1D'), header=F)
rownames(ctDataRaw) <- substr(ctIDList$V1,0, nchar(ctIDList$V1)-14)
ctROILabel <- read.csv(paste0(baseDir,'/anatomical/CTROILabels.txt'), header=F)
colnames(ctDataRaw) <- ctROILabel$V1
ctData <- ctDataRaw[,c(1:70,79:90)]
ctData$GlobalCT <- t(read.csv(paste0(baseDir,'/anatomical/CT_Global_Thickness.csv'), header=F))


# Analyze trajectories across CT ROI
prepCTAnalysis <- function(i) {
  ctVolumeROI <- ctData[,i, drop=FALSE]
  colnames(ctVolumeROI) <- 'CT'
  ctVolumeROI$CT <- as.numeric(ctVolumeROI$CT)
  ctVolumeROI$GlobalCT <- ctData$GlobalCT

  tempAnalyses <- merge(fetalMatched[,c('shortID','longID','scanAge', 'mchat_result')], ctVolumeROI, by="row.names")
  tempAnalyses <- tempAnalyses[complete.cases(tempAnalyses),]
  tempAnalyses$mchat_result <- factor(tempAnalyses$mchat_result)
  tempAnalyses$scanAge <- as.numeric(tempAnalyses$scanAge)

  tempAnalyses$CTNormed <- tempAnalyses$CT/tempAnalyses$GlobalCT
  tempAnalyses
}

ctModels <- cached('structural_ct_models',
                   inputs = list(ctData, fetalMatched, deparse(prepCTAnalysis), deparse(fitRoiModels), deparse(runRoiModels)),
                   compute = function() runRoiModels(colnames(ctData)[1:(ncol(ctData)-1)], prepCTAnalysis, rawVar='CT', normVar='CTNormed'))

glmCTTable <- ctModels$Raw
glmCTTable.Itx <- ctModels$RawItx
glmCTNormedTable <- ctModels$Normed
glmCTNormedTable.Itx <- ctModels$NormedItx


# FDR correction and significant ROIs

glmCTTable <- addFDR(glmCTTable)
reportSignificant(glmCTTable)
saveResultTable(glmCTTable, 'ct_raw_main')

glmCTTable.Itx <- addFDR(glmCTTable.Itx)
reportSignificant(glmCTTable.Itx)
saveResultTable(glmCTTable.Itx, 'ct_raw_intx')

glmCTNormedTable <- addFDR(glmCTNormedTable)
reportSignificant(glmCTNormedTable)
saveResultTable(glmCTNormedTable, 'ct_normed_main')

glmCTNormedTable.Itx <- addFDR(glmCTNormedTable.Itx)
reportSignificant(glmCTNormedTable.Itx)
saveResultTable(glmCTNormedTable.Itx, 'ct_normed_intx')


# 3. FUNCTIONAL ANALYSES (run after Setup + Shared Data) ----


## 3a. Subject Selection and Matching ----

### Functional Matching ----


#| echo: false

asdDemo <- read.csv(paste0(baseDir,'/functional/asdDemo_expand.csv'), stringsAsFactors = F)
tdDemo <- read.csv(paste0(baseDir,'/functional/tdDemo.csv'), stringsAsFactors = F)

asdDemo <- asdDemo[asdDemo$ID!='rescue_0187_01',] #Why is this excluded? Why are 715,748 included?
rownames(asdDemo) <- 1:nrow(asdDemo)


asdDemo$Grp <- 'ASD'
tdDemo$Grp <- 'TD'
combDemo <- rbind(asdDemo, tdDemo)

matchFrame <- MatchIt::matchit(as.factor(Grp) ~ ScanGA, data=combDemo)

tdDemoMatch <- tdDemo[matchFrame$match.matrix,]
tdDemoMatch <- tdDemoMatch[complete.cases(tdDemoMatch),]
tdDemoMatch <- tdDemoMatch[order(as.numeric(row.names(tdDemoMatch))), ]

dropPtList <- c('none')

keepRows <- as.numeric(rownames(asdDemo[!(asdDemo$ID %in% dropPtList), ]))
tdDemoFinal <- tdDemoMatch[keepRows,]
asdDemoFinal <- asdDemo[keepRows,]


mean(tdDemoFinal$ScanGA)
sd(tdDemoFinal$ScanGA)
mean(asdDemoFinal$ScanGA)
sd(asdDemoFinal$ScanGA)

t.test(tdDemoFinal$ScanGA, asdDemoFinal$ScanGA)


xiSqTable <- matrix(ncol=2, nrow=2)
xiSqTable[1,1] <- nrow(subset(tdDemoFinal, Sex=='m'))
xiSqTable[1,2] <- nrow(subset(tdDemoFinal, Sex=='f'))
xiSqTable[2,1] <- nrow(subset(asdDemoFinal, Sex=='m'))
xiSqTable[2,2] <- nrow(subset(asdDemoFinal, Sex=='f'))
chisq.test(xiSqTable)


tdBirthAge <- c(40.00,39.00,40.29,39.86,39.71,38.71,40.43,39.86,38.88,40.57,38.43,39.71,39.71,38.43,39.43,39.00,40.14,38.14,38.71)
asdBirthAge <- c(37.71,39.57,38.86,38.86,39.29,38.29,38.14,39.71,38.00,39.00,39.00,41.00,40.86,38.86,38.86,38.00,38.00,39.43,37.29)

t.test(tdBirthAge, asdBirthAge)

mean(tdBirthAge)
sd(tdBirthAge)
mean(asdBirthAge)
sd(asdBirthAge)

# tdDemoFinal <- tdDemoFinal[tdDemoFinal$ID != 'fetus_00715_02',]
tdDemoFinal <- tdDemoFinal[tdDemoFinal$ID != 'fetus_00748_02',]
tdDemoFinal <- tdDemoFinal[tdDemoFinal$ID != 'fetus_00556_01',]
asdDemoFinal <- asdDemoFinal[asdDemoFinal$ID != 'pl_0028_02',]

saveTable(rbind(asdDemoFinal, tdDemoFinal), 'functional_matchedSample')


#DO demos here!
fcMatchedDataTemp <- rbind(asdDemoFinal, tdDemoFinal)

fcMatchedDataTemp$study_id <- formatFCStudyID(fcMatchedDataTemp$ID)


testDataFC <- read.csv(paste0(baseDir,'/functional/combinedASDData.csv'),stringsAsFactors=F)
demoDataPaperFC <- merge(fcMatchedDataTemp, testDataFC[,c('study_id','Sex','BirthGA','X18MLTFAge','birthWeight', 'Race', 'mchat_result')], all.x=T, all.y=F, by='study_id')


xiSqTable <- matrix(ncol=2, nrow=2)
xiSqTable[1,1] <- nrow(subset(demoDataPaperFC, Sex.x=='m' & mchat_result=='Pass'))
xiSqTable[1,2] <- nrow(subset(demoDataPaperFC, Sex.x=='f' & mchat_result=='Pass'))
xiSqTable[2,1] <- nrow(subset(demoDataPaperFC, Sex.x=='m' & mchat_result=='Fail'))
xiSqTable[2,2] <- nrow(subset(demoDataPaperFC, Sex.x=='f' & mchat_result=='Fail'))
chisq.test(xiSqTable)


t.test(subset(demoDataPaperFC, mchat_result=='Pass')$birthWeight,subset(demoDataPaperFC, mchat_result=='Fail')$birthWeight, paired=F)
mean(subset(demoDataPaperFC, mchat_result=='Pass')$birthWeight, na.rm=T)
sd(subset(demoDataPaperFC, mchat_result=='Pass')$birthWeight, na.rm=T)
mean(subset(demoDataPaperFC, mchat_result=='Fail')$birthWeight, na.rm=T)
sd(subset(demoDataPaperFC, mchat_result=='Fail')$birthWeight, na.rm=T)


t.test(subset(demoDataPaperFC, mchat_result=='Pass')$BirthGA,subset(demoDataPaperFC, mchat_result=='Fail')$BirthGA, paired=F)
mean(subset(demoDataPaperFC, mchat_result=='Pass')$BirthGA, na.rm=T)
sd(subset(demoDataPaperFC, mchat_result=='Pass')$BirthGA, na.rm=T)
mean(subset(demoDataPaperFC, mchat_result=='Fail')$BirthGA, na.rm=T)
sd(subset(demoDataPaperFC, mchat_result=='Fail')$BirthGA, na.rm=T)


contingency_table <- table(demoDataPaperFC$mchat_result, demoDataPaperFC$Race)
chisq.test(contingency_table)


### Functional Data ----


#| echo: false
#| 
allTDFC <- readWideCSV(paste0(baseDir,'/functional/UseableData/tdCorrSummary.csv'))
allASDFC <- readWideCSV(paste0(baseDir,'/functional/UseableData/asdCorrSummary_Expand.csv'))

rownames(allTDFC) <- substr(rownames(allTDFC),0,nchar(rownames(allTDFC))-28)
#rownames(allASDFC) <- substr(rownames(allASDFC),0,nchar(rownames(allASDFC))-16)

tdFC <- allTDFC[rownames(allTDFC) %in% tdDemoFinal$ID,]
asdFC <- allASDFC[rownames(allASDFC) %in% asdDemoFinal$ID,]

rm(list=c('allTDFC', 'allASDFC'))


## 3b. Whole Brain Approach ----

# First step is to just do a whole brain search

# Here, let's run all 19,900 ROI-ROI connections

# It's more exploratory and FDR will kill us for sure


##### WHOLE BRAIN LINEAR MODEL APPROACH #####

# Fit both models for a chunk of connections (columns). Self-contained (no global objects) so it can run on a worker.
# tdBase / asdBase: each group's demographics merged once, with rowIdx = the row of that participant in the FC matrix.
fitWholeBrainChunk <- function(tdChunk, asdChunk, tdBase, asdBase) {

  lmeCtrl <- nlme::lmeControl(msMaxIter=200, tolerance=1e-4, msTol=1e-5, opt='optim')
  nConn <- ncol(tdChunk)
  mainMat <- matrix(NA_real_, nConn, 7)   # TD.Mean, ASD.Mean, Grp.B, Grp.p, GA.B, GA.p, nIncluded
  intxMat <- matrix(NA_real_, nConn, 11)  # TD.Mean, ASD.Mean, TD.SE, ASD.SE, Grp.B, Grp.p, GA.B, GA.p, Intx.B, Intx.p, nIncluded

  for(j in seq_len(nConn)){

    tempTD <- tdBase
    tempTD$FC <- tdChunk[tdBase$rowIdx, j]

    tempASD <- asdBase
    tempASD$FC <- asdChunk[asdBase$rowIdx, j]

    glmFrame <- rbind(tempTD, tempASD)
    glmFrame$shortID <- substr(glmFrame$ID, 0, nchar(glmFrame$ID)-3)

    # As an alternative, include the ones that have a missing here or there
    # Lose three significant connections, but gain two more that are BS to Cun and to FFG
    # glmFrame <- glmFrame[complete.cases(glmFrame),]

    # connections with any missing value are left as NA
    if(is.na(sum(glmFrame$FC))) { next }

    # whole-brain group means are raw means (the network tables use model-adjusted means)
    tdMean <- mean(as.numeric(tempTD$FC), na.rm = T)
    asdMean <- mean(as.numeric(tempASD$FC), na.rm = T)

    tempGLM <- nlme::lme(data=glmFrame, FC ~ as.factor(Grp) + ScanGA, random= ~1|factor(shortID), control=lmeCtrl)
    tTable <- summary(tempGLM)$tTable
    mainMat[j,] <- c(tdMean, asdMean, tTable[2,1], tTable[2,5], tTable[3,1], tTable[3,5], nrow(glmFrame))

    tempGLM <- nlme::lme(data=glmFrame, FC ~ as.factor(Grp) * ScanGA, random= ~1|factor(shortID), control=lmeCtrl)
    tTable <- summary(tempGLM)$tTable
    tdSE <- sd(as.numeric(tempTD$FC), na.rm = T)/sqrt(length(as.numeric(tempTD$FC)))
    asdSE <- sd(as.numeric(tempASD$FC), na.rm = T)/sqrt(length(as.numeric(tempASD$FC)))
    intxMat[j,] <- c(tdMean, asdMean, tdSE, asdSE, tTable[2,1], tTable[2,5], tTable[3,1], tTable[3,5],
                     tTable[4,1], tTable[4,5], nrow(glmFrame))
  }

  list(main = mainMat, intx = intxMat)
}

# Run fitWholeBrainChunk over every connection, split across nCores workers (nCores = 1 runs in this session).
# Returns the main-effects table (fcGLM) and the interaction table (fcIntxGLM), one row per connection.
runWholeBrainModels <- function(tdFC, asdFC, tdDemo, asdDemo, nCores) {

  # Merge each group's demographics once instead of once per connection
  tdBase <- merge(data.frame(ID=rownames(tdFC), rowIdx=seq_len(nrow(tdFC)), Grp='TD', stringsAsFactors=FALSE), tdDemo)
  asdBase <- merge(data.frame(ID=rownames(asdFC), rowIdx=seq_len(nrow(asdFC)), Grp='ASD', stringsAsFactors=FALSE), asdDemo)
  tdMat <- as.matrix(tdFC)
  asdMat <- as.matrix(asdFC)

  nConn <- ncol(tdMat)
  chunks <- parallel::splitIndices(nConn, max(1, min(nConn, nCores * 4)))
  # each worker gets only its own columns
  chunkData <- lapply(chunks, function(cols) list(td = tdMat[, cols, drop=FALSE], asd = asdMat[, cols, drop=FALSE]))

  startTime <- Sys.time()
  message('Fitting ', nConn, ' connections (', length(chunks), ' chunks, ', nCores, ' core(s))...')
  if (nCores > 1) {
    cl <- parallel::makeCluster(nCores)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterExport(cl, c('tdBase', 'asdBase'), envir = environment())
    parallel::clusterExport(cl, 'fitWholeBrainChunk', envir = globalenv())
    workerFun <- function(d) fitWholeBrainChunk(d$td, d$asd, tdBase, asdBase)
    environment(workerFun) <- globalenv()   # so only the exported objects (not this function's data) travel to the workers
    results <- parallel::parLapply(cl, chunkData, workerFun)
  } else {
    results <- lapply(seq_along(chunkData), function(k) {
      message('  chunk ', k, ' of ', length(chunkData))
      fitWholeBrainChunk(chunkData[[k]]$td, chunkData[[k]]$asd, tdBase, asdBase)
    })
  }
  message('Done in ', round(as.numeric(difftime(Sys.time(), startTime, units='mins')), 1), ' minutes')

  fcGLM <- data.frame(ROI = colnames(tdMat), do.call(rbind, lapply(results, `[[`, 'main')), stringsAsFactors = FALSE)
  colnames(fcGLM) <- c('ROI', 'TD.Mean', 'ASD.Mean','Grp.B','Grp.p','GA.B','GA.p', 'nIncluded')

  fcIntxGLM <- data.frame(ROI = colnames(tdMat), do.call(rbind, lapply(results, `[[`, 'intx')), stringsAsFactors = FALSE)
  colnames(fcIntxGLM) <- c('ROI', 'TD.Mean', 'ASD.Mean', 'TD.SE', 'ASD.SE','Grp.B','Grp.p','GA.B','GA.p', 'Intx.B','Intx.p', 'nIncluded')

  list(fcGLM = fcGLM, fcIntxGLM = fcIntxGLM)
}

# The fit is cached: it is recomputed only if the connectivity data, demographics or the model code change (or rerun is TRUE)
wholeBrainModels <- cached('wholebrain_models',
                           inputs = list(tdFC, asdFC, tdDemoFinal, asdDemoFinal, deparse(fitWholeBrainChunk)),
                           compute = function() runWholeBrainModels(tdFC, asdFC, tdDemoFinal, asdDemoFinal, nCores))

fcGLM <- addFDRfc(wholeBrainModels$fcGLM)
fcIntxGLM <- addFDRfc(wholeBrainModels$fcIntxGLM)

print(subset(fcGLM, Grp.p.adjust <.05))

print(subset(fcIntxGLM, Grp.p.adjust <.05))
print(subset(fcIntxGLM, Intx.p.adjust <.05))

saveResultTable(fcGLM, 'wholebrain_main')
saveResultTable(fcIntxGLM, 'wholebrain_intx')


write.table(fcIntxGLM, paste0(baseDir,'/functional/output/WholeBrain_WholeTable.csv'), col.names=T, row.names=F, sep=',')
write.table(subset(fcIntxGLM, Grp.p.adjust <.05), paste0(baseDir,'/functional/output/WholeBrain_P05FDR.csv'), col.names=T, row.names=F, sep=',')


### Whole Brain Figures ----


wholeBrainTable <- read.csv(paste0(baseDir,'/functional/output/WholeBrain_WholeTable.csv'))
wholeBrainTableFDR <- read.csv(paste0(baseDir,'/functional/output/WholeBrain_P05FDR.csv'))

## Connectogram here
# Create a dataframe for plotting
fullFCList <- read.csv(paste0(baseDir,'/functional/100PerHemiLabels_Decimal.1D'), header=F)

# Lobe colors: the single palette used for the chord diagram sectors and the heatmap lobe strips (Okabe-Ito, colorblind-safe).
# Change a color here and every figure follows.
lobeColors <- c(Frontal = '#E69F00', Parietal = '#F0E442', Temporal = '#009E73', Occipital = '#56B4E9',
                Subcortical = '#CC79A7', Cerebellum = '#999999', Brainstem = '#4d4d4d')

NetList <- read.csv(paste0(baseDir,'/functional/r200_Net.csv'))
NetList$LobeColor <- lobeColors[NetList$Lobe]

NetListNewOrder <- NetList[order(NetList$OrderNew),] 
ROIOrder <- NetListNewOrder$Name_Abrv
ROIColor <- NetListNewOrder$LobeColor


ROILobe <- c(ACG.L = 'Frontal',
ACG.L = 'Frontal',
ACG.L = 'Frontal',
IFGtriang.L = 'Frontal',
IFGtriang.L = 'Frontal',
IFGtriang.L = 'Frontal',
MFG.L = 'Frontal',
MFG.L = 'Frontal',
MFG.L = 'Frontal',
MFG.L = 'Frontal',
ORBinf.L = 'Frontal',
ORBinf.L = 'Frontal',
ORBmid.L = 'Frontal',
ORBmid.L = 'Frontal',
ORBsupb.L = 'Frontal',
ORBsupb.L = 'Frontal',
PreCG.L = 'Frontal',
PreCG.L = 'Frontal',
PreCG.L = 'Frontal',
PreCG.L = 'Frontal',
PreCG.L = 'Frontal',
REC.L = 'Frontal',
SFGdor.L = 'Frontal',
SFGdor.L = 'Frontal',
SFGdor.L = 'Frontal',
SFGdor.L = 'Frontal',
SFGdor.L = 'Frontal',
SFGdor.L = 'Frontal',
SFGmed.L = 'Frontal',
SFGmed.L = 'Frontal',
INS.L = 'Parietal',
INS.L = 'Parietal',
INS.L = 'Parietal',
INS.L = 'Parietal',
IPL.L = 'Parietal',
IPL.L = 'Parietal',
IPL.L = 'Parietal',
MCG.L = 'Parietal',
MCG.L = 'Parietal',
MCG.L = 'Parietal',
MCG.L = 'Parietal',
PCG.L = 'Parietal',
PCUN.L = 'Parietal',
PCUN.L = 'Parietal',
PCUN.L = 'Parietal',
PoCG.L = 'Parietal',
PoCG.L = 'Parietal',
PoCG.L = 'Parietal',
SMG.L = 'Parietal',
SMG.L = 'Parietal',
SPG.L = 'Parietal',
SPG.L = 'Parietal',
SPG.L = 'Parietal',
SPG.L = 'Parietal',
FFG.L = 'Temporal',
ITG.L = 'Temporal',
ITG.L = 'Temporal',
ITG.L = 'Temporal',
MTG.L = 'Temporal',
MTG.L = 'Temporal',
MTG.L = 'Temporal',
MTG.L = 'Temporal',
MTG.L = 'Temporal',
STG.L = 'Temporal',
STG.L = 'Temporal',
STG.L = 'Temporal',
STG.L = 'Temporal',
TPOmid.L = 'Temporal',
CAL.L = 'Occipital',
CAL.L = 'Occipital',
CUN.L = 'Occipital',
CUN.L = 'Occipital',
CUN.L = 'Occipital',
LING.L = 'Occipital',
LING.L = 'Occipital',
IOG.L = 'Occipital',
IOG.L = 'Occipital',
MOG.L = 'Occipital',
MOG.L = 'Occipital',
MOG.L = 'Occipital',
SOG.L = 'Occipital',
AMY.L = 'Subcortical',
CAU.L = 'Subcortical',
CAU.L = 'Subcortical',
CAU.L = 'Subcortical',
PAL.L = 'Subcortical',
PHG.L = 'Subcortical',
PHG.L = 'Subcortical',
PUT.L = 'Subcortical',
THA.L = 'Subcortical',
THA.L = 'Subcortical',
CRB.L = 'Cerebellum',
CRB.L = 'Cerebellum',
CRB.L = 'Cerebellum',
CRB.L = 'Cerebellum',
CRB.L = 'Cerebellum',
BS = 'Brainstem',
BS = 'Brainstem',
BS = 'Brainstem',
BS = 'Brainstem',
BS = 'Brainstem',
BS = 'Brainstem',
BS = 'Brainstem',
BS = 'Brainstem',
CRB.R = 'Cerebellum',
CRB.R = 'Cerebellum',
CRB.R = 'Cerebellum',
CRB.R = 'Cerebellum',
CRB.R = 'Cerebellum',
THA.R = 'Subcortical',
THA.R = 'Subcortical',
THA.R = 'Subcortical',
PUT.R = 'Subcortical',
PUT.R = 'Subcortical',
PUT.R = 'Subcortical',
PHG.R = 'Subcortical',
PAL.R = 'Subcortical',
HIP.R = 'Subcortical',
AMY.R = 'Subcortical',
SOG.R = 'Occipital',
SOG.R = 'Occipital',
MOG.R = 'Occipital',
MOG.R = 'Occipital',
MOG.R = 'Occipital',
IOG.R = 'Occipital',
IOG.R = 'Occipital',
LING.R = 'Occipital',
LING.R = 'Occipital',
CUN.R = 'Occipital',
CUN.R = 'Occipital',
CAL.R = 'Occipital',
CAL.R = 'Occipital',
TPOsup.R = 'Temporal',
TPOmid.R = 'Temporal',
STG.R = 'Temporal',
STG.R = 'Temporal',
STG.R = 'Temporal',
MTG.R = 'Temporal',
MTG.R = 'Temporal',
MTG.R = 'Temporal',
MTG.R = 'Temporal',
ITG.R = 'Temporal',
ITG.R = 'Temporal',
ITG.R = 'Temporal',
FFG.R = 'Temporal',
FFG.R = 'Temporal',
SPG.R = 'Parietal',
SPG.R = 'Parietal',
SPG.R = 'Parietal',
SPG.R = 'Parietal',
SMG.R = 'Parietal',
PoCG.R = 'Parietal',
PoCG.R = 'Parietal',
PoCG.R = 'Parietal',
PoCG.R = 'Parietal',
PoCG.R = 'Parietal',
PoCG.R = 'Parietal',
PCUN.R = 'Parietal',
PCUN.R = 'Parietal',
PCUN.R = 'Parietal',
PCL.R = 'Parietal',
PCL.R = 'Parietal',
MCG.R = 'Parietal',
MCG.R = 'Parietal',
MCG.R = 'Parietal',
IPL.R = 'Parietal',
IPL.R = 'Parietal',
IPL.R = 'Parietal',
IPL.R = 'Parietal',
INS.R = 'Parietal',
INS.R = 'Parietal',
HES.R = 'Parietal',
SMA.R = 'Frontal',
SMA.R = 'Frontal',
SFGmed.R = 'Frontal',
SFGmed.R = 'Frontal',
SFGdor.R = 'Frontal',
SFGdor.R = 'Frontal',
SFGdor.R = 'Frontal',
SFGdor.R = 'Frontal',
PreCG.R = 'Frontal',
PreCG.R = 'Frontal',
PreCG.R = 'Frontal',
ORBsupb.R = 'Frontal',
ORBsupb.R = 'Frontal',
ORBsupb.R = 'Frontal',
ORBsupb.R = 'Frontal',
ORBmid.R = 'Frontal',
ORBmed.R = 'Frontal',
ORBinf.R = 'Frontal',
ORBinf.R = 'Frontal',
MFG.R = 'Frontal',
MFG.R = 'Frontal',
MFG.R = 'Frontal',
MFG.R = 'Frontal',
IFGtriang.R = 'Frontal',
IFGtriang.R = 'Frontal',
IFGoperc.R = 'Frontal',
ACG.R = 'Frontal',
ACG.R = 'Frontal')

# the color of each chord sector comes from its lobe
ROIColor <- setNames(lobeColors[ROILobe], names(ROILobe))


# One row per (to, from) region pair, 'from' varying fastest. A connection's MeanDiff is its TD - ASD mean difference
# from the whole-brain table, or 0 if the pair is not in the table exactly once.
chordLabels <- fullFCList[,1]
nChordLabels <- length(chordLabels)
chordPairs <- expand.grid(from = seq_len(nChordLabels), to = seq_len(nChordLabels))

chordDiagram.Base <- data.frame(
  ConnFull = paste0(chordLabels[chordPairs$to],'.',chordPairs$to,'_',chordLabels[chordPairs$from],'.',chordPairs$from),
  ConnAbbrev = paste0(chordLabels[chordPairs$to],'_',chordLabels[chordPairs$from]),
  ToFull = paste0(chordLabels[chordPairs$to],'.',chordPairs$to),
  ToAbbrev = chordLabels[chordPairs$to],
  FromFull = paste0(chordLabels[chordPairs$from],'.',chordPairs$from),
  FromAbbrev = chordLabels[chordPairs$from],
  stringsAsFactors = FALSE)

tableRow <- match(chordDiagram.Base$ConnFull, wholeBrainTable$ROI)
inTableOnce <- !is.na(tableRow) & !(chordDiagram.Base$ConnFull %in% wholeBrainTable$ROI[duplicated(wholeBrainTable$ROI)])
chordDiagram.Base$MeanDiff <- ifelse(inTableOnce, wholeBrainTable$TD.Mean[tableRow] - wholeBrainTable$ASD.Mean[tableRow], 0)

chordDiagram.Base$ASDConn <- as.numeric(chordDiagram.Base$MeanDiff < 0)
chordDiagram.Base$TDConn <- as.numeric(chordDiagram.Base$MeanDiff > 0)
chordDiagram.Base$MeanDiff[is.na(chordDiagram.Base$MeanDiff)] <- 0
chordDiagram.Base$ASDConn[is.na(chordDiagram.Base$ASDConn)] <- 0
chordDiagram.Base$TDConn[is.na(chordDiagram.Base$TDConn)] <- 0

chordColorASD <- paste0(groupColors[['high']], chordLinkAlpha)
chordColorTD <- paste0(groupColors[['low']], chordLinkAlpha)

chordDiagram.Base$ASDCol <- '#00000000'
chordDiagram.Base$ASDCol[chordDiagram.Base$ASDConn==1] <- chordColorASD
chordDiagram.Base$ASDCol[!(chordDiagram.Base$ConnFull %in% wholeBrainTableFDR$ROI)] <- '#00000000'

chordDiagram.Base$TDCol <- '#00000000'
chordDiagram.Base$TDCol[chordDiagram.Base$TDConn==1] <- chordColorTD
chordDiagram.Base$TDCol[!(chordDiagram.Base$ConnFull %in% wholeBrainTableFDR$ROI)] <- '#00000000'


gapSizes <- c(rep(2:2, 9), #L Frontal
              5,
              rep(2:2, 7), #L Parietal
              5,
              rep(2:2, 4), #L Temporal
              5,
              rep(2:2, 5), #L Occipital
              5,
              rep(2:2, 5), #L Subcortical
              5,
              5,
              5,
              5,
              rep(2:2, 5), #R Subcortical
              5,
              rep(2:2, 5), #R Occipital
              5,
              rep(2:2, 5), #R Temporal
              5,
              rep(2:2, 8), #R Parietal
              5,
              rep(2:2, 11), #R Frontal)
              13) 


# circlize lays out all 40,000 region pairs, which takes minutes, but only the significant links are visible; the rest
# are transparent and only take up space. Within a pair of regions, a run of consecutive transparent links takes up
# exactly the same space as one transparent link of their combined value, so each run is merged into one link.
# This gives the identical picture from about 6,000 links instead of 40,000.
compressChordLinks <- function(links, colorCol, valueCol, invisibleColor = '#00000000') {
  n <- nrow(links)
  pair <- paste(links$ToAbbrev, links$FromAbbrev)
  ord <- order(factor(pair), seq_len(n))              # group rows by region pair, keeping original order within a pair
  samePair <- pair[ord]
  hidden <- links[[colorCol]][ord] == invisibleColor
  newRun <- c(TRUE, samePair[-1] != samePair[-n] | !hidden[-1] | !hidden[-n])   # visible links stay on their own
  runId <- cumsum(newRun)
  firstRow <- ord[newRun]                              # first original row of every run
  out <- links[firstRow, , drop = FALSE]
  out[[valueCol]] <- as.vector(rowsum(links[[valueCol]][ord], runId))
  out[order(firstRow), , drop = FALSE]                 # back to the original row order
}

chordASDLinks <- compressChordLinks(chordDiagram.Base, 'ASDCol', 'ASDConn')
chordTDLinks <- compressChordLinks(chordDiagram.Base, 'TDCol', 'TDConn')

# Draw a chord diagram from compressed links. connCol = column with the link values, colorCol = column with the link colors.
# Transparent links are not drawn (link.visible), which leaves the layout unchanged but is much faster.
# The left hemisphere is on the left of the picture; sectors are labeled with full region names outside the ring, and
# 'Left' / 'Right' are written underneath. labelSize scales the region names; canvas = half-width and half-height of the drawing area.
drawChord <- function(links, connCol, colorCol, labelSize = .8, canvas = c(1.5, 1.5)) {
  chordArgs <- list(links[,c('ToAbbrev','FromAbbrev',connCol)],
                    order = ROIOrder, grid.col=ROIColor, grid.border = '#000000',
                    col=links[[colorCol]], symmetric = TRUE, scale=TRUE,
                    link.border=links[[colorCol]], link.lwd=chordLinkBorderWidth,
                    annotationTrack = c("grid"), annotationTrackHeight = c(0.03, 0.03))
  if ('link.visible' %in% names(formals(circlize::chordDiagram))) {
    chordArgs$link.visible <- links[[colorCol]] != '#00000000'
  }

  circos.clear()
  # counter-clockwise from 99 degrees is the mirror image of the old layout (clockwise from 81), putting the left hemisphere on the left
  circos.par(start.degree=99, clock.wise=FALSE, gap.after=gapSizes, canvas.xlim=c(-canvas[1], canvas[1]), canvas.ylim=c(-canvas[2], canvas[2]),
             points.overflow.warning=FALSE)   # the names sit just outside the ring on purpose
  do.call(chordDiagram, chordArgs)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    # names start a fixed distance outside the ring (the same gap for every region, on both sides)
    circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(1.5), roiFullName(CELL_META$sector.index),
                facing = "clockwise", niceFacing = TRUE, adj = c(0, .5), cex = labelSize)
  }, bg.border = NA) # here set bg.border to NA is important
  text(c(-.8, .8), -(canvas[2] - .08), c('Left', 'Right'), cex = labelSize * 2, xpd = NA)
  circos.clear()
}

showBaseFigure('wholebrain_chord_ASD', width = 9, height = 9, drawFun = function() drawChord(chordASDLinks, 'ASDConn', 'ASDCol'), subfolder = componentFolder)

showBaseFigure('wholebrain_chord_TD', width = 9, height = 9, drawFun = function() drawChord(chordTDLinks, 'TDConn', 'TDCol'), subfolder = componentFolder)


asdConnections <- c(chordDiagram.Base[chordDiagram.Base$ASDCol == chordColorASD,]$ToAbbrev, chordDiagram.Base[chordDiagram.Base$ASDCol == chordColorASD,]$FromAbbrev)

tdConnections <- c(chordDiagram.Base[chordDiagram.Base$TDCol == chordColorTD,]$ToAbbrev, chordDiagram.Base[chordDiagram.Base$TDCol == chordColorTD,]$FromAbbrev)


sort(table(asdConnections))
sort(table(tdConnections))


## DESCRIPTIVE FOR REGIONAL CONNECTIVITY

lobeNames <- c('Frontal','Parietal','Temporal','Occipital','Subcortical','Cerebellum','Brainstem')

# Row of lobeNames for an ROI such as 'CAU.L.39': its first two '.'-separated fields ('CAU.L') are looked up in the
# network table. The brainstem has no L/R, so any name starting 'BS' is just 'BS'. NA if the ROI is not in the table.
roiLobeIndex <- function(roiPart) {
  fields <- strsplit(roiPart, '.', fixed=T)
  roiName <- paste0(vapply(fields, function(x) x[1], character(1)), '.', vapply(fields, function(x) x[2], character(1)))
  roiName[substr(roiName,0,2)=='BS'] <- 'BS'
  match(NetList$Lobe[match(roiName, NetList$Name_Abrv)], lobeNames)
}

# 7 x 7 lobe-by-lobe count of connections. A connection adds 1 to both (lobe1, lobe2) and (lobe2, lobe1), so a
# connection within one lobe adds 2 to that lobe's diagonal cell.
countLobePairs <- function(index1, index2) {
  keep <- !is.na(index1) & !is.na(index2)
  pairCounts <- table(factor(index1[keep], levels=seq_along(lobeNames)), factor(index2[keep], levels=seq_along(lobeNames)))
  counts <- matrix(as.numeric(pairCounts + t(pairCounts)), nrow=7, ncol=7, dimnames=list(lobeNames, lobeNames))
  as.data.frame(counts)
}

connParts <- strsplit(wholeBrainTable$ROI, '_', useBytes=T)
lobeIndex1 <- roiLobeIndex(vapply(connParts, function(x) x[1], character(1)))
lobeIndex2 <- roiLobeIndex(vapply(connParts, function(x) x[2], character(1)))

# Significant connections (FDR) and whether TD or ASD has the higher mean
sigRow <- match(wholeBrainTable$ROI, wholeBrainTableFDR$ROI)
isSig <- !is.na(sigRow)
tdHigher <- isSig & ((wholeBrainTableFDR$TD.Mean[sigRow] > wholeBrainTableFDR$ASD.Mean[sigRow]) %in% TRUE)

allConnDesc <- countLobePairs(lobeIndex1, lobeIndex2)
allConnDescSig <- countLobePairs(lobeIndex1[isSig], lobeIndex2[isSig])
allConnDescSig.TD <- countLobePairs(lobeIndex1[tdHigher], lobeIndex2[tdHigher])
allConnDescSig.ASD <- countLobePairs(lobeIndex1[isSig & !tdHigher], lobeIndex2[isSig & !tdHigher])

# 7 x 7 lobe-by-lobe connection counts (the lobe names are the row names, so they become a column)
lobeCountNames <- c(allConnDesc = 'all', allConnDescSig = 'sig', allConnDescSig.TD = 'sigTD', allConnDescSig.ASD = 'sigASD')
for (lobeTable in names(lobeCountNames)) {
  saveTable(cbind(Lobe = rownames(get(lobeTable)), get(lobeTable)), paste0('wholebrain_lobeCounts_', lobeCountNames[[lobeTable]]))
}


# Lower-triangle heatmap of the lobe-by-lobe connection counts: counts in the tiles, a colored strip for each lobe beside
# the axis labels, and a small colorbar in the empty upper-right corner.
#   counts = 7 x 7 table (rows and columns in lobeNames order);  low/high = ends of the color scale
#   limits = fill scale limits (NULL = 0 to the largest count)
lobeHeatmapPlot <- function(counts, high, limits = NULL, low = '#ffffff') {
  nLobes <- length(lobeNames)
  tiles <- data.frame(row = rep(seq_len(nLobes), each = nLobes), col = rep(seq_len(nLobes), times = nLobes),
                      num = as.vector(t(as.matrix(counts))))
  tiles <- tiles[tiles$col <= tiles$row,]                 # lower triangle, including the diagonal
  tiles$y <- nLobes - tiles$row + 1                       # row 1 at the top
  if (is.null(limits)) { limits <- c(0, max(tiles$num)) }

  # colorbar: gradient (same interpolation as the tiles), tick marks, and a title
  gradientColors <- scales::seq_gradient_pal(low, high)(seq(1, 0, length.out = 100))
  barX <- c(5.55, 5.85); barY <- c(5.15, 6.75)
  ticks <- pretty(limits, n = 5); ticks <- ticks[ticks >= limits[1] & ticks <= limits[2]]
  tickY <- barY[1] + (ticks - limits[1]) / diff(limits) * diff(barY)

  ggplot(tiles, aes(x = col, y = y)) +
    geom_tile(aes(fill = num), color = '#bfbfbf', linewidth = .4) +
    geom_text(aes(label = num), color = 'black', size = 5.5, fontface = 'bold') +
    scale_fill_gradient(low = low, high = high, limits = limits, guide = 'none') +
    # row labels and strips (left)
    annotate('rect', xmin = .22, xmax = .42, ymin = nLobes - seq_len(nLobes) + .5, ymax = nLobes - seq_len(nLobes) + 1.5,
             fill = lobeColors[lobeNames], colour = '#808080', linewidth = .3) +
    annotate('text', x = .12, y = nLobes - seq_len(nLobes) + 1, label = lobeNames, hjust = 1, size = 4.6) +
    # column labels and strips (bottom)
    annotate('rect', xmin = seq_len(nLobes) - .5, xmax = seq_len(nLobes) + .5, ymin = .08, ymax = .28,
             fill = lobeColors[lobeNames], colour = '#808080', linewidth = .3) +
    annotate('text', x = seq_len(nLobes), y = .02, label = lobeNames, angle = 90, hjust = 1, size = 4.6) +
    # colorbar
    annotation_raster(as.raster(matrix(gradientColors, ncol = 1)), xmin = barX[1], xmax = barX[2], ymin = barY[1], ymax = barY[2]) +
    annotate('rect', xmin = barX[1], xmax = barX[2], ymin = barY[1], ymax = barY[2], fill = NA, colour = 'black', linewidth = .3) +
    annotate('segment', x = barX[2], xend = barX[2] + .08, y = tickY, yend = tickY, linewidth = .3) +
    annotate('text', x = barX[2] + .12, y = tickY, label = ticks, hjust = 0, size = 3.4) +
    annotate('text', x = barX[2] + .6, y = mean(barY) + .15, label = 'Significant\nConnections', hjust = 0, size = 3.8, lineheight = .95) +
    coord_fixed(xlim = c(.5, nLobes + .5), ylim = c(.5, nLobes + .5), expand = FALSE, clip = 'off') +
    theme_void() +
    theme(plot.margin = margin(t = 8, r = 40, b = 100, l = 105, unit = 'pt'))
}

# Both heatmaps use the same color scale, 0 to 25, so their colors and colorbars can be compared directly
lobeHeatmapLimits <- c(0, 25)

# Low Risk (TD) is blue and High Risk (ASD) is red, as in the other figures
wholebrain_lobeHeatmap_TD <- lobeHeatmapPlot(allConnDescSig.TD, high = groupColors[['low']], limits = lobeHeatmapLimits)
showFigure(wholebrain_lobeHeatmap_TD, 'wholebrain_lobeHeatmap_TD', width = 7, height = 7, subfolder = componentFolder)

wholebrain_lobeHeatmap_ASD <- lobeHeatmapPlot(allConnDescSig.ASD, high = groupColors[['high']], limits = lobeHeatmapLimits)
showFigure(wholebrain_lobeHeatmap_ASD, 'wholebrain_lobeHeatmap_ASD', width = 7, height = 7, subfolder = componentFolder)


### Whole Brain Glass Brains ----

# Significant whole-brain connections (FDR) drawn on a translucent brain: lateral view above, posterior view below.
# Low Risk (blue) = connections where Low Risk is higher; High Risk (red) = connections where High Risk is higher.

loadGlassBrain()

wholeBrainEdges <- connectionRegions(wholeBrainTableFDR$ROI)
wholeBrainEdges$dir <- ifelse(wholeBrainTableFDR$TD.Mean > wholeBrainTableFDR$ASD.Mean, 'low', 'high')
allRegionNodes <- data.frame(id = seq_len(nrow(roiCoordinates)), fill = '#505050')
message('Whole-brain glass brains: ', sum(wholeBrainEdges$dir == 'low', na.rm = TRUE), ' Low > High Risk edges and ',
        sum(wholeBrainEdges$dir == 'high', na.rm = TRUE), ' High > Low Risk edges')

wholeBrainGlass <- function(direction) {
  edges <- wholeBrainEdges[which(wholeBrainEdges$dir == direction),]
  ggpubr::ggarrange(glassBrainView('lateral', glassBrainBackgrounds, roiCoordinates, allRegionNodes, edges, nodeSize = 1.6, edgeWidth = .9),
                    glassBrainView('posterior', glassBrainBackgrounds, roiCoordinates, allRegionNodes, edges, nodeSize = 1.6, edgeWidth = .9),
                    ncol = 1, nrow = 2)
}

wholebrain_glass_TD <- wholeBrainGlass('low')
showFigure(wholebrain_glass_TD, 'wholebrain_glass_TD', width = 4.5, height = 9, subfolder = componentFolder)

wholebrain_glass_ASD <- wholeBrainGlass('high')
showFigure(wholebrain_glass_ASD, 'wholebrain_glass_ASD', width = 4.5, height = 9, subfolder = componentFolder)


### Whole Brain Composite Figure ----

# Low Risk (left half: A, B, C) and High Risk (right half: D, E, F). In each half: glass brains (A / D),
# the lobe heatmap (B / E) and the chord diagram (C / F).

# Dashed placeholder box for a glass-brain image. x, y = lower-left corner and width, height as fractions of the figure.
glassBrainPlaceholder <- function(x, y, width, height, label) {
  grid::grid.roundrect(x = grid::unit(x, 'npc'), y = grid::unit(y, 'npc'), width = grid::unit(width, 'npc'), height = grid::unit(height, 'npc'),
                       just = c('left', 'bottom'), r = grid::unit(.02, 'snpc'),
                       gp = grid::gpar(col = '#999999', fill = '#f4f4f4', lty = 'dashed', lwd = 2))
  grid::grid.text(label, x = grid::unit(x + width/2, 'npc'), y = grid::unit(y + height/2, 'npc'), gp = grid::gpar(col = '#888888', fontsize = 14))
}

drawWholeBrainComposite <- function() {
  # start a base-graphics page first, so that base graphics (the chord diagrams) draw onto it instead of starting a new page
  # and erasing what grid (the heatmaps and boxes) has already drawn
  plot.new()

  halves <- list(list(title = 'Low Risk', x0 = 0, heatmap = wholebrain_lobeHeatmap_TD, links = chordTDLinks,
                      conn = 'TDConn', color = 'TDCol', letters = c('A.', 'B.', 'C.'), glass = wholebrain_glass_TD),
                 list(title = 'High Risk', x0 = .5, heatmap = wholebrain_lobeHeatmap_ASD, links = chordASDLinks,
                      conn = 'ASDConn', color = 'ASDCol', letters = c('D.', 'E.', 'F.'), glass = wholebrain_glass_ASD))

  for (half in halves) {
    # A / D: glass brains (lateral view above, posterior view below)
    print(half$glass, newpage = FALSE,
          vp = grid::viewport(x = grid::unit(half$x0 + .005, 'npc'), y = grid::unit(.49, 'npc'), width = grid::unit(.178, 'npc'),
                              height = grid::unit(.45, 'npc'), just = c('left', 'bottom')))

    # B / E: lobe heatmap
    print(half$heatmap, newpage = FALSE,
          vp = grid::viewport(x = grid::unit(half$x0 + .185, 'npc'), y = grid::unit(.50, 'npc'), width = grid::unit(.315, 'npc'),
                              height = grid::unit(.455, 'npc'), just = c('left', 'bottom')))

    # C / F: chord diagram in the lower half
    par(fig = c(half$x0, half$x0 + .5, 0, .5), new = TRUE, mar = c(0, 0, 0, 0))
    drawChord(half$links, half$conn, half$color, labelSize = .6, canvas = c(1.75, 1.4))

    # heading and panel letters
    grid::grid.text(half$title, x = grid::unit(half$x0 + .25, 'npc'), y = grid::unit(.972, 'npc'), gp = grid::gpar(fontsize = 44))
    letterGp <- grid::gpar(fontsize = 30, fontface = 'bold')
    grid::grid.text(half$letters[1], x = grid::unit(half$x0 + .008, 'npc'), y = grid::unit(.945, 'npc'), just = c('left', 'top'), gp = letterGp)
    grid::grid.text(half$letters[2], x = grid::unit(half$x0 + .185, 'npc'), y = grid::unit(.945, 'npc'), just = c('left', 'top'), gp = letterGp)
    grid::grid.text(half$letters[3], x = grid::unit(half$x0 + .008, 'npc'), y = grid::unit(.485, 'npc'), just = c('left', 'top'), gp = letterGp)
  }
}

showBaseFigure('wholebrain_composite', width = 16, height = 13, drawFun = drawWholeBrainComposite)


## 3c. Smaller Network Approach ----

# Here we are subdividing out into smaller networks

# This will (A) help with the crazy FDR correction and (B) more pointedly target questions about connectivity within networks such as the CRB

# Biggest job, however, is to define the networks in a rigorous way


## First -- divide out the connections into an ROI1 and ROI2 table

connROIs <- as.data.frame(matrix(ncol=2, nrow=ncol(tdFC)))
colnames(connROIs) <- c('ROI1','ROI2')

for(i in 1:nrow(connROIs)){
  tempROI <- colnames(tdFC)[i]
  
  connROIs$ROI1[i] <- strsplit(tempROI,'_')[[1]][1]
  connROIs$ROI2[i] <- strsplit(tempROI,'_')[[1]][2] 
  
}


### Network Analyses ----

# For each network, take the connections where BOTH regions belong to the network, then for every connection fit
#   (1) FC ~ Grp + ScanGA   and   (2) FC ~ Grp * ScanGA   (random intercept per subject).
# Each network writes its tables to functional/output and gets a bar chart (network_<name>_bars).


# ROI number = the third '.'-separated part of an ROI name (e.g. 'CAU.L.39' -> 39)
roiNumber <- function(roi) as.numeric(vapply(strsplit(roi, '.', fixed = TRUE), function(x) x[3], character(1)))

# One entry per network:
#   match         - function(roi names) -> TRUE for regions in the network
#   xLabel        - x-axis label of the bar chart
#   zeroLine      - draw a line at 0 on the bar chart
networkSpecs <- list(

  CRB = list(match = function(roi) grepl('CRB', roi),
             xLabel = 'Cerebellar Connections', zeroLine = TRUE),

  # Originally c(5,9,20,29,56,61,67,68,81,86,94,108,120,133,141,143,155,160,188,196,198,200); parahippocampal gyrus added (43,86,125)
  Limbic = list(match = function(roi) roiNumber(roi) %in% c(5,9,20,29,56,61,67,68,81,86,94,108,120,133,141,143,155,160,188,196,198,200,
                                                            43,86,125),
                xLabel = 'Limbic Connections', zeroLine = FALSE),

  # Originally also included 43,86,125 (hippocampus and parahippocampal gyrus); removed
  DMN = list(match = function(roi) roiNumber(roi) %in% c(4,11,18,22,23,34,40,42,48,53,70,85,88,96,100,102,104,106,123,124,128,137,150,152,156,163,168,171,176,178),
             xLabel = 'Default Mode Connections', zeroLine = FALSE),

  SomMot = list(match = function(roi) grepl('PreCG', roi) | grepl('PoCG', roi) | grepl('SMA', roi),
                xLabel = 'Sensorimotor Connections', zeroLine = FALSE),

  SubCort = list(match = function(roi) roiNumber(roi) %in% c(31,32,39,51,80,92,97,107,111,134,139,145,151),
                 xLabel = 'Subcortical Connections', zeroLine = FALSE),

  Auditory = list(match = function(roi) roiNumber(roi) %in% c(88,94,2,12,48,171,107,139,152,185),
                  xLabel = 'Auditory Connections', zeroLine = FALSE),

  Visual = list(match = function(roi) roiNumber(roi) %in% c(181,174,153,114,149,195,167,105,116,194,197,162,131,186,135,
                                                            72,60,19,52,35,75,25,15,42,30,91,87,98,57),
                xLabel = 'Visual Connections', zeroLine = FALSE)
)

# Fit both models for every connection (column) of the network.
# Returns the main-effects table (fcGLM) and the interaction table (fcIntxGLM).
fitNetworkModels <- function(tdFCNetwork, asdFCNetwork, tdDemo, asdDemo) {

  lmeCtrl <- lmeControl(msMaxIter=200, tolerance=1e-4, msTol=1e-5, opt='optim')

  fcGLM <- as.data.frame(matrix(ncol=10, nrow=ncol(tdFCNetwork)))
  colnames(fcGLM) <- c('ROI', 'TD.Mean','TD.SE', 'ASD.Mean','ASD.SE','Grp.B','Grp.p','GA.B','GA.p','pEta2')

  fcIntxGLM <- as.data.frame(matrix(ncol=12, nrow=ncol(tdFCNetwork)))
  colnames(fcIntxGLM) <- c('ROI', 'TD.Mean', 'TD.SE', 'ASD.Mean','ASD.SE','Grp.B','Grp.p','GA.B','GA.p','Intx.B','Intx.p','pEta2')

  for(i in 1:nrow(fcGLM)){

    tempTD <- as.data.frame(cbind(rownames(tdFCNetwork),tdFCNetwork[,i]))
    colnames(tempTD) <- c('ID','FC')
    tempTD$Grp <- 'TD'
    tempTD <- merge(tempTD, tdDemo)

    tempASD <- as.data.frame(cbind(rownames(asdFCNetwork),asdFCNetwork[,i]))
    colnames(tempASD) <- c('ID','FC')
    tempASD$Grp <- 'ASD'
    tempASD <- merge(tempASD, asdDemo)

    glmFrame <- rbind(tempTD, tempASD)
    glmFrame$FC <- as.numeric(glmFrame$FC)
    glmFrame$shortID <- substr(glmFrame$ID, 0, nchar(glmFrame$ID)-3)

    tempGLM <- lme(data=glmFrame, FC ~ Grp + ScanGA, random= ~1|factor(shortID), control=lmeCtrl)

    fcGLM$ROI[i] <- colnames(tdFCNetwork)[i]
    fcGLM$TD.Mean[i] <- mean(as.numeric(tempTD$FC), na.rm = T)
    fcGLM$TD.SE[i] <- sd(as.numeric(tempTD$FC), na.rm = T)/sqrt(length(as.numeric(tempTD$FC)))
    fcGLM$ASD.Mean[i] <- mean(as.numeric(tempASD$FC), na.rm = T)
    fcGLM$ASD.SE[i] <- sd(as.numeric(tempASD$FC), na.rm = T)/sqrt(length(as.numeric(tempASD$FC)))
    fcGLM$Grp.B[i] <- summary(tempGLM)$tTable[2,1]
    fcGLM$Grp.p[i] <- summary(tempGLM)$tTable[2,5]
    fcGLM$GA.B[i] <- summary(tempGLM)$tTable[3,1]
    fcGLM$GA.p[i] <- summary(tempGLM)$tTable[3,5]

    tempGLM <- lme(data=glmFrame, FC ~ Grp * ScanGA, random= ~1|factor(shortID), control=lmeCtrl)

    fcIntxGLM$ROI[i] <- colnames(tdFCNetwork)[i]
    # group means in the interaction table are model-adjusted (emmeans)
    adjustMeans <- as.data.frame(emmeans(tempGLM, 'Grp', data=glmFrame))
    fcIntxGLM$TD.Mean[i] <- adjustMeans$emmean[2]
    fcIntxGLM$TD.SE[i] <- adjustMeans$SE[2]
    fcIntxGLM$ASD.Mean[i] <- adjustMeans$emmean[1]
    fcIntxGLM$ASD.SE[i] <- adjustMeans$SE[1]
    fcIntxGLM$Grp.B[i] <- summary(tempGLM)$tTable[2,1]
    fcIntxGLM$Grp.p[i] <- summary(tempGLM)$tTable[2,5]
    fcIntxGLM$GA.B[i] <- summary(tempGLM)$tTable[3,1]
    fcIntxGLM$GA.p[i] <- summary(tempGLM)$tTable[3,5]
    fcIntxGLM$Intx.B[i] <- summary(tempGLM)$tTable[4,1]
    fcIntxGLM$Intx.p[i] <- summary(tempGLM)$tTable[4,5]
  }

  list(fcGLM = fcGLM, fcIntxGLM = fcIntxGLM)
}

# Write a network's tables to functional/output (these files are read back for the connectivity figures below)
writeNetworkTables <- function(name, fcGLM, fcIntxGLM) {
  out <- function(tbl, suffix) {
    write.table(tbl, paste0(baseDir,'/functional/output/', name, '_', suffix, '.csv'), col.names=T, row.names=F, sep=',')
  }
  out(fcGLM, 'WholeTable')
  out(subset(fcGLM, Grp.p <.05), 'P05')
  out(subset(fcGLM, Grp.p <.01), 'P01')
  out(subset(fcGLM, Grp.p <.001), 'P001')
  out(subset(fcGLM, Grp.p.adjust <.05), 'P05FDR')

  out(subset(fcIntxGLM, Grp.p.adjust <.05 | Intx.p.adjust < .05), 'Intx_P05FDR')
  out(subset(fcIntxGLM, Grp.p.adjust <.1 | Intx.p.adjust < .1), 'Intx_P10FDR')
  out(fcIntxGLM, 'Intx_Wholebrain')
}

# Run one network: select its connections, fit the models (cached), FDR-correct, print what survives, and write the tables
runNetworkAnalysis <- function(name, spec) {
  bothROIMatch <- spec$match(connROIs$ROI1) & spec$match(connROIs$ROI2)

  # The fit is cached per network: it is recomputed only if this network's connections, the demographics or the model code change
  tdNetwork <- tdFC[,bothROIMatch]
  asdNetwork <- asdFC[,bothROIMatch]
  models <- cached(paste0('network_', name),
                   inputs = list(tdNetwork, asdNetwork, tdDemoFinal, asdDemoFinal, deparse(fitNetworkModels)),
                   compute = function() fitNetworkModels(tdNetwork, asdNetwork, tdDemoFinal, asdDemoFinal))
  fcGLM <- addFDRfc(models$fcGLM)
  fcIntxGLM <- addFDRfc(models$fcIntxGLM)

  message(name, ': ', sum(bothROIMatch), ' connections')
  print(subset(fcGLM, Grp.p.adjust <.05))
  print(subset(fcIntxGLM, Grp.p.adjust <.05))
  print(subset(fcIntxGLM, Intx.p.adjust <.05))

  writeNetworkTables(name, fcGLM, fcIntxGLM)
  saveResultTable(fcGLM, paste0('network_', name, '_main'))
  saveResultTable(fcIntxGLM, paste0('network_', name, '_intx'))
  list(fcGLM = fcGLM, fcIntxGLM = fcIntxGLM)
}

# Bar chart of the group means for every connection in a network, ordered by group p-value, with a boxed legend to its right
plotNetworkBars <- function(fcGLM, spec) {
  longFCGLM <- melt(fcGLM, id.vars=c('ROI'), measure.vars = c('TD.Mean','ASD.Mean'), variable.name='Group', value.name='FC')
  pValOrder <- fcGLM$ROI[order(fcGLM$Grp.p, decreasing = T)]

  # y axis fitted to this network's data: from the lowest bar (or 0) to the highest, on round tick marks.
  # A network with few negative connections gets only a small negative range, not a symmetrical axis.
  yBreaks <- pretty(c(0, longFCGLM$FC), n = 6)

  barPlot <- ggplot(data=longFCGLM, aes(x=factor(ROI, levels=pValOrder, ordered=is.ordered(T)), y=FC, group=factor(Group))) +
    scale_y_continuous(limits=range(yBreaks), breaks=yBreaks, labels=fixedDecimalLabels, expand=c(0,0)) +
    scale_fill_manual(values=c(networkBarFill[['low']], networkBarFill[['high']])) +
    scale_color_manual(values=c(groupColors[['low']], groupColors[['high']]))
  if (spec$zeroLine) { barPlot <- barPlot + geom_hline(yintercept = 0, linewidth=1.5) }

  barPlot <- barPlot +
    geom_bar(stat='identity', position='dodge', aes(color=Group,fill=Group)) +
    ylab('Connectivity Strength') + xlab(spec$xLabel) +
    theme(axis.title.x = element_text(color = 'black', size = 16),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(color = 'black', size = 16),
          axis.text.y = element_text(color = 'black', size = 12),
          axis.ticks.y = element_blank(),
          axis.line.x = element_line(color='black', linewidth = 1),
          axis.line.y = element_line(color='black', linewidth = 1),
          plot.title = element_text(color = 'black', size = 24, hjust=0.5),
          panel.grid.major.y = element_line(color = '#cccccc', linewidth = 1),
          panel.grid.major.x = element_blank(),   # no vertical lines: with hundreds of bars they fill the panel with gray
          panel.background = element_rect(fill='white'),
          legend.position = 'none',
          legend.text = element_text(color = 'black', size = 12),
          legend.title=element_blank())

  addRiskLegend(barPlot, 'network')
}


### Run Network Analyses ----


networkResults <- list()
for (net in names(networkSpecs)) {
  networkResults[[net]] <- runNetworkAnalysis(net, networkSpecs[[net]])
  showFigure(plotNetworkBars(networkResults[[net]]$fcGLM, networkSpecs[[net]]),
             paste0('network_', net, '_bars'), width = 12, height = 6)
}

#### Network Connectivity Figures ----


# Read in the output of the above analyses for plotting
FCTable.TrendOnly.DMN <- read.csv(paste0(baseDir,'/functional/output/DMN_Intx_P10FDR.csv'), stringsAsFactors=F)
FCTable.TrendOnly.Limbic <- read.csv(paste0(baseDir,'/functional/output/Limbic_Intx_P10FDR.csv'), stringsAsFactors=F)
FCTable.TrendOnly.Subcort <- read.csv(paste0(baseDir,'/functional/output/SubCort_Intx_P10FDR.csv'), stringsAsFactors=F)
FCTable.TrendOnly.CRB <- read.csv(paste0(baseDir,'/functional/output/CRB_Intx_P10FDR.csv'), stringsAsFactors=F)

FCTable.TrendOnly.DMN$Net <- 'DMN'
FCTable.TrendOnly.Limbic$Net <- 'Limbic'
FCTable.TrendOnly.Subcort$Net <- 'Subcort'
FCTable.TrendOnly.CRB$Net <- 'CRB'

FCTable.TrendOnly.All <- rbind(rbind(rbind(FCTable.TrendOnly.DMN, FCTable.TrendOnly.Limbic),
                                    FCTable.TrendOnly.Subcort),
                              FCTable.TrendOnly.CRB)


# Plots

plotOrder <- c('CAU.L.39_THA.R.139',
  'PAL.L.97_PUT.R.111',
  'CAU.L.39_THA.L.80',
  'THA.L.51_PUT.L.92',
  'PAL.L.97_PAL.R.107',
  'PAL.L.97_THA.R.139',
  'PAL.R.107_PUT.R.145',
  'TPOmid.L.9_ITG.R.188',
  'REC.L.20_TPOmid.R.108',
  'REC.L.20_ORBsupb.R.155',
  'ITG.L.56_ITG.R.188',
  'SFGdor.L.67_AMY.R.133',
  'SFGdor.L.67_ORBsupb.R.196',
  'ORBsupb.R.143_ITG.R.188',
  'ITG.R.188_ORBsupb.R.196',
  'ITG.R.188_ORBsupb.R.198',
  'SFGmed.L.18_ACG.R.123',
  'MTG.L.22_ACG.L.34',
  'MTG.L.22_ACG.R.156',
  'PCUN.L.23_STG.L.88',
  'MTG.L.96_STG.R.102',
  'SFGmed.R.128_MTG.R.150',
  'PHG.L.43_PCUN.R.178',
  'CRB.L.13_CRB.R.142',
  'CRB.L.13_CRB.R.192',
  'CRB.L.17_CRB.L.99',
  'CRB.L.17_CRB.R.192',
  'CRB.L.54_CRB.R.192',
  'CRB.L.99_CRB.R.142',
  'CRB.L.99_CRB.R.192',
  'CRB.R.129_CRB.R.192')


plotFrame.TD <- FCTable.TrendOnly.All[,c('ROI','TD.Mean', 'TD.SE')]
colnames(plotFrame.TD) <- c('Connection','Mean','SE')
plotFrame.ASD <- FCTable.TrendOnly.All[,c('ROI','ASD.Mean', 'ASD.SE')]
colnames(plotFrame.ASD) <- c('Connection','Mean','SE')
plotFrame.TD$Grp <- 'Low Risk'
plotFrame.ASD$Grp <- 'High Risk'

plotFrame <- as.data.frame(rbind(plotFrame.TD, plotFrame.ASD))
plotFrame$Connection <- factor(plotFrame$Connection, levels=plotOrder)
plotFrame$Grp <- factor(plotFrame$Grp, levels=c('Low Risk', 'High Risk'))
plotFrame$Mean <- as.numeric(plotFrame$Mean)
plotFrame$SE <- as.numeric(plotFrame$SE)

##plot ggeffects maybe ggprdict?
##need marginal means? -- emmeans? Give the SE from the model itself

network_combined_bars_p10 <- ggplot(data=plotFrame, aes(x=Connection, y=Mean, Group=Grp)) +
  scale_y_continuous(labels=fixedDecimalLabels) +
  scale_x_discrete(labels=setNames(connectionLabel(plotOrder), plotOrder)) +
  scale_fill_manual(values=c(networkBarFill[['low']], networkBarFill[['high']])) +
  scale_color_manual(values=c(groupColors[['low']], groupColors[['high']])) +
  geom_hline(yintercept = 0, lwd=1) +
  geom_bar(stat='identity', position=position_dodge2(0.2), aes(color=Grp,fill=Grp)) +
  # geom_errorbar(position=position_dodge2(0.2), aes(ymin=Mean-SE, ymax=Mean+SE), size=.2) +
  ylab('Connectivity Strength') + xlab('Network Connections') +
  theme(axis.title.x = element_text(color = 'black', size = 16),
        axis.text.x = element_text(color = 'black', size = 6, angle=90, hjust=1, vjust=.5),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(color = 'black', size = 16),
        axis.text.y = element_text(color = 'black', size = 12),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(color='black', linewidth = 1),
        axis.line.y = element_line(color='black', linewidth = 1),
        plot.title = element_text(color = 'black', size = 24, hjust=0.5),
        panel.grid.major.y = element_line(color = '#aaaaaa'),
        panel.grid.major.x = element_line(color = '#aaaaaa'),
        panel.background = element_rect(fill='white'),
        legend.position = 'none',
        legend.text = element_text(color = 'black', size = 12),
        legend.title=element_blank())

showFigure(addRiskLegend(network_combined_bars_p10, 'network'), 'network_combined_bars_p10', width = 12, height = 6)


#### Network Glass Brain ----

# Nodes of the four networks colored by network; edges are the connections with p.adj < .10 (the same ones as the bar chart above),
# blue where Low Risk is higher and red where High Risk is higher. Needs the network tables written by the network analyses.

loadGlassBrain()

regionLabels <- read.csv(paste0(baseDir,'/functional/100PerHemiLabels_Decimal.1D'), header=F)
roiNames <- paste0(regionLabels$V1, '.', seq_len(nrow(regionLabels)))

networkGlassNodes <- do.call(rbind, lapply(names(networkNodeColors), function(net) {
  ids <- which(networkSpecs[[net]]$match(roiNames))
  data.frame(id = ids, fill = networkNodeColors[[net]], network = net, stringsAsFactors = FALSE)
}))

networkGlassEdges <- do.call(rbind, lapply(names(networkNodeColors), function(net) {
  tbl <- read.csv(paste0(baseDir,'/functional/output/', net, '_Intx_P10FDR.csv'), stringsAsFactors = FALSE)
  if (nrow(tbl) == 0) { return(NULL) }
  edges <- connectionRegions(tbl$ROI)
  edges$dir <- ifelse(tbl$TD.Mean > tbl$ASD.Mean, 'low', 'high')
  edges
}))

message('Network glass brain: ', sum(networkGlassEdges$dir == 'low', na.rm = TRUE), ' Low > High Risk edges and ',
        sum(networkGlassEdges$dir == 'high', na.rm = TRUE), ' High > Low Risk edges')

# Legend: one chip per network (node color with a lighter tint), then what the blue and red edges mean
networkLegendPanel <- function() {
  nNet <- length(networkNodeColors); yChip <- seq(.80, by = -.14, length.out = nNet)
  tints <- vapply(networkNodeColors, tintColor, character(1), amount = .65)
  ggplot() + xlim(0, 1) + ylim(0, 1) + coord_fixed() + theme_void() +
    annotate('rect', xmin = .02, xmax = .98, ymin = .04, ymax = .96, fill = 'white', colour = 'black', linewidth = 1) +
    annotate('rect', xmin = .08, xmax = .14, ymin = yChip - .05, ymax = yChip + .05, fill = unname(networkNodeColors), colour = '#444444', linewidth = .4) +
    annotate('rect', xmin = .14, xmax = .92, ymin = yChip - .05, ymax = yChip + .05, fill = unname(tints), colour = '#444444', linewidth = .4) +
    annotate('text', x = .17, y = yChip, label = paste(networkNodeLabels[names(networkNodeColors)], 'Network'), hjust = 0, size = 4.2) +
    annotate('segment', x = .04, xend = .96, y = .27, yend = .27, linewidth = .5) +
    annotate('text', x = c(.27, .73), y = .20, label = c('Low > High Risk', 'High > Low Risk'), size = 3.6) +
    annotate('segment', x = c(.12, .58), xend = c(.42, .88), y = .12, yend = .12, linewidth = 2.2,
             colour = c(groupLineColor[['low']], groupLineColor[['high']]), lineend = 'round')
}

network_glass <- ggpubr::ggarrange(glassBrainView('axial', glassBrainBackgrounds, roiCoordinates, networkGlassNodes, networkGlassEdges, nodeSize = 3.2, edgeWidth = 1.4),
                                   glassBrainView('posterior', glassBrainBackgrounds, roiCoordinates, networkGlassNodes, networkGlassEdges, nodeSize = 3.2, edgeWidth = 1.4),
                                   glassBrainView('lateral', glassBrainBackgrounds, roiCoordinates, networkGlassNodes, networkGlassEdges, nodeSize = 3.2, edgeWidth = 1.4),
                                   networkLegendPanel(), nrow = 2, ncol = 2)
showFigure(network_glass, 'network_glass', width = 9, height = 9)


# 4. Anatomical-Functional Overlap (optional; needs both 2a and 3a) ----


fetalMatched

asdDemoFinal$ID[!c(asdDemoFinal$ID %in% fetalMatched$longID)]

subset(fetalMatched, mchat_result=='Fail')$longID[!c(subset(fetalMatched, mchat_result=='Fail')$longID %in% asdDemoFinal$ID)]

