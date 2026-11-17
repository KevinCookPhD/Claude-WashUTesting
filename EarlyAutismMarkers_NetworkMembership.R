# NETWORK MEMBERSHIP MOCK-UP (fetal parcels, Yeo-style networks) ----
#
# Proposes a network for every one of the 200 fetal parcels and draws it next to the Yeo 7-network (Schaefer 400) parcellation.
#
#   Cerebellar    every parcel in the cerebellum
#   Subcortical   caudate, putamen, thalamus, pallidum (amygdala, hippocampus and parahippocampal: see limbicStructures)
#   Brainstem     brainstem parcels (kept apart)
#   Cortical      one of the 7 Yeo networks, taken from where the parcel falls on the Schaefer/Yeo atlas after an APPROXIMATE
#                 registration (the fetal cortex is scaled onto the adult MNI cortex; two simple versions are run and the parcel is
#                 marked 'high' confidence when both agree). This is a mock-up to compare with your own hand classification, not a
#                 proper registration: it is most reliable toward the back of the brain and least reliable frontally and centrally.
#
# Needs only base R plus ggplot2 and ggpubr; the NIfTI files are read with base R. It also compares with your current network lists
# (networkSpecs) when EarlyAutismMarkers_edit.R has been run in the same session.
#
# OUTPUT: a table with one row per parcel (tableDir/networkMembership/networkMembership_proposed.csv), a printed summary, and the figure
# figDir/networkMembership/networkMembership_slices (a slice through the Yeo atlas, the fetal parcels as proposed here, and the fetal parcels as
# in your current lists).

# Settings ----

# The three input files
fetalAtlasFile  <- 'tcorr_cluster_100ROIs_per_hemi.nii'                         # fetal parcellation (labels 1-200)
fetalLabelFile  <- 'r200_Net.csv'                                               # one row per parcel: ROINum, Name_Abrv, Lobe, X, Y, Z, ...
schaeferFile    <- 'Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm.nii'  # Yeo 7-network cortical atlas in MNI space

# Amygdala, hippocampus and parahippocampal parcels are subcortical structures but are in your Limbic network. Which network
# should they get? 'Limbic' (as in your current lists) or 'Subcortical'.
limbicStructures <- c('AMY', 'HIP', 'PHG')
limbicStructureNetwork <- 'Limbic'

# How far (mm, in MNI space) from a parcel voxel to look for the nearest Yeo parcel
searchRadius <- 6

# Slices shown in the figure, as the fraction of the way through the fetal cortex (axial: bottom to top, sagittal: left to right, coronal: back to front)
sliceFractions <- c(axial = .45, sagittal = .28, coronal = .55)

# Outlines: draw a black outline around the fetal parcels that are in dispute. 'In dispute' = the parcel is in one of your current lists and
# its proposed network is different. outlineLowConfidence also outlines (in orange) the cortical parcels where the two registrations disagree.
# Without your current lists (networkSpecs) only the low-confidence parcels can be outlined.
outlineDisagreements <- TRUE
outlineLowConfidence <- FALSE

# Parcel numbers: write the fetal parcel number (ROINum, the same number as in the tables) on each fetal parcel in the figures, so a decision can be
# referred to by number ("parcel 27 ...").
labelParcels <- TRUE
labelSize2D <- 2.3             # text size in the slice figure
labelSize3D <- 0.9             # text size in the 3D figure

# 3D rendering (needs the rgl package and a display; set FALSE to skip). Drawn as the left hemisphere, blocky parcels (voxels), lateral and medial views.
make3D <- TRUE
yeo3DStep <- 2                 # the Yeo atlas is 1 mm; 2 draws every second voxel (faster, still smooth enough)
interactive3D <- FALSE         # TRUE also opens an interactive rgl window of the fetal parcels as proposed (rotate with the mouse)

# Folder name (inside the figure and table folders) for the outputs
membershipFolder <- 'networkMembership'
figureTypes <- c('pdf')

# Reading the NIfTI files (base R) ----

readNifti <- function(path) {
  raw <- if (grepl('\\.gz$', path)) { con <- gzfile(path, 'rb'); on.exit(close(con)); readBin(con, 'raw', n = 5e8) } else readBin(path, 'raw', n = file.size(path))
  endian <- if (readBin(raw[1:4], 'integer', 1, size = 4, endian = 'little') == 348) 'little' else 'big'
  dims <- readBin(raw[41:56], 'integer', 8, size = 2, endian = endian)
  datatype <- readBin(raw[71:72], 'integer', 1, size = 2, endian = endian)
  voxOffset <- readBin(raw[109:112], 'numeric', 1, size = 4, endian = endian)
  srow <- matrix(readBin(raw[281:328], 'numeric', 12, size = 4, endian = endian), 3, 4, byrow = TRUE)
  n <- prod(dims[2:4])
  body <- raw[(voxOffset + 1):length(raw)]
  values <- switch(as.character(datatype),
                   '2' = as.integer(body[1:n]), '4' = readBin(body, 'integer', n, size = 2, endian = endian),
                   '8' = readBin(body, 'integer', n, size = 4, endian = endian), '16' = readBin(body, 'numeric', n, size = 4, endian = endian),
                   '64' = readBin(body, 'numeric', n, size = 8, endian = endian), stop('Unsupported NIfTI data type ', datatype))
  list(data = array(values, dim = dims[2:4]), affine = rbind(srow, c(0, 0, 0, 1)))
}
voxelToWorld <- function(index, affine) (cbind(index - 1, 1) %*% t(affine))[, 1:3, drop = FALSE]      # index: n x 3, 1-based
worldToVoxel <- function(xyz, affine) round((cbind(xyz, 1) %*% t(solve(affine)))[, 1:3, drop = FALSE]) + 1

# Yeo networks ----

yeoNetworks <- c('Vis', 'SomMot', 'DorsAttn', 'SalVentAttn', 'Limbic', 'Cont', 'Default')
yeoNames <- c(Vis = 'Visual', SomMot = 'Somatomotor', DorsAttn = 'Dorsal Attention', SalVentAttn = 'Salience / Ventral Attention',
              Limbic = 'Limbic', Cont = 'Control', Default = 'Default Mode')
# The Schaefer 400-parcel labels run in blocks by network: 1-200 left hemisphere, 201-400 right hemisphere
schaeferBlocks <- data.frame(network = rep(yeoNetworks, 2),
                             first = c(1, 32, 69, 92, 114, 127, 149, 201, 231, 271, 294, 319, 332, 362),
                             last  = c(31, 68, 91, 113, 126, 148, 200, 230, 270, 293, 318, 331, 361, 400), stringsAsFactors = FALSE)
networkColors <- c('Visual' = '#781286', 'Somatomotor' = '#4682B4', 'Dorsal Attention' = '#00760E', 'Salience / Ventral Attention' = '#C43AFA',
                   'Limbic' = '#DCF8A4', 'Control' = '#E69422', 'Default Mode' = '#CD3E4E',
                   'Cerebellar' = '#00B3B3', 'Subcortical' = '#8C510A', 'Brainstem' = '#BDBDBD', 'Auditory (your list)' = '#E377C2', 'Not in a list' = '#F0F0F0')

# Load the data ----

for (f in c(fetalAtlasFile, fetalLabelFile, schaeferFile)) { if (!file.exists(f)) stop('File not found: ', f, '. Set the file paths at the top of the script.') }
message('Reading the atlases...')
fetal <- readNifti(fetalAtlasFile)
schaefer <- readNifti(schaeferFile)
labels <- read.csv(fetalLabelFile, stringsAsFactors = FALSE)
labels <- labels[order(labels$ROINum), ]
if (!identical(sort(unique(as.integer(fetal$data[fetal$data > 0]))), 1:200)) { stop('The fetal atlas must contain the labels 1-200') }

schaeferNetworkArray <- array(0L, dim = dim(schaefer$data))
for (b in seq_len(nrow(schaeferBlocks))) {
  schaeferNetworkArray[schaefer$data >= schaeferBlocks$first[b] & schaefer$data <= schaeferBlocks$last[b]] <- match(schaeferBlocks$network[b], yeoNetworks)
}

# Which fetal parcels are cortical, and where their voxels are
fetalIndex <- which(fetal$data > 0, arr.ind = TRUE)
fetalLabel <- fetal$data[fetalIndex]
fetalWorld <- voxelToWorld(fetalIndex, fetal$affine)
isCortical <- labels$Lobe %in% c('Frontal', 'Parietal', 'Temporal', 'Occipital') & !sub('\\..*', '', labels$Name_Abrv) %in% limbicStructures
corticalParcels <- labels$ROINum[isCortical]
corticalVoxels <- fetalLabel %in% corticalParcels

# Approximate registration: fetal cortex onto the adult cortex ----

schaeferWorld <- voxelToWorld(which(schaeferNetworkArray > 0, arr.ind = TRUE), schaefer$affine)
box <- function(xyz) rbind(lo = apply(xyz, 2, quantile, .01), hi = apply(xyz, 2, quantile, .99))
fetalBox <- box(fetalWorld[corticalVoxels, ]); adultBox <- box(schaeferWorld)
scaleA <- (adultBox['hi', ] - adultBox['lo', ]) / (fetalBox['hi', ] - fetalBox['lo', ])
offsetA <- adultBox['lo', ] - fetalBox['lo', ] * scaleA
registerA <- function(xyz) sweep(sweep(xyz, 2, scaleA, '*'), 2, offsetA, '+')
fetalCentre <- colMeans(fetalWorld[corticalVoxels, ]); adultCentre <- colMeans(schaeferWorld)
scaleB <- sqrt(mean(rowSums(sweep(schaeferWorld, 2, adultCentre)^2)) / mean(rowSums(sweep(fetalWorld[corticalVoxels, ], 2, fetalCentre)^2)))
registerB <- function(xyz) sweep(sweep(xyz, 2, fetalCentre) * scaleB, 2, adultCentre, '+')

# Nearest Yeo network within searchRadius of each mapped voxel
offsets <- expand.grid(dx = -searchRadius:searchRadius, dy = -searchRadius:searchRadius, dz = -searchRadius:searchRadius)
offsets <- offsets[sqrt(rowSums(offsets^2)) <= searchRadius, ]
offsets <- offsets[order(rowSums(offsets^2)), ]
nearestYeo <- function(mniXyz) {
  index <- worldToVoxel(mniXyz, schaefer$affine)
  dimS <- dim(schaeferNetworkArray)
  vapply(seq_len(nrow(index)), function(i) {
    cand <- cbind(index[i, 1] + offsets$dx, index[i, 2] + offsets$dy, index[i, 3] + offsets$dz)
    ok <- cand[, 1] >= 1 & cand[, 1] <= dimS[1] & cand[, 2] >= 1 & cand[, 2] <= dimS[2] & cand[, 3] >= 1 & cand[, 3] <= dimS[3]
    values <- schaeferNetworkArray[cand[ok, , drop = FALSE]]
    hit <- which(values > 0)
    if (length(hit) == 0) 0L else values[hit[1]]
  }, integer(1))
}
parcelNetwork <- function(register) {
  voxels <- which(corticalVoxels)
  codes <- nearestYeo(register(fetalWorld[voxels, , drop = FALSE]))
  out <- lapply(corticalParcels, function(p) {
    tab <- table(codes[fetalLabel[voxels] == p & codes > 0])
    if (length(tab) == 0) return(c(NA, NA))
    c(as.integer(names(tab)[which.max(tab)]), max(tab) / sum(tab))
  })
  data.frame(ROINum = corticalParcels, code = vapply(out, `[`, numeric(1), 1), purity = vapply(out, `[`, numeric(1), 2))
}
message('Placing the fetal parcels on the Yeo atlas (two approximate registrations)...')
placeA <- parcelNetwork(registerA); placeB <- parcelNetwork(registerB)

# The proposed network for every parcel ----

membership <- labels[, c('ROINum', 'Name_Abrv', 'ROIName', 'Hemisphere', 'Lobe', 'X', 'Y', 'Z')]
abbreviation <- sub('\\..*', '', membership$Name_Abrv)
membership$yeoA <- yeoNames[yeoNetworks[placeA$code[match(membership$ROINum, placeA$ROINum)]]]
membership$yeoB <- yeoNames[yeoNetworks[placeB$code[match(membership$ROINum, placeB$ROINum)]]]
membership$purityA <- round(placeA$purity[match(membership$ROINum, placeA$ROINum)], 2)
membership$proposed <- ifelse(membership$Lobe == 'Cerebellum', 'Cerebellar',
                       ifelse(abbreviation == 'BS', 'Brainstem',
                       ifelse(abbreviation %in% limbicStructures, limbicStructureNetwork,
                       ifelse(membership$Lobe == 'Subcortical', 'Subcortical', membership$yeoA))))
membership$basis <- ifelse(membership$proposed %in% c('Cerebellar', 'Brainstem', 'Subcortical') | abbreviation %in% limbicStructures, 'anatomy',
                    ifelse(!is.na(membership$yeoA) & !is.na(membership$yeoB) & membership$yeoA == membership$yeoB & membership$purityA >= .6, 'Yeo, high confidence', 'Yeo, low confidence'))

# Your current network lists, if they are available in this session
currentOrder <- c('CRB', 'SubCort', 'Limbic', 'DMN', 'Visual', 'SomMot', 'Auditory')
currentLabels <- c(CRB = 'Cerebellar', SubCort = 'Subcortical', Limbic = 'Limbic', DMN = 'Default Mode', Visual = 'Visual', SomMot = 'Somatomotor', Auditory = 'Auditory (your list)')
haveCurrent <- exists('networkSpecs') && exists('roiNumber')
if (haveCurrent) {
  roiNames <- ifelse(grepl('\\.[LR]$', membership$Name_Abrv), paste0(membership$Name_Abrv, '.', membership$ROINum), paste0(membership$Name_Abrv, '.', membership$ROINum))
  inList <- sapply(currentOrder, function(n) if (n %in% names(networkSpecs)) networkSpecs[[n]]$match(roiNames) else rep(FALSE, length(roiNames)))
  membership$currentLists <- apply(inList, 1, function(r) paste(currentOrder[r], collapse = '/'))
  membership$current <- apply(inList, 1, function(r) if (!any(r)) 'Not in a list' else currentLabels[[currentOrder[which(r)[1]]]])
}

# Summary ----

cat('\n=== Proposed network for the 200 fetal parcels ===\n')
print(table(membership$proposed, membership$basis))
cat('\nParcels where the two registrations disagree (low confidence): ', sum(membership$basis == 'Yeo, low confidence'), ' of ', sum(grepl('^Yeo', membership$basis)), ' cortical parcels\n', sep = '')
if (haveCurrent) {
  cat('\nYour current lists (rows) against the proposed networks (columns):\n')
  print(table(`current list` = membership$current, proposed = membership$proposed))
  changed <- membership[membership$current != 'Not in a list' & membership$current != membership$proposed, c('ROINum', 'Name_Abrv', 'currentLists', 'current', 'proposed', 'basis')]
  cat('\nParcels in one of your lists whose proposed network is different:', nrow(changed), '\n')
  oldWidth <- options(width = 200); print(changed, row.names = FALSE); options(oldWidth)
}
if (exists('tableDir')) {
  outDir <- file.path(tableDir, membershipFolder); dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
  write.csv(membership, file.path(outDir, 'networkMembership_proposed.csv'), row.names = FALSE)
}

# Parcels to outline ----

disputedIds <- if (haveCurrent) membership$ROINum[membership$current != 'Not in a list' & membership$current != membership$proposed] else integer(0)
lowConfidenceIds <- membership$ROINum[membership$basis == 'Yeo, low confidence']
outlineGroups <- list()
if (outlineDisagreements) {
  if (haveCurrent) outlineGroups[['disputed']] <- list(ids = disputedIds, color = '#000000')
  if (!haveCurrent || outlineLowConfidence) outlineGroups[['lowConfidence']] <- list(ids = lowConfidenceIds, color = if (haveCurrent) '#E65100' else '#000000')
}
cat('\nOutlined parcels:', paste(vapply(names(outlineGroups), function(g) paste0(g, ' = ', length(outlineGroups[[g]]$ids)), character(1)), collapse = '; '), '\n')

# Figure: slices through the Yeo atlas, the proposed membership and your current lists ----

networkOfParcel <- function(column) { v <- membership[[column]][match(seq_len(200), membership$ROINum)]; v }
fetalNetworkArray <- function(column) {
  v <- networkOfParcel(column); out <- array(NA_character_, dim = dim(fetal$data)); pos <- fetal$data > 0
  out[pos] <- v[fetal$data[pos]]; out
}
yeoNetworkArray <- array(NA_character_, dim = dim(schaeferNetworkArray))
yeoNetworkArray[schaeferNetworkArray > 0] <- yeoNames[yeoNetworks[schaeferNetworkArray[schaeferNetworkArray > 0]]]

# One slice of a labelled volume: view = 'axial', 'sagittal' or 'coronal'; coord = the world coordinate of the slice
sliceOf <- function(volume, affine, view, coord, parcelVolume = NULL) {
  axis <- c(axial = 3, sagittal = 1, coronal = 2)[[view]]
  scaleAxis <- diag(affine[1:3, 1:3])[axis]; offsetAxis <- affine[axis, 4]
  idx <- round((coord - offsetAxis) / scaleAxis) + 1
  idx <- max(1, min(dim(volume)[axis], idx))
  plane <- switch(view, axial = volume[, , idx], sagittal = volume[idx, , ], coronal = volume[, idx, ])
  parcelPlane <- if (is.null(parcelVolume)) NULL else switch(view, axial = parcelVolume[, , idx], sagittal = parcelVolume[idx, , ], coronal = parcelVolume[, idx, ])
  axes <- switch(view, axial = c(1, 2), sagittal = c(2, 3), coronal = c(1, 3))
  pos <- which(!is.na(plane), arr.ind = TRUE)
  index3 <- matrix(idx, nrow(pos), 3); index3[, axes[1]] <- pos[, 1]; index3[, axes[2]] <- pos[, 2]
  world <- voxelToWorld(index3, affine)
  data.frame(h = world[, axes[1]], v = world[, axes[2]], network = plane[pos], parcel = if (is.null(parcelPlane)) NA_integer_ else parcelPlane[pos],
             step = abs(diag(affine[1:3, 1:3])[axes[1]]), stringsAsFactors = FALSE)
}

limitsF <- fetalBox
sliceCoordF <- c(axial = unname(limitsF['lo', 3] + sliceFractions[['axial']] * (limitsF['hi', 3] - limitsF['lo', 3])),
                 sagittal = unname(limitsF['lo', 1] + sliceFractions[['sagittal']] * (limitsF['hi', 1] - limitsF['lo', 1])),
                 coronal = unname(limitsF['lo', 2] + sliceFractions[['coronal']] * (limitsF['hi', 2] - limitsF['lo', 2])))
axisOf <- c(axial = 3, sagittal = 1, coronal = 2)
sliceCoordM <- vapply(names(sliceCoordF), function(v) unname(sliceCoordF[[v]] * scaleA[axisOf[[v]]] + offsetA[axisOf[[v]]]), numeric(1))

# Segments along the edges of the listed parcels (edges that touch a different parcel or empty space)
outlineSegments <- function(df, ids) {
  if (length(ids) == 0 || all(is.na(df$parcel))) return(NULL)
  step <- df$step[1]; h0 <- min(df$h); v0 <- min(df$v)
  ih <- round((df$h - h0) / step); iv <- round((df$v - v0) / step)
  lookup <- setNames(df$parcel, paste(ih, iv))
  segs <- list()
  for (d in list(c(1, 0), c(-1, 0), c(0, 1), c(0, -1))) {
    inside <- which(df$parcel %in% ids)
    neighbour <- lookup[paste(ih[inside] + d[1], iv[inside] + d[2])]
    edge <- inside[is.na(neighbour) | neighbour != df$parcel[inside]]
    if (length(edge) == 0) next
    if (d[1] != 0) segs[[length(segs) + 1]] <- data.frame(x = df$h[edge] + d[1] * step / 2, xend = df$h[edge] + d[1] * step / 2, y = df$v[edge] - step / 2, yend = df$v[edge] + step / 2)
    else segs[[length(segs) + 1]] <- data.frame(x = df$h[edge] - step / 2, xend = df$h[edge] + step / 2, y = df$v[edge] + d[2] * step / 2, yend = df$v[edge] + d[2] * step / 2)
  }
  if (length(segs) == 0) NULL else do.call(rbind, segs)
}

textColorFor <- function(hex) { rgbValues <- grDevices::col2rgb(hex) / 255; ifelse(.299 * rgbValues[1, ] + .587 * rgbValues[2, ] + .114 * rgbValues[3, ] < .5, 'white', 'black') }

panelPlot <- function(df, title, view, tile, outlines = list(), label = FALSE) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = h, y = v, fill = network)) +
    ggplot2::geom_tile(width = tile, height = tile) +
    ggplot2::scale_fill_manual(values = networkColors, limits = names(networkColors), drop = FALSE, name = NULL) +
    ggplot2::coord_fixed() + ggplot2::labs(title = title) + ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10, hjust = .5))
  for (g in outlines) {
    seg <- outlineSegments(df, g$ids)
    if (!is.null(seg)) p <- p + ggplot2::geom_segment(data = seg, ggplot2::aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE, color = g$color, linewidth = .7, lineend = 'square')
  }
  if (label && any(!is.na(df$parcel))) {
    centres <- do.call(rbind, lapply(split(df[!is.na(df$parcel), ], df$parcel[!is.na(df$parcel)]), function(d)
      data.frame(h = mean(d$h), v = mean(d$v), parcel = d$parcel[1], network = d$network[1], n = nrow(d))))
    centres$textColor <- textColorFor(networkColors[centres$network])
    p <- p + ggplot2::geom_text(data = centres, ggplot2::aes(x = h, y = v, label = parcel, color = textColor), inherit.aes = FALSE, size = labelSize2D) +
      ggplot2::scale_color_identity()
  }
  if (view == 'sagittal') p <- p + ggplot2::scale_x_reverse()      # front on the left
  p
}
fetalProposed <- fetalNetworkArray('proposed')
fetalCurrent <- if (haveCurrent) fetalNetworkArray('current') else NULL
columns <- c('Yeo 7 (Schaefer 400, MNI)', 'Fetal parcels: proposed', if (haveCurrent) 'Fetal parcels: your current lists')
plots <- list()
for (view in c('axial', 'sagittal', 'coronal')) {
  yeoSlice <- sliceOf(yeoNetworkArray, schaefer$affine, view, sliceCoordM[[view]])
  plots[[length(plots) + 1]] <- panelPlot(yeoSlice, if (view == 'axial') columns[1] else '', view, yeoSlice$step[1])
  fs <- sliceOf(fetalProposed, fetal$affine, view, sliceCoordF[[view]], fetal$data)
  plots[[length(plots) + 1]] <- panelPlot(fs, if (view == 'axial') columns[2] else '', view, fs$step[1], outlineGroups, labelParcels)
  if (haveCurrent) {
    fc <- sliceOf(fetalCurrent, fetal$affine, view, sliceCoordF[[view]], fetal$data)
    plots[[length(plots) + 1]] <- panelPlot(fc, if (view == 'axial') columns[3] else '', view, fc$step[1], outlineGroups, labelParcels)
  }
}
networkMembership_slices <- ggpubr::ggarrange(plotlist = plots, ncol = length(columns), nrow = 3, common.legend = TRUE, legend = 'bottom')
networkMembership_slices <- ggpubr::annotate_figure(networkMembership_slices,
  bottom = ggpubr::text_grob(paste0('Rows: axial, sagittal (front on the left), coronal slices; the fetal slices are at the same relative positions as the Yeo slices (approximate registration).', if (length(outlineGroups) > 0) ' Outlined = parcels in dispute.' else '', if (labelParcels) ' Numbers = fetal parcel numbers (ROINum).' else ''), size = 8, color = '#444444'))

if (exists('showFigure') && exists('saveFigure')) {
  oldSave <- saveFigures; oldTypes <- figureFileTypes
  saveFigures <- TRUE; figureFileTypes <- figureTypes
  showFigure(networkMembership_slices, 'networkMembership_slices', width = 4 * length(columns) + 1, height = 11, subfolder = membershipFolder)
  saveFigures <- oldSave; figureFileTypes <- oldTypes
} else {
  print(networkMembership_slices)
  outFig <- if (exists('figDir')) file.path(figDir, membershipFolder) else file.path(getwd(), membershipFolder)
  dir.create(outFig, showWarnings = FALSE, recursive = TRUE)
  for (type in figureTypes) ggplot2::ggsave(file.path(outFig, paste0('networkMembership_slices.', type)), networkMembership_slices, width = 4 * length(columns) + 1, height = 11)
}

# 3D rendering ----
#
# The left hemisphere drawn as voxel blocks: the Yeo atlas, the fetal parcels as proposed here, and (when available) the fetal parcels as in
# your current lists; lateral view (top row) and medial view (bottom row). Outlines are drawn along the borders of the parcels in dispute.

# Surface of a labelled volume: the faces of every occupied voxel that touch empty space.
#   parcelVol: integer array, 0 = empty, otherwise a parcel number. Returns the quads (4 vertices each), the parcel of each quad and its edges.
cross3 <- function(a, b) c(a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
voxelSurface <- function(parcelVol, affine, step = 1) {
  if (step > 1) {
    parcelVol <- parcelVol[seq(1, dim(parcelVol)[1], by = step), seq(1, dim(parcelVol)[2], by = step), seq(1, dim(parcelVol)[3], by = step)]
    affine[, 1:3] <- affine[, 1:3] * step
  }
  d <- dim(parcelVol)
  padded <- array(0L, dim = d + 2); padded[2:(d[1] + 1), 2:(d[2] + 1), 2:(d[3] + 1)] <- parcelVol
  idx <- which(parcelVol > 0, arr.ind = TRUE); id <- parcelVol[idx]
  corners <- list(c(-.5, -.5), c(.5, -.5), c(.5, .5), c(-.5, .5))
  toWorld <- function(m) (cbind(m, 1) %*% t(affine))[, 1:3, drop = FALSE]
  quadParts <- list(); idParts <- list(); normalParts <- list(); edgeParts <- list()
  for (axis in 1:3) for (side in c(-1, 1)) {
    nb <- idx; nb[, axis] <- nb[, axis] + side
    exposed <- padded[nb + 1] == 0
    if (!any(exposed)) next
    base <- idx[exposed, , drop = FALSE] - 1; fid <- id[exposed]; others <- setdiff(1:3, axis)
    cornerIndex <- lapply(corners, function(cc) { m <- base; m[, axis] <- m[, axis] + side * .5; m[, others[1]] <- m[, others[1]] + cc[1]; m[, others[2]] <- m[, others[2]] + cc[2]; m })
    cornerWorld <- lapply(cornerIndex, toWorld)
    n <- nrow(base)
    # make the corner order run counter-clockwise seen from outside (the voxel-to-world transform can mirror it), so lighting is right
    outward <- affine[1:3, axis] * side
    cornerOrder <- 1:4
    if (sum(cross3(cornerWorld[[2]][1, ] - cornerWorld[[1]][1, ], cornerWorld[[3]][1, ] - cornerWorld[[1]][1, ]) * outward) < 0) {
      cornerIndex <- cornerIndex[c(1, 4, 3, 2)]; cornerWorld <- cornerWorld[c(1, 4, 3, 2)]
    }
    quadParts[[length(quadParts) + 1]] <- do.call(rbind, cornerWorld)[order(rep(seq_len(n), 4)), , drop = FALSE]
    idParts[[length(idParts) + 1]] <- fid
    normalDirection <- affine[1:3, axis] * side; normalDirection <- normalDirection / sqrt(sum(normalDirection^2))
    normalParts[[length(normalParts) + 1]] <- matrix(normalDirection, n, 3, byrow = TRUE)
    # the four edges of every face, as integer corner codes (index coordinates times two)
    code <- function(m) { m2 <- round(m * 2) + 2; (m2[, 1] * 300 + m2[, 2]) * 300 + m2[, 3] }
    for (e in 1:4) {
      a <- e; b <- e %% 4 + 1
      ca <- code(cornerIndex[[a]]); cb <- code(cornerIndex[[b]])
      edgeParts[[length(edgeParts) + 1]] <- data.frame(key = pmin(ca, cb) * 27000000 + pmax(ca, cb), id = fid,
                                                       x1 = cornerWorld[[a]][, 1], y1 = cornerWorld[[a]][, 2], z1 = cornerWorld[[a]][, 3],
                                                       x2 = cornerWorld[[b]][, 1], y2 = cornerWorld[[b]][, 2], z2 = cornerWorld[[b]][, 3],
                                                       nx = normalDirection[1], ny = normalDirection[2], nz = normalDirection[3])
    }
  }
  list(vertices = do.call(rbind, quadParts), id = unlist(idParts), normals = do.call(rbind, normalParts), edges = do.call(rbind, edgeParts),
       size = mean(sqrt(colSums(affine[1:3, 1:3]^2))))
}

# Border segments around the parcels in `ids`: edges shared by faces of two different parcels where one of them is in ids
borderEdges <- function(surface, ids) {
  e <- surface$edges
  if (length(ids) == 0) return(NULL)
  nParcels <- tapply(e$id, e$key, function(v) length(unique(v)))
  touches <- tapply(e$id %in% ids, e$key, any)
  keep <- names(nParcels)[nParcels > 1 & touches]
  inKeep <- as.character(e$key) %in% keep
  if (!any(inKeep)) return(NULL)
  # lift each border segment slightly off the surface, along the average direction the faces meeting there point, so it is not hidden by them
  meanNormal <- aggregate(e[inKeep, c('nx', 'ny', 'nz')], by = list(key = e$key[inKeep]), FUN = mean)
  sel <- e[inKeep & !duplicated(e$key), ]
  m <- meanNormal[match(sel$key, meanNormal$key), c('nx', 'ny', 'nz')]
  len <- sqrt(rowSums(m^2)); len[len == 0] <- 1
  lift <- liftFraction * surface$size / len
  for (axisName in c('x', 'y', 'z')) {
    shift <- m[[paste0('n', axisName)]] * lift
    sel[[paste0(axisName, '1')]] <- sel[[paste0(axisName, '1')]] + shift
    sel[[paste0(axisName, '2')]] <- sel[[paste0(axisName, '2')]] + shift
  }
  sel
}
liftFraction <- 0.25

render3D <- function(parcelVol, affine, step, colorOfParcel, view, outlines = list(), pixels = 900, labelIds = NULL) {
  surface <- voxelSurface(parcelVol, affine, step)
  rgl::open3d(useNULL = FALSE, windowRect = c(0, 0, pixels, pixels))
  on.exit(rgl::close3d())
  rgl::bg3d('white'); rgl::clear3d('lights')
  rgl::light3d(theta = 0, phi = 0, viewpoint.rel = TRUE, ambient = '#808080', diffuse = '#E0E0E0', specular = '#000000')
  rgl::quads3d(surface$vertices, color = rep(colorOfParcel[surface$id], each = 4), normals = surface$normals[rep(seq_len(nrow(surface$normals)), each = 4), ],
               specular = 'black', lit = TRUE)
  for (g in outlines) {
    border <- borderEdges(surface, g$ids)
    if (!is.null(border)) rgl::segments3d(rbind(as.matrix(border[, c('x1', 'y1', 'z1')]), as.matrix(border[, c('x2', 'y2', 'z2')]))[order(rep(seq_len(nrow(border)), 2)), ],
                                          color = g$color, lwd = 4, lit = FALSE)
  }
  if (!is.null(labelIds)) {
    # a number on each parcel that can be seen from this view: the faces that point at the camera and have nothing in front of them
    facing <- if (view == 'lateral') surface$normals[, 1] < -.5 else surface$normals[, 1] > .5
    nQuad <- nrow(surface$vertices) / 4
    centres <- rowsum(surface$vertices, rep(seq_len(nQuad), each = 4)) / 4
    key <- paste(round(centres[, 2] / surface$size), round(centres[, 3] / surface$size))
    depth <- if (view == 'lateral') centres[, 1] else -centres[, 1]
    front <- ave(depth, key, FUN = function(v) v == min(v))
    visible <- facing & front == 1 & surface$id %in% labelIds
    if (any(visible)) {
      pos <- do.call(rbind, lapply(split(which(visible), surface$id[visible]), function(i) colMeans(centres[i, , drop = FALSE])))
      pos[, 1] <- pos[, 1] + (if (view == 'lateral') -1 else 1) * surface$size * .6
      ids <- as.integer(rownames(pos))
      rgl::text3d(pos, texts = as.character(ids), color = textColorFor(colorOfParcel[ids]), cex = labelSize3D, font = 2, adj = c(.5, .5))
    }
  }
  R <- rgl::rotationMatrix
  rgl::par3d(userMatrix = switch(view, lateral = R(pi / 2, 0, 1, 0) %*% R(-pi / 2, 1, 0, 0), medial = R(-pi / 2, 0, 1, 0) %*% R(-pi / 2, 1, 0, 0)), zoom = .75)
  file <- tempfile(fileext = '.png'); rgl::snapshot3d(file, webshot = FALSE)
  img <- png::readPNG(file)[, , 1:3, drop = FALSE]; unlink(file)
  used <- which(apply(img < .98, c(1, 2), any), arr.ind = TRUE); pad <- 10
  img <- img[max(1, min(used[, 1]) - pad):min(dim(img)[1], max(used[, 1]) + pad), max(1, min(used[, 2]) - pad):min(dim(img)[2], max(used[, 2]) + pad), , drop = FALSE]
  ggplot2::ggplot() + ggplot2::annotation_raster(as.raster(img), xmin = 0, xmax = ncol(img), ymin = 0, ymax = nrow(img), interpolate = TRUE) +
    ggplot2::xlim(0, ncol(img)) + ggplot2::ylim(0, nrow(img)) + ggplot2::coord_fixed(expand = FALSE) + ggplot2::theme_void()
}

canDraw3D <- make3D && requireNamespace('rgl', quietly = TRUE) && requireNamespace('png', quietly = TRUE) &&
  tryCatch({ rgl::open3d(useNULL = FALSE); ok <- !rgl::rgl.useNULL(); rgl::close3d(); ok }, error = function(e) FALSE)
if (make3D && !canDraw3D) { message('3D rendering skipped: the rgl and png packages and a display are needed') }
if (canDraw3D) {
  message('Drawing the 3D views...')
  # left hemisphere only
  fetalLeft <- fetal$data; fetalLeft[!(fetal$data %in% membership$ROINum[membership$Hemisphere == 'Left'])] <- 0L
  fetalLeft <- array(as.integer(fetalLeft), dim = dim(fetal$data))
  yeoLeft <- array(as.integer(ifelse(schaefer$data >= 1 & schaefer$data <= 200, schaefer$data, 0)), dim = dim(schaefer$data))
  colorFetalProposed <- networkColors[membership$proposed[match(seq_len(200), membership$ROINum)]]
  colorFetalCurrent <- if (haveCurrent) networkColors[membership$current[match(seq_len(200), membership$ROINum)]] else NULL
  yeoBlock <- schaeferBlocks[schaeferBlocks$first <= 200, ]
  colorYeo <- networkColors[yeoNames[vapply(seq_len(200), function(l) yeoBlock$network[yeoBlock$first <= l & yeoBlock$last >= l], character(1))]]
  panels3D <- list()
  for (view in c('lateral', 'medial')) {
    panels3D[[length(panels3D) + 1]] <- render3D(yeoLeft, schaefer$affine, yeo3DStep, colorYeo, view)
    leftIds <- if (labelParcels) membership$ROINum[membership$Hemisphere == 'Left'] else NULL
    panels3D[[length(panels3D) + 1]] <- render3D(fetalLeft, fetal$affine, 1, colorFetalProposed, view, outlineGroups, labelIds = leftIds)
    if (haveCurrent) panels3D[[length(panels3D) + 1]] <- render3D(fetalLeft, fetal$affine, 1, colorFetalCurrent, view, outlineGroups, labelIds = leftIds)
  }
  legendPlot <- ggplot2::ggplot(data.frame(network = factor(names(networkColors), levels = names(networkColors)), x = 1, y = 1), ggplot2::aes(x, y, fill = network)) +
    ggplot2::geom_tile() + ggplot2::scale_fill_manual(values = networkColors, name = NULL) + ggplot2::theme(legend.position = 'bottom')
  headers <- lapply(columns, function(txt) ggplot2::ggplot() + ggplot2::annotate('text', x = 0, y = 0, label = txt, size = 4) + ggplot2::theme_void())
  grid3D <- ggpubr::ggarrange(plotlist = c(headers, panels3D), ncol = length(columns), nrow = 3, heights = c(.08, 1, 1), legend.grob = ggpubr::get_legend(legendPlot), legend = 'bottom')
  networkMembership_3d <- ggpubr::annotate_figure(grid3D, bottom = ggpubr::text_grob(
    paste0('Left hemisphere as voxel blocks: lateral view (top) and medial view (bottom). ', if (length(outlineGroups) > 0) 'Outlined = parcels in dispute. ' else '', if (labelParcels) 'Numbers = fetal parcel numbers (ROINum).' else ''), size = 8, color = '#444444'))
  if (exists('showFigure') && exists('saveFigure')) {
    oldSave <- saveFigures; oldTypes <- figureFileTypes
    saveFigures <- TRUE; figureFileTypes <- figureTypes
    showFigure(networkMembership_3d, 'networkMembership_3d', width = 4 * length(columns) + 1, height = 9, subfolder = membershipFolder)
    saveFigures <- oldSave; figureFileTypes <- oldTypes
  } else {
    print(networkMembership_3d)
    for (type in figureTypes) ggplot2::ggsave(file.path(outFig, paste0('networkMembership_3d.', type)), networkMembership_3d, width = 4 * length(columns) + 1, height = 9)
  }
  if (interactive3D) {
    surface <- voxelSurface(fetalLeft, fetal$affine, 1)
    rgl::open3d(); rgl::bg3d('white')
    rgl::quads3d(surface$vertices, color = rep(colorFetalProposed[surface$id], each = 4), normals = surface$normals[rep(seq_len(nrow(surface$normals)), each = 4), ], specular = 'black')
    for (g in outlineGroups) { b <- borderEdges(surface, g$ids); if (!is.null(b)) rgl::segments3d(rbind(as.matrix(b[, c('x1', 'y1', 'z1')]), as.matrix(b[, c('x2', 'y2', 'z2')]))[order(rep(seq_len(nrow(b)), 2)), ], color = g$color, lwd = 4, lit = FALSE) }
  }
}
