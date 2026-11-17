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

# Figure: slices through the Yeo atlas, the proposed membership and your current lists ----

networkOfParcel <- function(column) { v <- membership[[column]][match(seq_len(200), membership$ROINum)]; v }
fetalNetworkArray <- function(column) {
  v <- networkOfParcel(column); out <- array(NA_character_, dim = dim(fetal$data)); pos <- fetal$data > 0
  out[pos] <- v[fetal$data[pos]]; out
}
yeoNetworkArray <- array(NA_character_, dim = dim(schaeferNetworkArray))
yeoNetworkArray[schaeferNetworkArray > 0] <- yeoNames[yeoNetworks[schaeferNetworkArray[schaeferNetworkArray > 0]]]

# One slice of a labelled volume: view = 'axial', 'sagittal' or 'coronal'; coord = the world coordinate of the slice
sliceOf <- function(volume, affine, view, coord) {
  axis <- c(axial = 3, sagittal = 1, coronal = 2)[[view]]
  scaleAxis <- diag(affine[1:3, 1:3])[axis]; offsetAxis <- affine[axis, 4]
  idx <- round((coord - offsetAxis) / scaleAxis) + 1
  idx <- max(1, min(dim(volume)[axis], idx))
  plane <- switch(view, axial = volume[, , idx], sagittal = volume[idx, , ], coronal = volume[, idx, ])
  axes <- switch(view, axial = c(1, 2), sagittal = c(2, 3), coronal = c(1, 3))
  pos <- which(!is.na(plane), arr.ind = TRUE)
  index3 <- matrix(idx, nrow(pos), 3); index3[, axes[1]] <- pos[, 1]; index3[, axes[2]] <- pos[, 2]
  world <- voxelToWorld(index3, affine)
  data.frame(h = world[, axes[1]], v = world[, axes[2]], network = plane[pos], step = abs(diag(affine[1:3, 1:3])[axes[1]]), stringsAsFactors = FALSE)
}

limitsF <- fetalBox
sliceCoordF <- c(axial = unname(limitsF['lo', 3] + sliceFractions[['axial']] * (limitsF['hi', 3] - limitsF['lo', 3])),
                 sagittal = unname(limitsF['lo', 1] + sliceFractions[['sagittal']] * (limitsF['hi', 1] - limitsF['lo', 1])),
                 coronal = unname(limitsF['lo', 2] + sliceFractions[['coronal']] * (limitsF['hi', 2] - limitsF['lo', 2])))
axisOf <- c(axial = 3, sagittal = 1, coronal = 2)
sliceCoordM <- vapply(names(sliceCoordF), function(v) unname(sliceCoordF[[v]] * scaleA[axisOf[[v]]] + offsetA[axisOf[[v]]]), numeric(1))

panelPlot <- function(df, title, view, tile) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = h, y = v, fill = network)) +
    ggplot2::geom_tile(width = tile, height = tile) +
    ggplot2::scale_fill_manual(values = networkColors, limits = names(networkColors), drop = FALSE, name = NULL) +
    ggplot2::coord_fixed() + ggplot2::labs(title = title) + ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10, hjust = .5))
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
  fs <- sliceOf(fetalProposed, fetal$affine, view, sliceCoordF[[view]])
  plots[[length(plots) + 1]] <- panelPlot(fs, if (view == 'axial') columns[2] else '', view, fs$step[1])
  if (haveCurrent) {
    fc <- sliceOf(fetalCurrent, fetal$affine, view, sliceCoordF[[view]])
    plots[[length(plots) + 1]] <- panelPlot(fc, if (view == 'axial') columns[3] else '', view, fc$step[1])
  }
}
networkMembership_slices <- ggpubr::ggarrange(plotlist = plots, ncol = length(columns), nrow = 3, common.legend = TRUE, legend = 'bottom')
networkMembership_slices <- ggpubr::annotate_figure(networkMembership_slices,
  bottom = ggpubr::text_grob('Rows: axial, sagittal (front on the left), coronal slices. The fetal slices are at the same relative positions as the Yeo slices (approximate registration).', size = 8, color = '#444444'))

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
