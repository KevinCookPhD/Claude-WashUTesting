# REVISED NETWORK MEMBERSHIP (Cerebellar, Limbic, Default Mode, Subcortical) ----
#
# Overrides the definitions of the four a priori networks in networkSpecs with the membership decided parcel by parcel (see
# membership_decisions.csv for each decision and its reason, and membership_revised.csv for the resulting table).
#
# HOW TO USE: run the Setup of EarlyAutismMarkers_edit.R and the code that defines networkSpecs and roiNumber, then source this file, then run the
# analyses (for example EarlyAutismMarkers_ConnectivityTiers.R). The old definitions are kept in networkSpecsOriginal.
#
# Summary of the changes
#   Limbic       removed 61, 67, 200 (moved to Default Mode) and 94 (insula); added orbitofrontal 50 and 148
#   Default Mode removed 4, 42, 102, 150, 152, 156, 171, 176; added 21, 24, 28, 36, 46, 55, 61, 63, 67, 69, 74, 84, 110, 159, 183, 195, 200
#   Subcortical  added 173 (a putamen parcel that was missing)
#   Cerebellar   unchanged (every parcel named CRB)
# Parcels in more than one of these four networks: none.

if (!exists('networkSpecs')) stop('networkSpecs is not defined; run the Setup and network definitions of EarlyAutismMarkers_edit.R first')
if (!exists('networkSpecsOriginal')) networkSpecsOriginal <- networkSpecs

networkSpecs$Limbic <- list(match = function(roi) roiNumber(roi) %in% c(5,9,20,29,43,50,56,68,81,86,108,120,125,133,141,143,148,155,160,188,196,198),
                            xLabel = 'Limbic Connections', zeroLine = FALSE)

networkSpecs$DMN <- list(match = function(roi) roiNumber(roi) %in% c(11,18,21,22,23,24,28,34,36,40,46,48,53,55,61,63,67,69,70,74,84,85,88,96,100,104,106,110,123,124,128,137,159,163,168,178,183,195,200),
                         xLabel = 'Default Mode Connections', zeroLine = FALSE)

networkSpecs$SubCort <- list(match = function(roi) roiNumber(roi) %in% c(31,32,39,51,80,92,97,107,111,134,139,145,151,173),
                             xLabel = 'Subcortical Connections', zeroLine = FALSE)

# The cerebellar network (every parcel named CRB) is unchanged
networkSpecs$CRB <- list(match = function(roi) grepl('CRB', roi), xLabel = 'Cerebellar Connections', zeroLine = TRUE)

message('Revised network membership loaded: Limbic ', length(c(5,9,20,29,43,50,56,68,81,86,108,120,125,133,141,143,148,155,160,188,196,198)), ', Default Mode ', length(c(11,18,21,22,23,24,28,34,36,40,46,48,53,55,61,63,67,69,70,74,84,85,88,96,100,104,106,110,123,124,128,137,159,163,168,178,183,195,200)),
        ', Subcortical ', length(c(31,32,39,51,80,92,97,107,111,134,139,145,151,173)), ' parcels (plus the 10 cerebellar parcels)')
