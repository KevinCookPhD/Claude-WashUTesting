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



## Figure Output Settings ----

# Set saveFigures to TRUE to write every figure to figDir (existing files with the same name are overwritten).
saveFigures <- FALSE

# File types to save, in this order. Choose any of: 'tiff', 'pdf', 'svg', 'eps'
figureFileTypes <- c('tiff', 'pdf', 'svg', 'eps')

# Every ggplot figure is also kept here by name (e.g. figures[['regional_group_p05']]) so it can be re-shown or re-saved
figures <- list()


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

# Save a figure to figDir in each of figureFileTypes (only if saveFigures is TRUE).
# drawFun is a function with no arguments that draws the figure on the current device.
# eps files are written without transparency (the eps format cannot hold it); svg and pdf keep it.
saveFigure <- function(name, width, height, drawFun) {
  if (!saveFigures) { return(invisible(NULL)) }

  badTypes <- setdiff(figureFileTypes, c('tiff', 'pdf', 'svg', 'eps'))
  if (length(badTypes) > 0) { stop('Unknown figure file type(s): ', paste(badTypes, collapse=', '), ". Use 'tiff', 'pdf', 'svg', 'eps'.") }

  dir.create(figDir, showWarnings = FALSE, recursive = TRUE)
  for (type in figureFileTypes) {
    file <- file.path(figDir, paste0(name, '.', type))
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

# Show a ggplot/ggarrange figure, keep it in figures[[name]], and save it. Does nothing for a NULL figure.
showFigure <- function(plotObj, name, width, height) {
  if (is.null(plotObj)) { return(invisible(NULL)) }
  figures[[name]] <<- plotObj
  print(plotObj)
  saveFigure(name, width, height, function() print(plotObj))
  invisible(plotObj)
}

# Same for a base-graphics figure (e.g. circlize chord diagram): drawFun draws it.
showBaseFigure <- function(name, width, height, drawFun) {
  drawFun()
  saveFigure(name, width, height, drawFun)
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

## Behavioral Data (MCHAT) ----


# ----

mchatData <- read.csv(paste0(baseDir,'/behavioral/MCHAT.csv'))




# 2. STRUCTURAL ANALYSES (run after Setup + Shared Data) ----


## 2a. Subject Selection and Matching ----

### Anatomical Data ----


# ----

testData <- read.csv(paste0(baseDir,'/anatomical/combinedASDData.csv'),stringsAsFactors=F)
testPreemie <- subset(testData, PretermStatus=='Preterm' & (Kidokoro <8 | is.na(Kidokoro)))
testFetal <- subset(testData, PretermStatus=='Term' & CHDStatus !='CHD' & AnyFetal==T)

testFetalPass <- subset(testFetal, mchat_result=='Pass')
testFetalFail <- subset(testFetal, mchat_result=='Fail')



fetalPassID <- formatID(testFetalPass$study_id, testFetalPass$study)

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


# ----

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


# ----

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


# ----

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


# ----

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

# Axis label and y-axis settings for one ROI: the hand-tuned values if it is in regionalFigureROIs, otherwise automatic.
# type is 'group' or 'intx'; maxValue is the largest value that must fit on the y axis (used only for automatic limits).
roiAxis <- function(roi, type, maxValue) {
  spec <- regionalFigureROIs[regionalFigureROIs$roi == roi,]
  yMaxCol <- paste0(type, 'YMax'); breakMaxCol <- paste0(type, 'BreakMax'); stepCol <- paste0(type, 'BreakStep')

  if (nrow(spec) == 1) {
    list(label = spec$label, limits = c(0, spec[[yMaxCol]]), breaks = seq(0, spec[[breakMaxCol]], spec[[stepCol]]))
  } else {
    breaks <- pretty(c(0, maxValue * 1.1), n = 6)
    list(label = paste(gsub('_', ' ', roi), 'Volume'), limits = c(0, max(breaks)), breaks = breaks)
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
    scale_y_continuous(limits=axisInfo$limits, breaks=axisInfo$breaks, expand=c(0,0)) +
    scale_fill_manual(values=c('#a3c0e9','#e09999')) +
    geom_bar(stat='identity',aes(fill=`Autism Risk`), color='#000000', linewidth=2) +
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
    scale_y_continuous(limits=axisInfo$limits, breaks=axisInfo$breaks, expand=c(0,0)) +
    scale_x_continuous(limits=c(19.5,40.5), breaks=seq(20,40,5), expand=c(0,0)) +
    scale_shape_manual(values=c(23,21)) +   # Low Risk = diamond, High Risk = circle (as before)
    scale_fill_manual(values=c('#a3c0e955','#e0999955')) +
    geom_point(aes(fill=`Autism Risk`, shape=`Autism Risk`), color='#000000', size=4.5, stroke=.75) +
    scale_color_manual(values=c('#72a3e8','#e37575')) +
    geom_smooth(method='lm', formula='y~x', aes(color=`Autism Risk`), se=F, linewidth=2) +
    ylab(axisInfo$label) + xlab('Gestional Age (Weeks)') +
    themeRegionalFigure(showXTitle = TRUE)
}

# Arrange one plot per ROI in a grid (4 columns wide). Prints a message and draws nothing if no ROI qualifies.
regionalFigure <- function(rois, plotFun, label) {
  if (length(rois) == 0) {
    message('No regions to plot for: ', label)
    return(invisible(NULL))
  }
  message(label, ': ', length(rois), ' region(s): ', paste(rois, collapse=', '))
  plots <- lapply(rois, plotFun)
  names(plots) <- rois
  ggpubr::ggarrange(plotlist = plots, nrow = ceiling(length(plots)/4), ncol = min(4, length(plots)), legend = 'none')
}


## Significant regions

# Group effect and interaction p-values both come from the normalized interaction model table
groupROIs <- lapply(figureThresholds, function(thr) significantROIs(glmRegionalNormedVolTable.Itx, 'Grp.p.Adj', thr))
intxROIs <- lapply(figureThresholds, function(thr) significantROIs(glmRegionalNormedVolTable.Itx, 'Itx.p.Adj', thr))
names(groupROIs) <- names(intxROIs) <- paste0('p', figureThresholds)


## Group main effect figures

for (k in seq_along(figureThresholds)) {
  thrTag <- paste0('p', sprintf('%02d', round(figureThresholds[k]*100)))
  showFigure(regionalFigure(groupROIs[[k]], plotRoiGroup, paste('Group effect, p.adj <', figureThresholds[k])),
             paste0('regional_group_', thrTag), width = 14, height = 7)
}


## Interaction figures

for (k in seq_along(figureThresholds)) {
  thrTag <- paste0('p', sprintf('%02d', round(figureThresholds[k]*100)))
  showFigure(regionalFigure(intxROIs[[k]], plotRoiInteraction, paste('Interaction, p.adj <', figureThresholds[k])),
             paste0('regional_intx_', thrTag), width = 14, height = 7)
}



## 2d. Cortical Thickness ----


# ----


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


# ----

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


# ----


#| echo: false
#| 
allTDFC <- read.csv(paste0(baseDir,'/functional/UseableData/tdCorrSummary.csv'), row.names = 1)
allASDFC <- read.csv(paste0(baseDir,'/functional/UseableData/asdCorrSummary_Expand.csv'), row.names = 1)

rownames(allTDFC) <- substr(rownames(allTDFC),0,nchar(rownames(allTDFC))-28)
#rownames(allASDFC) <- substr(rownames(allASDFC),0,nchar(rownames(allASDFC))-16)

tdFC <- allTDFC[rownames(allTDFC) %in% tdDemoFinal$ID,]
asdFC <- allASDFC[rownames(allASDFC) %in% asdDemoFinal$ID,]

rm(list=c('allTDFC', 'allASDFC'))



## 3b. Whole Brain Approach ----

# First step is to just do a whole brain search

# Here, let's run all 19,900 ROI-ROI connections

# It's more exploratory and FDR will kill us for sure


# ----

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


SigROIList <- subset(fcIntxGLM, Grp.p.adjust <.05)$ROI
ROIList1 <- vapply(strsplit(SigROIList, '_'), function(x) x[1], character(1))
ROIList2 <- vapply(strsplit(SigROIList, '_'), function(x) x[2], character(1))

ROIListFull <- c(ROIList1, ROIList2)
# sort(table(ROIListFull))


write.table(fcIntxGLM, paste0(baseDir,'/functional/output/WholeBrain_WholeTable.csv'), col.names=T, row.names=F, sep=',')
write.table(subset(fcIntxGLM, Grp.p.adjust <.05), paste0(baseDir,'/functional/output/WholeBrain_P05FDR.csv'), col.names=T, row.names=F, sep=',')



### Whole Brain Figures ----


# ----

wholeBrainTable <- read.csv(paste0(baseDir,'/functional/output/WholeBrain_WholeTable.csv'))
wholeBrainTableFDR <- read.csv(paste0(baseDir,'/functional/output/WholeBrain_P05FDR.csv'))

## Connectogram here
# Create a dataframe for plotting
fullFCList <- read.csv(paste0(baseDir,'/functional/100PerHemiLabels_Decimal.1D'), header=F)

NetList <- read.csv(paste0(baseDir,'/functional/r200_Net.csv'))
NetList$LobeColor <- '#000000'
NetList$LobeColor[NetList$Lobe == "Temporal"]  <- '#a9dd8c'
NetList$LobeColor[NetList$Lobe == "Frontal"]  <- '#efa973'
NetList$LobeColor[NetList$Lobe == "Cerebellum"]  <- '#96a5d7'
NetList$LobeColor[NetList$Lobe == "Parietal"]  <- '#d9c97e'
NetList$LobeColor[NetList$Lobe == "Occipital"]  <- '#79f09f'
NetList$LobeColor[NetList$Lobe == "Brainstem"]  <- '#a67de8'
NetList$LobeColor[NetList$Lobe == "Subcortical"] <- '#89d4c1'

NetListNewOrder <- NetList[order(NetList$OrderNew),] 
ROIOrder <- NetListNewOrder$Name_Abrv
ROIColor <- NetListNewOrder$LobeColor


ROIColor <- c(ACG.L = '#efa973',
ACG.L = '#efa973',
ACG.L = '#efa973',
IFGtriang.L = '#efa973',
IFGtriang.L = '#efa973',
IFGtriang.L = '#efa973',
MFG.L = '#efa973',
MFG.L = '#efa973',
MFG.L = '#efa973',
MFG.L = '#efa973',
ORBinf.L = '#efa973',
ORBinf.L = '#efa973',
ORBmid.L = '#efa973',
ORBmid.L = '#efa973',
ORBsupb.L = '#efa973',
ORBsupb.L = '#efa973',
PreCG.L = '#efa973',
PreCG.L = '#efa973',
PreCG.L = '#efa973',
PreCG.L = '#efa973',
PreCG.L = '#efa973',
REC.L = '#efa973',
SFGdor.L = '#efa973',
SFGdor.L = '#efa973',
SFGdor.L = '#efa973',
SFGdor.L = '#efa973',
SFGdor.L = '#efa973',
SFGdor.L = '#efa973',
SFGmed.L = '#efa973',
SFGmed.L = '#efa973',
INS.L = '#d9c97e',
INS.L = '#d9c97e',
INS.L = '#d9c97e',
INS.L = '#d9c97e',
IPL.L = '#d9c97e',
IPL.L = '#d9c97e',
IPL.L = '#d9c97e',
MCG.L = '#d9c97e',
MCG.L = '#d9c97e',
MCG.L = '#d9c97e',
MCG.L = '#d9c97e',
PCG.L = '#d9c97e',
PCUN.L = '#d9c97e',
PCUN.L = '#d9c97e',
PCUN.L = '#d9c97e',
PoCG.L = '#d9c97e',
PoCG.L = '#d9c97e',
PoCG.L = '#d9c97e',
SMG.L = '#d9c97e',
SMG.L = '#d9c97e',
SPG.L = '#d9c97e',
SPG.L = '#d9c97e',
SPG.L = '#d9c97e',
SPG.L = '#d9c97e',
FFG.L = '#a9dd8c',
ITG.L = '#a9dd8c',
ITG.L = '#a9dd8c',
ITG.L = '#a9dd8c',
MTG.L = '#a9dd8c',
MTG.L = '#a9dd8c',
MTG.L = '#a9dd8c',
MTG.L = '#a9dd8c',
MTG.L = '#a9dd8c',
STG.L = '#a9dd8c',
STG.L = '#a9dd8c',
STG.L = '#a9dd8c',
STG.L = '#a9dd8c',
TPOmid.L = '#a9dd8c',
CAL.L = '#79f09f',
CAL.L = '#79f09f',
CUN.L = '#79f09f',
CUN.L = '#79f09f',
CUN.L = '#79f09f',
LING.L = '#79f09f',
LING.L = '#79f09f',
IOG.L = '#79f09f',
IOG.L = '#79f09f',
MOG.L = '#79f09f',
MOG.L = '#79f09f',
MOG.L = '#79f09f',
SOG.L = '#79f09f',
AMY.L = '#89d4c1',
CAU.L = '#89d4c1',
CAU.L = '#89d4c1',
CAU.L = '#89d4c1',
PAL.L = '#89d4c1',
PHG.L = '#89d4c1',
PHG.L = '#89d4c1',
PUT.L = '#89d4c1',
THA.L = '#89d4c1',
THA.L = '#89d4c1',
CRB.L = '#96a5d7',
CRB.L = '#96a5d7',
CRB.L = '#96a5d7',
CRB.L = '#96a5d7',
CRB.L = '#96a5d7',
BS = '#a67de8',
BS = '#a67de8',
BS = '#a67de8',
BS = '#a67de8',
BS = '#a67de8',
BS = '#a67de8',
BS = '#a67de8',
BS = '#a67de8',
CRB.R = '#96a5d7',
CRB.R = '#96a5d7',
CRB.R = '#96a5d7',
CRB.R = '#96a5d7',
CRB.R = '#96a5d7',
THA.R = '#89d4c1',
THA.R = '#89d4c1',
THA.R = '#89d4c1',
PUT.R = '#89d4c1',
PUT.R = '#89d4c1',
PUT.R = '#89d4c1',
PHG.R = '#89d4c1',
PAL.R = '#89d4c1',
HIP.R = '#89d4c1',
AMY.R = '#89d4c1',
SOG.R = '#79f09f',
SOG.R = '#79f09f',
MOG.R = '#79f09f',
MOG.R = '#79f09f',
MOG.R = '#79f09f',
IOG.R = '#79f09f',
IOG.R = '#79f09f',
LING.R = '#79f09f',
LING.R = '#79f09f',
CUN.R = '#79f09f',
CUN.R = '#79f09f',
CAL.R = '#79f09f',
CAL.R = '#79f09f',
TPOsup.R = '#a9dd8c',
TPOmid.R = '#a9dd8c',
STG.R = '#a9dd8c',
STG.R = '#a9dd8c',
STG.R = '#a9dd8c',
MTG.R = '#a9dd8c',
MTG.R = '#a9dd8c',
MTG.R = '#a9dd8c',
MTG.R = '#a9dd8c',
ITG.R = '#a9dd8c',
ITG.R = '#a9dd8c',
ITG.R = '#a9dd8c',
FFG.R = '#a9dd8c',
FFG.R = '#a9dd8c',
SPG.R = '#d9c97e',
SPG.R = '#d9c97e',
SPG.R = '#d9c97e',
SPG.R = '#d9c97e',
SMG.R = '#d9c97e',
PoCG.R = '#d9c97e',
PoCG.R = '#d9c97e',
PoCG.R = '#d9c97e',
PoCG.R = '#d9c97e',
PoCG.R = '#d9c97e',
PoCG.R = '#d9c97e',
PCUN.R = '#d9c97e',
PCUN.R = '#d9c97e',
PCUN.R = '#d9c97e',
PCL.R = '#d9c97e',
PCL.R = '#d9c97e',
MCG.R = '#d9c97e',
MCG.R = '#d9c97e',
MCG.R = '#d9c97e',
IPL.R = '#d9c97e',
IPL.R = '#d9c97e',
IPL.R = '#d9c97e',
IPL.R = '#d9c97e',
INS.R = '#d9c97e',
INS.R = '#d9c97e',
HES.R = '#d9c97e',
SMA.R = '#efa973',
SMA.R = '#efa973',
SFGmed.R = '#efa973',
SFGmed.R = '#efa973',
SFGdor.R = '#efa973',
SFGdor.R = '#efa973',
SFGdor.R = '#efa973',
SFGdor.R = '#efa973',
PreCG.R = '#efa973',
PreCG.R = '#efa973',
PreCG.R = '#efa973',
ORBsupb.R = '#efa973',
ORBsupb.R = '#efa973',
ORBsupb.R = '#efa973',
ORBsupb.R = '#efa973',
ORBmid.R = '#efa973',
ORBmed.R = '#efa973',
ORBinf.R = '#efa973',
ORBinf.R = '#efa973',
MFG.R = '#efa973',
MFG.R = '#efa973',
MFG.R = '#efa973',
MFG.R = '#efa973',
IFGtriang.R = '#efa973',
IFGtriang.R = '#efa973',
IFGoperc.R = '#efa973',
ACG.R = '#efa973',
ACG.R = '#efa973')


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

chordDiagram.Base$ASDCol <- '#00000000'
chordDiagram.Base$ASDCol[chordDiagram.Base$ASDConn==1] <-'#a60000aa'
chordDiagram.Base$ASDCol[!(chordDiagram.Base$ConnFull %in% wholeBrainTableFDR$ROI)] <- '#00000000'

chordDiagram.Base$TDCol <- '#00000000'
chordDiagram.Base$TDCol[chordDiagram.Base$TDConn==1] <- '#0050bfaa'
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


drawASDChord <- function() {
  circos.clear()
  circos.par(start.degree=81, gap.after=gapSizes)
  chordDiagram(chordDiagram.Base[,c('ToAbbrev','FromAbbrev','ASDConn')],
               order = ROIOrder, grid.col=ROIColor, grid.border = '#000000',
               col=chordDiagram.Base$ASDCol, symmetric = TRUE, scale=TRUE, 
               link.border=chordDiagram.Base$ASDCol, link.lwd=3,
               annotationTrack = c("grid"), annotationTrackHeight = c(0.03, 0.03))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(-.25, .5))
  }, bg.border = NA) # here set bg.border to NA is important
}

showBaseFigure('wholebrain_chord_ASD', width = 8, height = 8, drawFun = drawASDChord)


drawTDChord <- function() {
  circos.clear()
  circos.par(start.degree=81, gap.after=gapSizes)
  chordDiagram(chordDiagram.Base[,c('ToAbbrev','FromAbbrev','TDConn')],
               order = ROIOrder, grid.col=ROIColor, grid.border = '#000000',
               col=chordDiagram.Base$TDCol, symmetric = TRUE, scale=TRUE, 
               link.border=chordDiagram.Base$TDCol, link.lwd=3,
               annotationTrack = c("grid"), annotationTrackHeight = c(0.03, 0.03))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(-.25, .5))
  }, bg.border = NA) # here set bg.border to NA is important
}

showBaseFigure('wholebrain_chord_TD', width = 8, height = 8, drawFun = drawTDChord)


asdConnections <- c(chordDiagram.Base[chordDiagram.Base$ASDCol == '#a60000aa',]$ToAbbrev, chordDiagram.Base[chordDiagram.Base$ASDCol == '#a60000aa',]$TFromAbbrev)

tdConnections <- c(chordDiagram.Base[chordDiagram.Base$TDCol == '#0050bfaa',]$ToAbbrev, chordDiagram.Base[chordDiagram.Base$TDCol == '#0050bfaa',]$TFromAbbrev)


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


meltTDConn <- allConnDescSig.TD

meltTDConn$ROWNAMES <- rownames(allConnDescSig.TD)
tdConnPlot <- melt(meltTDConn, id.vars='ROWNAMES')
colnames(tdConnPlot) <- c('Row','Col','Num')

tdConnPlot$Row <- factor(tdConnPlot$Row, levels=lobeNames)
tdConnPlot$Col <- factor(tdConnPlot$Col, levels=rev(lobeNames))

wholebrain_lobeHeatmap_TD <- ggplot(tdConnPlot, aes(x=Row, y=Col)) +
  scale_fill_gradient(low='#ffffff', high='#0050bf') +
  geom_tile(aes(fill=Num), color='#000000') +
  geom_text(aes(label=Num), color='#000000')  +
  theme(axis.title.x = element_text(color = 'black', size = 16),
        axis.text.x = element_text(color = 'black', size = 12),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(color = 'black', size = 16),
        axis.text.y = element_text(color = 'black', size = 12),
        axis.ticks.y = element_blank(),
        axis.line.x = element_blank(),
        plot.title = element_text(color = 'black', size = 24, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill='white'),
        legend.position = 'Right',
        legend.text = element_text(color = 'black', size = 12),
        legend.title=element_blank())
showFigure(wholebrain_lobeHeatmap_TD, 'wholebrain_lobeHeatmap_TD', width = 8.95, height = 8)


meltASDConn <- allConnDescSig.ASD
meltASDConn$ROWNAMES <- rownames(allConnDescSig.ASD)
asdConnPlot <- melt(meltASDConn, id.vars='ROWNAMES')
colnames(asdConnPlot) <- c('Row','Col','Num')

asdConnPlot$Row <- factor(asdConnPlot$Row, levels=lobeNames)
asdConnPlot$Col <- factor(asdConnPlot$Col, levels=rev(lobeNames))

wholebrain_lobeHeatmap_ASD <- ggplot(asdConnPlot, aes(x=Row, y=Col)) +
  scale_fill_gradient(limits=c(0,25), low='#ffffff', high='#a60000') +
  geom_tile(aes(fill=Num), color='#000000') +
  geom_text(aes(label=Num), color='#000000')  +
  theme(axis.title.x = element_text(color = 'black', size = 16),
        axis.text.x = element_text(color = 'black', size = 12),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(color = 'black', size = 16),
        axis.text.y = element_text(color = 'black', size = 12),
        axis.ticks.y = element_blank(),
        axis.line.x = element_blank(),
        plot.title = element_text(color = 'black', size = 24, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill='white'),
        legend.position = 'Right',
        legend.text = element_text(color = 'black', size = 12),
        legend.title=element_blank())
showFigure(wholebrain_lobeHeatmap_ASD, 'wholebrain_lobeHeatmap_ASD', width = 8.95, height = 8)


## 3c. Smaller Network Approach ----

# Here we are subdividing out into smaller networks

# This will (A) help with the crazy FDR correction and (B) more pointedly target questions about connectivity within networks such as the CRB

# Biggest job, however, is to define the networks in a rigorous way


# ----
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


# ----

# ROI number = the third '.'-separated part of an ROI name (e.g. 'CAU.L.39' -> 39)
roiNumber <- function(roi) as.numeric(vapply(strsplit(roi, '.', fixed = TRUE), function(x) x[3], character(1)))

# One entry per network:
#   match         - function(roi names) -> TRUE for regions in the network
#   xLabel        - x-axis label of the bar chart
#   yLimits/yBreaks - y-axis of the bar chart (yLimits NULL = automatic)
#   zeroLine      - draw a line at 0 on the bar chart
networkSpecs <- list(

  CRB = list(match = function(roi) grepl('CRB', roi),
             xLabel = 'Cerebellar Connections', yLimits = c(0,1.2), yBreaks = seq(0,2,.2), zeroLine = TRUE),

  # Originally c(5,9,20,29,56,61,67,68,81,86,94,108,120,133,141,143,155,160,188,196,198,200); parahippocampal gyrus added (43,86,125)
  Limbic = list(match = function(roi) roiNumber(roi) %in% c(5,9,20,29,56,61,67,68,81,86,94,108,120,133,141,143,155,160,188,196,198,200,
                                                            43,86,125),
                xLabel = 'Limbic Connections', yLimits = c(-.4,1.4), yBreaks = seq(-2,2,.2), zeroLine = FALSE),

  # Originally also included 43,86,125 (hippocampus and parahippocampal gyrus); removed
  DMN = list(match = function(roi) roiNumber(roi) %in% c(4,11,18,22,23,34,40,42,48,53,70,85,88,96,100,102,104,106,123,124,128,137,150,152,156,163,168,171,176,178),
             xLabel = 'Default Mode Connections', yLimits = c(-.6,1), yBreaks = seq(-2,2,.2), zeroLine = FALSE),

  SomMot = list(match = function(roi) grepl('PreCG', roi) | grepl('PoCG', roi) | grepl('SMA', roi),
                xLabel = 'Sensorimotor Connections', yLimits = c(-.3,1.1), yBreaks = seq(-2,2,.2), zeroLine = FALSE),

  SubCort = list(match = function(roi) roiNumber(roi) %in% c(31,32,39,51,80,92,97,107,111,134,139,145,151),
                 xLabel = 'Subcortical Connections', yLimits = c(-.2,1.2), yBreaks = seq(-2,2,.2), zeroLine = FALSE),

  Auditory = list(match = function(roi) roiNumber(roi) %in% c(88,94,2,12,48,171,107,139,152,185),
                  xLabel = 'Auditory Connections', yLimits = NULL, yBreaks = NULL, zeroLine = FALSE),

  Visual = list(match = function(roi) roiNumber(roi) %in% c(181,174,153,114,149,195,167,105,116,194,197,162,131,186,135,
                                                            72,60,19,52,35,75,25,15,42,30,91,87,98,57),
                xLabel = 'Visual Connections', yLimits = NULL, yBreaks = NULL, zeroLine = FALSE)
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

# Bar chart of the group means for every connection in a network, ordered by group p-value
plotNetworkBars <- function(fcGLM, spec) {
  longFCGLM <- melt(fcGLM, id.vars=c('ROI'), measure.vars = c('TD.Mean','ASD.Mean'), variable.name='Group', value.name='FC')
  pValOrder <- fcGLM$ROI[order(fcGLM$Grp.p, decreasing = T)]

  barPlot <- ggplot(data=longFCGLM, aes(x=factor(ROI, levels=pValOrder, ordered=is.ordered(T)), y=FC, group=factor(Group)))
  if (!is.null(spec$yLimits)) {
    barPlot <- barPlot + scale_y_continuous(limits=spec$yLimits, breaks=spec$yBreaks, expand=c(0,0))
  }
  barPlot <- barPlot +
    scale_fill_manual(values=c('#d88080cc','#8cb0e3cc')) +
    scale_color_manual(values=c('#bb0000','#0000bb'))
  if (spec$zeroLine) { barPlot <- barPlot + geom_hline(yintercept = 0, linewidth=1.5) }
  barPlot +
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
          panel.grid.major.y = element_line(color = 'dark gray'),
          panel.grid.major.x = element_line(color = 'dark gray'),
          panel.background = element_rect(fill='white'),
          legend.position = 'none',
          legend.text = element_text(color = 'black', size = 12),
          legend.title=element_blank())
}


### Run Network Analyses ----


# ----

networkResults <- list()
for (net in names(networkSpecs)) {
  networkResults[[net]] <- runNetworkAnalysis(net, networkSpecs[[net]])
  showFigure(plotNetworkBars(networkResults[[net]]$fcGLM, networkSpecs[[net]]),
             paste0('network_', net, '_bars'), width = 10, height = 6)
}

# Effect size conversion and group-difference summary for the cerebellar network
effectsize::eta2_to_f2(0.2953587)
effectsize::eta2_to_f(0.2953587)

deltaSS <- mean(networkResults$CRB$fcGLM$TD.Mean) - mean(networkResults$CRB$fcGLM$ASD.Mean)
sigmaSS <- sd(networkResults$CRB$fcGLM$TD.Mean) - sd(networkResults$CRB$fcGLM$ASD.Mean)



#### Network Connectivity Figures ----


# ----

# Read in the output of the above analyses for plotting
FCTable.AllConn.DMN <- read.csv(paste0(baseDir,'/functional/output/DMN_Intx_Wholebrain.csv'), stringsAsFactors=F)
FCTable.AllConn.Limbic <- read.csv(paste0(baseDir,'/functional/output/Limbic_Intx_Wholebrain.csv'), stringsAsFactors=F)
FCTable.AllConn.Subcort <- read.csv(paste0(baseDir,'/functional/output/SubCort_Intx_Wholebrain.csv'), stringsAsFactors=F)
FCTable.AllConn.CRB <- read.csv(paste0(baseDir,'/functional/output/CRB_Intx_Wholebrain.csv'), stringsAsFactors=F)

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
  #scale_y_continuous(limits=c(-.2,1.2), breaks=seq(-2,2,.2), expand=c(0,0)) +
  scale_fill_manual(values=c('#d88080cc','#8cb0e3cc')) +
  scale_color_manual(values=c('#bb0000','#0000bb')) +
  geom_hline(yintercept = 0, lwd=1) +
  geom_bar(stat='identity', position=position_dodge2(0.2), aes(color=Grp,fill=Grp)) +
  # geom_errorbar(position=position_dodge2(0.2), aes(ymin=Mean-SE, ymax=Mean+SE), size=.2) +
  ylab('Connectivity Strength') + xlab('Subcortical Connections') +
  theme(axis.title.x = element_text(color = 'black', size = 16),
        axis.text.x = element_text(color = 'black', size = 4, angle=90),
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
        legend.position = 'right',
        legend.text = element_text(color = 'black', size = 12),
        legend.title=element_blank())

showFigure(network_combined_bars_p10, 'network_combined_bars_p10', width = 10, height = 6)




# 4. Anatomical-Functional Overlap (optional; needs both 2a and 3a) ----


# ----

fetalMatched

asdDemoFinal$ID[!c(asdDemoFinal$ID %in% fetalMatched$longID)]

subset(fetalMatched, mchat_result=='Fail')$longID[!c(subset(fetalMatched, mchat_result=='Fail')$longID %in% asdDemoFinal$ID)]

