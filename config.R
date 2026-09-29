# Project paths. Sourced by every script; contains no data.
# Add a new computer by adding its project folder to baseDirs.
baseDirs <- c(
  'C:/Users/Kevin/OneDrive - Childrens National Hospital/ComputerSyncing/Projects/Active Projects/WashUActigraphy',
  'C:/Users/kmcook/OneDrive - Childrens National Hospital/ComputerSyncing/Projects/Active Projects/WashUActigraphy'
)
baseDir <- baseDirs[dir.exists(baseDirs)][1]
if (is.na(baseDir)) stop('None of the project folders in baseDirs exist on this computer. Add this computer\'s folder to baseDirs.', call. = FALSE)
dataDir <- paste0(baseDir,'/data')
figureDir <- paste0(dataDir,'/figures')

# Input files (inside dataDir). Edit the filename to match yours.
actigraphyFile <- paste0(dataDir,'/actigraphy.csv')

# Any count below this is shown as "<minCell" in printed diagnostics,
# so output is safe to paste back into the chat.
minCell <- 5
