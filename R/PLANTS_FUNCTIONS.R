## All data is downloaded from the NRCS,
## available at http://websoilsurvey.sc.egov.usda.gov

## Author: R. Kyle Bocinsky
## Date: 02/14/2014

loadPLANTS <- function(label, raw.dir, extraction.dir=NULL, force.redo=F){  
  if(is.null(extraction.dir)){
    extraction.dir <- paste(raw.dir,"/EXTRACTIONS",sep='')
  }
  dir.create(paste(extraction.dir,"/",label,"/",sep=''), showWarnings = FALSE, recursive=T)
  dsn.vectors <- paste(extraction.dir,"/",label,"/vectors",sep='')
  dir.create(dsn.vectors, showWarnings = FALSE, recursive = TRUE)
  
  tables.dir <- paste(extraction.dir,"/",label,"/tables",sep='')
  dir.create(tables.dir, showWarnings = FALSE, recursive = TRUE)
  
  if(file.exists(paste(tables.dir,"/all_plants.csv",sep='')) & !force.redo){
    plants.all <- read.csv(paste(tables.dir,"/all_plants.csv",sep=''))
    return(plants.all)
  }
  
  
  # Load the NRCS PLANTS database.
  # This search command returns a CSV of all plants in the PLANTS database that are native to the lower 48
  # US States and their "Growth Habit" (forb, tree, shrub, etc.). These data will be used to weed out (ha!)
  # invasive species from the Soil Survey plant data, and to build rasters of plant communities by Growth Habit.
  #   cat("Loading NRCS PLANTS database.\n\n")
  if(!file.exists(paste(tables.dir,"/L48_native_plants.csv",sep=''))){
    L48_native_plants <- read.csv(url("http://plants.usda.gov/java/AdvancedSearchServlet?pfa=l48&nativestatuscode=nus48&dsp_grwhabt=on&dsp_vernacular=on&Synonyms=all&viewby=sciname&download=on"))
    write.csv(L48_native_plants, paste(tables.dir,"/L48_native_plants.csv",sep=''), row.names=FALSE)
  }else{
    L48_native_plants <- read.csv(paste(tables.dir,"/L48_native_plants.csv",sep=''))
  }
  
  # Load the NRCS PLANTS database of unknown plant symbols.
  # This downloads the text file, assigns each symbol a "Growth Habit" code, and appends the unknown plants 
  # to the lower 48 native plants data.
  #   cat("Loading NRCS PLANTS database unknown plant codes.\n\n")
  if(!file.exists(paste(tables.dir,"/unknown_plants.csv", sep=''))){
    unknown_plants <- read.csv(url("http://plants.usda.gov/Data/unknown_plants.txt"))
    write.csv(unknown_plants, paste(tables.dir,"/unknown_plants.csv", sep=''), row.names=FALSE)
  }else{
    unknown_plants <- read.csv(paste(tables.dir,"/unknown_plants.csv", sep=''))
  }
  # Scan each row, and assign "Growth Habit" accordingly.
  unknown_plants$Growth.Habit <- NA
  unknown_plants$Synonym.Symbol <- NA
  unknown_plants$Scientific.Name <- NA
  for(i in 1:nrow(unknown_plants)){
    if(!is.na(charmatch("Forb",unknown_plants$Common.Name[i]))) unknown_plants$Growth.Habit[i] <- "Forb/herb"
    else if(!is.na(charmatch("Tree",unknown_plants$Common.Name[i]))) unknown_plants$Growth.Habit[i] <- "Tree"
    else if(!is.na(charmatch("Shrub",unknown_plants$Common.Name[i]))) unknown_plants$Growth.Habit[i] <- "Shrub"
    else if(!is.na(charmatch("Subshrub",unknown_plants$Common.Name[i]))) unknown_plants$Growth.Habit[i] <- "Shrub"
    else if(!is.na(charmatch("Grass",unknown_plants$Common.Name[i]))) unknown_plants$Growth.Habit[i] <- "Grass"
    else unknown_plants$Growth.Habit[i] <- "Other"
  }
  unknown_plants <- unknown_plants[,c(1,4,5,2,3)]
  names(unknown_plants) <- c("Accepted.Symbol", "Synonym.Symbol", "Scientific.Name", "Common.Name", "Growth.Habit")
  # Append the unknown plants to the lower 48 plants
  plants.all <- rbind(L48_native_plants,unknown_plants)
  # A function that takes only the first (primary) Growth Habit of the plant
  select.first <- function(x) strsplit(x,",")[[1]][1]
  # Select only first (primary) Growth Habit of the plant
  plants.all$Growth.Habit <- as.vector(sapply(as.character(plants.all$Growth.Habit),FUN=select.first))
  # A function that renames Growth Habits
  simplify <- function(x){
    if(!is.na(charmatch("Forb",x))) y <- "Forb"
    else if(!is.na(charmatch("Tree",x))) y <- "Tree"
    else if(!is.na(charmatch("Shrub",x))) y <- "Shrub"
    else if(!is.na(charmatch("Subshrub",x))) y <- "Shrub"
    else if(!is.na(charmatch("Grass",x))) y <- "Grass"
    else if(!is.na(charmatch("Graminoid",x))) y <- "Grass"
    else y <- "Other"
    return(y)
  }
  # Simplify the Growth Habits
  plants.all$Growth.Habit <- as.vector(sapply(as.character(plants.all$Growth.Habit),FUN=simplify))
  
  write.csv(plants.all, paste(tables.dir,"/all_plants.csv", sep=''), row.names=FALSE)
  
  return(plants.all)
  
}
