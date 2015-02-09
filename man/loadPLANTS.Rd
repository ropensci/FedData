\name{loadPLANTS}
\alias{loadPLANTS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
loadPLANTS(label, raw.dir, extraction.dir = NULL, force.redo = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{label}{
%%     ~~Describe \code{label} here~~
}
  \item{raw.dir}{
%%     ~~Describe \code{raw.dir} here~~
}
  \item{extraction.dir}{
%%     ~~Describe \code{extraction.dir} here~~
}
  \item{force.redo}{
%%     ~~Describe \code{force.redo} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (label, raw.dir, extraction.dir = NULL, force.redo = F) 
{
    if (is.null(extraction.dir)) {
        extraction.dir <- paste(raw.dir, "/EXTRACTIONS", sep = "")
    }
    dir.create(paste(extraction.dir, "/", label, "/", sep = ""), 
        showWarnings = FALSE, recursive = T)
    dsn.vectors <- paste(extraction.dir, "/", label, "/vectors", 
        sep = "")
    dir.create(dsn.vectors, showWarnings = FALSE, recursive = TRUE)
    tables.dir <- paste(extraction.dir, "/", label, "/tables", 
        sep = "")
    dir.create(tables.dir, showWarnings = FALSE, recursive = TRUE)
    if (file.exists(paste(tables.dir, "/all_plants.csv", sep = "")) & 
        !force.redo) {
        plants.all <- read.csv(paste(tables.dir, "/all_plants.csv", 
            sep = ""))
        return(plants.all)
    }
    if (!file.exists(paste(tables.dir, "/L48_native_plants.csv", 
        sep = ""))) {
        L48_native_plants <- read.csv(url("http://plants.usda.gov/java/AdvancedSearchServlet?pfa=l48&nativestatuscode=nus48&dsp_grwhabt=on&dsp_vernacular=on&Synonyms=all&viewby=sciname&download=on"))
        write.csv(L48_native_plants, paste(tables.dir, "/L48_native_plants.csv", 
            sep = ""), row.names = FALSE)
    }
    else {
        L48_native_plants <- read.csv(paste(tables.dir, "/L48_native_plants.csv", 
            sep = ""))
    }
    if (!file.exists(paste(tables.dir, "/unknown_plants.csv", 
        sep = ""))) {
        unknown_plants <- read.csv(url("http://plants.usda.gov/Data/unknown_plants.txt"))
        write.csv(unknown_plants, paste(tables.dir, "/unknown_plants.csv", 
            sep = ""), row.names = FALSE)
    }
    else {
        unknown_plants <- read.csv(paste(tables.dir, "/unknown_plants.csv", 
            sep = ""))
    }
    unknown_plants$Growth.Habit <- NA
    unknown_plants$Synonym.Symbol <- NA
    unknown_plants$Scientific.Name <- NA
    for (i in 1:nrow(unknown_plants)) {
        if (!is.na(charmatch("Forb", unknown_plants$Common.Name[i]))) 
            unknown_plants$Growth.Habit[i] <- "Forb/herb"
        else if (!is.na(charmatch("Tree", unknown_plants$Common.Name[i]))) 
            unknown_plants$Growth.Habit[i] <- "Tree"
        else if (!is.na(charmatch("Shrub", unknown_plants$Common.Name[i]))) 
            unknown_plants$Growth.Habit[i] <- "Shrub"
        else if (!is.na(charmatch("Subshrub", unknown_plants$Common.Name[i]))) 
            unknown_plants$Growth.Habit[i] <- "Shrub"
        else if (!is.na(charmatch("Grass", unknown_plants$Common.Name[i]))) 
            unknown_plants$Growth.Habit[i] <- "Grass"
        else unknown_plants$Growth.Habit[i] <- "Other"
    }
    unknown_plants <- unknown_plants[, c(1, 4, 5, 2, 3)]
    names(unknown_plants) <- c("Accepted.Symbol", "Synonym.Symbol", 
        "Scientific.Name", "Common.Name", "Growth.Habit")
    plants.all <- rbind(L48_native_plants, unknown_plants)
    select.first <- function(x) strsplit(x, ",")[[1]][1]
    plants.all$Growth.Habit <- as.vector(sapply(as.character(plants.all$Growth.Habit), 
        FUN = select.first))
    simplify <- function(x) {
        if (!is.na(charmatch("Forb", x))) 
            y <- "Forb"
        else if (!is.na(charmatch("Tree", x))) 
            y <- "Tree"
        else if (!is.na(charmatch("Shrub", x))) 
            y <- "Shrub"
        else if (!is.na(charmatch("Subshrub", x))) 
            y <- "Shrub"
        else if (!is.na(charmatch("Grass", x))) 
            y <- "Grass"
        else if (!is.na(charmatch("Graminoid", x))) 
            y <- "Grass"
        else y <- "Other"
        return(y)
    }
    plants.all$Growth.Habit <- as.vector(sapply(as.character(plants.all$Growth.Habit), 
        FUN = simplify))
    write.csv(plants.all, paste(tables.dir, "/all_plants.csv", 
        sep = ""), row.names = FALSE)
    return(plants.all)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
