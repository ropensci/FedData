\name{getITRDB}
\alias{getITRDB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getITRDB(raw.dir, output.dir, type = "standard", download = F, force.redo = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{raw.dir}{
%%     ~~Describe \code{raw.dir} here~~
}
  \item{output.dir}{
%%     ~~Describe \code{output.dir} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{download}{
%%     ~~Describe \code{download} here~~
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
function (raw.dir, output.dir, type = "standard", download = F, 
    force.redo = F) 
{
    raw.dir <- paste(raw.dir, type, "/", sep = "")
    dir.create(raw.dir, recursive = T, showWarnings = F)
    dir.create(output.dir, recursive = T, showWarnings = F)
    if (download) {
        if (type == "standard") {
            system(paste("wget -ignore-case -np -nH -nd -nc -r -A \"*[0-9][0-9].crn\" --directory-prefix=", 
                raw.dir, " ftp://ftp.ncdc.noaa.gov/pub/data/paleo/treering/chronologies/", 
                sep = ""))
        }
        else if (type == "whitened") {
            system(paste("wget -ignore-case -np -nH -nd -nc -r -A \"*[0-9][0-9]r.crn\" --directory-prefix=", 
                raw.dir, " ftp://ftp.ncdc.noaa.gov/pub/data/paleo/treering/chronologies/", 
                sep = ""))
        }
        else if (type == "arstan") {
            system(paste("wget -ignore-case -np -nH -nd -nc -r -A \"*[0-9][0-9]a.crn\" --directory-prefix=", 
                raw.dir, " ftp://ftp.ncdc.noaa.gov/pub/data/paleo/treering/chronologies/", 
                sep = ""))
        }
        else {
            stop("ERROR! Unknown chronology type. Please correct and retry.\n")
        }
    }
    dendros <- paste(raw.dir, list.files(raw.dir), sep = "")
    if (file.exists(paste(output.dir, "/ITRDB_DATA.csv", sep = "")) & 
        !force.redo) {
        existing.data <- read.csv(paste(output.dir, "/ITRDB_DATA.csv", 
            sep = ""))
        if ((ncol(existing.data) - 1) == length(dendros)) {
            return(existing.data)
        }
    }
    series.list <- lapply(dendros, function(x, ...) {
        seriesToMatrix(x, raw.dir)
    })
    series.list.error <- dendros[which(grepl("ERROR", series.list))]
    if (length(series.list.error) > 0) {
        cat("ERROR! One or more improperly formatted chronologies. Please correct and retry.\n")
        print(series.list.error)
        stop()
    }
    meta <- as.data.frame(do.call("rbind", lapply(series.list, 
        FUN = function(x) {
            x$meta
        })))
    names(meta) <- c("SERIES", "NAME", "SPECIES", "ELEVATION", 
        "LAT", "LON", "START", "END", "CONTRIBUTOR")
    meta$ELEVATION <- gsub("M", "", meta$ELEVATION)
    meta$ELEVATION <- gsub("m", "", meta$ELEVATION)
    meta$CONTRIBUTOR <- gsub("  ", " ", meta$CONTRIBUTOR)
    SHemi <- grepl("-", meta$LAT)
    WHemi <- grepl("-", meta$LON)
    meta$LAT <- gsub("+", "", meta$LAT, fixed = TRUE)
    meta$LON <- gsub("-", "", meta$LON, fixed = TRUE)
    meta$CONTRIBUTOR <- toupper(meta$CONTRIBUTOR)
    meta$NAME <- toupper(meta$NAME)
    meta$SERIES <- toupper(meta$SERIES)
    meta$LAT <- format(as.numeric(meta$LAT)/100, nsmall = 2)
    meta$LON <- format(as.numeric(meta$LON)/100, nsmall = 2)
    locations <- meta[c("LAT", "LON")]
    locations$LAT <- within(data.frame(lats = locations$LAT), 
        {
            dms <- do.call(rbind, strsplit(as.character(lats), 
                "\\."))
            LAT <- as.numeric(dms[, 1]) + (as.numeric(dms[, 2]))/60
            rm(dms)
            rm(lats)
        })
    locations$LON <- within(data.frame(lons = locations$LON), 
        {
            dms <- do.call(rbind, strsplit(as.character(lons), 
                "\\."))
            LON <- as.numeric(dms[, 1]) + (as.numeric(dms[, 2]))/60
            rm(dms)
            rm(lons)
        })
    meta$LAT <- locations$LAT
    meta$LON <- locations$LON
    if (any(SHemi)) {
        meta[SHemi, ]$LAT <- 0 - meta[SHemi, ]$LAT
    }
    if (any(WHemi)) {
        meta[WHemi, ]$LON <- 0 - meta[WHemi, ]$LON
    }
    meta <- as.data.frame(lapply(meta, function(x) if (is.character(x) | 
        is.factor(x)) 
        gsub(", ", "-", x)
    else x))
    data <- lapply(series.list, FUN = function(x) {
        data.table::data.table(x$data, key = names(x$data))
    })
    series.merged <- as.data.frame(Reduce(function(...) {
        merge(..., by = "YEAR", all = T, sort = F)
    }, data))
    series.merged$YEAR <- as.numeric(series.merged$YEAR)
    series.merged <- series.merged[order(series.merged$YEAR), 
        ]
    YEAR <- series.merged$YEAR
    series.merged <- series.merged[, c(-1)]
    non.nas <- as.numeric(sapply(series.merged, FUN = function(x) {
        match(FALSE, is.na(x))
    }))
    series.merged <- series.merged[, order(non.nas)]
    names(series.merged) <- toupper(names(series.merged))
    series.merged <- cbind(YEAR, series.merged)
    meta$SERIES <- as.character(meta$SERIES)
    ordering <- unlist(lapply(2:length(names(series.merged)), 
        FUN = function(x, ...) {
            which(meta$SERIES == names(series.merged)[x])
        }))
    meta <- meta[ordering, ]
    write.csv(meta, paste(output.dir, "/ITRDB_METADATA.csv", 
        sep = ""), row.names = F)
    write.csv(series.merged, paste(output.dir, "/ITRDB_DATA.csv", 
        sep = ""), row.names = F)
    return(series.merged)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
