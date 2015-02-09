\name{getNRCSStudyAreaData}
\alias{getNRCSStudyAreaData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getNRCSStudyAreaData(area, date, raw.dir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{area}{
%%     ~~Describe \code{area} here~~
}
  \item{date}{
%%     ~~Describe \code{date} here~~
}
  \item{raw.dir}{
%%     ~~Describe \code{raw.dir} here~~
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
function (area, date, raw.dir) 
{
    cat("(Down)Loading soils data for", area, "\n")
    state <- substring(area, 1, 2)
    url <- paste("http://websoilsurvey.sc.egov.usda.gov/DSD/Download/Cache/SSA/wss_SSA_", 
        area, "_soildb_", state, "_2003_[", date, "].zip", sep = "")
    destdir <- raw.dir
    tryCatch(wgetDownload(url = url, destdir = destdir, nc = T), 
        warning = function(w) {
            url <- paste("http://websoilsurvey.sc.egov.usda.gov/DSD/Download/Cache/SSA/wss_SSA_", 
                area, "_soildb_US_2003_[", date, "].zip", sep = "")
            wgetDownload(url = url, destdir = destdir, nc = T)
            state <<- "US"
        })
    unzip(paste(raw.dir, "/wss_SSA_", area, "_soildb_", state, 
        "_2003_[", date, "].zip", sep = ""), exdir = raw.dir)
    files <- list.files(paste(raw.dir, area, "/tabular/", sep = ""))
    tablesData <- lapply(files, function(file) {
        tryCatch(return(read.delim(paste(raw.dir, area, "/tabular/", 
            file, sep = ""), header = F, sep = "|")), error = function(e) {
            return(NULL)
        })
    })
    names(tablesData) <- files
    tablesData <- tablesData[!sapply(tablesData, is.null)]
    tablesHeaders <- Hmisc::mdb.get(paste(raw.dir, area, "/soildb_", 
        state, "_2003.mdb", sep = ""))
    NRCSTableMapping <- tablesData[["mstab.txt"]][, c(1, 5)]
    names(NRCSTableMapping) <- c("TABLE", "FILE")
    NRCSTableMapping[, "FILE"] <- paste(NRCSTableMapping[, "FILE"], 
        ".txt", sep = "")
    tablesData <- tablesData[as.character(NRCSTableMapping[, 
        "FILE"])]
    tablesHeaders <- tablesHeaders[as.character(NRCSTableMapping[, 
        "TABLE"])]
    notNull <- (!sapply(tablesData, is.null) & !sapply(tablesHeaders, 
        is.null))
    tablesData <- tablesData[notNull]
    tablesHeaders <- tablesHeaders[notNull]
    tables <- mapply(tablesData, tablesHeaders, FUN = function(theData, 
        theHeader) {
        names(theData) <- names(theHeader)
        return(theData)
    })
    names(tables) <- names(tablesHeaders)
    unlink(paste(raw.dir, "/", area, sep = ""), recursive = T, 
        force = T)
    return(tables)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
