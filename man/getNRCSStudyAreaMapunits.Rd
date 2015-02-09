\name{getNRCSStudyAreaMapunits}
\alias{getNRCSStudyAreaMapunits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getNRCSStudyAreaMapunits(area, date, raw.dir)
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
    cat("(Down)Loading soils polygons for", area, "\n")
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
    mapunits <- rgdal::readOGR(paste(raw.dir, "/", area, "/spatial", 
        sep = ""), layer = paste("soilmu_a_", tolower(area), 
        sep = ""), verbose = F)
    unlink(paste(raw.dir, "/", area, sep = ""), recursive = T, 
        force = T)
    return(mapunits)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
