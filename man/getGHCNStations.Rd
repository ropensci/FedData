\name{getGHCNStations}
\alias{getGHCNStations}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getGHCNStations(template = NULL, elements = NULL, standardize = F, raw.dir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{template}{
%%     ~~Describe \code{template} here~~
}
  \item{elements}{
%%     ~~Describe \code{elements} here~~
}
  \item{standardize}{
%%     ~~Describe \code{standardize} here~~
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
function (template = NULL, elements = NULL, standardize = F, 
    raw.dir) 
{
    if (!is.null(template) & !is(template, "SpatialPolygonsDataFrame")) {
        template <- polygonFromExtent(template)
    }
    url <- "ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/ghcnd-inventory.txt"
    destdir <- raw.dir
    wgetDownload(url = url, destdir = destdir)
    system(paste("sed -i -E 's/#/ /' ", paste(raw.dir, "ghcnd-inventory.txt", 
        sep = ""), sep = ""))
    system(paste("rm ", paste(raw.dir, "ghcnd-inventory.txt-E", 
        sep = ""), sep = ""))
    station.inventory <- utils::read.fwf(paste(raw.dir, "ghcnd-inventory.txt", 
        sep = ""), c(11, 1, 8, 1, 9, 1, 4, 1, 4, 1, 4))[, seq(1, 
        11, 2)]
    names(station.inventory) <- c("ID", "LATITUDE", "LONGITUDE", 
        "ELEMENT", "YEAR_START", "YEAR_END")
    stations.sp <- SpatialPointsDataFrame(coords = station.inventory[, 
        c("LONGITUDE", "LATITUDE")], station.inventory, proj4string = sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
    if (!is.null(template)) {
        stations.sp <- stations.sp[!is.na(sp::over(stations.sp, 
            sp::spTransform(template, sp::CRS(projection(stations.sp))))[, 
            1]), ]
    }
    return(stations.sp)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
