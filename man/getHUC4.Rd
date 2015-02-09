\name{getHUC4}
\alias{getHUC4}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getHUC4(template, raw.dir, force.redo = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{template}{
%%     ~~Describe \code{template} here~~
}
  \item{raw.dir}{
%%     ~~Describe \code{raw.dir} here~~
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
function (template, raw.dir, force.redo = F) 
{
    url <- "ftp://ftp.igsb.uiowa.edu/gis_library/USA/huc_04.zip"
    destdir <- raw.dir
    wgetDownload(url = url, destdir = destdir)
    cat("Unzipping the NHD HUC4 dataset.\n")
    unzip(paste(raw.dir, "/huc_04.zip", sep = ""), exdir = paste(raw.dir, 
        "/huc_04", sep = ""))
    cat("Loading the NHD HUC4 dataset.\n")
    HUC4 <- rgdal::readOGR(normalizePath(paste(raw.dir, "/huc_04/", 
        sep = "")), layer = "huc_04", verbose = FALSE)
    HUC4@proj4string <- sp::CRS("+proj=utm +zone=15 +datum=NAD83 +ellps=WGS84")
    HUC4 <- raster::crop(HUC4, spTransform(template, CRS(projection(HUC4))))
    unlink(paste(raw.dir, "/huc_04", sep = ""), recursive = TRUE)
    return(HUC4)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
