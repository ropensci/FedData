\name{getNRCSStudyAreas}
\alias{getNRCSStudyAreas}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getNRCSStudyAreas(template = NULL, raw.dir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{template}{
%%     ~~Describe \code{template} here~~
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
function (template = NULL, raw.dir) 
{
    url <- "http://websoilsurvey.sc.egov.usda.gov/DataAvailability/SoilDataAvailabilityShapefile.zip"
    destdir <- raw.dir
    wgetDownload(url = url, destdir = destdir)
    cat("Unzipping the NRCS study areas.\n")
    unzip(paste(raw.dir, "/SoilDataAvailabilityShapefile.zip", 
        sep = ""), exdir = paste(raw.dir, "/SoilDataAvailabilityShapefile", 
        sep = ""))
    cat("Loading the NRCS study areas.\n")
    NRCSAreas <- rgdal::readOGR(normalizePath(paste(raw.dir, 
        "/SoilDataAvailabilityShapefile/", sep = "")), layer = "soilsa_a_nrcs", 
        verbose = FALSE)
    if (is.null(template)) {
        return(NRCSAreas)
    }
    NRCSAreas <- raster::crop(NRCSAreas, sp::spTransform(template, 
        sp::CRS(projection(NRCSAreas))))
    unlink(paste(raw.dir, "/SoilDataAvailabilityShapefile", sep = ""), 
        recursive = TRUE)
    if (0 \%in\% NRCSAreas@data$iscomplete) {
        cat("WARNING! Some of the soil surveys in your area are unavailable.\n")
        cat("Soils and productivity data will have holes.\n")
        cat("Missing areas:\n")
        cat(as.vector(NRCSAreas@data[NRCSAreas@data$iscomplete == 
            0, ]$areasymbol))
        cat("\n\n")
        cat("Continuing with processing available soils.\n\n")
    }
    return(NRCSAreas)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
