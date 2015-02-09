\name{getNRCS}
\alias{getNRCS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getNRCS(template, label, raw.dir = "./RAW/NRCS/", extraction.dir = "./EXTRACTIONS/NRCS/", force.redo = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{template}{
%%     ~~Describe \code{template} here~~
}
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
function (template, label, raw.dir = "./RAW/NRCS/", extraction.dir = "./EXTRACTIONS/NRCS/", 
    force.redo = FALSE) 
{
    vectors.dir <- paste(extraction.dir, "/", label, "/spatial", 
        sep = "")
    tables.dir <- paste(extraction.dir, "/", label, "/tabular", 
        sep = "")
    dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(extraction.dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(vectors.dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(tables.dir, showWarnings = FALSE, recursive = TRUE)
    if (!force.redo & length(list.files(vectors.dir)) > 0 & length(list.files(tables.dir)) > 
        0 & file.exists(paste(vectors.dir, "/NRCSMapunits.shp", 
        sep = ""))) {
        NRCSMapunits <- rgdal::readOGR(normalizePath(vectors.dir), 
            "NRCSMapunits", verbose = F)
        files <- list.files(tables.dir)
        files <- files[grepl("csv", files)]
        files <- files[order(files)]
        tables <- lapply(files, function(file) {
            read.csv(paste(normalizePath(tables.dir), "/", file, 
                sep = ""))
        })
        names(tables) <- files
        return(list(spatial = NRCSMapunits, tabular = tables))
    }
    if (class(template) \%in\% c("RasterLayer", "RasterStack", 
        "RasterBrick")) {
        template <- SPDFfromPolygon(sp::spTransform(polygonFromExtent(template), 
            sp::CRS("+proj=longlat +ellps=GRS80")))
    }
    NRCSAreas <- getNRCSStudyAreas(template = template, raw.dir = raw.dir)
    NRCSAreas <- NRCSAreas[NRCSAreas@data$iscomplete != 0, ]
    NRCSMapunits <- getNRCSMapunits(template = template, areas = NRCSAreas, 
        raw.dir = raw.dir)
    suppressWarnings(rgdal::writeOGR(NRCSMapunits, vectors.dir, 
        "NRCSMapunits", "ESRI Shapefile", overwrite_layer = TRUE))
    NRCSData <- getNRCSData(areas = NRCSAreas, raw.dir = raw.dir)
    NRCSData <- extractNRCSData(tables = NRCSData, mapunits = NRCSMapunits)
    junk <- lapply(names(NRCSData), function(tab) {
        write.csv(NRCSData[[tab]], file = paste(tables.dir, "/", 
            tab, ".csv", sep = ""), row.names = F)
    })
    return(list(spatial = NRCSMapunits, tabular = NRCSData))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
