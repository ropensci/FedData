\name{getNHD}
\alias{getNHD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getNHD(template, label, raw.dir = "./RAW/NHD/", extraction.dir = "./EXTRACTIONS/NHD/", force.redo = FALSE)
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
function (template, label, raw.dir = "./RAW/NHD/", extraction.dir = "./EXTRACTIONS/NHD/", 
    force.redo = FALSE) 
{
    vectors.dir <- paste(extraction.dir, "/", label, "/spatial", 
        sep = "")
    dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(extraction.dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(vectors.dir, showWarnings = FALSE, recursive = TRUE)
    if (!force.redo & length(list.files(vectors.dir)) > 0) {
        files <- list.files(vectors.dir)
        files <- files[grepl("shp", files)]
        files <- files[!grepl("template", files)]
        files <- gsub(".shp", "", files)
        files <- files[order(files)]
        shapes <- lapply(files, function(file) {
            rgdal::readOGR(normalizePath(vectors.dir), file, 
                verbose = F)
        })
        names(shapes) <- files
        return(shapes)
    }
    if (class(template) \%in\% c("RasterLayer", "RasterStack", 
        "RasterBrick")) {
        template <- SPDFfromPolygon(sp::spTransform(polygonFromExtent(template), 
            sp::CRS("+proj=longlat +ellps=GRS80")))
    }
    HUC4 <- getHUC4(template = template, raw.dir = raw.dir, force.redo = force.redo)
    area.list <- getAreaList(HUC4)
    shapes <- loadNHDSubregions(template = template, area.list = area.list, 
        raw.dir = raw.dir, vectors.dir = vectors.dir, force.redo = force.redo)
    return(shapes)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
