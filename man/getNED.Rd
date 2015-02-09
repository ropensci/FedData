\name{getNED}
\alias{getNED}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getNED(template, label, res, raw.dir = "./RAW/NED/", extraction.dir = "./EXTRACTIONS/NED/", force.redo = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{template}{
%%     ~~Describe \code{template} here~~
}
  \item{label}{
%%     ~~Describe \code{label} here~~
}
  \item{res}{
%%     ~~Describe \code{res} here~~
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
function (template, label, res, raw.dir = "./RAW/NED/", extraction.dir = "./EXTRACTIONS/NED/", 
    force.redo = F) 
{
    rasters.dir <- paste(extraction.dir, "/", label, "/rasters", 
        sep = "")
    dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(rasters.dir, showWarnings = FALSE, recursive = TRUE)
    if (file.exists(paste(rasters.dir, "/NED_", res, ".tif", 
        sep = "")) & !force.redo) {
        extracted.DEM <- raster::raster(paste(rasters.dir, "/NED_", 
            res, ".tif", sep = ""))
        return(extracted.DEM)
    }
    extent.latlon <- raster::extent(raster::projectExtent(template, 
        sp::CRS("+proj=longlat +ellps=WGS84")))
    if (is.null(res) & (class(template) \%in\% c("RasterLayer", 
        "RasterStack", "RasterBrick"))) {
        y.dim <- extent.latlon@ymax - extent.latlon@ymin
        y.res <- y.dim/nrow(template)
        if (y.res < (1/60)/60) {
            res <- "13"
        }
        else {
            res <- "1"
        }
    }
    else if (is.null(res)) {
        warning("If template is a sp* or extent object, res should be provided! \n Defaulting to 1 arc-second NED.\n", 
            immediate. = T)
        res <- "1"
    }
    wests <- c(ceiling(abs(extent.latlon@xmax)), ceiling(abs(extent.latlon@xmin)))
    wests <- unique(append(wests, seq(wests[1], wests[2])))
    norths <- c(ceiling(abs(extent.latlon@ymin)), ceiling(abs(extent.latlon@ymax)))
    norths <- unique(append(norths, seq(norths[1], norths[2])))
    wests <- formatC(wests, width = 3, format = "d", flag = "0")
    norths <- formatC(norths, width = 2, format = "d", flag = "0")
    tiles <- vector("list", length(wests) * length(norths))
    cat("\nArea of interest includes", length(wests) * length(norths), 
        "NED tiles.")
    cat("\nChecking availability of 1x1 degree tiles from the USGS for the study area.\n")
    dir.create(paste(raw.dir, "/", res, sep = ""), showWarnings = FALSE, 
        recursive = TRUE)
    t <- 1
    for (w in wests) {
        for (n in norths) {
            url <- paste("ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/NED/", 
                res, "/ArcGrid/n", n, "w", w, ".zip", sep = "")
            destdir <- paste(raw.dir, "/", res, "/", sep = "")
            wgetDownload(url = url, destdir = destdir)
            unzip(paste(raw.dir, "/", res, "/n", n, "w", w, ".zip", 
                sep = ""), exdir = paste(raw.dir, "/", res, "/n", 
                n, "w", w, sep = ""))
            tiles[t] <- raster::raster(rgdal::readGDAL(paste(raw.dir, 
                "/", res, "/n", n, "w", w, "/grdn", n, "w", w, 
                "_", res, sep = "")))
            unlink(paste(raw.dir, "/", res, "/n", n, "w", w, 
                sep = ""), recursive = TRUE)
            t <- t + 1
        }
    }
    if (length(tiles) > 1) {
        cat("Mosaic-ing NED tiles.\n\n")
        flush.console()
        tiles$fun <- mean
        mosaic.all <- do.call(raster::mosaic, tiles)
        rm(tiles)
        gc()
    }
    else {
        mosaic.all <- tiles[[1]]
    }
    mosaic.all <- raster::crop(mosaic.all, sp::spTransform(template, 
        sp::CRS(raster::projection(mosaic.all))), snap = "out")
    rgdal::writeGDAL(as(mosaic.all, "SpatialGridDataFrame"), 
        paste(rasters.dir, "/DEM_", res, ".tif", sep = ""), drivername = "GTiff", 
        type = "Float32", options = c("INTERLEAVE=PIXEL", "COMPRESS=DEFLATE", 
            "ZLEVEL=9"))
    return(mosaic.all)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
