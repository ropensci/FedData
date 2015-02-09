\name{loadNHDLayers}
\alias{loadNHDLayers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
loadNHDLayers(template, area.list, raw.dir, vectors.dir, force.redo = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{template}{
%%     ~~Describe \code{template} here~~
}
  \item{area.list}{
%%     ~~Describe \code{area.list} here~~
}
  \item{raw.dir}{
%%     ~~Describe \code{raw.dir} here~~
}
  \item{vectors.dir}{
%%     ~~Describe \code{vectors.dir} here~~
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
function (template, area.list, raw.dir, vectors.dir, force.redo = FALSE) 
{
    all.data <- lapply(area.list, function(area) {
        unzip(paste(raw.dir, "/NHDM", area, "_92v200.zip", sep = ""), 
            exdir = paste(raw.dir, sep = ""))
        dsn <- normalizePath(paste(raw.dir, "/NHDM", area, ".gdb", 
            sep = ""))
        layers <- rgdal::ogrListLayers(dsn)
        shapes <- lapply(layers, function(layer) {
            tryCatch(suppressWarnings(rgdal::readOGR(dsn = dsn, 
                layer = layer, verbose = F)), error = function(e) NULL)
        })
        names(shapes) <- layers
        shapes <- lapply(shapes, function(shape) {
            tryCatch(sp::spChFIDs(shape, as.character(paste(area, 
                "_", row.names(shape@data), sep = ""))), error = function(e) NULL)
        })
        unlink(paste(raw.dir, "/NHDM", area, ".gdb", sep = ""), 
            recursive = TRUE)
        return(shapes)
    })
    layers <- unique(unlist(lapply(all.data, names)))
    merged.data <- lapply(layers, function(layer) {
        shapes <- sapply(all.data, "[[", layer)
        null.shapes <- sapply(shapes, is.null)
        shapes <- do.call("rbind", shapes[!null.shapes])
        if (is.null(shapes)) 
            return(shapes)
        shapes <- raster::crop(shapes, sp::spTransform(template, 
            CRS(projection(shapes))))
        if (is.null(shapes)) 
            return(shapes)
        suppressWarnings(rgdal::writeOGR(shapes, vectors.dir, 
            layer, "ESRI Shapefile", overwrite_layer = TRUE))
        return(shapes)
    })
    names(merged.data) <- layers
    null.layers <- sapply(merged.data, is.null)
    merged.data <- merged.data[!null.layers]
    return(merged.data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
