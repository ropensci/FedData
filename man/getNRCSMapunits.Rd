\name{getNRCSMapunits}
\alias{getNRCSMapunits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getNRCSMapunits(template, areas, raw.dir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{template}{
%%     ~~Describe \code{template} here~~
}
  \item{areas}{
%%     ~~Describe \code{areas} here~~
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
function (template, areas, raw.dir) 
{
    cat("Loading the NRCS soil survey polygons for each region.\n")
    NRCSPolys <- vector("list", length(areas))
    template <- sp::spTransform(template, CRS(projection(areas)))
    for (i in 1:length(areas)) {
        area <- as.character(areas$areasymbol[i])
        date <- as.Date(areas$saverest[i])
        poly <- getNRCSStudyAreaMapunits(area = area, date = date, 
            raw.dir = raw.dir)
        poly <- poly[!is.na(poly \%over\% template), ]
        poly <- sp::spChFIDs(poly, as.character(paste(area, "_", 
            row.names(poly@data), sep = "")))
        NRCSPolys[[i]] <- poly
    }
    cat("Merging all NRCS Map Unit polygons\n")
    NRCSPolys <- do.call("rbind", NRCSPolys)
    cat("Cropping NRCS Map Unit polygons to the extent of the template\n")
    NRCSPolys <- raster::crop(NRCSPolys, template)
    return(NRCSPolys)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
