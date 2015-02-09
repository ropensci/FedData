\name{getITRDBMetadata}
\alias{getITRDBMetadata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getITRDBMetadata(series.names, data.dir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{series.names}{
%%     ~~Describe \code{series.names} here~~
}
  \item{data.dir}{
%%     ~~Describe \code{data.dir} here~~
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
function (series.names, data.dir) 
{
    dendros.study <- read.csv(paste(data.dir, "ITRDB_METADATA.csv", 
        sep = ""))
    dendros.study <- dendros.study[!is.na(dendros.study$LON) & 
        (dendros.study$SERIES \%in\% series.names), ]
    dendros.study <- dendros.study[order(dendros.study$SERIES), 
        ]
    dendros.study <- SpatialPointsDataFrame(dendros.study[, c("LON", 
        "LAT")], dendros.study, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
    suppressWarnings(writeOGR(dendros.study, paste(data.dir, 
        "ITRDB_METADATA.kml", sep = ""), "ITRDB_METADATA", driver = "KML", 
        overwrite_layer = TRUE))
    return(dendros.study)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
