\name{createArea}
\alias{createArea}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
createArea(North, South, East, West, projection.string)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{North}{
%%     ~~Describe \code{North} here~~
}
  \item{South}{
%%     ~~Describe \code{South} here~~
}
  \item{East}{
%%     ~~Describe \code{East} here~~
}
  \item{West}{
%%     ~~Describe \code{West} here~~
}
  \item{projection.string}{
%%     ~~Describe \code{projection.string} here~~
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
function (North, South, East, West, projection.string) 
{
    datainUTM <- matrix(c(East, West, West, East, East, North, 
        North, South, South, North), nrow = 5)
    master.proj <- sp::CRS(projection.string)
    sim.poly <- sp::Polygons(list(Polygon(datainUTM, hole = FALSE)), 
        ID = "A")
    sim.poly <- sp::SpatialPolygons(list(sim.poly), proj4string = master.proj)
    IDs <- sapply(slot(sim.poly, "polygons"), function(x) slot(x, 
        "ID"))
    df <- data.frame(rep(0, length(IDs)), row.names = IDs)
    sim.poly <- sp::SpatialPolygonsDataFrame(sim.poly, df)
    return(sim.poly)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
