\name{getNHDSubregions}
\alias{getNHDSubregions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getNHDSubregions(area.list, raw.dir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{area.list}{
%%     ~~Describe \code{area.list} here~~
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
function (area.list, raw.dir) 
{
    for (area in area.list) {
        url <- paste("ftp://nhdftp.usgs.gov/DataSets/Staged/SubRegions/FileGDB/MediumResolution/NHDM", 
            area, "_92v200.zip", sep = "")
        destdir <- raw.dir
        wgetDownload(url = url, destdir = destdir)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
