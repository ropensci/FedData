\name{downloadGHCNDaily}
\alias{downloadGHCNDaily}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
downloadGHCNDaily(ID, raw.dir = "../DATA/", force.redo = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ID}{
%%     ~~Describe \code{ID} here~~
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
function (ID, raw.dir = "../DATA/", force.redo = F) 
{
    if (!file.exists(raw.dir)) {
        dir.create(raw.dir, recursive = T)
    }
    urls <- paste("ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/all/", 
        ID, ".dly", sep = "")
    out.files <- paste(raw.dir, ID, ".dly", sep = "")
    if (!force.redo) {
        existing.files <- sapply(out.files, file.exists)
        urls <- urls[!existing.files]
        out.files <- out.files[!existing.files]
    }
    out <- mapply(function(x, y) download.file(url = x, destfile = y, 
        mode = "wb"), urls, out.files)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
