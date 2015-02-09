\name{wgetDownload}
\alias{wgetDownload}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
wgetDownload(url, destdir, timestamping = T, nc = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{url}{
%%     ~~Describe \code{url} here~~
}
  \item{destdir}{
%%     ~~Describe \code{destdir} here~~
}
  \item{timestamping}{
%%     ~~Describe \code{timestamping} here~~
}
  \item{nc}{
%%     ~~Describe \code{nc} here~~
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
function (url, destdir, timestamping = T, nc = T) 
{
    if (nc) {
        status <- system(paste("wget -nc -nd --quiet --directory-prefix=", 
            destdir, " ", url, sep = ""))
    }
    else if (timestamping) {
        status <- system(paste("wget -N -nd --quiet --directory-prefix=", 
            destdir, " ", url, sep = ""))
    }
    else {
        status <- system(paste("wget -nd --quiet --directory-prefix=", 
            destdir, " ", url, sep = ""))
    }
    if (status != 0) 
        warning("Download of ", url, " had nonzero exit status")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
