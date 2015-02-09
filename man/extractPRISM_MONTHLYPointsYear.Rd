\name{extractPRISM_MONTHLYPointsYear}
\alias{extractPRISM_MONTHLYPointsYear}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extractPRISM_MONTHLYPointsYear(type, out.dir, LT81.dir, sp.points, year)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{out.dir}{
%%     ~~Describe \code{out.dir} here~~
}
  \item{LT81.dir}{
%%     ~~Describe \code{LT81.dir} here~~
}
  \item{sp.points}{
%%     ~~Describe \code{sp.points} here~~
}
  \item{year}{
%%     ~~Describe \code{year} here~~
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
function (type, out.dir, LT81.dir, sp.points, year) 
{
    unzip(paste(LT81.dir, type, "/", type, "_", year, ".zip", 
        sep = ""), exdir = "TEMP/")
    months <- 1:12
    year.data <- mclapply(months, function(month) {
        extractPRISM_MONTHLYPointsMonth(type = type, sp.points = sp.points, 
            year = year, month = month)
    }, mc.cores = detectCores())
    year.data <- do.call("rbind", year.data)
    unlink("TEMP/", recursive = T)
    return(year.data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
