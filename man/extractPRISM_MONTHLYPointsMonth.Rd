\name{extractPRISM_MONTHLYPointsMonth}
\alias{extractPRISM_MONTHLYPointsMonth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extractPRISM_MONTHLYPointsMonth(type, sp.points, year, month)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{sp.points}{
%%     ~~Describe \code{sp.points} here~~
}
  \item{year}{
%%     ~~Describe \code{year} here~~
}
  \item{month}{
%%     ~~Describe \code{month} here~~
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
function (type, sp.points, year, month) 
{
    month.char <- formatC(month, width = 2, format = "d", flag = "0")
    file.path <- paste("TEMP/", year, "/cai_", type, "_us_us_30s_", 
        year, month.char, ".bil", sep = "")
    rast <- extractPRISMPointsFromFile(file.path = file.path, 
        sp.points = sp.points)
    return(rast)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
