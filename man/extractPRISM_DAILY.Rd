\name{extractPRISM_DAILY}
\alias{extractPRISM_DAILY}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extractPRISM_DAILY(template, types = c("ppt", "tmin", "tmax"), out.dir, daily.dir, label, year.range, force.redo = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{template}{
%%     ~~Describe \code{template} here~~
}
  \item{types}{
%%     ~~Describe \code{types} here~~
}
  \item{out.dir}{
%%     ~~Describe \code{out.dir} here~~
}
  \item{daily.dir}{
%%     ~~Describe \code{daily.dir} here~~
}
  \item{label}{
%%     ~~Describe \code{label} here~~
}
  \item{year.range}{
%%     ~~Describe \code{year.range} here~~
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
function (template, types = c("ppt", "tmin", "tmax"), out.dir, 
    daily.dir, label, year.range, force.redo = FALSE) 
{
    poly <- polygonFromExtent(extent(template), projection(template))
    junk <- lapply(types, function(type) {
        getPRISM_DAILYData(type = type, label = label, template = template, 
            out.dir = out.dir, daily.dir = daily.dir, year.range = year.range, 
            force.redo = force.redo)
    })
    rm(junk)
    return(cat("PRISM Extracted\n"))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
