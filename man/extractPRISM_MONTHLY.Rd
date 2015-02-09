\name{extractPRISM_MONTHLY}
\alias{extractPRISM_MONTHLY}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extractPRISM_MONTHLY(template, out.dir, LT81.dir, label, force.redo = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{template}{
%%     ~~Describe \code{template} here~~
}
  \item{out.dir}{
%%     ~~Describe \code{out.dir} here~~
}
  \item{LT81.dir}{
%%     ~~Describe \code{LT81.dir} here~~
}
  \item{label}{
%%     ~~Describe \code{label} here~~
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
function (template, out.dir, LT81.dir, label, force.redo = FALSE) 
{
    types <- c("ppt", "tmin", "tmax")
    lapply(types, function(type) {
        getPRISM_MONTHLYData(type = type, label = label, template = template, 
            out.dir = out.dir, monthly.dir = LT81.dir, force.redo = force.redo)
    })
    return(cat("PRISM Extracted\n"))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
