\name{annualizePRISM_MONTHLY}
\alias{annualizePRISM_MONTHLY}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
annualizePRISM_MONTHLY(prism.brick, months = c(1:12), fun)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prism.brick}{
%%     ~~Describe \code{prism.brick} here~~
}
  \item{months}{
%%     ~~Describe \code{months} here~~
}
  \item{fun}{
%%     ~~Describe \code{fun} here~~
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
function (prism.brick, months = c(1:12), fun) 
{
    brick.names <- names(prism.brick)
    brick.years <- gsub("([.].*)", "", brick.names)
    brick.years <- as.numeric(gsub("Y", "", brick.years))
    brick.months <- as.numeric(gsub("Y\\d{4}[.]M", "", brick.names))
    if (length(months) > 12) {
        stop("ERROR! Too many months.")
    }
    previous.year <- 12 + months[months < 1]
    current.year <- months[months >= 1 & months <= 12]
    next.year <- months[months > 12] - 12
    no.year <- (1:12)[!(1:12 \%in\% c(previous.year, current.year, 
        next.year))]
    brick.years[brick.months \%in\% previous.year] <- brick.years[brick.months \%in\% 
        previous.year] + 1
    brick.years[brick.months \%in\% next.year] <- brick.years[brick.months \%in\% 
        next.year] - 1
    brick.years[brick.months \%in\% no.year] <- 0
    signal.brick.temp <- setZ(prism.brick, brick.years, "year")
    if (fun == "sum") {
        signal.brick.temp <- zApply(signal.brick.temp, by = getZ(signal.brick.temp), 
            fun = sum)
    }
    else if (fun == "mean") {
        signal.brick.temp <- zApply(signal.brick.temp, by = getZ(signal.brick.temp), 
            fun = mean)
    }
    signal.brick.temp <- signal.brick.temp[[(1:nlayers(signal.brick.temp))[names(signal.brick.temp) != 
        "X0"]]]
    return(signal.brick.temp)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
