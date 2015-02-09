\name{getPRISM_MONTHLYPointsData}
\alias{getPRISM_MONTHLYPointsData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getPRISM_MONTHLYPointsData(sp.points, type, LT81.dir, year.range, months = c(1:12), fun = "mean")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sp.points}{
%%     ~~Describe \code{sp.points} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{LT81.dir}{
%%     ~~Describe \code{LT81.dir} here~~
}
  \item{year.range}{
%%     ~~Describe \code{year.range} here~~
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
function (sp.points, type, LT81.dir, year.range, months = c(1:12), 
    fun = "mean") 
{
    brick.years <- rep(year.range, each = 12)
    brick.months <- rep(1:12, length(year.range))
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
    all.data <- lapply(year.range, function(year) {
        extractPRISM_MONTHLYPointsYear(type = type, sp.points = sp.points, 
            LT81.dir = LT81.dir, year = year)
    })
    all.data <- do.call("rbind", all.data)
    if (fun == "mean") {
        all.data <- aggregate(all.data, by = list(brick.years), 
            FUN = mean)
    }
    else if (fun == "sum") {
        all.data <- aggregate(all.data, by = list(brick.years), 
            FUN = sum)
    }
    all.data <- all.data[all.data$Group.1 != 0, -1]
    if (fun == "mean") {
        all.data <- as.numeric(colMeans(all.data))
    }
    else if (fun == "sum") {
        all.data <- as.numeric(colSums(all.data))
    }
    return(all.data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
