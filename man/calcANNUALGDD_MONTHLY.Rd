\name{calcANNUALGDD_MONTHLY}
\alias{calcANNUALGDD_MONTHLY}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calcANNUALGDD_MONTHLY(extraction.dir, months, t.base, t.cap = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{extraction.dir}{
%%     ~~Describe \code{extraction.dir} here~~
}
  \item{months}{
%%     ~~Describe \code{months} here~~
}
  \item{t.base}{
%%     ~~Describe \code{t.base} here~~
}
  \item{t.cap}{
%%     ~~Describe \code{t.cap} here~~
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
function (extraction.dir, months, t.base, t.cap = NULL) 
{
    tmin.brick <- brick(paste(extraction.dir, "/tmin_MONTHLY_FINAL_800.tif", 
        sep = ""))
    names(tmin.brick) <- read.csv(paste(extraction.dir, "/tmin_MONTHLY_FINAL_800_BANDS.csv", 
        sep = ""), stringsAsFactors = F)[, 1]
    tmax.brick <- brick(paste(extraction.dir, "/tmax_MONTHLY_FINAL_800.tif", 
        sep = ""))
    names(tmax.brick) <- read.csv(paste(extraction.dir, "/tmax_MONTHLY_FINAL_800_BANDS.csv", 
        sep = ""), stringsAsFactors = F)[, 1]
    if (nlayers(tmin.brick) != nlayers(tmax.brick) | any(names(tmin.brick) != 
        names(tmax.brick))) {
        stop("tmin and tmax bricks must have same number of layers!")
    }
    tmin.brick <- calc(tmin.brick, function(x) {
        x[x < t.base] <- t.base
        return(x)
    })
    tmax.brick <- calc(tmax.brick, function(x) {
        x[x < t.base] <- t.base
        return(x)
    })
    if (!is.null(t.cap)) {
        tmin.brick <- calc(tmin.brick, function(x) {
            x[x > t.cap] <- t.cap
            return(x)
        })
        tmax.brick <- calc(tmax.brick, function(x) {
            x[x > t.cap] <- t.cap
            return(x)
        })
    }
    GDD.brick <- ((tmin.brick + tmax.brick)/2) - t.base
    GDD.brick.months <- as.numeric(gsub("Y\\d{4}[.]M", "", names(GDD.brick)))
    year.months <- 1:12
    days.per.month <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 
        30, 31)
    GDD.brick.days <- as.numeric(mapply(gsub, year.months, days.per.month, 
        GDD.brick.months))
    GDD.brick <- GDD.brick * GDD.brick.days * 1.8
    GDD.brick.total <- annualizePRISM_MONTHLY(GDD.brick, months = months, 
        fun = "sum")
    return(GDD.brick.total)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
