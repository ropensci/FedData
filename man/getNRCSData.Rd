\name{getNRCSData}
\alias{getNRCSData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getNRCSData(areas, raw.dir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{areas}{
%%     ~~Describe \code{areas} here~~
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
function (areas, raw.dir) 
{
    cat("Loading the NRCS soil survey data for each region.\n")
    NRCSTables <- vector("list", length(areas))
    for (i in 1:length(areas)) {
        area <- as.character(areas$areasymbol[i])
        date <- as.Date(areas$saverest[i])
        tables <- getNRCSStudyAreaData(area = area, date = date, 
            raw.dir = raw.dir)
        NRCSTables[[i]] <- tables
    }
    cat("Merging all NRCS data tables\n")
    tableNames <- unique(unlist(sapply(NRCSTables, names)))
    tableNames <- tableNames[order(tableNames)]
    NRCSTables <- lapply(tableNames, function(name) {
        tables <- lapply(NRCSTables, "[[", name)
        tables <- do.call("rbind", tables)
        tables <- unique(tables)
        return(tables)
    })
    names(NRCSTables) <- tableNames
    return(NRCSTables)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
