\name{extractPRISM_DAILYYear}
\alias{extractPRISM_DAILYYear}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extractPRISM_DAILYYear(type, out.dir, daily.dir, template = NULL, year)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{out.dir}{
%%     ~~Describe \code{out.dir} here~~
}
  \item{daily.dir}{
%%     ~~Describe \code{daily.dir} here~~
}
  \item{template}{
%%     ~~Describe \code{template} here~~
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
function (type, out.dir, daily.dir, template = NULL, year) 
{
    unlink(paste(daily.dir, "/", type, "/", year, "/TEMP/", sep = ""), 
        recursive = T, force = T)
    daily.zips <- list.files(paste(daily.dir, "/", type, "/", 
        year, "/", sep = ""), full.names = T)
    junk <- sapply(daily.zips, function(i) unzip(i, exdir = paste(daily.dir, 
        "/", type, "/", year, "/TEMP/", sep = "")))
    rm(junk)
    daily.files <- list.files(paste(daily.dir, "/", type, "/", 
        year, "/TEMP", sep = ""), full.names = T)
    daily.files <- grep("*\\.bil$", daily.files, value = TRUE)
    daily.files <- grep("early", daily.files, value = TRUE, invert = T)
    year.data <- suppressWarnings(stack(daily.files, quick = T, 
        native = T))
    projection(year.data) <- CRS("+proj=longlat +datum=NAD83")
    year.data <- crop(year.data, template, snap = "out")
    names.data <- names(year.data)
    names.data <- gsub(paste("PRISM_", type, "_stable_4kmD1_", 
        sep = ""), "", names.data, )
    names.data <- gsub(paste("PRISM_", type, "_provisional_4kmD1_", 
        sep = ""), "", names.data, )
    names.data <- gsub("_bil", "", names.data, )
    names.data <- as.Date(names.data, "\%Y\%m\%d")
    names.data <- as.character(names.data, "Y\%Y.M\%m.D\%d")
    names(year.data) <- names.data
    year.data <- year.data[[order(names(year.data))]]
    unlink(paste(daily.dir, "/", type, "/", year, "/TEMP/", sep = ""), 
        recursive = T)
    return(year.data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
