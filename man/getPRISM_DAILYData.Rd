\name{getPRISM_DAILYData}
\alias{getPRISM_DAILYData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getPRISM_DAILYData(type, label, out.dir, daily.dir, template, year.range, date.format = "X%Y.%m.%d", force.redo = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{label}{
%%     ~~Describe \code{label} here~~
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
  \item{year.range}{
%%     ~~Describe \code{year.range} here~~
}
  \item{date.format}{
%%     ~~Describe \code{date.format} here~~
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
function (type, label, out.dir, daily.dir, template, year.range, 
    date.format = "X\%Y.\%m.\%d", force.redo = FALSE) 
{
    if (!file.exists(paste(out.dir, label, "/", sep = ""))) {
        dir.create(paste(out.dir, label, "/", sep = ""))
    }
    if (!force.redo & file.exists(paste(out.dir, label, "/", 
        type, "_DAILY_FINAL_800.tif", sep = ""))) {
        all.data <- brick(paste(out.dir, label, "/", type, "_DAILY_FINAL_800.tif", 
            sep = "")) * 1
        names(all.data) <- base::as.Date(read.csv(paste(out.dir, 
            label, "/", type, "_DAILY_FINAL_800_BANDS.csv", sep = ""), 
            stringsAsFactors = F)$x, format = date.format)
        return(all.data)
    }
    sim.poly.latlon <- spTransform(template, CRS("+proj=longlat +datum=NAD83"))
    all.data <- mclapply(year.range, function(year) {
        extractPRISM_DAILYYear(type = type, template = sim.poly.latlon, 
            out.dir = out.dir, daily.dir = daily.dir, year = year)
    }, mc.cores = detectCores())
    all.data <- stack(unlist(all.data), quick = T)
    if (!file.exists(paste(out.dir, label, "/", sep = ""))) {
        dir.create(paste(out.dir, label, "/", sep = ""))
    }
    writeRaster(all.data, paste(out.dir, label, "/", type, "_DAILY_FINAL_800.tif", 
        sep = ""), datatype = "FLT4S", options = c("COMPRESS=DEFLATE", 
        "ZLEVEL=9", "INTERLEAVE=BAND"), overwrite = T, setStatistics = FALSE)
    write.csv(names(all.data), paste(out.dir, label, "/", type, 
        "_DAILY_FINAL_800_BANDS.csv", sep = ""), row.names = F, 
        col.names = F)
    gc()
    return(all.data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
