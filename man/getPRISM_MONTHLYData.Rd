\name{getPRISM_MONTHLYData}
\alias{getPRISM_MONTHLYData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getPRISM_MONTHLYData(template, type, label, out.dir, monthly.dir, force.redo = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{template}{
%%     ~~Describe \code{template} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{label}{
%%     ~~Describe \code{label} here~~
}
  \item{out.dir}{
%%     ~~Describe \code{out.dir} here~~
}
  \item{monthly.dir}{
%%     ~~Describe \code{monthly.dir} here~~
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
function (template, type, label, out.dir, monthly.dir, force.redo = FALSE) 
{
    if (!file.exists(paste(out.dir, label, "/", sep = ""))) {
        dir.create(paste(out.dir, label, "/", sep = ""))
    }
    if (!force.redo & file.exists(paste(out.dir, label, "/", 
        type, "_MONTHLY_FINAL_800.tif", sep = ""))) {
        all.data <- brick(paste(out.dir, label, "/", type, "_MONTHLY_FINAL_800.tif", 
            sep = "")) * 1
        names(all.data) <- read.csv(paste(out.dir, label, "/", 
            type, "_MONTHLY_FINAL_800_BANDS.csv", sep = ""))$x
        return(all.data)
    }
    sim.poly.latlon <- spTransform(template, CRS("+proj=longlat +datum=NAD83"))
    monthly.files <- list.files(paste(monthly.dir, type, sep = ""), 
        recursive = T, full.names = T)
    if (!length(grep("*\\.bil$", monthly.files, value = TRUE)) > 
        0) {
        system(paste("ls ", paste(monthly.dir, type, sep = ""), 
            "/* | parallel unzip -qqo {} -d ", paste(monthly.dir, 
                type, sep = ""), sep = ""))
        monthly.files <- list.files(paste(monthly.dir, type, 
            sep = ""), recursive = T, full.names = T)
    }
    monthly.files <- grep("*\\.bil$", monthly.files, value = TRUE)
    monthly.files <- grep("spqc", monthly.files, value = TRUE, 
        invert = T)
    monthly.files <- grep("/cai", monthly.files, value = TRUE)
    all.data <- mclapply(as.list(monthly.files), function(monthly.file, 
        ...) {
        extractPRISM_MONTHLYMonth(template = sim.poly.latlon, 
            file = monthly.file)
    }, mc.cores = detectCores())
    all.data <- brick(all.data)
    all.data <- all.data[[order(names(all.data))]]
    if (!file.exists(paste(out.dir, label, "/", sep = ""))) {
        dir.create(paste(out.dir, label, "/", sep = ""))
    }
    writeRaster(all.data, paste(out.dir, label, "/", type, "_MONTHLY_FINAL_800.tif", 
        sep = ""), datatype = "FLT4S", options = c("COMPRESS=DEFLATE", 
        "ZLEVEL=9", "INTERLEAVE=BAND"), overwrite = T, setStatistics = FALSE)
    write.csv(names(all.data), paste(out.dir, label, "/", type, 
        "_MONTHLY_FINAL_800_BANDS.csv", sep = ""), row.names = F, 
        col.names = F)
    gc()
    gc()
    return(NULL)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
