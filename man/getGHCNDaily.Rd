\name{getGHCNDaily}
\alias{getGHCNDaily}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getGHCNDaily(template = NULL, elements = NULL, label = NULL, raw.dir = "./RAW/GHCN/", extraction.dir = "./EXTRACTIONS/GHCN/", standardize = F, force.redo = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{template}{
%%     ~~Describe \code{template} here~~
}
  \item{elements}{
%%     ~~Describe \code{elements} here~~
}
  \item{label}{
%%     ~~Describe \code{label} here~~
}
  \item{raw.dir}{
%%     ~~Describe \code{raw.dir} here~~
}
  \item{extraction.dir}{
%%     ~~Describe \code{extraction.dir} here~~
}
  \item{standardize}{
%%     ~~Describe \code{standardize} here~~
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
function (template = NULL, elements = NULL, label = NULL, raw.dir = "./RAW/GHCN/", 
    extraction.dir = "./EXTRACTIONS/GHCN/", standardize = F, 
    force.redo = F) 
{
    dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
    if (is.null(template)) {
        label <- "allStations"
    }
    vectors.dir <- paste(extraction.dir, "/", label, "/spatial", 
        sep = "")
    tables.dir <- paste(extraction.dir, "/", label, "/tabular", 
        sep = "")
    dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(extraction.dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(vectors.dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(tables.dir, showWarnings = FALSE, recursive = TRUE)
    cat("\nGetting spatial data of GHCN stations")
    if (!force.redo & file.exists(paste(vectors.dir, "/stations.shp", 
        sep = ""))) {
        stations.sp <- rgdal::readOGR(dsn = vectors.dir, layer = "stations")
    }
    else {
        stations.sp <- getGHCNStations(template = template, raw.dir = raw.dir)
        suppressWarnings(rgdal::writeOGR(stations.sp, vectors.dir, 
            "stations", "ESRI Shapefile", overwrite_layer = TRUE))
    }
    if (is.null(elements)) {
        elements <- unique(stations.sp$ELEMENT)
    }
    stations.sp <- stations.sp[stations.sp@data[, "ELEMENT"] \%in\% 
        toupper(elements), ]
    if (standardize & !is.null(elements)) {
        stations.sp.splits <- split(as.character(stations.sp$ELEMENT), 
            f = stations.sp$ID, drop = T)
        stations.sp.splits.all <- sapply(stations.sp.splits, 
            function(x) {
                all(sapply(toupper(elements), function(y) {
                  y \%in\% x
                }))
            })
        stations.sp <- stations.sp[stations.sp$ID \%in\% names(stations.sp.splits.all)[stations.sp.splits.all], 
            ]
    }
    stations.out <- stations.sp[, c("ID", "ELEMENT", "YEAR_START", 
        "YEAR_END")]
    stations.sp <- stations.sp[!duplicated(stations.sp@data[, 
        c("ID", "LATITUDE", "LONGITUDE")]), c("ID", "LATITUDE", 
        "LONGITUDE")]
    if (!force.redo) {
        daily <- tryCatch(lapply(elements, function(element) {
            readRDS(paste(tables.dir, "/", element, ".Rds", sep = ""))
        }), warning = function(w) {
            return(NULL)
        })
        if (!is.null(daily)) {
            names(daily) <- elements
            daily <- lapply(as.character(stations.sp$ID), function(station) {
                stationDaily <- tryCatch(lapply(daily, "[[", 
                  station), error = function(e) {
                  return(NULL)
                })
                stationDaily <- stationDaily[!sapply(stationDaily, 
                  is.null)]
                return(stationDaily)
            })
            names(daily) <- as.character(stations.sp$ID)
            daily <- daily[!sapply(daily, is.null)]
            if (setequal(names(daily), stations.sp$ID) & all(sapply(daily, 
                function(dat) {
                  setequal(names(dat), elements)
                }))) {
                return(list(spatial = stations.out, tabular = daily))
            }
        }
    }
    cat("\nDownloading daily GHCN data")
    downloadGHCNDaily(stations.sp$ID, raw.dir = raw.dir, force.redo = force.redo)
    daily <- lapply(1:length(stations.sp$ID), function(i) {
        id <- stations.sp$ID[i]
        cat("\nProcessing daily GHCH station ", i, "of", length(stations.sp$ID))
        daily <- utils::read.fwf(paste(raw.dir, id, ".dly", sep = ""), 
            c(11, 4, 2, 4, rep(c(5, 1, 1, 1), 31)))
        names(daily)[1:4] <- c("STATION", "YEAR", "MONTH", "ELEMENT")
        daily <- daily[daily$ELEMENT \%in\% toupper(elements), 
            c(2:4, seq(5, 125, 4))]
        daily[daily == -9999] <- NA
        names(daily) <- c("YEAR", "MONTH", "ELEMENT", paste("D", 
            1:31, sep = ""))
        out.list <- lapply(elements, function(element) {
            return(daily[daily$ELEMENT == toupper(element), -3])
        })
        if (standardize) {
            yearMonths <- lapply(out.list, function(element) {
                element <- element[order(element$YEAR, element$MONTH), 
                  ]
                return(paste("Y", element[, c("YEAR")], "M", 
                  element[, c("MONTH")], sep = ""))
            })
            all.yearMonths <- Reduce(intersect, yearMonths)
            out.list <- lapply(out.list, function(element) {
                element.yearMonths <- paste("Y", element[, c("YEAR")], 
                  "M", element[, c("MONTH")], sep = "")
                return(element[match(all.yearMonths, element.yearMonths), 
                  ])
            })
        }
        names(out.list) <- elements
        return(out.list)
    })
    names(daily) <- stations.sp$ID
    daily.split <- lapply(elements, function(element) {
        lapply(daily, "[[", element)
    })
    names(daily.split) <- elements
    junk <- lapply(elements, function(element) {
        saveRDS(daily.split[[element]], paste(tables.dir, "/", 
            element, ".Rds", sep = ""), compress = "xz")
    })
    return(list(spatial = stations.out, tabular = daily))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
