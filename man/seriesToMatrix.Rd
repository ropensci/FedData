\name{seriesToMatrix}
\alias{seriesToMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
seriesToMatrix(series, output.dir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{series}{
%%     ~~Describe \code{series} here~~
}
  \item{output.dir}{
%%     ~~Describe \code{output.dir} here~~
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
function (series, output.dir) 
{
    cat("Processing ", series, "\n")
    if (suppressWarnings(tail(scan(series, what = "character", 
        nlines = 1, quiet = T), n = 1) != "RAW")) {
        meta.1 <- as.character(read.fwf(series, c(6, 3, 52, 4), 
            skip = 0, n = 1, colClasses = "character", strip.white = T))
        meta.2 <- as.character(read.fwf(series, c(6, 3, 13, 18, 
            6, 5, 6, 9, 6, 5), skip = 1, n = 1, colClasses = "character", 
            strip.white = T))
        meta.3 <- as.character(read.fwf(series, c(6, 3, 60, 12), 
            skip = 2, n = 1, colClasses = "character", strip.white = T))
        name <- gsub(".crn", "", series)
        name <- gsub(paste(output.dir, sep = ""), "", name)
        meta <- c(name, meta.1[3:4], meta.2[c(5:7, 9:10)], meta.3[3])
        digits.year <- max(nchar(meta[7]), nchar(meta[8]), 4)
        fname <- series
        header = T
        encoding = getOption("encoding")
        con <- file(fname, encoding = encoding)
        on.exit(close(con))
        if (is.null(header)) {
            hdr1 <- readLines(con, n = 1)
            if (length(hdr1) == 0) 
                stop("file is empty")
            if (nchar(hdr1) < 10) 
                stop("first line in the crn file ends before col 10")
            yrcheck <- suppressWarnings(as.numeric(substr(hdr1, 
                7, 10)))
            if (is.null(yrcheck) || length(yrcheck) != 1 || is.na(yrcheck) | 
                yrcheck < -10000 || yrcheck > 10000) {
                is.head <- TRUE
            }
            else {
                is.head <- FALSE
            }
        }
        else if (!is.logical(header)) {
            stop("'header' must be NULL or logical")
        }
        else {
            is.head <- header
        }
        if (is.head) {
            dat1 <- readLines(con, n = 4)
            if (length(dat1) < 4) 
                stop("file has under 4 lines")
            dat1 <- dat1[4]
        }
        else {
            dat1 <- readLines(con, n = 1)
            if (length(dat1) == 0) 
                stop("file is empty")
        }
        if (nchar(dat1) < 10) 
            stop("first data line ends before col 10")
        nlines <- length(readLines(con, n = -1))
        skip.lines <- ifelse(is.head, 3, 0)
        on.exit()
        suppressWarnings(chron.stats <- read.fwf(con, c(6, digits.year, 
            6, 6, 6, 7, 9, 9, 10), skip = nlines - 1, strip.white = TRUE, 
            colClasses = "character"))
        con <- file(fname, encoding = encoding)
        if ((!is.int(chron.stats[[3]]) | (chron.stats[[3]] == 
            0)) & !grepl(" ", chron.stats[[3]])) {
            names(chron.stats) <- c("SiteID", "nYears", "AC[1]", 
                "StdDev", "MeanSens", "MeanRWI", "IndicesSum", 
                "IndicesSS", "MaxSeries")
            dat <- read.fwf(con, c(6, digits.year, rep(c(4, 3), 
                10)), skip = skip.lines, n = nlines - skip.lines - 
                1, strip.white = TRUE, colClasses = "character")
        }
        else {
            dat <- read.fwf(con, c(6, digits.year, rep(c(4, 3), 
                10)), skip = skip.lines, n = nlines - skip.lines, 
                strip.white = TRUE, colClasses = "character")
        }
        dat[[1]] <- as.character(dat[[1]])
        for (i in 2:22) {
            dat[[i]] <- as.numeric(as.character(dat[[i]]))
        }
        dat <- dat[!is.na(dat[[2]]), , drop = FALSE]
        series.label <- dat[[1]]
        series.ids <- unique(series.label)
        decade.yr <- as.numeric(dat[[2]])
        nseries <- length(series.ids)
        series.index <- match(series.label, series.ids)
        min.year <- (min(decade.yr)\%/\%10) * 10
        max.year <- ((max(decade.yr) + 10)\%/\%10) * 10
        span <- max.year - min.year + 1
        ncol.crn.mat <- nseries + 1
        crn.mat <- matrix(NA_real_, ncol = ncol.crn.mat, nrow = span)
        colnames(crn.mat) <- c(as.character(series.ids), "samp.depth")
        rownames(crn.mat) <- min.year:max.year
        x <- as.matrix(dat[seq(from = 3, to = 21, by = 2)])
        y <- as.matrix(dat[seq(from = 4, to = 22, by = 2)])
        for (i in seq_len(nseries)) {
            idx <- which(series.index == i)
            for (j in idx) {
                yr <- (decade.yr[j]\%/\%10) * 10
                row.seq <- seq(from = yr - min.year + 1, by = 1, 
                  length.out = 10)
                crn.mat[row.seq, i] <- x[j, ]
                if (i == 1) {
                  crn.mat[row.seq, ncol.crn.mat] <- y[j, ]
                }
            }
        }
        crn.mat[which(crn.mat[, -ncol.crn.mat] == 9990)] <- NA
        crn.mat <- crn.mat[!apply(is.na(crn.mat[, -ncol.crn.mat, 
            drop = FALSE]), 1, all), , drop = FALSE]
        sd.one <- all(crn.mat[, ncol.crn.mat] == 1)
        if (!is.na(sd.one) & sd.one) {
            save.names <- colnames(crn.mat)[-ncol.crn.mat]
            crn.mat <- crn.mat[, -ncol.crn.mat, drop = FALSE]
            crn.mat <- crn.mat/1000
            crn.df <- as.data.frame(crn.mat)
            names(crn.df) <- save.names
        }
        else {
            seq.series <- seq_len(nseries)
            crn.mat[, seq.series] <- crn.mat[, seq.series]/1000
            crn.df <- as.data.frame(crn.mat)
        }
        output <- crn.df
        output$YEAR <- row.names(output)
        if (ncol(output) == 2) {
            output <- output[, c(2, 1)]
        }
        else if (ncol(output) == 3) {
            output <- output[, c(3, 1)]
        }
        else {
            return("ERROR: IMPROPERLY FORMATTED CRN")
        }
        names(output) <- c("YEAR", name)
        out.mean <- mean(output[, 2], na.rm = T)
        output[, 2] <- output[, 2]/out.mean
        final.out <- list(meta = meta, data = output)
        return(final.out)
    }
    else if (tail(scan(series, what = "character", nlines = 1, 
        quiet = T), n = 1) == "RAW") {
        raw.data <- scan(series, what = "character", multi.line = F, 
            fill = T, sep = "\n", quiet = T)
        raw.data <- sub("s\237dl.", "     ", raw.data)
        raw.data <- raw.data[grepl("*STD", raw.data)]
        meta <- raw.data[1:3]
        data <- raw.data[-c(1:3)]
        data <- sub("  STD", "", data)
        meta.1 <- unlist(strsplit(meta[1], "\\s+"))[c(-1, -2)]
        place <- paste(meta.1[1:which(meta.1 == "WIDTH_RING") - 
            1], collapse = " ")
        type <- meta.1[which(meta.1 == "WIDTH_RING") + 1]
        meta.2 <- unlist(strsplit(meta[2], "\\s+"))[c(-1, -2)]
        meta.2 <- meta.2[-length(meta.2)]
        if (length(which(meta.2 == "-")) > 0) {
            meta.2 <- meta.2[-which(meta.2 == "-")]
        }
        title <- paste(meta.2[1:(length(meta.2) - 4)], collapse = " ")
        elev <- meta.2[(length(meta.2) - 3)]
        location <- meta.2[(length(meta.2) - 2)]
        location.split <- unlist(strsplit(location, ""))
        splits <- which(location.split == "-" | location.split == 
            " ") - 1
        if (splits[1] == 0) {
            splits <- splits[-1]
        }
        lat <- paste(location.split[1:splits[1]], collapse = "")
        lon <- paste(location.split[(splits[1] + 1):length(location.split)], 
            collapse = "")
        begin <- meta.2[(length(meta.2) - 1)]
        end <- meta.2[(length(meta.2))]
        widths <- c(9, 71)
        starts <- c(1, cumsum(widths) + 1)
        stops <- cumsum(widths)
        starts <- head(starts, n = length(stops))
        meta.3 <- sapply(1:length(starts), function(ii) {
            substr(meta[3], starts[ii], stops[ii])
        })[2]
        meta.3 <- sub("-", "", meta.3)
        meta.3 <- sub("STD", "", meta.3)
        meta.3 <- gsub("^\\s+|\\s+$", "", meta.3)
        name <- gsub(".crn", "", series)
        name <- gsub(paste(output.dir, sep = ""), "", name)
        meta <- c(name, place, type, elev, lat, lon, begin, end, 
            meta.3)
        digits.year <- max(nchar(meta[7]), nchar(meta[8]), 4)
        widths <- c(6, digits.year, rep(c(4, 3), 10))
        starts <- c(1, cumsum(widths) + 1)
        stops <- cumsum(widths)
        starts <- head(starts, n = length(stops))
        dat <- as.data.frame(matrix(unlist(lapply(data, FUN = function(x) {
            sapply(1:length(starts), function(ii) {
                substr(x, starts[ii], stops[ii])
            })
        })), ncol = 22, byrow = T))
        dat[[1]] <- as.character(dat[[1]])
        for (i in 2:22) {
            dat[[i]] <- as.numeric(as.character(dat[[i]]))
        }
        series.name <- dat[[1]]
        series.ids <- unique(series.name)
        decade.yr <- dat[[2]]
        nseries <- length(series.ids)
        series.index <- match(series.name, series.ids)
        min.year <- (min(decade.yr)\%/\%10) * 10
        max.year <- ((max(decade.yr) + 10)\%/\%10) * 10
        span <- max.year - min.year + 1
        ncol.crn.mat <- nseries + 1
        crn.mat <- matrix(NA_real_, ncol = ncol.crn.mat, nrow = span)
        colnames(crn.mat) <- c(as.character(series.ids), "samp.depth")
        rownames(crn.mat) <- min.year:max.year
        x <- as.matrix(dat[seq(from = 3, to = 21, by = 2)])
        y <- as.matrix(dat[seq(from = 4, to = 22, by = 2)])
        for (i in seq_len(nseries)) {
            idx <- which(series.index == i)
            for (j in idx) {
                yr <- (decade.yr[j]\%/\%10) * 10
                row.seq <- seq(from = yr - min.year + 1, by = 1, 
                  length.out = 10)
                crn.mat[row.seq, i] <- x[j, ]
                if (i == 1) {
                  crn.mat[row.seq, ncol.crn.mat] <- y[j, ]
                }
            }
        }
        crn.mat[which(crn.mat[, -ncol.crn.mat] == 9990)] <- NA
        crn.mat <- crn.mat[!apply(is.na(crn.mat[, -ncol.crn.mat, 
            drop = FALSE]), 1, all), , drop = FALSE]
        seq.series <- seq_len(nseries)
        crn.mat[, seq.series] <- crn.mat[, seq.series]/1000
        crn.df <- as.data.frame(crn.mat)
        output <- crn.df
        output$YEAR <- row.names(output)
        if (ncol(output) == 2) {
            output <- output[, c(2, 1)]
        }
        else if (ncol(output) == 3) {
            output <- output[, c(3, 1)]
        }
        else {
            return("ERROR: IMPROPERLY FORMATTED CRN")
        }
        names(output) <- c("YEAR", name)
        out.mean <- mean(output[, 2], na.rm = T)
        output[, 2] <- output[, 2]/out.mean
        final.out <- list(meta = meta, data = output)
        return(final.out)
    }
    else {
        return(paste("ERROR! IMPROPERLY FORMATTED CRN: ", series, 
            sep = ""))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
