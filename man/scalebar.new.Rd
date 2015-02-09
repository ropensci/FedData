\name{scalebar.new}
\alias{scalebar.new}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scalebar.new(d, xy = NULL, height = NULL, line.offset = c(0, 0), side = "right", lab.side = "top", lonlat = NULL, label, adj = c(0.5, -0.5), lwd = 2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
%%     ~~Describe \code{d} here~~
}
  \item{xy}{
%%     ~~Describe \code{xy} here~~
}
  \item{height}{
%%     ~~Describe \code{height} here~~
}
  \item{line.offset}{
%%     ~~Describe \code{line.offset} here~~
}
  \item{side}{
%%     ~~Describe \code{side} here~~
}
  \item{lab.side}{
%%     ~~Describe \code{lab.side} here~~
}
  \item{lonlat}{
%%     ~~Describe \code{lonlat} here~~
}
  \item{label}{
%%     ~~Describe \code{label} here~~
}
  \item{adj}{
%%     ~~Describe \code{adj} here~~
}
  \item{lwd}{
%%     ~~Describe \code{lwd} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (d, xy = NULL, height = NULL, line.offset = c(0, 0), 
    side = "right", lab.side = "top", lonlat = NULL, label, adj = c(0.5, 
        -0.5), lwd = 2, ...) 
{
    pr <- par()
    if (is.null(lonlat)) {
        if (pr$usr[1] > -181 & pr$usr[2] < 181 & pr$yaxp[1] > 
            -200 & pr$yaxp[2] < 200) {
            lonlat <- TRUE
        }
        else {
            lonlat <- FALSE
        }
    }
    if (lonlat) {
        lat <- mean(pr$yaxp[1:2])
        if (missing(d)) {
            dx <- (pr$usr[2] - pr$usr[1])/10
            d <- raster::pointDistance(cbind(0, lat), cbind(dx, 
                lat), TRUE)
            d <- signif(d/1000, 2)
            label <- NULL
        }
        p <- cbind(0, lat)
        dd <- raster:::.destPoint(p, d * 1000)
        dd <- dd[1, 1]
    }
    else {
        if (missing(d)) {
            d <- round(10 * (pr$usr[2] - pr$usr[1])/10)/10
            label <- NULL
        }
        dd <- d
    }
    if (is.null(xy)) {
        padding = c(5, 5)/100
        parrange <- c(pr$usr[2] - pr$usr[1], pr$usr[4] - pr$usr[3])
        xy <- c(pr$usr[1] + (padding[1] * parrange[1]), pr$usr[3] + 
            (padding[2] * parrange[2]))
    }
    xy <- xy + line.offset
    if (side == "right") {
        xstart = xy[1]
        xend = xy[1] + dd
    }
    else {
        xstart = xy[1] - dd
        xend = xy[1]
    }
    if (is.null(height)) {
        height <- dd * 0.1
    }
    rect(xleft = xstart, ybottom = xy[2], xright = xend, ytop = xy[2] + 
        height, col = "black", border = NA, lend = 1, xpd = T)
    if (missing(label)) {
        label <- paste(d)
    }
    if (is.null(label)) {
        label <- paste(d)
    }
    if (missing(adj)) {
        adj <- c(0.5, -0.2 - lwd/20)
    }
    if (lab.side == "top") {
        text(mean(c(xstart, xend)), xy[2], labels = label, adj = c(0.5, 
            -0.5), ...)
    }
    else if (lab.side == "right") {
        text(xend, xy[2], labels = label, adj = c(-0.1, 0), ...)
    }
    else if (lab.side == "left") {
        text(xstart, xy[2], labels = label, adj = c(1.1, 0), 
            ...)
    }
    else if (lab.side == "bottom") {
        text(mean(c(xstart, xend)), xy[2], labels = label, adj = c(0.5, 
            1.5), ...)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
