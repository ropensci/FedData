\name{extractNRCSData}
\alias{extractNRCSData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extractNRCSData(tables, mapunits)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tables}{
%%     ~~Describe \code{tables} here~~
}
  \item{mapunits}{
%%     ~~Describe \code{mapunits} here~~
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
function (tables, mapunits) 
{
    mapunits <- as.character(unique(mapunits$MUKEY))
    mapping <- tables[["mdstatrshipdet"]]
    mappingGraph <- igraph::graph.edgelist(as.matrix(mapping[, 
        c("ltabphyname", "rtabphyname")]))
    igraph::E(mappingGraph)$mapVar <- as.character(mapping[, 
        "ltabcolphyname"])
    mappingGraph <- igraph::graph.neighborhood(mappingGraph, 
        order = max(sapply(igraph::decompose.graph(mappingGraph), 
            diameter)) + 1, nodes = "mapunit", mode = "out")[[1]]
    mapHierarchy <- igraph::shortest.paths(mappingGraph, "mapunit")
    mapHierarchy <- colnames(mapHierarchy)[order(mapHierarchy)]
    mapHierarchy <- mapHierarchy[-1]
    mapEdges <- cbind(igraph::get.edgelist(mappingGraph), igraph::E(mappingGraph)$mapVar)
    mapEdges <- mapEdges[match(mapHierarchy, mapEdges[, 2]), 
        ]
    tables[["mapunit"]] <- tables[["mapunit"]][tables[["mapunit"]][, 
        "mukey"] \%in\% mapunits, ]
    newTables <- apply(mapEdges, 1, function(X) {
        return(tables[[X[2]]][tables[[X[2]]][, X[3]] \%in\% tables[[X[1]]][, 
            X[3]], ])
    })
    names(newTables) <- mapEdges[, 2]
    tables[names(newTables)] <- newTables
    tables <- tables[!sapply(tables, is.null)]
    return(tables)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
