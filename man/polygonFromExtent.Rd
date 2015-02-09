\name{polygonFromExtent}
\alias{polygonFromExtent}

\title{
Convert an Extent object to a SpatialPolygons.
}
\description{
This function takes an Extent object (or the extent of any spatial object passed to it) and creates a rectangular SpatialPolygons object of the Extent. If projection information isn't provided, it will attempt to get projection information from the object passed to it.
}
\usage{
polygonFromExtent(x, proj4string = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Extent object, or any object from which an Extent object can be extracted.
}
  \item{proj4string}{
PROJ.4 description of the coordinate reference system.
}
}
\value{
A SpatialPolygons object.
}