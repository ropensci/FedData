% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/SSURGO_FUNCTIONS.R
\name{get_ssurgo}
\alias{get_ssurgo}
\title{Download and crop data from the NRCS SSURGO soils database.}
\usage{
get_ssurgo(template, label, raw.dir = "./RAW/SSURGO/",
  extraction.dir = "./EXTRACTIONS/SSURGO/", force.redo = FALSE)
}
\arguments{
\item{template}{A Raster* or Spatial* object to serve
as a template for cropping; optionally, a vector of area names [e.g., c("IN087","IN088")] may be provided.}

\item{label}{A character string naming the study area.}

\item{raw.dir}{A character string indicating where raw downloaded files should be put.
The directory will be created if missing. Defaults to "./RAW/SSURGO/".}

\item{extraction.dir}{A character string indicating where the extracted and cropped SSURGO shapefiles should be put.
The directory will be created if missing. Defaults to "./EXTRACTIONS/SSURGO/".}

\item{force.redo}{If an extraction for this template and label already exists, should a new one be created? Defaults to FALSE.}
}
\value{
A named list containing the "spatial" and "tabular" data.
}
\description{
This is an efficient method for spatially merging several different soil survey areas
as well as merging their tabular data.
}
\details{
\code{get_ssurgo} returns a named list of length 2:
\enumerate{
\item "spatial": A \code{SpatialPolygonsDataFrame} of soil mapunits
in the template, and
\item "tabular": A named list of \code{\link{data.frame}s} with the SSURGO tabular data.
}
}
