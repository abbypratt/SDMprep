% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SDMprep_fxns.R
\name{bg.bb}
\alias{bg.bb}
\title{Background Bounding Box Formation}
\usage{
bg.bb(x)
}
\arguments{
\item{param1}{x = output from the occ.prep function, or a data frame with occurrence points and column names 'name', 'latitude', 'longitude'.}
}
\value{
The output is a Formal Class SpatialPolygon.
}
\description{
This function creates the spatial extent from which background points will be selected.
}
\details{
This function finds the minimum and maximum lat and long values from the occurrence points and uses that to create the extent from which background points will be selected.
}
