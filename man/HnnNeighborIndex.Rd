% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hnn.R
\name{HnnNeighborIndex}
\alias{HnnNeighborIndex}
\title{Hexagonal nearest neighbor index}
\usage{
HnnNeighborIndex(dist.hnn, k, include.self = TRUE)
}
\arguments{
\item{dist.hnn}{A hexagonal nearest neighbor distance matrix.}

\item{k}{Number of neighbors.}

\item{include.self}{Whether to include self as 1st neighbor, default is \code{TRUE}.}
}
\value{
An n x k matrix for the nearest neighbor indice.
}
\description{
Get hexagonal nearest neighbor indices.
}
\keyword{internal}
