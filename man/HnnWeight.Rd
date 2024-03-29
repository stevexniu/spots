% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hnn.R
\name{HnnWeight}
\alias{HnnWeight}
\title{Hexagonal nearest neighbor weight}
\usage{
HnnWeight(dist.hnn, dist.k = NULL, mu = 0, sigma = 1)
}
\arguments{
\item{dist.hnn}{A hexagonal nearest neighbor distance matrix.}

\item{dist.k}{The maximum distance used to calculate the weight. Default is \code{NULL} and all neighbor weights are calculated.}

\item{mu}{The mean of Gaussian filter, default is 0.}

\item{sigma}{The standard deviation of Gaussian filter, default is 1.}
}
\value{
A weight matrix.
}
\description{
Calculate hexagonal nearest neighbor weights using Gaussian filter.
}
\examples{
{
data.use <- quakes[1:100,]
dist.use <- as.matrix(dist(data.use[,1:2]))
res <- HnnWeight(dist.use, mu = 0, sigma = 0.5)
}

}
\concept{hnn}
