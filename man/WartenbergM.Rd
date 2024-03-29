% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sca.R
\name{WartenbergM}
\alias{WartenbergM}
\title{Bivariate Moran's I using Wartenberg's M matrix}
\usage{
WartenbergM(X, W)
}
\arguments{
\item{X}{A matrix with observations as rows and features as columns.}

\item{W}{A weight matrix across all observations, i.e inverse of a pairwise distance matrix.}
}
\value{
A spatial cross-correlation matrix.
}
\description{
Calculate bivariate Moran's I using Wartenberg's method.
}
\examples{
{
data.use <- quakes[1:100,]
W <- 1/as.matrix(dist(data.use[,1:2]))
diag(W) <- 0
M <- WartenbergM(data.use[,3:4], W)
}
}
\references{
Wartenberg, D. Multivariate spatial correlation:
A method for exploratory geographical analysis. Geogr. Anal. 17, 263–283 (1985)
}
\concept{sca}
