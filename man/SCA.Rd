% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sca.R
\name{SCA}
\alias{SCA}
\title{Spatial Component Analysis}
\usage{
SCA(X, W, n.eigen = 20, method = c("L", "M"), scaled.data = NULL, ...)
}
\arguments{
\item{X}{A matrix with observations as rows and features as columns.}

\item{W}{A weight matrix across all observations, i.e inverse of a pairwise distance matrix.}

\item{n.eigen}{Number of spatial components (eigenvectors) to compute. Default is 20.}

\item{method}{Method used to calculate spatial cross-correlation. See \code{\link[spots]{SpatialXCorr}}.
\itemize{
\item M, using the Wartenburg's M (Default).
\item L, using the Lee's L.
}}

\item{scaled.data}{Centered and scaled data used for SVD. Default is \code{NULL}.}

\item{...}{Additional arguments passed for eigenvalue decomposition. See \code{\link[RSpectra]{eigs_sym}}.}
}
\value{
A list of Spatial Component Analysis results.
\itemize{
\item X, raw or scaled input data.
\item rotation, computed eigenvectors.
\item eigenvalues, computed eigenvalues.
\item xcor, spatial cross-correlation matrix calculated using \code{\link[spots]{SpatialXCorr}}.
}
}
\description{
Performs spatial component analysis (SCA) on the given data and weight matrices.
}
\examples{
{
data.use <- quakes[1:100,]
W <- 1/as.matrix(dist(data.use[,1:2]))
diag(W) <- 0
sca.res <- SCA(data.use[,3:5], W, n.eigen = 2)
}
}
\references{
Wartenberg, D. Multivariate spatial correlation:
A method for exploratory geographical analysis. Geogr. Anal. 17, 263–283 (1985)

Lee, S.-I. Developing a bivariate spatial association measure:
An integration of Pearson’s r and Moran's I. J. Geogr. Syst. 3, 369–385 (2001)
}
\concept{sca}
