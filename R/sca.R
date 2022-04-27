#' @import RSpectra Matrix
#'
NULL

#' Bivariate Moran's I using Wartenberg's M matrix
#'
#' Calculate bivariate Moran's I using Wartenberg's method.
#' @param X A matrix with observations as rows and features as columns.
#' @param W A weight matrix across all observations, i.e inverse of a pairwise distance matrix.
#' @return A spatial cross-correlation matrix.
#' @export
#' @concept sca
#' @examples \dontrun{
#' M <- WartenbergM(X, W)
#' }
#' @references
#' Wartenberg, D. Multivariate spatial correlation:
#' A method for exploratory geographical analysis. Geogr. Anal. 17, 263–283 (1985)
#'
WartenbergM <- function(X, W){
  N <- nrow(X)
  X <- scale(X, center = TRUE, scale = FALSE)
  W.sum <- sum(W)
  W.sum[W.sum==0] <- 1
  W <- W/W.sum
  S <- sqrt(colSums(X*X)/N)
  M <- t(X)%*%W%*%X/tcrossprod(S)
  if(is(M,"denseMatrix")) M <- as.matrix(M)
  return(M)
}

#' Bivariate Moran's I using Lee's L measurement
#'
#' Calculate bivariate Moran's I using Lee's L measurement.
#'
#' @param X A matrix with observations as rows and features as columns.
#' @param W A weight matrix across all observations, i.e inverse of a pairwise distance matrix.
#' @return A spatial cross-correlation matrix.
#' @export
#' @concept sca
#' @examples \dontrun{
#' L <- LeeL(X, W)
#' }
#' @references
#' Lee, S.-I. Developing a bivariate spatial association measure:
#' An integration of Pearson’s r and Moran's I. J. Geogr. Syst. 3, 369–385 (2001)
#'
LeeL <- function(X, W){
  N <- nrow(X)
  X <- scale(X, center = TRUE, scale = TRUE)
  W.rowsum <- rowSums(W)
  W.rowsum[W.rowsum==0] <- 1
  W <- W/W.rowsum
  L <- t(X)%*%t(W)%*%W%*%X/N
  if(is(L,"denseMatrix")) L <- as.matrix(L)
  return(L)
}

#' Spatial cross-correlation
#'
#' Calculate spatial cross-correlation matrix.
#'
#' @param X A matrix with observations as rows and features as columns.
#' @param W A weight matrix across all observations, i.e inverse of a pairwise distance matrix.
#' @param method Method used to calculate spatial cross-correlation.
#' \itemize{
#'   \item M, using the Wartenburg's M (Default).
#'   \item L, using the Lee's L.
#' }
#' @return A spatial cross-correlation matrix.
#' @export
#' @concept sca
#' @seealso \code{\link[rspca]{WartenbergM}} \code{\link[rspca]{LeeL}}.
#' @examples \dontrun{
#' M <- SpatialXCorr(X, W, method = "M")
#' L <- SpatialXCorr(X, W, method = "L")
#' }
#' @references
#' Wartenberg, D. Multivariate spatial correlation:
#' A method for exploratory geographical analysis. Geogr. Anal. 17, 263–283 (1985)
#'
#' Lee, S.-I. Developing a bivariate spatial association measure:
#' An integration of Pearson’s r and Moran's I. J. Geogr. Syst. 3, 369–385 (2001)
#'
SpatialXCorr <- function(X, W, method = c("M", "L")){
  xcor <- switch(match.arg(method),
                 M = WartenbergM,
                 L = LeeL
                 )
  xcor(X, W)
}

#' Spatial Component Analysis
#'
#' Performs spatial component analysis (SCA) on the given data and weight matrices.
#'
#' @param X A matrix with observations as rows and features as columns.
#' @param W A weight matrix across all observations, i.e inverse of a pairwise distance matrix.
#' @param n.eigen Number of spatial components (eigenvectors) to compute. Default is 20.
#' @param method Method used to calculate spatial cross-correlation. See \code{\link[rspca]{SpatialXCorr}}.
#' \itemize{
#'   \item M, using the Wartenburg's M (Default).
#'   \item L, using the Lee's L.
#' }
#' @param scaled.data Centered and scaled data used for SVD. Default is \code{NULL}.
#' @param ... Additional arguments passed for eigenvalue decomposition. See \code{\link[RSpectra]{eigs_sym}}.
#' @return A list of Spatial Component Analysis results.
#' \itemize{
#'   \item X, raw or scaled input data.
#'   \item rotation, computed eigenvectors.
#'   \item eigenvalues, computed eigenvalues.
#'   \item xcor, spatial cross-correlation matrix calculated using \code{\link[rspca]{SpatialXCorr}}.
#' }
#' @export
#' @concept sca
#' @importFrom RSpectra eigs_sym
#' @examples \dontrun{
#' sca.res <- SCA(X, W, n.eigen = 30)
#' }
#' @references
#' Wartenberg, D. Multivariate spatial correlation:
#' A method for exploratory geographical analysis. Geogr. Anal. 17, 263–283 (1985)
#'
#' Lee, S.-I. Developing a bivariate spatial association measure:
#' An integration of Pearson’s r and Moran's I. J. Geogr. Syst. 3, 369–385 (2001)
#'
SCA <- function(X, W, n.eigen = 20, method = c("L","M"), scaled.data = NULL, ...){
  COV <- SpatialXCorr(X, W, match.arg(method))
  Eig <- eigs_sym(COV, n.eigen, ...)
  names(Eig$values) <- paste0("SC_", 1:n.eigen)
  colnames(Eig$vectors) <- paste0("SC_", 1:n.eigen)
  rownames(Eig$vectors) <- colnames(X)
  if(!is.null(scaled.data)){
    X <- scaled.data %*% Eig$vectors
  } else {
    X <- scale(X) %*% Eig$vectors
  }
  return(list(X = X, rotation = Eig$vectors, eigenvalues = Eig$values, xcor = COV))
}
