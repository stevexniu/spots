#' @import Matrix
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
#' @examples {
#' data.use <- quakes[1:100,]
#' W <- 1/as.matrix(dist(data.use[,1:2]))
#' diag(W) <- 0
#' M <- WartenbergM(data.use[,3:4], W)
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
#' @examples {
#' data.use <- quakes[1:100,]
#' W <- 1/as.matrix(dist(data.use[,1:2]))
#' diag(W) <- 0
#' L <- LeeL(data.use[,3:4], W)
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
#' @seealso \code{\link[spots]{WartenbergM}} \code{\link[spots]{LeeL}}.
#' @examples {
#' data.use <- quakes[1:100,]
#' W <- 1/as.matrix(dist(data.use[,1:2]))
#' diag(W) <- 0
#' M <- SpatialXCorr(data.use[,3:4], W, method = "M")
#' L <- SpatialXCorr(data.use[,3:4], W, method = "L")
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
