#' @include stats.R sca.R
#'
NULL

#' Univariate Moran's I
#'
#' Calculate univariate (canonical) Moran's I.
#'
#' @param X A matrix with observations as rows and features as columns.
#' @param W A weight matrix across all observations, i.e inverse of a pairwise distance matrix.
#' @param normalize Whether to normalize the weight matrix such that each row adds up to one. Default is \code{TRUE}.
#' @param alternative Alternative hypothesis used, default is \code{two.sided}.
#' @param p.adjust.method Method used for multiple comparisons correction, default is \code{BH}. See \code{\link[stats]{p.adjust}}.
#' @return A list containing the following:
#' \itemize{
#'   \item Morans.I, the Moran's I.
#'   \item Z.I, the Z score of Moran's I.
#'   \item X, data matrix used for calculating Moran's I.
#'   \item Y, a matrix of spatial lags.
#'   \item Expected.I, the expectation of Moran's I under the null hypothesis.
#'   \item SD.I, the standard deviation of Moran's I under the null hypothesis.
#'   \item p.val, p-values.
#'   \item p.adj, adjusted p-values.
#'   \item normalize, whether to normalize the weight matrix.
#'   \item alternative, alternative hypothesis used.
#'   \item p.adjust.method, method used for multiple comparisons correction.
#' }
#' @export
#' @concept moransi
#' @examples \dontrun{
#' UnivariateMoransI(X, W)
#' }
#' @references
#' Moran, P. A. P. Notes on continuous stochastic phenomena. Biometrika 37, 17–23 (1950)
#'
UnivariateMoransI <- function(X, W, normalize = TRUE, alternative = c("two.sided", "less", "greater"), p.adjust.method = "BH"){
  alternative <- match.arg(arg = alternative)
  X <- scale(X, center = TRUE, scale = FALSE)
  N <- nrow(X)
  if(normalize){
    W.rowsum <- rowSums(W)
    W.rowsum[W.rowsum==0] <- 1
    W <- W/W.rowsum
    W2 <- N^2
    S2 <- sum((1+colSums(W))^2)
  } else {
    W.sum <- sum(W)
    W <- N/W.sum*W
    W2 <- W.sum^2
    S2 <- sum((rowSums(W)+colSums(W))^2)
  }
  Y <- W%*%X
  if(is(Y,"denseMatrix")) Y <- as.matrix(Y)
  X2 <- colSums(X*X)
  I <- colSums(X*Y)/X2
  EI <- -1/(N-1)
  S1 <- 0.5*sum((W+t(W))^2)
  S3 <- (colSums(X^4)*N)/(X2^2)
  S4 <- (N^2-3*N+3)*S1-N*S2+3*W2
  S5 <- N*(N-1)*S1-2*N*S2+6*W2
  SDI <- sqrt((N*S4-S3*S5) / ((N-1)*(N-2)*(N-3)*W2) - EI^2)
  Z <- (I-EI)/SDI
  P.list <- ZPvalue(as.matrix(Z), alternative, p.adjust.method)
  return(list(Morans.I = I, Z.I = Z, X = X, Y = Y, Expected.I = EI, SD.I = SDI,  p.val = P.list[[1]][,1], p.adj = P.list[[2]][,1],
             normalize = normalize,alternative = alternative, p.adjust.method = p.adjust.method))
}

#' Permutation based Moran's I
#'
#' Calculate permutation based p-value for Moran's I.
#'
#' @param x A numerical vector.
#' @param W A weight matrix across all observations, i.e inverse of a pairwise distance matrix.
#' @param n The number of permutations to be conducted, set to 999 by default.
#' @param seed Random seed used. Default is 1.
#' @param normalize Whether to normalize the weight matrix such that each row adds up to one. Default is \code{TRUE}.
#' @param return.permutation Return permutations. Default is \code{FALSE}.
#' @return A list containing the following:
#' \itemize{
#'   \item Morans.I, the Moran's I.
#'   \item p.val, permutation based p-value.
#'   \item return.permutation, permutation used if returned.
#' }
#' @export
#' @concept moransi
#' @importFrom stats sd
#' @examples \dontrun{
#' res <- PermutationMoransI(x, W)
#' }
#'
PermutationMoransI <- function(x, W, n = 999, seed = 1, normalize = TRUE, return.permutation = FALSE){
  x <- x - mean(x)
  N <- length(x)
  if(normalize){
    W.rowsum <- rowSums(W)
    W.rowsum[W.rowsum==0] <- 1
    W <- W/W.rowsum
  } else {
    W.sum <- sum(W)
    W.sum[W.sum==0] <- 1
    W <- N/W.sum*W
  }
  y <- (W%*%x)[,1]
  sd.x <- sd(x)
  sd.y <- sd(y)
  r <- PermutationCorr(x, y, n, seed, return.permutation)
  if(return.permutation) return.permutation = r[[3]]*sd.y/sd.x
  return(list(Morans.I = r[[1]]*sd.y/sd.x, p.val = r[[2]], return.permutation = return.permutation))
}

#' Bivariate Moran's I
#'
#' Calculate bivariate (multivariate) Moran's I using Wartenberg's method.
#'
#' @param X A matrix with observations as rows and features as columns.
#' @param W A weight matrix across all observations, i.e inverse of a pairwise distance matrix.
#' @param alternative Alternative hypothesis used, default is \code{two.sided}.
#' @param p.adjust.method Method used for multiple comparisons correction, default is \code{BH}. See \code{\link[stats]{p.adjust}}.
#' @return A list containing the following:
#' \itemize{
#'   \item Morans.I, the Moran's I.
#'   \item Z.I, the Z score of Moran's I.
#'   \item Expected.I, the expectation of Moran's I under the null hypothesis.
#'   \item SD.I, the standard deviation of Moran's I under the null hypothesis.
#'   \item alternative, alternative hypothesis used.
#'   \item p.adjust.method, method used for multiple comparisons correction.
#' }
#' @export
#' @concept moransi
#' @examples \dontrun{
#' BivariateMoransI(X, W)
#' }
#' @references
#' Wartenberg, D. Multivariate spatial correlation:
#' A method for exploratory geographical analysis. Geogr. Anal. 17, 263–283 (1985)
#'
#' Czaplewski, R. L. Expected Value and Variance of Moran’s Bivariate Spatial Autocorrelation Statistic for a Permutation Test.
#' (U.S. Department of Agriculture, Forest Service, Rocky Mountain Forest and Range Experiment Station, 1993).
#'
BivariateMoransI <- function(X, W, alternative = c("two.sided", "less", "greater"), p.adjust.method = "BH"){
  alternative <- match.arg(arg = alternative)
  I <- WartenbergM(X, W)
  I.stats <- BivariateMoransIStats(X, W)
  Z <- (I-I.stats[[1]])/I.stats[[2]]
  P.list <- ZPvalue(Z, alternative, p.adjust.method)
  return(list(Morans.I = I, Z.I = Z, Expected.I = I.stats[[1]], SD.I = I.stats[[2]],
              p.val = P.list[[1]], p.adj = P.list[[2]], alternative = alternative, p.adjust.method = p.adjust.method))
}

#' Moran's I using ordinary least squares (OLS)
#'
#' Calculate Moran's I as linear regression using ordinary least squares (OLS).
#'
#' @param X A matrix with observations as rows and features as columns.
#' @param W A weight matrix across all observations, i.e inverse of a pairwise distance matrix.
#' @param normalize Whether to normalize the weight matrix such that each row adds up to one. Default is \code{TRUE}.
#' @param alternative Alternative hypothesis used, default is \code{two.sided}.
#' @param p.adjust.method Method used for multiple comparisons correction, default is \code{BH}. See \code{\link[stats]{p.adjust}}.
#' @return A list containing the following:
#' \itemize{
#'   \item Morans.I, the Moran's I.
#'   \item Z.I, the Z score of Moran's I.
#'   \item X, data matrix used for calculating Moran's I.
#'   \item Y, a matrix of spatial lags.
#'   \item Expected.I, the expectation of Moran's I under the null hypothesis.
#'   \item SD.I, the standard deviation of Moran's I under the null hypothesis.
#'   \item p.val, p-values.
#'   \item p.adj, adjusted p-values.
#'   \item alternative, alternative hypothesis used.
#'   \item p.adjust.method, method used for multiple comparisons correction.
#' }
#' @export
#' @concept moransi
#' @examples \dontrun{
#' OLSMoransI(X, W)
#' }
#' @references
#' Anselin, L. Local indicators of spatial association-LISA. Geogr. Anal. 27, 93–115 (2010)
#'
OLSMoransI <- function(X, W, normalize = TRUE, alternative = c("two.sided", "less", "greater"), p.adjust.method = "BH"){
  alternative <- match.arg(arg = alternative)
  X <- scale(X, center = TRUE, scale = FALSE)
  N <- nrow(X)
  if(normalize){
    W.rowsum <- rowSums(W)
    W.rowsum[W.rowsum==0] <- 1
    W <- W/W.rowsum
  } else {
    W.sum <- sum(W)
    W.sum[W.sum==0] <- 1
    W <- N/W.sum*W
  }
  Y <- W%*%X
  if(is(Y,"denseMatrix")) Y <- as.matrix(Y)
  XY <- crossprod(X,Y)
  X2 <- colSums(X*X)
  I <- XY/X2
  S <- sqrt(tcrossprod(1/X2,X2))
  I.stats <- BivariateMoransIStats(X,W)
  I.stats[[1]] <- I.stats[[1]]*S
  I.stats[[2]] <- I.stats[[2]]*S
  Z <- (I-I.stats[[1]])/I.stats[[2]]
  P.list <- ZPvalue(Z, alternative, p.adjust.method)
  return(list(Morans.I = I, Z.I = Z, X = X, Y = Y, Expected.I = I.stats[[1]], SD.I = I.stats[[2]],
              p.val = P.list[[1]], p.adj = P.list[[2]], alternative = alternative, p.adjust.method = p.adjust.method))
}

#' Fitted values of Moran's I
#'
#' Calculate fitted values of Moran's I as linear regression.
#'
#' @param X A matrix with observations as rows and features as columns.
#' @param Y A matrix of spatial lags.
#' @param I A matrix of Moran's I, output by \code{\link[rspca]{OLSMoransI}}..
#' @param all Whether to use multivariate Moran's I, default is \code{FALSE}.
#' @return A list containing the following:
#' \itemize{
#'   \item residuals, residuals of fitted Moran's I.
#'   \item fitted.values, fitted values of Moran's I.
#'   \item all.fitted.values, a three-dimensional tensor of all fitted multivariate Moran's I, if calculated.
#'   \item all.residuals, a three-dimensional tensor of all residuals using fitted multivariate Moran's I, if calculated.
#' }
#' @export
#' @concept moransi
#' @examples \dontrun{
#' res <- OLSMoransI(X, W)
#' FitMoransI(X, res$Y, res$Morans.I)
#' }
#' @references
#' Anselin, L. Local indicators of spatial association-LISA. Geogr. Anal. 27, 93–115 (2010)
#'
FitMoransI<- function(X, Y, I, all = FALSE){
  N <- nrow(X)
  M <- ncol(X)
  names.dim <- dimnames(X)
  if(all){
    All.X <- array(rep(X, M), dim = c(N, M, M))
    All.Y.hat <- array(0, dim = c(N, M, M), dimnames = list(names.dim[[1]], names.dim[[2]], paste0(names.dim[[2]],".y")))
    All.R <- array(0, dim = c(N, M, M), dimnames = list(names.dim[[1]], names.dim[[2]], paste0(names.dim[[2]],".y")))
    Y.hat <- matrix(0, nrow = N, ncol = M, dimnames = names.dim)
    R <- matrix(0, nrow = N, ncol = M, dimnames = names.dim)
    for(i in 1:M){
      All.Y.hat[,,i] <- t(t(All.X[,,i])*I[,i])
      All.R[,,i] <- Y[,i]-All.Y.hat[,,i]
      Y.hat[,i] <- All.Y.hat[,i,i]
      R[,i] <- All.R[,i,i]
    }
    list(residuals = R, fitted.values = Y.hat, all.fitted.values = All.Y.hat, all.residuals = All.R)
  } else {
    Y.hat1 <- t(t(X)*diag(I))
    R1 <- Y-Y.hat
    list(residuals = R, fitted.values = Y.hat)
  }
}

#' Process irregular data shape (between two groups) for Moran's I
#'
#' Process irregular data shape for Moran's I between two groups.
#'
#' @param X A matrix with observations as rows and features as columns.
#' @param W A weight matrix across all observations, i.e inverse of a pairwise distance matrix.
#' @param group1 The indices or names for the first group of observations.
#' @param group2 The indices or names for the second group of observations.
#' @return A list containing the following:
#' \itemize{
#'   \item X, new data matrix of group1 and group2.
#'   \item W, new weight matrix of group1 and group2.
#' }
#' @export
#' @concept moransi
#' @examples \dontrun{
#' IrregularData(X, W, 1:10, 101:110)
#' }
#'
IrregularData <- function(X, W, group1, group2){
  X <- X[c(group1, group2),]
  X <- X[,colSums(X)!=0]
  W[group1, group1] <- 0
  W[group2, group2] <- 0
  W <- W[c(group1, group2), c(group1, group2)]
  return(list(X = X, W = W))
}

#' Moran's I between two groups (irregular data shape)
#'
#' Calculate Moran's I between two groups with irregular data shape.
#'
#' @param X A matrix with observations as rows and features as columns.
#' @param W A weight matrix across all observations, i.e inverse of a pairwise distance matrix.
#' @param group1 The indices or names for the first group of observations.
#' @param group2 The indices or names for the second group of observations.
#' @param OLS Whether to use \code{\link[rspca]{OLSMoransI}}, default is \code{FALSE} and \code{\link[rspca]{BivariateMoransI}} is used.
#' @param normalize Whether to normalize the weight matrix such that each row adds up to one. Default is \code{TRUE}.
#' @param alternative Alternative hypothesis used, default is \code{two.sided}.
#' @param p.adjust.method Method used for multiple comparisons correction, default is \code{BH}. See \code{\link[stats]{p.adjust}}.
#' @return A list containing the output from \code{\link[rspca]{BivariateMoransI}} or \code{\link[rspca]{OLSMoransI}}.
#' @export
#' @concept moransi
#' @examples \dontrun{
#' IrregularMoransI(X, W, 1:10, 101:110)
#' }
#'
IrregularMoransI <- function(X, W, group1, group2, OLS = FALSE, normalize = TRUE, alternative = c("two.sided", "less", "greater"), p.adjust.method = "BH"){
  alternative <- match.arg(alternative)
  new.data <- IrregularData(X, W, group1, group2)
  if(OLS){
    OLSMoransI(new.data[[1]], new.data[[2]], normalize, alternative, p.adjust.method)
  } else {
    BivariateMoransI(new.data[[1]], new.data[[2]], alternative, p.adjust.method)
  }
}
