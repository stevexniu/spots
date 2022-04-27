#' @include stats.R
#'
NULL

#' Local Moran's I
#'
#' Calculate Local Moran's I.
#'
#' @param X A matrix with observations as rows and features as columns.
#' @param W A weight matrix across all observations, i.e inverse of a pairwise distance matrix.
#' @param normalize Whether to normalize the weight matrix such that each row adds up to one. Default is \code{TRUE}.
#' @param alternative Alternative hypothesis used, default is \code{two.sided}.
#' @param p.adjust.method Method used for multiple comparisons correction, default is \code{BH}. See \code{\link[stats]{p.adjust}}.
#' @param scale Whether to scale the data. Default is \code{TRUE}.
#' @return A list containing the following:
#' \itemize{
#'   \item Local.Morans.I, the local Moran's I.
#'   \item Z.I, the Z score of local Moran's I.
#'   \item X, data matrix used for calculating local Moran's I.
#'   \item Y, a matrix of spatial lags.
#'   \item Expected.I, the expectation of Moran's I under the null hypothesis.
#'   \item SD.I, the standard deviation of Moran's I under the null hypothesis.
#'   \item p.val, p-values.
#'   \item p.adj, adjusted p-values.
#'   \item normalize, whether to normalize the weight matrix.
#'   \item scale, whether to scale the data.
#'   \item scale.factor, number of observations.
#'   \item alternative, alternative hypothesis used.
#'   \item p.adjust.method, method used for multiple comparisons correction.
#' }
#' @export
#' @concept local
#' @examples \dontrun{
#' LocalMoransI(X, W)
#' }
#' @references
#' Anselin, L. Local indicators of spatial association-LISA. Geogr. Anal. 27, 93–115 (1995)
#' Sokal, R. R., Oden, N. L. & Thomson, B. A. Local spatial autocorrelation in a biological model. Geogr. Anal. 30, 331–354 (1998)
#'
LocalMoransI <- function(X, W, normalize = TRUE, alternative = c("two.sided", "less", "greater"), p.adjust.method = "BH", scale = TRUE){
  alternative <- match.arg(arg = alternative)
  X <- scale(X, center = TRUE, scale = FALSE)
  N <- nrow(X)
  if(normalize){
    W.rowsum <- rowSums(W)
    W.rowsum[W.rowsum==0] <- 1
    W <- W/W.rowsum
  }
  Y <- W%*%X
  if(is(Y,"denseMatrix")) Y <- as.matrix(Y)
  X2 <- colSums(X*X)
  I <- t(t(X*Y)/X2)
  if(!scale) I <- I*N
  I.stats <- LocalMoransIStats(X, W, scale)
  Z <- (I-I.stats[[1]])/I.stats[[2]]
  P.list <- ZPvalue(Z, alternative, p.adjust.method)
  return(list(Local.Morans.I = I, Z.I = Z, X = X, Y = Y, Expected.I = I.stats[[1]], SD.I = I.stats[[2]],
              p.val = P.list[[1]], p.adj = P.list[[2]], normalize = normalize, scale = scale, scale.factor = N,
              alternative = alternative, p.adjust.method = p.adjust.method))
}

#' Permutation based local Moran's I
#'
#' Calculate permutation based p-value for local Moran's I.
#'
#' @param x A numerical vector.
#' @param W A weight matrix across all observations, i.e inverse of a pairwise distance matrix.
#' @param n The number of permutations to be conducted, set to 999 by default.
#' @param seed Random seed used. Default is 1.
#' @param alternative Alternative hypothesis used, default is \code{two.sided}.
#' @param p.adjust.method Method used for multiple comparisons correction, default is \code{BH}. See \code{\link[stats]{p.adjust}}.
#' @param normalize Whether to normalize the weight matrix such that each row adds up to one. Default is \code{TRUE}.
#' @param scale Whether to scale the data. Default is \code{TRUE}.
#' @param return.permutation Return permutations. Default is \code{FALSE}.
#' @param condition Value under null hypothesis to compare with, default is 0.
#' @return A list containing the following:
#' \itemize{
#'   \item Local.Morans.I, local Moran's I.
#'   \item p.val, permutation based p-value.
#'   \item p.adj, adjusted p-values.
#'   \item scale.factor, number of observations.
#'   \item return.permutation, permutation used if returned.
#'   \item params, parameters used to calculate local Moran's I.
#' }
#' @export
#' @concept local
#' @examples \dontrun{
#' res <- PermutationLocalI(x, W)
#' }
#'
PermutationLocalI <- function(x, W, n = 999, seed = 1, alternative = c("two.sided", "less", "greater"), p.adjust.method = "BH",
                              normalize = TRUE, scale = TRUE, return.permutation = FALSE, condition = 0){
  alternative <- match.arg(arg = alternative)
  nx <- length(x)
  x <- x - mean(x)
  x2 <- sum(x^2)
  if(normalize){
    W.rowsum <- rowSums(W)
    W.rowsum[W.rowsum==0] <- 1
    W <- W/W.rowsum
  }
  i.list <- PermutationLag(x, W, n, seed, alternative, p.adjust.method, condition, return.permutation)
  i.list[[1]] <- x*i.list[[1]]/x2
  if(!scale) i.list[[1]] <- i.list[[1]]*nx
  return(list(Local.Morans.I = i.list[[1]], p.val = i.list[[2]], p.adj = i.list[[3]], return.permutation = i.list[[4]], scale.factor = nx,
              params = list(n = n, seed = seed, alternative = alternative, p.adjust.method = p.adjust.method, normalize = normalize, scale = scale)))
}

#' Local Getis-Ord Gi or Gi* statistics for spatial hot spot analysis
#'
#' Calculate local Getis-Ord Gi or Gi* statistics.
#'
#' @param X A matrix with observations as rows and features as columns.
#' @param W A weight matrix across all observations, i.e inverse of a pairwise distance matrix.
#' @param gstar Whether to calculate the Gi* statistics, default is \code{TRUE}.
#' @param alternative Alternative hypothesis used, default is \code{two.sided}.
#' @param p.adjust.method Method used for multiple comparisons correction, default is \code{BH}. See \code{\link[stats]{p.adjust}}.
#' @return A list containing the following:
#' \itemize{
#'   \item Gi, Gi or Gi* statistics.
#'   \item p.val, permutation based p-value.
#'   \item p.adj, adjusted p-values.
#'   \item gstar, permutation used if returned.
#'   \item alternative, alternative hypothesis used.
#'   \item p.adjust.method, method used for multiple comparisons correction.
#' }
#' @export
#' @concept local
#' @examples \dontrun{
#' res <- LocalGi(X, W)
#' }
#' @references
#' Getis, A. & Ord, J. K. The analysis of spatial association by use of distance statistics. Geogr. Anal. 24, 189–206 (1992)
#'
LocalGi <- function(X, W, gstar = TRUE, alternative = c("two.sided", "less", "greater"), p.adjust.method = "BH"){
  alternative <- match.arg(arg = alternative)
  if(!gstar) diag(W) <- 0
  N <- ifelse(gstar, nrow(X), nrow(X)-1)
  X <- scale(X, center = TRUE, scale = FALSE)
  W.rowsum <- rowSums(W)
  W.rowsum[W.rowsum==0] <- 1
  W <- W/W.rowsum
  G <- W%*%X
  if(is(G,"denseMatrix")) G <- as.matrix(G)
  W2 <- (N*rowSums(W^2)-1) / (N-1)
  V <- colSums(X^2)/N
  SD <- sqrt(tcrossprod(W2,V))
  G <- G/SD
  P.list <- ZPvalue(G, alternative, p.adjust.method)
  return(list(Gi = G, p.val = P.list[[1]], p.adj = P.list[[2]],
              gstar = gstar, alternative = alternative, p.adjust.method = p.adjust.method))
}

#' Permutation based Getis-Ord Gi or Gi* statistics for spatial hots pot analysis
#'
#' Calculate permutation based p-value for Getis-Ord Gi or Gi* statistics.
#'
#' @param x A numerical vector.
#' @param W A weight matrix across all observations, i.e inverse of a pairwise distance matrix.
#' @param n The number of permutations to be conducted, set to 999 by default.
#' @param gstar Whether to calculate the Gi* statistics, default is \code{TRUE}.
#' @param seed Random seed used. Default is 1.
#' @param alternative Alternative hypothesis used, default is \code{two.sided}.
#' @param p.adjust.method Method used for multiple comparisons correction, default is \code{BH}. See \code{\link[stats]{p.adjust}}.
#' @param return.permutation Return permutations. Default is \code{FALSE}.
#' @param condition Value under null hypothesis to compare with, default is 0.
#' @return A list containing the following:
#' \itemize{
#'   \item Gi, Gi or Gi* statistics.
#'   \item p.val, permutation based p-value.
#'   \item p.adj, adjusted p-values.
#'   \item return.permutation, permutation used if returned.
#'   \item params, parameters used to calculate local Moran's I.
#' }
#' @export
#' @concept local
#' @examples \dontrun{
#' res <- PermutationGi(x, W)
#' }
#'
PermutationGi <- function(x, W, gstar = TRUE, n = 999, seed = 1, alternative = c("two.sided", "less", "greater"), p.adjust.method = "BH",
                          return.permutation = FALSE, condition = 0){
  alternative <- match.arg(arg = alternative)
  if(!gstar) diag(W) <- 0
  N <- ifelse(gstar, length(x), length(x)-1)
  x <- x - mean(x)
  W.rowsum <- rowSums(W)
  W.rowsum[W.rowsum==0] <- 1
  W <- W/W.rowsum
  W2 <- (N*rowSums(W^2)-1) / (N-1)
  V <- sum(x^2)/N
  SD <- sqrt(V*W2)
  gi.list <- PermutationLag(x, W, n, seed, alternative, p.adjust.method, condition, return.permutation)
  gi.list[[1]] <- gi.list[[1]]/SD
  return(list(Gi = gi.list[[1]], p.val = gi.list[[2]], p.adj = gi.list[[3]], return.permutation = gi.list[[4]],
              params = list(gstar = gstar, n = n, seed = seed, alternative = alternative, p.adjust.method = p.adjust.method, condition = condition)))
}
