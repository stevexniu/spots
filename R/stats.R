#' @import Matrix
#' @importFrom stats cor p.adjust pnorm
#' @importFrom methods is
#'
NULL

#' Bivariate Moran's I statistics
#'
#' Calculate theoretical expectation and standard deviation of bivariate Moran's I under the null hypothesis.
#'
#' @param X A matrix with observations as rows and features as columns.
#' @param W A weight matrix across all observations, i.e inverse of a pairwise distance matrix.
#' @return A list containing the following:
#' \itemize{
#'   \item E.I, the expectation of Moran's I under the null hypothesis.
#'   \item SD.I, the standard deviation of Moran's I under the null hypothesis.
#' }
#' @concept stats
#' @references
#' Czaplewski, R. L. Expected Value and Variance of Moran’s Bivariate Spatial Autocorrelation Statistic for a Permutation Test.
#' (U.S. Department of Agriculture, Forest Service, Rocky Mountain Forest and Range Experiment Station, 1993).
#'
BivariateMoransIStats <- function(X, W){
  N <- nrow(X)
  O <- matrix(1, nrow = N)
  OT <- matrix(1, ncol = N)
  WT <- t(W)
  W2 <- (OT%*%(W%*%O))[1]^2
  if(isSymmetric(W)){
    S1 <- (2*OT%*%((W*W)%*%O))[1]
    S2 <- (4*OT%*%W%*%W%*%O)[1]
    S3 <- S4 <- S1/2
    S5 <- S2/4
    S6 <- S2/2
  } else {
    S3 <- (OT%*%((W*WT)%*%O))[1]
    S4 <- (OT%*%((W*W)%*%O))[1]
    S5 <- (OT%*%(W%*%(W%*%O)))[1]
    S6 <- (t(W%*%O)%*%(W%*%O)+t(WT%*%O)%*%(WT%*%O))[1]
    S1 <- S4+S3
    S2 <- 2*S5+S6
  }
  X <- scale(X, center = T,scale = F)
  EI <- -cor(X)/(N-1)
  COV <- crossprod(X)/N
  VAR <- tcrossprod(diag(COV))
  X2 <- crossprod(X*X)/N
  COV <- COV*COV
  VI <- (((COV*N) * (2*(W2-S2+S1)+(2*S3-2*S5)*(N-3)+S3*(N-2)*(N-3)) -
            (X2 * (6*(W2-S2+S1)+(4*S1-2*S2)*(N-3)+S1*(N-2)*(N-3))))/VAR +
           N * (W2-S2+S1+(2*S4-S6)*(N-3)+S4*(N-2)*(N-3))) /
    ((N-1)*(N-2)*(N-3)*(W2)) - (COV/VAR/(N-1)^2)
  return(list(E.I = EI, SD.I = sqrt(VI)))
}

#' Local Moran's I statistics
#'
#' Calculate theoretical expectation and standard deviation of local Moran's I under the null hypothesis.
#'
#' @param X A matrix with observations as rows and features as columns.
#' @param W A weight matrix across all observations, i.e inverse of a pairwise distance matrix.
#' @param scale Whether to scale the data. Default is \code{FALSE}.
#' @return A list containing the following:
#' \itemize{
#'   \item E.I, the expectation of local Moran's I under the null hypothesis.
#'   \item SD.I, the standard deviation of local Moran's I under the null hypothesis.
#' }
#' @concept stats
#' @references
#' Anselin, L. Local indicators of spatial association-LISA. Geogr. Anal. 27, 93–115 (1995)
#'
#' Sokal, R. R., Oden, N. L. & Thomson, B. A. Local spatial autocorrelation in a biological model. Geogr. Anal. 30, 331–354 (1998)
#'
LocalMoransIStats <- function(X, W, scale = FALSE){
  N <- nrow(X)
  WI <- rowSums(W)
  EI <- -WI/(N-1)
  WI2 <- rowSums(W^2)
  M2 <- colSums(X^2)/N
  M4 <- colSums(X^4)/N
  B2 <- M4/M2^2
  A <- tcrossprod(WI2, (N-B2)/(N-1))
  B <- tcrossprod((WI^2-WI2), (2*B2-N)/((N-1)*(N-2)))
  SDI <- sqrt(A+B-EI^2)
  dimnames(SDI) <- dimnames(X)
  if(scale){
    list(E.I = EI/N, SD.I = SDI/N)
  } else {
    list(E.I = EI, SD.I = SDI)
  }
}

#' P-value of Z scores
#'
#' Calculate p-value of Z scores.
#'
#' @param Z A matrix of Z scores.
#' @param alternative Alternative hypothesis used, default is \code{two.sided}.
#' @param p.adjust.method Method used for multiple comparisons correction, default is \code{BH}. See \code{\link[stats]{p.adjust}}.
#' @return A list containing the following:
#' \itemize{
#'   \item p.val, a matrix of p-values.
#'   \item p.adj, a matrix of adjusted p-values.
#' }
#' @concept stats
#'
ZPvalue <- function(Z, alternative = c("two.sided", "less", "greater"), p.adjust.method = "BH"){
  if (alternative == "two.sided") {
    p.val <- 2*pnorm(abs(Z), lower.tail = FALSE)
    p.adj <- p.adjust(p = p.val, method = p.adjust.method)
  } else if (alternative == "greater") {
    p.val <- pnorm(Z, lower.tail = FALSE)
    p.adj <- p.val
    p.adj[Z>=0] <- p.adjust(p.adj[Z>=0], method = p.adjust.method)
  } else {
    p.val <- pnorm(Z)
    p.adj <- p.val
    p.adj[Z<=0] <- p.adjust(p.adj[Z<=0], method = p.adjust.method)
  }
  p.adj <- matrix(p.adj, nrow = nrow(p.val), dimnames = dimnames(p.val))
  return(list(p.val = p.val, p.adj = p.adj))
}

#' Permutation based p-value for Pearson's correlation coefficient
#'
#' Calculate permutation based p-value for Pearson's correlation coefficient.
#'
#' It's adapted from \code{permcor} function in \href{https://cran.r-project.org/package=Rfast}{Rfast} with corrections in calculating the number of permutations and the p-values.
#'
#' @param x A numerical vector with the first variable.
#' @param y A numerical vector with the second variable.
#' @param R The number of permutations to be conducted, set to 999 by default.
#' @param seed Random seed used. Default is 1.
#' @param return.permutation Return permutations. Default is \code{FALSE}.
#' @return A list containing the following:
#' \itemize{
#'   \item correlation, Pearson's correlation coefficient.
#'   \item p.val, permutation based p-value.
#'   \item return.permutation, permutation used if returned.
#' }
#' @export
#' @concept stats
#' @examples \dontrun{
#' x <- iris[, 1]
#' y <- iris[, 2]
#' res <- PermutationCorr(x, y, R = 9999)
#' }
#'
PermutationCorr <- function(x, y, R = 999, seed = 1, return.permutation = FALSE){
  n <- length(x)
  m1 <- sum(x)
  m12 <- sum(x^2)
  m2 <- sum(y)
  m22 <- sum(y^2)
  up <-  m1*m2/n
  down <- sqrt( (m12-m1^2/n) * (m22-m2^2/n))
  r <- (sum(x*y)-up) / down
  test <- log((1+r) / (1-r))
  B <- ceiling(sqrt(R)) # change from round to ceiling
  xp <- matrix(0, n, B)
  yp <- matrix(0, n, B)
  set.seed(seed)
  for(i in 1:B){
    xp[,i] <- sample(x, n)
    yp[,i] <- sample(y, n)
  }
  sxy <- crossprod(xp, yp)
  rb <- (sxy-up)/down
  tb <- log((1+rb) / (1-rb))
  pvalue <- (sum(abs(tb)>=abs(test))+1) / (B^2+1) # change the '>' to '>=' for p-value calculation
  if(return.permutation) return.permutation <- c(rb)
  return(list(correlation = r, p.val = pvalue, return.permutation = return.permutation))
}

#' Permutation based p-value
#'
#' Calculate permutation based p-values.
#'
#' @param x A numerical vector.
#' @param X.perm A matrix contains permutation result of x.
#' @param alternative Alternative hypothesis used, default is \code{two.sided}.
#' @param p.adjust.method Method used for multiple comparisons correction, default is \code{BH}. See \code{\link[stats]{p.adjust}}.
#' @param condition Value under null hypothesis to compare with, default is 0.
#' @return A list containing the following:
#' \itemize{
#'   \item p.val, permutation based p-value.
#'   \item p.adj, adjusted p-values.
#' }
#' @keywords internal
#'
PermutationPval <- function(x, X.perm, alternative = c("two.sided", "less", "greater"), p.adjust.method = "BH", condition = 0){
  `%>=%` <- function(x, y, e=1e-10) x + e > y
  `%<=%` <- function(x, y, e=1e-10) x - e < y
  compare.func <- switch(alternative,
                         two.sided = `%>=%`,
                         less = `%<=%`,
                         greater = `%>=%`)
  n <- ncol(X.perm)
  if(alternative == "two.sided"){
    X.perm <- abs(X.perm)
    x <- abs(x)
    x.count <- sweep(X.perm, 1, x, compare.func)
    p.val <- (rowSums(x.count)+1) / (n+1)
    p.adj <- p.adjust(p = p.val, method = p.adjust.method)
  } else {
    x.count <- sweep(X.perm, 1, x, compare.func)
    p.val <- (rowSums(x.count)+1) / (n+1)
    p.adj <- p.val
    p.adj[compare.func(x, condition)] <- p.adjust(p.adj[compare.func(x, condition)], method = p.adjust.method)
  }
  return(list(p.val = p.val, p.adj = p.adj))
}

#' Permutation based p-value for spatial lag
#'
#' Calculate permutation based p-value for spatial lag.
#'
#' @param x A numerical vector.
#' @param W A weight matrix across all observations, i.e inverse of a pairwise distance matrix.
#' @param n The number of permutations to be conducted, set to 999 by default.
#' @param seed Random seed used. Default is 1.
#' @param alternative Alternative hypothesis used, default is \code{two.sided}.
#' @param p.adjust.method Method used for multiple comparisons correction, default is \code{BH}. See \code{\link[stats]{p.adjust}}.
#' @param condition Value under null hypothesis to compare with, default is 0.
#' @param return.permutation Return permutations. Default is \code{FALSE}.
#' @return A list containing the following:
#' \itemize{
#'   \item lag, spatial lag.
#'   \item p.val, permutation based p-value.
#'   \item p.adj, adjusted p-values.
#'   \item return.permutation, permutation used if returned.
#' }
#' @keywords internal
#'
PermutationLag <- function(x, W, n = 999, seed = 1, alternative = c("two.sided", "less", "greater"), p.adjust.method = "BH",
                           condition = 0, return.permutation = FALSE){
  y <- (W%*%x)[,1]
  X.perm <- matrix(0, nrow = length(x), ncol = n)
  set.seed(seed)
  for(i in 1:n) X.perm[,i] <- sample(x)
  Y.perm <- W%*%X.perm
  if(is(Y.perm,"denseMatrix")) Y.perm <- as.matrix(Y.perm)
  pval.list <- PermutationPval(y, Y.perm, alternative, p.adjust.method, condition)
  if(return.permutation) return.permutation <- Y.perm
  return(list(lag = y, p.val = pval.list[[1]], p.adj = pval.list[[2]], return.permutation = return.permutation))
}
