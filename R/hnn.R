#' @include stats.R
#'
NULL

#' Hexagonal nearest neighbor
#'
#' Calculate hexagonal nearest neighbors.
#'
#' @param dist.hnn A hexagonal nearest neighbor distance matrix.
#' @param k Number of neighbors.
#' @param include.self Whether to include self as 1st neighbor, default is \code{TRUE}.
#' @return A list containing the following:
#' \itemize{
#'   \item knn.idx, an n x k matrix for the nearest neighbor indice.
#'   \item knn.dist, an n x k matrix for the nearest neighbor hexagonal distances.
#'   \item dist.mat, a connectivity-based distance matrix.
#' }
#' @export
#' @concept hnn
#' @examples {
#' data.use <- quakes[1:100,]
#' dist.use <- as.matrix(dist(data.use[,1:2]))
#' res <- HnnNeighbor(dist.use, k = 10)
#' }
#' @references
#' Middleton, L. & Sivaswamy, J.
#' Edge detection in a hexagonal-image processing framework. Image Vis. Comput. 19, 1071–1081 (2001)
#'
HnnNeighbor <- function(dist.hnn, k, include.self = TRUE){
  n <- nrow(x = dist.hnn)
  k.start <- ifelse(test = include.self, yes = 1, no = 2)
  k.end <- ifelse(test = include.self, yes = k, no = k+1)
  knn.idx <- matrix(data = 0, ncol = k, nrow = n)
  knn.dist <- matrix(data = 0, ncol = k, nrow = n)
  dist.mat <- dist.hnn
  for(i in 1:n){
    knn.idx[i,] <- order(dist.hnn[i,])[k.start:k.end]
    knn.dist[i,] <- dist.hnn[i,knn.idx[i,]]
    dist.mat[i,-knn.idx[i,]] <- Inf
  }
  rownames(x = knn.idx) <- rownames(x = dist.hnn)
  return(list(knn.idx = knn.idx, knn.dist = knn.dist, dist.mat = dist.mat))
}

#' Hexagonal nearest neighbor index
#'
#' Get hexagonal nearest neighbor indices.
#'
#' @param dist.hnn A hexagonal nearest neighbor distance matrix.
#' @param k Number of neighbors.
#' @param include.self Whether to include self as 1st neighbor, default is \code{TRUE}.
#' @return An n x k matrix for the nearest neighbor indice.
#' @keywords internal
#'
HnnNeighborIndex <- function(dist.hnn, k, include.self = TRUE){
  k.start <- ifelse(test = include.self, yes = 1, no = 2)
  k.end <- ifelse(test = include.self, yes = k, no = k+1)
  knn.idx <- t(x = apply(X = dist.hnn, MARGIN = 1, FUN = function(x){
    order(x)[k.start:k.end]
  }))
  rownames(x = knn.idx) <- rownames(x = dist.hnn)
  return(knn.idx)
}

#' Hexagonal nearest neighbor weight
#'
#' Calculate hexagonal nearest neighbor weights using Gaussian filter.
#'
#' @param dist.hnn A hexagonal nearest neighbor distance matrix.
#' @param dist.k The maximum distance used to calculate the weight. Default is \code{NULL} and all neighbor weights are calculated.
#' @param mu The mean of Gaussian filter, default is 0.
#' @param sigma The standard deviation of Gaussian filter, default is 1.
#' @return A weight matrix.
#' @importFrom stats dnorm
#' @export
#' @concept hnn
#' @examples {
#' data.use <- quakes[1:100,]
#' dist.use <- as.matrix(dist(data.use[,1:2]))
#' res <- HnnWeight(dist.use, mu = 0, sigma = 0.5)
#' }
#'
HnnWeight <- function(dist.hnn, dist.k = NULL, mu = 0, sigma=1){
  if(!is.null(x = dist.k)) dist.hnn[dist.hnn>dist.k] <- Inf
  hnn.weight <- dnorm(x = dist.hnn, mean = mu, sd = sigma)
  return(hnn.weight)
}

#' Hexagonal nearest neighbor based imputation
#'
#' Data imputation and smoothing using hexagonal nearest neighbor.
#'
#' @param data A data matrix with features as rows and observations as columns.
#' @param dist.hnn A hexagonal nearest neighbor distance matrix.
#' @param dist.k The maximum distance used to calculate the weight. Default is \code{NULL} and all neighbor weights are calculated.
#' @param mu The mean of Gaussian filter, default is 0.
#' @param sigma The standard deviation of Gaussian filter, default is 1.
#' @return Imputed data.
#' @seealso \code{\link[spots]{HnnWeight}}
#' @export
#' @concept hnn
#' @examples {
#' data.use <- quakes[1:100,]
#' dist.use <- as.matrix(dist(data.use[,1:2]))
#' # transpose the data to have features in rows and observations in columns
#' res <- HnnImpute(t(data.use[,3:4]), dist.use)
#' }
#'
HnnImpute <- function(data, dist.hnn, dist.k = NULL, mu = 0, sigma=1){
  if(!is.null(x = dist.k)) dist.hnn[dist.hnn>dist.k] <- Inf
  hnn.weight <- HnnWeight(dist.hnn = dist.hnn, dist.k = dist.k, mu = mu, sigma = sigma)
  hnn.weight <- sweep(x = hnn.weight, MARGIN = 2, STATS = colSums(hnn.weight), FUN = "/")
  data <- data %*% hnn.weight
  return(data)
}

#' 10x Genomics Visium hexagonal nearest neighbor
#'
#' Get hexagonal nearest neighbor distance between 10x Visium spatial barcodes.
#'
#' @param path Path to downloaded HNN distance for Visium barcodes. See \code{\link[spots]{available.data}} for download url.
#' @param barcodes 10x Visium whitelisted spatial barcodes.
#' @return A hexagonal distance matrix for selected barcodes.
#' @seealso \code{\link[spots]{available.data}}.
#' @export
#' @concept hnn
#' @examples \dontrun{
#' barcodes <- c("AAACAACGAATAGTTC-1", "AAACTTAATTGCACGC-1")
#' VisiumHnn("~/Downloads/", barcodes)
#' }
#'
VisiumHnn <- function(path, barcodes){
  Visium.hnn.dist <- readRDS(path)
  dist.use <- as.matrix(Visium.hnn.dist[barcodes, barcodes])
  return(dist.use + t(dist.use))
}
