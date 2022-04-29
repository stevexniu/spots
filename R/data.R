#' @importFrom tools md5sum
#' @importFrom utils download.file
#'
NULL

#' Available data
#'
#' Available data can be loaded using \code{\link[spots]{LoadData}} function.
#'
#' @format A data frame with 2 rows (datasets) and 4 columns (attributes):
#' \describe{
#'   \item{filename}{Dataset}
#'   \item{md5}{MD5 checksum}
#'   \item{url}{URL link}
#'   \item{date}{Date updated}
#' }
#' @concept data
#'
"available.data"
