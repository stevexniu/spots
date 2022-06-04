#' Greater than 
#'
#' Greater than for floating point errors.
#'
#' @param x A numerical vector.
#' @param y A numerical vector to compare to.
#' @param e Rounding error allowed, default is 1e-10.
#' @return TRUE or FALSE.
#' @keywords internal
#'
`%>=%` <- function(x, y, e=1e-10) x + e > y

#' Smaller than 
#'
#' Smaller than for floating point errors.
#'
#' @param x A numerical vector.
#' @param y A numerical vector to compare to.
#' @param e Rounding error allowed, default is 1e-10.
#' @return TRUE or FALSE.
#' @keywords internal
#'
`%<=%` <- function(x, y, e=1e-10) x - e < y
