#' @importFrom tools md5sum
#' @importFrom utils download.file
#'
NULL

#' Available data
#'
#' Available data can be loaded using \code{\link[rspca]{LoadData}} function.
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

#' Load selected data
#'
#' If not previously run, it will download the selected data.
#'
#' @param data.use Selected data to load. See \code{\link[rspca]{available.data}}.
#' \itemize{
#'   \item Visium.HNN, hexagonal nearest neighbor distance matrix of
#'   \href{https://www.10xgenomics.com/spatial-transcriptomics}{10x Visium} 4992 whitelist spatial barcodes.
#'   \item Visium.Brain, sample \href{https://support.10xgenomics.com/spatial-gene-expression/datasets}{10x Visium Sagittal-Posterior Mouse Brain}
#'   data with 3353 spatial barcodes and 2000 highly variable genes, load as a \href{Seurat object}{https://satijalab.org/seurat/articles/spatial_vignette.html}.
#' }
#' @return Selected data:
#' @export
#' @concept data
#' @examples \dontrun{
#' Visium.hnn.dist <- LoadData("Visium.HNN")
#' }
#'
LoadData <- function(data.use = c("Visium.HNN", "Visium.Brain")){
  meta.data <- rspca::available.data[match.arg(data.use, choices = rownames(rspca::available.data)), ]
  file.path <- paste0(system.file("extdata", package = "rspca"), "/", meta.data[['filename']])
  file.md5 <- md5sum(file.path)
  if(is.na(file.md5) | file.md5 != meta.data[['md5']]){
    message("Downloading ", data.use)
    download.file(url = meta.data[['url']], destfile = file.path)
  }
  data.use <- readRDS(file.path)
  return(data.use)
}
