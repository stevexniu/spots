#' @importFrom tools md5sum
#' @importFrom utils download.file
#'
NULL

#' Available data
#'
#' Available data can be loaded using \code{\link[spots]{LoadData}} function.
#'
#' @format A data frame with 2 rows (datasets) and 5 columns (attributes):
#' \describe{
#'   \item{filename}{Dataset}
#'   \item{md5}{MD5 checksum}
#'   \item{url}{URL link}
#'   \item{date}{Date updated}
#'   \item{size}{Data size}
#' }
#' @concept data
#'
"available.data"

#' Load selected data
#'
#' If not previously run, it will download the selected data.
#'
#' @param file.path Path to the directory contains or to download the data object, i.e. "~/Downloads".
#' @param data.use Selected data to load. See \code{\link[spots]{available.data}}.
#' \itemize{
#'   \item Visium.HNN, hexagonal nearest neighbor distance matrix of
#'   \href{https://www.10xgenomics.com/spatial-transcriptomics}{10x Visium} 4992 whitelist spatial barcodes.
#'   \item Visium.Brain, sample \href{https://www.10xgenomics.com/resources/datasets/mouse-brain-serial-section-1-sagittal-posterior-1-standard-1-0-0}{10x Visium Sagittal-Posterior Mouse Brain}
#'   data with 3353 spatial barcodes and 2000 highly variable genes, load as a \href{https://satijalab.org/seurat/articles/spatial_vignette.html}{Seurat object}.
#' }
#' @return Selected data:
#' @export
#' @concept data
#' @examples \dontrun{
#' Visium.hnn.dist <- LoadData("~/Downloads", "Visium.HNN")
#' }
#'
LoadData <- function(file.path, data.use = c("Visium.HNN", "Visium.Brain")){
  meta.data <- spots::available.data[match.arg(data.use, choices = rownames(spots::available.data)), ]
  file.path <- paste0(file.path, "/", meta.data[['filename']])
  file.md5 <- md5sum(file.path)
  if(is.na(file.md5) | file.md5 != meta.data[['md5']]){
    message("Downloading ", data.use)
    download.file(url = meta.data[['url']], destfile = file.path)
  }
  data.use <- readRDS(file.path)
  return(data.use)
}


