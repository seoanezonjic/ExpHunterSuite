#' Table of counts dummy dataset for DEgeneHunter examples
#'
#' @docType data
#'
#' @usage data(toc)
#'
#' @format A truncated data.frame of aligned read counts for genes (rows) and 
#' samples (columns)
#'
#' @keywords datasets
#'
#' @examples
#' data(toc)
"toc"

#' Target data containing sample information for the toc data
#'
#' @docType data
#'
#' @usage data(target)
#'
#' @format A data.frame of sample information. Rows are samples, columns are 
#' different experimental variables/factors
#'
#' @keywords datasets
#'
#' @examples
#' data(target)
"target"

#' Output data from the main_degenes_Hunter function
#'
#' @docType data
#'
#' @usage data(degh_output)
#'
#' @format An output object from the main_degenes_Hunter function. 
#'
#' @keywords datasets
#'
#' @examples
#' data(degh_output)
"degh_output"

#' Subset of the pbmc dataset distributed in package SeuratObject
#'
#' @docType data
#'
#' @usage data(pbmc_tiny)
#'
#' @format A seurat object containing the 15 first samples of the
#' pbmc_small dataset, subset to 10 features.
#'
#' @keywords datasets
#'
#' @source Created with data-raw/pbmc_tiny.R. pbmc_small source:
#' https://search.r-project.org/CRAN/refmans/SeuratObject/html/pbmc_small.html
#'
"pbmc_tiny"