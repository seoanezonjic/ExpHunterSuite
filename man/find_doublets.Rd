% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc_library.R
\name{find_doublets}
\alias{find_doublets}
\title{find_doublets
\code{find_doublets} is a wrapper for the recommended steps for doublet
calculation in package DoubletFinder}
\usage{
find_doublets(seu)
}
\arguments{
\item{seu}{Seurat object}
}
\value{
A list. Item "seu" contains seurat object with tagged doublets
in metadata slot. Item "barcodes" contains a vector of cell barcodes
corresponding to barcodes.
}
\description{
find_doublets
\code{find_doublets} is a wrapper for the recommended steps for doublet
calculation in package DoubletFinder
}
\examples{
 \dontrun{
   find_doublets(seu = pbmc_tiny)
 }
}
