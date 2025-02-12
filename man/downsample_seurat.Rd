% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc_library.R
\name{downsample_seurat}
\alias{downsample_seurat}
\title{downsample_seurat
\code{downsample_seurat} takes a seurat object as input, and downsamples it to
specified number of cells and features. You can also input specific lists
if you know which cells and/or features you want to retrieve.}
\usage{
downsample_seurat(
  seu,
  cells = NULL,
  features = NULL,
  keep = "",
  assay = "RNA",
  layer = "counts"
)
}
\arguments{
\item{seu}{Seurat object.}

\item{cells, features}{Lists or integers. An integer will trigger random
downsampling by its respective variable.}

\item{keep}{A vector of genes to keep when downsampling features.}
}
\value{
A downsampled seurat object.
}
\description{
downsample_seurat
\code{downsample_seurat} takes a seurat object as input, and downsamples it to
specified number of cells and features. You can also input specific lists
if you know which cells and/or features you want to retrieve.
}
\examples{
data(pbmc_tiny)
print(pbmc_tiny)
print(downsample_seurat(seu = pbmc_tiny, cells = 2, features = 5))
}
