% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc_library.R
\name{subset_seurat}
\alias{subset_seurat}
\title{subset_seurat}
\usage{
subset_seurat(seu, column, value, expr = FALSE, layer = "data")
}
\arguments{
\item{seu}{Seurat object.}

\item{column}{Column with value by which to subset input.}

\item{value}{Value to search within column.}

\item{expr}{Directly subset expression matrix instead of seurat object.}

\item{layer}{Seurat object layer to subset.}
}
\value{
A subset of the seurat object, which itself is a seurat object or an
expression matrix if expr is set to TRUE.
}
\description{
\code{subset_seurat} subsets a seurat object by a specified value of provided
column.
}
\examples{
data(pbmc_tiny)
head(subset_seurat(seu = pbmc_tiny, column = "groups", value = "g2",
expr = FALSE))
}
