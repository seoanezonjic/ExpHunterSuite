% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc_library.R
\name{process_doublets}
\alias{process_doublets}
\title{process_doublets}
\usage{
process_doublets(
  seu,
  name = NULL,
  doublet_path = getwd(),
  method = "scDblFinder",
  assay = "counts",
  nfeatures = 1352,
  includePCs = 19,
  BPPARAM = NULL
)
}
\arguments{
\item{seu}{Seurat object to analyze.}

\item{name}{Sample name.}

\item{doublet_path}{A string. Path where cell IDs identified as doublets
will be written.}

\item{method}{A string specifying which doublet processing Library to use.
Values: "scDblFinder" (the default) or "DoubletFinder"}

\item{assay}{A string. Seurat assay to process. Default "counts"
(recommended).}

\item{nfeatures}{The number of top features to use. Alternatively, a 
character vectors of feature names (e.g. highly-variable genes) to use.}

\item{includePCs}{The index of principal components to include in the
predictors (e.g. `includePCs=1:2`), or the number of top components to use
(e.g. `includePCs=10`, equivalent to 1:10).}

\item{BPPARAM}{Used for multithreading when splitting by samples (i.e. when
`samples!=NULL`); otherwise passed to eventual PCA and K/SNN calculations.}
}
\value{
A list containing the updated seurat object and tagged qc object.
}
\description{
\code{process_doublets} simplifies doublet detection and tagging process, as well
as removes doublets from seurat object.
}
\examples{
 \dontrun{
   process_doublets(seu = sample_seu, name = "sample1", qc = qc)
 }
}
