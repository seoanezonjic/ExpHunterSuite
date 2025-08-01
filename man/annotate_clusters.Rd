% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc_library.R
\name{annotate_clusters}
\alias{annotate_clusters}
\title{annotate_clusters}
\usage{
annotate_clusters(
  seu,
  subset_by = NULL,
  cell_annotation,
  assay = "RNA",
  integrate = FALSE,
  verbose = FALSE,
  idents = "seurat_clusters",
  p_adj_cutoff = 1e-05
)
}
\arguments{
\item{seu}{Seurat object to analyze.}

\item{subset_by}{Metadata column by which seurat object will be subset for
marker calculation. Default NULL.}

\item{cell_annotation}{Table of cell types and their associated markers}

\item{assay}{A string. Assay whose markers will be calculated. Default "RNA",
as per usual workflow.}

\item{integrate}{Whether or not integrative analysis is active. Does not
bother with subsetting if FALSE.}

\item{verbose}{A boolean. Will be passed to Seurat function calls.}

\item{idents}{Identity class to which to set seurat object before calculating
markers in case conserved mode cannot be triggered.}

\item{p_adj_cutoff}{Minimum adjusted p-value of markers to consider.}
}
\value{
Annotated seurat object.
}
\description{
\code{annotate_clusters} is a function which implements our cell type annotation
algorithm.
}
\examples{
 \dontrun{
   annotate_clusters(seu = seu, subset_by = "genotype", cell_annotation =
     cell_types)
 }
}
