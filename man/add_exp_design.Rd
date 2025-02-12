% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc_library.R
\name{add_exp_design}
\alias{add_exp_design}
\title{add_exp_design
Add experimental condition to single-sample Seurat metadata}
\usage{
add_exp_design(seu, name, exp_design)
}
\arguments{
\item{seu}{Seurat object}

\item{name}{Sample name}

\item{exp_design}{Data frame containing experiment design}
}
\value{
Seurat object with the experimental conditions added as metadata
}
\description{
add_exp_design
Add experimental condition to single-sample Seurat metadata
}
\examples{
data(pbmc_tiny)
pbmc_tiny$sample <- "sampleA"
exp_design <- data.frame(sample = c("sampleA", "sampleB", "sampleC"),
                         new_field = c("info_A", "info_B", "info_C"))
add_exp_design(seu = pbmc_tiny, name = "sampleA", exp_design = exp_design)
}
\keyword{integration}
\keyword{preprocessing,}
\keyword{subsetting,}
