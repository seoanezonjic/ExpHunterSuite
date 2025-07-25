% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{get_genes}
\alias{get_genes}
\title{get_genes}
\usage{
get_genes(enrich_obj, showCategory = 30)
}
\arguments{
\item{enrich_obj}{Enrichment object.}

\item{showCategory}{showCategory.}
}
\value{
An integer. Number of traits in object.
}
\description{
\code{get_genes} Extracts the vector of unique genes from an enrichment object.
}
\examples{
 \dontrun{
   get_genes(enrich_obj, showCategory = 30)
 }
}
