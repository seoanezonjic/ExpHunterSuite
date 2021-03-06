% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functional_analysis_library.R
\name{enrichment_ORA}
\alias{enrichment_ORA}
\title{Performs Over Representation Analysis (ORA)
enrichment of specified ontology}
\usage{
enrichment_ORA(
  genes,
  organism,
  keyType = "ENTREZID",
  pvalueCutoff,
  pAdjustMethod = "BH",
  ont,
  useInternal = FALSE,
  qvalueCutoff = 0.2,
  ENRICH_DATA = NULL,
  semsim = TRUE
)
}
\arguments{
\item{genes}{significant genes to be used}

\item{organism}{target organism}

\item{keyType}{genes code type}

\item{pvalueCutoff}{p-value threshold}

\item{pAdjustMethod}{p-valued adjust method to be applied}

\item{ont}{ontology to be used. Allowed (GO_MF, GO_CC, GO_BP, KEGG, REACT)}

\item{useInternal}{used only for KEGG enrichment, activate internal data
usage mode}

\item{qvalueCutoff}{q-value threshold}

\item{ENRICH_DATA}{optional enrichment universe already loaded}

\item{semsim}{flag to indicate if semantic similitud must be calculated.
Necessary for emaplots}
}
\value{
enrichment table obtained
}
\description{
Performs Over Representation Analysis (ORA)
enrichment of specified ontology
}
\keyword{enrich}
