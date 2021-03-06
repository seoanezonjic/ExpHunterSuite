% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functional_analysis_library.R
\name{enrich_GSEA}
\alias{enrich_GSEA}
\title{Performs Gene Set Enrichment Analysis (GSEA) using a specified ontology}
\usage{
enrich_GSEA(
  geneList,
  organism,
  keyType = "ENTREZID",
  pvalueCutoff,
  pAdjustMethod = "BH",
  ont,
  useInternal = FALSE
)
}
\arguments{
\item{geneList}{significant genes to be used}

\item{organism}{target organism}

\item{keyType}{gene code type}

\item{pvalueCutoff}{p-value threshold}

\item{pAdjustMethod}{p-value adjust method to be applied}

\item{ont}{ontology to be used. Allowed (GO_MF, GO_CC, GO_BP, KEGG, REACT)}

\item{useInternal}{used only for KEGG enrichment, activate internal data
usage mode}
}
\value{
enrichment performed
}
\description{
Performs Gene Set Enrichment Analysis (GSEA) using a specified ontology
}
\keyword{enrich}
