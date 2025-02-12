% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main_sc_Hunter.R
\name{main_sc_Hunter}
\alias{main_sc_Hunter}
\title{main_sc_Hunter
\code{main_sc_Hunter} is the main seurat analysis function. Can be used
for integrative or non-integrative analysis.}
\usage{
main_sc_Hunter(
  seu,
  minqcfeats,
  percentmt,
  query,
  sigfig = 2,
  resolution,
  p_adj_cutoff = 0.005,
  name = NULL,
  integrate = FALSE,
  cluster_annotation = NULL,
  cell_annotation = NULL,
  DEG_columns = NULL,
  scalefactor = 10000,
  hvgs,
  int_columns = NULL,
  normalmethod = "LogNormalize",
  ndims,
  verbose = FALSE,
  output = getwd(),
  save_RDS = FALSE,
  reduce = FALSE,
  ref_label,
  SingleR_ref = NULL,
  ref_de_method = NULL,
  ref_n = NULL,
  BPPARAM = NULL,
  doublet_list = NULL
)
}
\arguments{
\item{seu}{A seurat object.}

\item{minqcfeats}{An integer. Minimum features to consider a cell valid}

\item{percentmt}{A float. Maximun MT percentage to consider a cell valid}

\item{query}{A string vector. List of genes to explore in dataset}

\item{sigfig}{An integer. Significant figures to output}

\item{resolution}{An integer. Controls clustering granularity}

\item{p_adj_cutoff}{A float. Adjusted p-value cutoff by which to consider a maker
valid for cell type annotation}

\item{name}{Project name. Default NULL (no project name)}

\item{integrate}{A boolean.
\itemize{
\item \code{TRUE}: Integrate seurat object using the \code{harmony} package.
\item \code{FALSE} (the default): Do not perform integration.
}}

\item{cluster_annotation}{A data frame. Table to use to rename clusters.}

\item{cell_annotation}{A data frame. Table of markers to use for cell type
annotation}

\item{DEG_columns}{A string vector. Categories by which DEG analysis will be
performed}

\item{scalefactor}{An integer. Factor by which to scale data in normalisation}

\item{hvgs}{An integer. Number of highly-variable features to select}

\item{int_columns}{A string vector. Categories to consider in integrative
analysis.}

\item{normalmethod}{A string. Method to use in normalisation. Default is
"LogNormalize", the Seurat default.}

\item{ndims}{An integer. Target dimensions in dimensionality reduction.}

\item{output}{A string. Path to output directory.}

\item{save_RDS}{A boolean.
\itemize{
\item \code{TRUE}: Save output as an rds file.
\item \code{FALSE} (the default): Do not save output as an rds file.
}}

\item{reduce}{A boolean.
\itemize{
\item \code{TRUE}: Skip QC filtering. Intended for development and testing.
\item \code{FALSE} (the default): QC filtering will be performed.
}}

\item{ref_label}{Reference column to use in annotation.}

\item{SingleR_ref}{SummarizedExperiment object to use as reference
for SingleR cell type annotation. If NULL (default value),
SingleR will not be used.}

\item{ref_de_method}{Method to use for marker calculation in single-cell
reference.}

\item{ref_n}{Top N reference markers to consider in annotation. Higher values
provide a more accurate annotation, but increase noise and computational
time. Will not be used if ref_de_method is empty.}

\item{doublet_list}{A vector containing barcodes to be marked as doublets.
NULL by default. Per-sample analysis finds this vector
for every sample, integrative mode requires this vector.}

\item{verbose.}{A boolean.
\itemize{
\item \code{TRUE}: Prints progress bars and messages to closely monitor progress.
\item \code{FALSE} (the default): Print minimal progress messages.
}}
}
\description{
main_sc_Hunter
\code{main_sc_Hunter} is the main seurat analysis function. Can be used
for integrative or non-integrative analysis.
}
