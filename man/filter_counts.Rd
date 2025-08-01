% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DROP_functions.R
\name{filter_counts}
\alias{filter_counts}
\title{Wrapper for DROP OUTRIDER filtering script.
\code{filter_counts} takes a counts table, a txdb object and a fpkm_cutoff to
build an OutriderDataSet object and filter it by fpkm.}
\usage{
filter_counts(counts, txdb, fpkm_cutoff)
}
\arguments{
\item{counts}{A counts table.}

\item{txdb}{Txdb database object.}

\item{fpkm_cutoff}{Min FPKM value for read filtering. Default 3.}
}
\value{
A filtered OutriderDataSet built from the merged table.
}
\description{
Wrapper for DROP OUTRIDER filtering script.
\code{filter_counts} takes a counts table, a txdb object and a fpkm_cutoff to
build an OutriderDataSet object and filter it by fpkm.
}
