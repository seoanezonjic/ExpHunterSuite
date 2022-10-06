% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_report.R
\name{write_clusters_to_enrichment}
\alias{write_clusters_to_enrichment}
\title{Write Main clusters to enrichment output
This function allows you to report the Functional analysis.}
\usage{
write_clusters_to_enrichment(
  output_path = "results",
  output_file = "results",
  mode = "PR",
  enrichments_ORA = NULL,
  enrichments_ORA_merged = NULL,
  task_size = 1,
  workers = 1,
  template_folder = template_folder,
  top_categories = 50,
  group_results = opt$group_results,
  n_category = 30,
  sim_thr = 0.7,
  summary_common_name = "ancestor",
  pvalcutoff = 0.1,
  gene_attributes = NULL,
  gene_attribute_name = NULL
)
}
\arguments{
\item{output_path}{output folder}

\item{output_file}{output file name for heatmaps}

\item{mode}{type of output to produce}

\item{enrichments_ORA}{list of enrich results for all clusters}

\item{enrichments_ORA_merged}{merged enrich results}

\item{task_size}{number of elements per packages used}

\item{workers}{(OPTIONAL) cores for parallel features}

\item{template_folder}{(OPTIONAL) RMD templates folder}

\item{top_categories}{numbers of categories from each cluster to use for merge}

\item{group_results}{experimental - whether to group results in the emap plot}

\item{n_category}{number of categories in the figures (per cluster)}

\item{sim_thr}{value to use when combining similar categories in summary mode}

\item{summary_common_name}{'significant' to use the most significant term to label each summarized group
'ancestor' to use the common ancestor of the group"}

\item{pvalcutoff}{used to select terms for summarizing}

\item{gene_attributes}{named list of attributes e.g. FCs for emap plot coloured nodes (genes)}

\item{gene_attribute_name}{name for the legend in the emap plot for the nodes (genes)}
}
\value{
void
}
\description{
Write Main clusters to enrichment output
This function allows you to report the Functional analysis.
}