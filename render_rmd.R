#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2)) 
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(DESeq2)) 
suppressPackageStartupMessages(library(NOISeq))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(FSA))
suppressPackageStartupMessages(require(rmarkdown))
suppressPackageStartupMessages(require(reshape2))
suppressPackageStartupMessages(require(PerformanceAnalytics))
suppressPackageStartupMessages(require(WGCNA))


script_path <- tryCatch(normalizePath(parent.frame(2)$ofile),  # works when using source
               error=function(e) # works when using R CMD
              normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', commandArgs())], '='))[2]))


# print(file.path(dirname(script_path), 'lib', 'plotting_functions.R'))
source(file.path(dirname(script_path), 'lib', 'plotting_functions.R'))


option_list <- list(
  make_option(c("-t", "--template"), type="character", default="main_report",
    help="Specify DEgenesHunter report template: 'main_report' for main DEGhunter execution, 'functional_report' for functional enrichment for differential expression packages and 'clusters_main_report' for functional enrichments of all clusters"),
  make_option(c("-i", "--input"), type="character", default=NULL,
    help="Specify DEGhunter sesion for generate report. ADVICE: --debug flag must be activated in previous DEGhunter execution"),
  make_option(c("-o", "--output"), type="character", default=NULL,
    help="Specify output path. Report html will be created with the name of the report template")
  )
actual_opt <- parse_args(OptionParser(option_list=option_list))

if(actual_opt$template == "functional_report"){
  suppressPackageStartupMessages(require(topGO))
  suppressPackageStartupMessages(require(biomaRt)) 
  suppressPackageStartupMessages(require(KEGGREST))
  suppressPackageStartupMessages(require(clusterProfiler))
}


if(is.null(actual_opt$output)){
  actual_opt$output <- paste0("./", actual_opt$template)
}
if(!grepl(".html", actual_opt$output)){
  actual_opt$output <- paste0(actual_opt$output, ".html")
}

load(actual_opt$input)
### load enviroment

### load template
DEGhunter_templates <- file.path(dirname(script_path), "/", "templates")


if(actual_opt$template == "functional_report"){
  results_path = dirname(actual_opt$output)
}

### Render
rmarkdown::render(
	normalizePath(file.path(DEGhunter_templates, paste0(actual_opt$template, ".Rmd"))), 
	output_file = file.path(actual_opt$output),
	intermediates_dir = file.path(normalizePath(dirname(actual_opt$output)), "tmp")
	)
