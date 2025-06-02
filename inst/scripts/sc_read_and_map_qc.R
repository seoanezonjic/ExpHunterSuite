#! /usr/bin/env Rscript


##########################################
## LOAD LIBRARIES
##########################################
# Obtain this script directory
full.fpath <- normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', 
                commandArgs())], '='))[2])

main_path_script <- dirname(full.fpath)
root_path <- file.path(main_path_script)
template_path <- file.path(root_path, "..", "templates")
# Load custom libraries
# devtools::load_all(file.path(root_path))

library(optparse)

##########################################
## OPTPARSE
##########################################

option_list <- list(
  optparse::make_option(c("-m", "--metrics"), type = "character",
              help="Metrics file in wide format"),
  optparse::make_option(c("--cellranger_metrics"), type = "character",
              help="Cell Ranger metrics file in wide format"),
  optparse::make_option(c("-o", "--output"), type = "character",
              help="Output folder"),
  optparse::make_option(c("-e", "--experiment_name"), type = "character",
              help="Experiment name"),
  optparse::make_option(c("-d", "--doublet_file"), type = "character",
              help="List of doublets")
)  

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))


##########################################
## MAIN
##########################################

write_qc_report(name = "All_samples",
                experiment = opt$experiment_name,
                template = file.path(template_path, "fastqc_report.Rmd"),
                outdir = opt$output,
                intermediate_files = "int_files",
                metrics = opt$metrics,
                long_metrics = opt$long_metrics,
                cellranger_metrics = opt$cellranger_metrics,
                cellranger_long_metrics = opt$cellranger_long_metrics)