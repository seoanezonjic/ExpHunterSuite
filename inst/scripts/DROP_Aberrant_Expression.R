#! /usr/bin/env Rscript


devtools::load_all("~aestebanm/dev_R/ExpHunterSuite") ## Temporary until it is installed
options(error = function() { 
  traceback(2)
  options(error = NULL)
  stop("exiting after script error") 
})
### ORIGIN:  MergeCounts.R
option_list <- list(
  optparse::make_option(c("-s", "--sample_annotation"), type="character", default=NULL,
    help="Sample annotation table in tsv format."),
  optparse::make_option(c("-a", "--anno_database"), type="character", default=NULL,
    help="TxDb annotation database in db format."),
  optparse::make_option(c("-r", "--count_ranges"), type="character", default=NULL,
    help="Count ranges R object in rds format."),
  optparse::make_option(c("-g", "--gene_mapping_file"), type="character", default=NULL,
    help="Gene mapping file in tsv format."),
  optparse::make_option(c("-c", "--count_files"), type="character", default=NULL,
    help="Count tables, can be rds objects or text files."),
  optparse::make_option(c("-d", "--dataset"), type="character", default=NULL,
    help="Dataset name."),
  optparse::make_option(c("-n", "--cpu"), type="integer", default=NULL,
    help="Number of CPUs provided to job."),
  optparse::make_option(c("-o", "--config_file"), type="character", default=NULL,
    help="Configuration file in yaml format. Will be loaded as a list."),
  optparse::make_option(c("-f", "--hpo_file"), type="character", default=NULL,
    help="Genes and associated HPO terms in compressed tsv format."),
  optparse::make_option(c("-b", "--sample_bam_stats"), type="character", default=NULL,
    help="Bam stats files in plain text format."),
  optparse::make_option(c("-t", "--top_N"), type="integer", default=10,
    help="Top N genes by adjusted p-value to be selected.")
  )

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

final_results <- main_aberrant_expression(
  sample_annotation = opt$sample_annotation,
  anno_database = opt$anno_database,
  count_ranges = opt$count_ranges,
  gene_mapping_file = opt$gene_mapping_file,
  count_files = opt$count_files,
  dataset = opt$dataset,
  cpu = opt$cpu,
  config_file = opt$config_file,
  hpo_file = opt$hpo_file,
  sample_bam_stats = opt$sample_bam_stats,
  top_N = opt$top_N)
