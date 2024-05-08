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
  optparse::make_option(c("--fpkm_cutoff"), type="numeric", default=1,
    help="FPKM cutoff to filter expressed genes."),
  optparse::make_option(c("--implementation"), type="character", default="autoencoder",
    help="Implementation for sample covariation removal."),
  optparse::make_option(c("--max_tested_dimension_proportion"), type="numeric", default=3,
    help="Quotient by which to divide the number of samples to calculate
          maximum value for encoding dimension."),
  optparse::make_option(c("-f", "--hpo_file"), type="character", default=NULL,
    help="Genes and associated HPO terms in compressed tsv format."),
  optparse::make_option(c("-b", "--sample_bam_stats"), type="character", default=NULL,
    help="Bam stats files in plain text format."),
  optparse::make_option(c("-t", "--top_N"), type="integer", default=10,
    help="Top N genes by adjusted p-value to be selected.")
  )

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

final_results <- main_abgenes_Hunter(
  sample_annotation = opt$sample_annotation,
  anno_database = opt$anno_database,
  count_ranges = opt$count_ranges,
  gene_mapping_file = opt$gene_mapping_file,
  count_files = opt$count_files,
  dataset = opt$dataset,
  cpu = opt$cpu,
  fpkm_cutoff = opt$fpkm_cutoff,
  implementation = opt$implementation,
  max_dim_proportion = opt$max_tested_dimension_proportion,
  hpo_file = opt$hpo_file,
  sample_bam_stats = opt$sample_bam_stats,
  top_N = opt$top_N)

data.table::fwrite(data.table::as.data.table(SummarizedExperiment::assay(
                   final_results$counts)$counts, keep.rownames = 'geneID'),
                   file = paste0(dataset, "_geneCounts.tsv.gz"),
                   quote = FALSE, row.names = FALSE, sep = '\t',
                   compress = 'gzip')
