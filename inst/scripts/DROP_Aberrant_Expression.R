#! /usr/bin/env Rscript


##########################################
## LOAD LIBRARIES
##########################################

if( Sys.getenv('DEGHUNTER_MODE') == 'DEVELOPMENT' ){
  # Obtain this script directory
  full.fpath <- normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', 
                  commandArgs())], '='))[2])
  main_path_script <- dirname(full.fpath)
  root_path <- file.path(main_path_script, '..', '..')
  # Load custom libraries
  custom_libraries <- c('main_abgenes_Hunter.R', 'DROP_functions.R')
  for (lib in custom_libraries){
    source(file.path(root_path, 'R', lib))
  }
  template_folder <- file.path(root_path, 'inst', 'templates')
} else {
  require('ExpHunterSuite')
  root_path <- find.package('ExpHunterSuite')
  template_folder <- file.path(root_path, 'templates')
}

##########################################
## OPTPARSE
##########################################

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
    help="Quotient by which to divide the number of samples to calculate maximum value for encoding dimension."),
  optparse::make_option(c("--p_adj_cutoff"), type="numeric", default=NULL,
    help="adjusted P-value to use as cutoff to mark a gene as aberrantly expressed."),
  optparse::make_option(c("--z_score_cutoff"), type="numeric", default=NULL,
    help="z score value to use as cutoff to mark a gene as aberrantly expressed."),
  optparse::make_option(c("-f", "--hpo_file"), type="character", default=NULL,
    help="Genes and associated HPO terms in compressed tsv format."),
  optparse::make_option(c("-b", "--stats_path"), type="character", default=NULL,
    help="Path to directory containing bam stats."),
  optparse::make_option(c("-t", "--top_N"), type="integer", default=10,
    help="Top N genes by adjusted p-value to be selected."),
  optparse::make_option(c("-o", "--report_dir"), type="character", default=getwd(),
    help="Directory where report will be generated.")
  )

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

##########################################
## MAIN
##########################################

BiocParallel::register(BiocParallel::MulticoreParam(opt$cpu))

final_results <- main_abgenes_Hunter(sample_annotation = opt$sample_annotation, anno_database = opt$anno_database,
                                     count_ranges = opt$count_ranges, gene_mapping_file = opt$gene_mapping_file,
                                     count_files = opt$count_files, dataset = opt$dataset, cpu = opt$cpu,
                                     fpkm_cutoff = opt$fpkm_cutoff, implementation = opt$implementation,
                                     max_dim_proportion = opt$max_tested_dimension_proportion, p_adj_cutoff = opt$p_adj_cutoff,
                                     z_score_cutoff = opt$z_score_cutoff, hpo_file = opt$hpo_file, stats_path = opt$stats_path,
                                     top_N = opt$top_N)

write_abgenes_report(final_results = final_results, template_folder = template_folder, output_dir = opt$report_dir,
                     p_adj_cutoff = opt$p_adj_cutoff, z_score_cutoff = opt$z_score_cutoff, top_N = opt$top_N)
