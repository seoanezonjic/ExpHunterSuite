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
  devtools::load_all(root_path)
  template_folder <- file.path(root_path, 'inst', 'templates')
} else {
  require('ExpHunterSuite')
  root_path <- find.package('ExpHunterSuite')
  template_folder <- file.path(root_path, 'templates')
}

option_list <- list(
  optparse::make_option(c("-g", "--gtf"), type="character", default=NULL,
    help="GTF file for txdb."),
  optparse::make_option(c("-o", "--output_dir"), type="character", default=NULL,
    help="Path to output directory.")
  )
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

preprocessing_results <- preprocess_gtf(opt$gtf)

if(is.null(opt$output_dir)) {
  output_dir <- getwd()
} else {
  output_dir <- opt$output_dir
}

gtf_name <- basename(tools::file_path_sans_ext(opt$gtf))
db_file <- file.path(output_dir, paste0(gtf_name, "_txdb.db"))
cr_file <- file.path(output_dir, paste0(gtf_name, "_count_ranges.rds"))
map_file <- file.path(output_dir, paste0(gtf_name, "_gene_name_mapping.tsv"))

AnnotationDbi::saveDb(preprocessing_results$txdb, db_file)
saveRDS(preprocessing_results$count_ranges, cr_file)
write.table(preprocessing_results$gene_name_mapping, file = map_file, sep = "\t", quote = FALSE,
            row.names = FALSE)
