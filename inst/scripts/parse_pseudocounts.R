#!/usr/bin/env Rscript
################################################################################
############################### Parse pseudocounts #################################
################################################################################


############################################################
##                      SETUP PROGRAM                     ##
############################################################

options(warn=1)
if( Sys.getenv('DEGHUNTER_MODE') == 'DEVELOPMENT' ){
  # Obtain this script directory
  full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile), 
                 error=function(e) # works when using R CMD
                normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', 
                  commandArgs())], '='))[2]))
  main_path_script <- dirname(full.fpath)
  root_path <- file.path(main_path_script, '..', '..')
  # Load custom libraries
  custom_libraries <- c('main_parse_pseudocounts.R',
    'general_functions.R')
  for (lib in custom_libraries){
    source(file.path(root_path, 'R', lib))
  }
}else{
  require('ExpHunterSuite')
  root_path <- find.package('ExpHunterSuite')
}

################### OPTIONS
option_list <- list(
    optparse::make_option(c("-i", "--input_folder"), type="character", 
        default="",
        help="Folder that include subfolders names as sequencing samples that contains quntification files"),
    optparse::make_option(c("-m", "--mapper"), type="character", 
        help="Program that were use to quantify mappings. salmon, kallisto, rsem and stringtie are available.",
        default="salmon"),
     optparse::make_option(c("-a", "--annotation_file"), type="character", 
        help="GFF or GTF annotation file.",
        default=NULL),
    optparse::make_option(c("--output_type"), type="character", 
        help="Indicate if the output must include [G] gene counts (gene_counts.txt), [T] transcript counts (transcript_counts.txt) or both",
        default="GT"),
   optparse::make_option(c("-o", "--output_folder"), type="character",
     	default=".", 
        help = "Output folder"))
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

quant_data <- parse_pseudocounts(opt$input_folder, opt$mapper, opt$output_type, opt$annotation_file)

if(grepl("G", opt$output_type))
	write.table(quant_data$genes, file = file.path(opt$output_folder, "genes_counts.txt"), quote = FALSE, sep = "\t")

if(grepl("T", opt$output_type))
	write.table(quant_data$transcripts, file = file.path(opt$output_folder, "transcipts_counts.txt"), quote = FALSE, sep = "\t")