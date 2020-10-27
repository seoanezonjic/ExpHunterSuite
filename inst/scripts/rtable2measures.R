#! /usr/bin/env Rscript

#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>

#############################################
### CONFIGURE 
#############################################

if( Sys.getenv('DEGHUNTER_MODE') == 'DEVELOPMENT' ){
	full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile),  # works when using source
               error=function(e) # works when using R CMD
              normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', commandArgs())], '='))[2]))
	main_path_script <- dirname(full.fpath)
	root_path <- file.path(main_path_script, '..', '..')

	custom_libraries <- c('statistics_functions.R')
	for (lib in custom_libraries){
		source(file.path(root_path, 'R', lib))
	}

	#############################################
	### LOAD LIBRARIES
	#############################################
	suppressPackageStartupMessages(require(optparse))
	suppressPackageStartupMessages(library(ROCR)) 
}else{
	require('DEgenesHunter')
	root_path <- find.package('DEgenesHunter')
}



option_list <- list(
  optparse::make_option(c("-i", "--inputfile"), type="character", 
    help="DEGenesHunter expression resulta table file"),
  optparse::make_option(c("-r", "--realprediction"), type="character",
    help="Real prediction values file"),
  optparse::make_option(c("-e", "--experiment"), type="character", default = NULL,
    help="[OPTIONAL] Experiment metrics file to be included to output."),
  optparse::make_option(c("-f", "--logFC_threshold"), type="double", default=1,
    help="Log Fold Change threshold. Default : %default"),
  optparse::make_option(c("-p", "--pval_threshold"), type="double", default=0.01,
    help="P-value threshold. Default : %default"),
  optparse::make_option(c("-a", "--aucfile"), type="character", default=NULL,
    help="If provided, AUC values for each method will be stored on this file"),
  optparse::make_option(c("-o", "--outfile"), type="character",
    help="Output file")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))



#############################################
### LOAD & PREPARE
#############################################

df_cuts <- rtable2measures(htfile = opt$inputfile,
						  realprediction = opt$realprediction,
						  aucfile = opt$aucfile,
						  experiment = opt$experiment,
						  logFC_threshold = opt$logFC_threshold,
						  pval_threshold = opt$pval_threshold)

# Store results table
write.table(df_cuts, file = opt$outfile, quote = FALSE, col.names = TRUE, sep = "\t", row.names = FALSE)