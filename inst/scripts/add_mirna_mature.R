#!/usr/bin/env Rscript

option_list <- list(
	optparse::make_option(c("-i", "--input_file"), type="character", default=NULL,
                        help="DEgenesHunter results table"),
	optparse::make_option(c("-o", "--output_file"), type="character", default='./',
                        help="Define the output path. Default = './'")
)

if( Sys.getenv('DEGHUNTER_MODE') == 'DEVELOPMENT' ){ 

    full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile),  
                   error=function(e) # works when using R CMD
                  normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', 
                    commandArgs())], '='))[2]))
    main_path_script <- dirname(full.fpath)
    root_path <- file.path(main_path_script, '..', '..')
    custom_libraries <- c("miRNA_RNA_functions.R")
    for (lib in custom_libraries){
        source(file.path(root_path, 'R', lib))
      }
    template_folder <- file.path(root_path, 'inst', 'templates')
    organism_table_path <- file.path(root_path,"inst","external_data", 
        "organism_table.txt")

} else {
    require('ExpHunterSuite')
    root_path <- find.package('ExpHunterSuite')
    template_folder <- file.path(root_path, 'templates')
    organism_table_path <- file.path(root_path, "inst","external_data",
        "organism_table.txt")
}



opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

original_table <- read.table(opt$input_file,header=TRUE, row.names=1, sep="\t")

annotated_table <- original_table

annotated_table$miRNA_name <- translate_miRNA_ids(rownames(annotated_table))

write.table(annotated_table, file.path(opt$output_file, "hunter_results_table_translated.txt"), quote=FALSE, row.names=TRUE, sep="\t")
