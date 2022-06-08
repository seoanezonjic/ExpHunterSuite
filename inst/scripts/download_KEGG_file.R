#!/usr/bin/env Rscript

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
  custom_libraries <- c('io_handling.R', "functional_analysis_library.R", "functional_analysis_library_new.R")
  for (lib in custom_libraries){
    source(file.path(root_path, 'R', lib))
  }
  template_folder <- file.path(root_path, 'inst', 'templates')
  organisms_table_file <- file.path(root_path, "inst", "external_data", 
        "organism_table.txt")
} else {
  require('ExpHunterSuite')
  root_path <- find.package('ExpHunterSuite')
  template_folder <- file.path(root_path, 'templates')
  organisms_table_file <- file.path(root_path, "external_data", 
        "organism_table.txt")
}

########################## OPTIONS
option_list <- list(
	  optparse::make_option(c("-K", "--kegg_data_file"), ,type = "character", default=NULL,
                        help=paste0("path to download KEGG data file. If unspecified, defaults to inst/kegg_data_files")), 
      optparse::make_option(c("-O", "--model_organism"), type="character", default="Human", 
                        help="Model organism.")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))


kegg_files_path <- file.path(root_path, 'inst', 'kegg_data_files')
if(! file.exists(kegg_files_path)) dir.create(kegg_files_path)

organisms_table <- get_organism_table(organisms_table_file)
current_organism_info <- organisms_table[rownames(organisms_table) %in% opt$model_organism,]
print("current_organism_info:")
print(current_organism_info)
kegg_data_file <- get_kegg_db_path(opt$kegg_data_file, current_organism_info=current_organism_info)
print("kegg_data_file:")
print(kegg_data_file)
download_latest_kegg_db(current_organism_info, kegg_data_file)
