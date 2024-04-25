#!/usr/bin/env Rscript

options(warn=1)
if( Sys.getenv('DEGHUNTER_MODE') == 'DEVELOPMENT' ){
  # Obtain this script directory
  full.fpath <- normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', 
                  commandArgs())], '='))[2])

  main_path_script <- dirname(full.fpath)
  root_path <- file.path(main_path_script, '..', '..')
  organisms_table_file <- file.path(root_path, "inst","external_data", 
        "organism_table.txt")
  # Load custom libraries
   custom_libraries <- c('functional_analysis_library.R')
  for (lib in custom_libraries){
    source(file.path(root_path, 'R', lib))
  }

  source_folder <- file.path(root_path, 'inst')
}else{
  require('ExpHunterSuite')
  root_path <- find.package('ExpHunterSuite')
  source_folder <- file.path(root_path)
  organisms_table_file <- file.path(root_path, "external_data", 
        "organism_table.txt")

}



option_list <- list(
	optparse::make_option(c("-i", "--input_file"), type="character", default=NULL,
                        help="Table to annotate"),
	optparse::make_option(c("-o", "--output_file"), type="character", default='Annotated_table.txt',
                        help="Define the output path. Default = %default"),
	optparse::make_option(c("-c", "--column"), type="character", default=1,
                        help="Column name or index with ensembl_IDs. Default = %default"),
	optparse::make_option(c("-I", "--input_keytype"), type="character", default=NULL,
                        help="Set the input keytype. Default=%default : All possible keytypes will be printed"),
	optparse::make_option(c("-K", "--output_keytype"), type="character", default="SYMBOL",
                        help="Set the output keytype. Default=%default"),
	optparse::make_option(c("-O", "--organism"), type="character", default="Human",
                        help="Set the model organism. Default = %default")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

organisms_table <- get_organism_table(organisms_table_file)
current_organism_info <- organisms_table[rownames(organisms_table) %in% opt$organism,]
org_db <- get_org_db(current_organism_info)

if (grepl("[0-9]+", opt$column)) opt$column <- as.numeric(opt$column)

table_to_annot <- read.table(opt$input_file, sep = "\t", header = TRUE)

translated_keytypes <- translate_ids_orgdb(ids = table_to_annot[,opt$column], 
                    input_id=opt$input_keytype, output_id = opt$output_keytype, org_db=org_db)

translated_ids <- translated_keytypes[match(table_to_annot[,opt$column], 
                                        translated_keytypes[,opt$input_keytype]) ,opt$output_keytype]

table_to_annot[,opt$output_keytype] <- translated_ids

write.table(table_to_annot, sep = "\t", quote = FALSE, row.names = FALSE, file = opt$output_file)