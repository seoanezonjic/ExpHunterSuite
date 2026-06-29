#! /usr/bin/env Rscript


##########################################
## OPTION PARSER
##########################################

# Prepare command line input 
option_list <- list(
  optparse::make_option(c("-i", "--input_file"), type = "character", default = NULL,
                        help = "Input file with gene table."),
  optparse::make_option("--gtf", type = "character", default = NULL,
                        help = "gtf file."),
  optparse::make_option("--column_name", type = "character", default = NULL,
                        help = "Name of the column in the input file that contains the IDs. By default, the rownames are used."),
  optparse::make_option(c("-o", "--output_file"), type="character", default='gene_table_feature_type.txt',
                        help="Define the output path.")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

##########################################
## LOAD LIBRARIES
##########################################

options(warn=1)
if (Sys.getenv('DEGHUNTER_MODE') == 'DEVELOPMENT'){
  # Loading libraries
  # Obtain this script directory
  full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile),  
                 error=function(e) # works when using R CMD
                normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', 
                  commandArgs())], '='))[2]))
  main_path_script <- dirname(full.fpath)
  root_path <- file.path(main_path_script, '..', '..')
  # Load custom libraries
  devtools::load_all(root_path)
} else {
  require('ExpHunterSuite')
  root_path <- find.package('ExpHunterSuite')
}

# Package to import the .gtf file
library(rtracklayer)

gene_table <- read.table(opt$input_file, header = TRUE, sep = "\t")
gtf_data <- as.data.frame(import(opt$gtf))

genes_gtf <- subset(gtf_data, type == "gene")
genes_gtf <- unique(genes_gtf[, c("gene_id", "gene_type", "gene_name")])

# In the .gtf file, the versions are found in the gene_id column (example ENSG00000290825.2), but not in the gene_table:
# the version is removed from .gtf
genes_gtf$gene_id <- sub("\\..*$", "", genes_gtf$gene_id)

# For the functional analysis results, the name of the column containing the IDs must be specified.
if (!is.null(opt$column_name)) {
  vector_to_match <- gene_table[, opt$column_name]
} else {
  vector_to_match <- rownames(gene_table)
}
gene_table$ref_gene_type <- genes_gtf$gene_type[match(vector_to_match, genes_gtf$gene_id)]
gene_table$ref_gene_name <- genes_gtf$gene_name[match(vector_to_match, genes_gtf$gene_id)]

write.table(gene_table, sep = "\t", quote = FALSE, row.names = TRUE, file = opt$output_file)
