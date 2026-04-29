#!/usr/bin/env Rscript

option_list <- list(
        optparse::make_option(c("-i", "--input"), type= "character", default = NULL,
                help = "Multimir RData to retrieve columns from"),
        optparse::make_option(c("-o","--output"), type = "character", default = "./gene_mirna",
                help = "Output file"),
        optparse::make_option(c("-c","--columns_of_interest"), type = "character", default = "rm",
                help = "Columns of interest to retrieve: r (RNA), m (miRNA), rm (RNA and miRNA) "),
	optparse::make_option(c("-d","--databases"), type = "character", default = "",
                help = "Comma-separated list of databases to consult. Call flag -l to list available databases"),
	optparse::make_option(c("-l","--list_databases"), type = "logical", default = FALSE, action = "store_true",
                help = "List all available databases")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
load(opt$input)
if(opt$list_databases){
	options(show.error.messages = FALSE)
	print(colnames(multimir_summary)[-which(colnames(multimir_summary) %in% c("target_ensembl", "mature_mirna_acc"))])
	stop()
}

columns_to_retrieve <- unlist(strsplit(opt$databases, ","))

if(grepl("r", opt$columns_of_interest)) columns_to_retrieve <- c(columns_to_retrieve, "target_ensembl") 
if(grepl("m", opt$columns_of_interest)) columns_to_retrieve <- c(columns_to_retrieve, "mature_mirna_acc") 

unique_cols <- unique(multimir_summary[, columns_to_retrieve])
write.table(unique_cols, opt$output, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

