#!/usr/bin/env Rscript

option_list <- list(
        optparse::make_option(c("-i", "--input"), type= "character", default = NULL,
                help = "Multimir RData to retrieve columns from"),
        optparse::make_option(c("-o","--output"), type = "character", default = ".",
                help = "Output file"),
        optparse::make_option(c("-c","--columns_of_interest"), type = "character", default = ".",
                help = "Columns of interest to retrieve: r (RNA), m (miRNA), rm (RNA and miRNA) ")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

load(opt$input)
dir.create(opt$output, recursive = TRUE)

columns_to_retrieve <- NULL

if(grepl("r", opt$columns_of_interest)) columns_to_retrieve <- c(columns_to_retrieve, "target_ensembl") 
if(grepl("m", opt$columns_of_interest)) columns_to_retrieve <- c(columns_to_retrieve, "mature_mirna_acc") 

unique_cols <- unique(multimir_summary[, columns_to_retrieve])
write.table(unique_cols, opt$output, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

