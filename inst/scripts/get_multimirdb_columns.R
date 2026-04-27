#!/usr/bin/env Rscript

option_list <- list(
        optparse::make_option(c("-i", "--input"), type= "character", default = NULL,
                help = "Multimir RData to retrieve columns from"),
        optparse::make_option(c("-o","--output"), type = "character", default = ".",
                help = "Set the output path."),
        optparse::make_option(c("-c","--columns_of_interest"), type = "character", default = ".",
                help = "Columns of interest to retrieve: r (RNA), m (miRNA), rm (RNA and miRNA) ")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

load(opt$input)
dir.create(opt$output, recursive = TRUE)

if (grepl("r", opt$columns_of_interest)){ 
        target_list <- unique(multimir_summary$target_ensembl)
        writeLines(target_list, file.path(opt$output, "mRNA.txt"))
}
if (grepl("m", opt$columns_of_interest)){
        miRNA_list <- unique(multimir_summary$mature_mirna_acc)
        writeLines(miRNA_list, file.path(opt$output, "miRNA.txt"))
}

