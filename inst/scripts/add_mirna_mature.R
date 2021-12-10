#!/usr/bin/env Rscript

option_list <- list(
	optparse::make_option(c("-i", "--input_file"), type="character", default=NULL,
                        help="DEgenesHunter results table"),
	optparse::make_option(c("-o", "--output_file"), type="character", default='./',
                        help="Define the output path. Default = './'")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

original_table <- read.table(opt$input_file,header=TRUE, row.names=1, sep="\t")

annotated_table <- original_table

annotated_table$miRNA_name <- miRBaseConverter::miRNA_AccessionToName(rownames(annotated_table))$TargetName

write.table(annotated_table, file.path(opt$output_file, "hunter_results_table_translated.txt"), quote=FALSE, row.names=TRUE, sep="\t")
