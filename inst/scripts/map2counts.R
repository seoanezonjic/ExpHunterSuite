#!/usr/bin/env Rscript


##########################################
## OPTION PARSER
##########################################

option_list <- list(
  optparse::make_option(c("-i", "--input_bam"), type = "character", default = NULL,
    help = "Path to a BAM file (must have a .bai index in the same folder)."),
  optparse::make_option("--exp_type", type = "character", default = "RNAseq_transcriptome",
    help = "Experiment type: RNAseq_transcriptome or RNAseq_genome"),
  optparse::make_option(c("-o", "--output_file"), type = "character", default = "selected_counts.txt",
    help = "Output file containing miRNA count matrix.")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

# Currently, it only works for RNAseq_transcriptome
counts <- Rsamtools::idxstatsBam(file = opt$input_bam)

counts <- counts[, c("seqnames", "mapped")]
sample_name <- basename(opt$input_bam)
colnames(counts) <- c("id", sample_name)
write.table(counts, file = opt$output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
