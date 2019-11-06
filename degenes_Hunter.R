#! /usr/bin/env Rscript

#' @author refactor by Fernando Moreno Jabato. Original authors Isabel Gonzalez Gayte

######################################################################################################
############################################# DEgenes Hunter #########################################
######################################################################################################


############################################################
##                      SETUP PROGRAM                     ##
############################################################

# Loading libraries
suppressPackageStartupMessages(library(ggplot2)) 
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(DESeq2)) 
suppressPackageStartupMessages(library(NOISeq))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(FSA))
suppressPackageStartupMessages(require(rmarkdown))
suppressPackageStartupMessages(require(reshape2))
suppressPackageStartupMessages(require(PerformanceAnalytics))

# Obtain this script directory
full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile),  # works when using source
               error=function(e) # works when using R CMD
              normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', commandArgs())], '='))[2]))
main_path_script <- dirname(full.fpath)

# Load custom libraries
source(file.path(main_path_script, 'lib', 'general_functions.R'))
source(file.path(main_path_script, 'lib', 'dif_expression_packages.R'))
source(file.path(main_path_script, 'lib', 'qc_and_benchmarking_functions.R'))

# Prepare command line input 
option_list <- list(
  make_option(c("-i", "--input_file"), type="character", default=NULL,
    help="Input file with read counts"),
  make_option(c("-C", "--Control_columns"), type="character", default=NULL,
    help="Control columns. Please indicate column names of control samples separated by commas"),
  make_option(c("-T", "--Treatment_columns"), type="character", default=NULL,
    help="Treatment columns. Please indicate column names of treatment samples separated by commas"),  
  make_option(c("-r", "--reads"), type="integer", default=2,
    help="Used in filtering. Reads (counts per million mapped reads) per gene per sample must be higher than this value to count towards the --minlibraries value. Default=%default.
    0 = No filtering"),
  make_option(c("-l", "--minlibraries"), type="integer", default=2,
    help="For each gene, the minimum number of samples that must have more than --reads. Default=%default"),
  make_option(c("-o", "--output_files"), type="character", default="hunter_DE_results",
    help="Output path. Default=%default"),
  make_option(c("-p", "--p_val_cutoff"), type="double", default=0.05,
    help="Adjusted p-value cutoff for the differential expression analysis. Default=%default"),
  make_option(c("-f", "--lfc"), type="double", default=1,
    help="Minimum log2 fold change in expression. Note this is on a log2 scale, so a value of 1 would mean a 2 fold change. Default=%default"),
  make_option(c("-q", "--q_value"), type="double", default=0.95,
    help="q value for NOISeq package. Default=%default"),
  make_option(c("-a", "--adjust_method"), type="character", default=c("BH"), #D = DESeq2, E = edgeR, L = limma, N = NOISeq
    help="Method selection to adjust the combined nominal p-values. By default method Default=%default is performed"),
  make_option(c("-n", "--name_exp"), type="character", default="experiment1",
    help="Type the name of your experiment."),
  make_option(c("-m", "--modules"), type="character", default=c("DELN"), #D = DESeq2, E = edgeR, L = limma, N = NOISeq
    help="Differential expression packages to able/disable (D = DESeq2, E = edgeR, L = limma, N= NOISeq.).
    By default the following modules Default=%default are performed"),
  make_option(c("-c", "--minpack_common"), type="integer", default=4,
    help="Number of minimum package to consider a gene as a 'PREVALENT' DEG"),
  make_option(c("-L", "--linked_samples"), type="logical", default=FALSE,
    help="If TRUE, samples will be linked (paired) by the order they appear in the arguments Treatment_columns and Control_columns Default=%default")
)
opt <- parse_args(OptionParser(option_list=option_list))

############################################################
##                       CHECK INPUTS                     ##
############################################################

# Check inputs not allowed values
if (is.null(opt$input_file)){
  stop(cat("No file with RNA-seq counts is provided.\nPlease use -i to submit file"))
}
if (is.null(opt$Control_columns)){
  stop(cat("No control samples are indicated.\nPlease select column names of the count table with -C"))
}
if (is.null(opt$Treatment_columns)){
  stop(cat("No treatment samples are indicated.\nPlease select column names of the count table with -T"))
}
if (opt$minlibraries < 1){
  stop(cat("Minimum library number to check minimum read counts cannot be less than 1.\nIf you want to avoid filtering, set --reads to 0."))
}
if (opt$reads == 0){
  cat("Minimum reads option is set to zero. Raw count table will not be filtered.")
}
if (opt$linked_samples == TRUE & (length(opt$Treatment_columns) != length(opt$Control_columns))) {
  cat("If a linked (paired) design is used there should be equivalent numbers of treatment and control samples.")
}
if (opt$linked_samples == TRUE & grepl("N", opt$modules) & nchar(opt$modules) == 1) {
  stop(cat("You cannot run a linked (paired) experimental design analysis with NOISeq only."))
  opt$modules <- gsub("N", "", opt$modules)
}
if (opt$linked_samples == TRUE & grepl("N", opt$modules) & nchar(opt$modules) > 1) {
  warning("As you are using a linked (paired) experiment design, NOISeq will not be run.")
  opt$modules <- gsub("N", "", opt$modules)
}

active_modules <- nchar(opt$modules)
if(opt$minpack_common > active_modules){
  opt$minpack_common <- active_modules
  warning("The number of active modules is lower than the thresold for tag PREVALENT DEG. The thresold is set to the number of active modules.")
}


# Control replicate number VS method selection
index_control_cols <- unlist(strsplit(opt$Control_columns, ",")) 
index_treatmn_cols <- unlist(strsplit(opt$Treatment_columns, ","))
replicatesC <- length(index_control_cols)
replicatesT <- length(index_treatmn_cols)

# Check if there are enough replicates for specified method
if((replicatesC < 2) | (replicatesT < 2)) stop('At least two replicates per class (i.e. treatment and control) are required\n')

############################################################
##                         I/O DATA                       ##
############################################################

# Create output folder 
dir.create(opt$output_files)

# Load raw count data
raw <- read.table(opt$input_file, header=TRUE, row.names=1, sep="\t")
raw <- raw[c(index_control_cols,index_treatmn_cols)] #Indexing selected columns from input count file

# Substitute NA values
raw[is.na(raw)] <- 0

############################################################
##                          FILTER                        ##
############################################################
# Prepare filtered set
if(opt$reads != 0){
  keep_cmp <- rowSums(edgeR::cpm(raw) > opt$reads) >= opt$minlibraries # genes with cpm greater than --reads value for at least --minlibrariess samples
  raw_filter <- raw[keep_cmp,] # Filter out count data frame
}else{ # NO FILTER
  raw_filter <- raw
}

# Defining contrast - Experimental design
design_vector <- c(rep("C", replicatesC), rep("T", replicatesT))

############################################################
##                       D.EXP ANALYSIS                   ##
############################################################

# Prepare results containers
all_data_normalized     <- list()
all_counts_for_plotting <- list()
all_package_results     <- list()
package_objects         <- list()
all_FDR_names           <- c()
all_LFC_names           <- c()
all_pvalue_names        <- c()
final_logFC_names       <- c()
final_FDR_names         <- c()
final_pvalue_names      <- c()
DEG_pack_columns        <- c()
lfc <- opt$lfc

#####
################## CASE D: DESeq2
#####
if(grepl("D",opt$modules)){
  if(replicatesC >= 2 & replicatesT >= 2) {
    # Verbose
    cat('Gene expression analysis is performed with DESeq2.\n')
    dir.create(file.path(opt$output_files, "Results_DESeq2"))
    # Calculate results
    results <- analysis_DESeq2(data   = raw_filter, 
                               p_val_cutoff  = opt$p_val_cutoff, 
                               lfc    = lfc, 
                               groups = design_vector,
                               linked_samples = opt$linked_samples)
    # Store results
    all_data_normalized[['DESeq2']] <- results[[1]]
    all_counts_for_plotting[['DESeq2']] <- results[[2]]
    package_objects[['DESeq2']] <- results[[3]]

    # Result Plot Visualization
    if (!is.null(all_counts_for_plotting[['DESeq2']])){
      all_FDR_names      <- c(all_FDR_names, 'padj')
      all_LFC_names      <- c(all_LFC_names, 'log2FoldChange')
      all_pvalue_names   <- c(all_pvalue_names, 'pvalue')
      final_pvalue_names <- c(final_pvalue_names, 'pvalue_DESeq2')
      final_logFC_names  <- c(final_logFC_names, 'logFC_DESeq2')
      final_FDR_names    <- c(final_FDR_names, 'FDR_DESeq2')
      DEG_pack_columns   <- c(DEG_pack_columns, 'DESeq2_DEG')
    } 
  } else {
  warning("DESeq2 will not be performed due to too few replicates")
  }
}

#####
################## CASE E: edgeR (AT LEAST 2 REPLICLATES)
#####
if(grepl("E", opt$modules)){ 
  if(replicatesC >= 2 & replicatesT >= 2) {
    # Verbose point
    cat('Gene expression analysis is performed with edgeR.\n')
    path <- file.path(opt$output_files, "Results_edgeR")
    dir.create(path)
    # Calculate results
    results <- analysis_edgeR(data   = raw_filter, 
                     p_val_cutoff = opt$p_val_cutoff, 
                     lfc   = lfc, 
                     path  = path,
                     groups = design_vector,
                     linked_samples = opt$linked_samples)
    # Store results
    all_data_normalized[['edgeR']] <- results[[1]]
    all_counts_for_plotting[['edgeR']] <- results[[2]]
    package_objects[['edgeR']] <- results[[3]]

    # Result Plot Visualization
    if (!is.null(all_counts_for_plotting[['edgeR']])){
      all_FDR_names      <- c(all_FDR_names, 'FDR')
      all_LFC_names      <- c(all_LFC_names, 'logFC')
      all_pvalue_names   <- c(all_pvalue_names, 'PValue')
      final_pvalue_names <- c(final_pvalue_names, 'pvalue_edgeR')
      final_logFC_names  <- c(final_logFC_names, 'logFC_edgeR')
      final_FDR_names    <- c(final_FDR_names, 'FDR_edgeR')
      DEG_pack_columns   <- c(DEG_pack_columns, 'edgeR_DEG')
    }
  } else {
  warning("edgeR will not be performed due to too few replicates")
  }
}

#####
################## CASE L: limma (AT LEAST 3 REPLICATES)
#####
if(grepl("L", opt$modules)){ 
  if(replicatesC >= 3 & replicatesT >= 3) {
    # Verbose
    cat('Gene expression analysis is performed with limma.\n')
    dir.create(file.path(opt$output_files, "Results_limma"))

    # Calculate results
    results <- analysis_limma(data = raw_filter, 
                              num_controls  = replicatesC, 
                              num_treatmnts = replicatesT, 
                              p_val_cutoff  = opt$p_val_cutoff, 
                              lfc  = lfc,
                              linked_samples = opt$linked_samples)
    # Store results
    all_data_normalized[['limma']]     <- results[[1]]
    all_counts_for_plotting[['limma']] <- results[[2]]

    # Result Plot Visualization
    if (!is.null(all_counts_for_plotting[['limma']])){
      all_FDR_names      <- c(all_FDR_names, 'adj.P.Val')
      all_LFC_names      <- c(all_LFC_names, 'logFC')
      all_pvalue_names   <- c(all_pvalue_names, 'P.Value')
      final_pvalue_names <- c(final_pvalue_names, 'pvalue_limma')
      final_logFC_names  <- c(final_logFC_names, 'logFC_limma')
      final_FDR_names    <- c(final_FDR_names, 'FDR_limma')
      DEG_pack_columns   <- c(DEG_pack_columns, 'limma_DEG')
    }
  } else {
    warning("limma will not be performed due to too few replicates")
  }
}

#####
################## CASE N: NOISeq (AT LEAST 3 REPLICLATES)
#####
if(grepl("N", opt$modules)){ 
    if(replicatesC >= 3 & replicatesT >= 3) {

    # Verbose
    cat("Gene expression analysis is performed with NOISeqBIO function within NOISeq.\n")
    path <- file.path(opt$output_files, "Results_NOISeq")
    dir.create(path)
    # Calculate results
    results <- analysis_NOISeq(data    = raw_filter, 
                               num_controls  = replicatesC, 
                               num_treatmnts = replicatesT, 
                               q_value = opt$q_value, 
                               groups  = design_vector, 
                               path = path,
                               lfc  = lfc)
    # Store results
    all_data_normalized[['NOISeq']]     <- results[[1]]
    all_counts_for_plotting[['NOISeq']] <- results[[2]]
    package_objects[['NOISeq']] <- results[[3]]

    #Result Plot Visualization
    if (!is.null(all_counts_for_plotting[['NOISeq']])){
      all_FDR_names      <- c(all_FDR_names, 'adj.p')
      all_LFC_names      <- c(all_LFC_names, 'log2FC')
      all_pvalue_names   <- c(all_pvalue_names, 'prob')
      final_pvalue_names <- c(final_pvalue_names, 'pvalue_NOISeq')
      final_logFC_names  <- c(final_logFC_names, 'logFC_NOISeq')
      final_FDR_names    <- c(final_FDR_names, 'FDR_NOISeq')
      DEG_pack_columns   <- c(DEG_pack_columns, 'NOISeq_DEG')
    }
  } else {
    warning("NOISeq will not be performed due to too few replicates")
  }
}

#################################################################################
##                       EXPORT FINAL RESULTS AND OTHER FILES                  ##
#################################################################################

# Write filtered count data
write.table(raw_filter, file=file.path(opt$output_files, "filtered_count_data.txt"), quote=FALSE, col.names=NA, sep="\t")

# Write Counts data and all DE for each  design to file as 2 column data frame
write.table(data.frame(class = design_vector, name = c(index_control_cols, index_treatmn_cols)), 
  file=file.path(opt$output_files, "control_treatment.txt"), row.names=FALSE, quote=FALSE, sep="\t")
write_df_list_as_tables(all_data_normalized, prefix = 'Normalized_counts_', root = opt$output_files)
write_df_list_as_tables(all_counts_for_plotting, prefix = 'allgenes_', root = opt$output_files)

# Write main results file, normalized counts per package and DE results per package
DE_all_genes <- unite_DEG_pack_results(all_counts_for_plotting, all_FDR_names, all_LFC_names, all_pvalue_names, final_pvalue_names, 
  final_logFC_names, final_FDR_names, raw, opt$p_val_cutoff, opt$lfc, opt$minpack_common)
# New structure - row names are now actually row names
dir.create(file.path(opt$output_files, "Common_results"))
write.table(DE_all_genes, file=file.path(opt$output_files, "Common_results", "hunter_results_table.txt"), quote=FALSE, row.names=TRUE, sep="\t")

############################################################
##                    GENERATE REPORT                     ##
############################################################
outf <- file.path(opt$output_files, "DEG_report.html")

rmarkdown::render(file.path(main_path_script, 'templates', 'main_report.Rmd'), 
                  output_file = outf, intermediates_dir = opt$output_files)