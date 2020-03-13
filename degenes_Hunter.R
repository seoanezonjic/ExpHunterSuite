#! /usr/bin/env Rscript

#' @author refactor by Fernando Moreno Jabato. Original authors Isabel Gonzalez Gayte

######################################################################################################
############################################# DEgenes Hunter #########################################
######################################################################################################


############################################################
##                      SETUP PROGRAM                     ##
############################################################
options(warn=1)
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
suppressPackageStartupMessages(require(WGCNA))
#suppressPackageStartupMessages(require(PCIT))


# Obtain this script directory
full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile),  # works when using source
               error=function(e) # works when using R CMD
              normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', commandArgs())], '='))[2]))
main_path_script <- dirname(full.fpath)

# Load custom libraries
source(file.path(main_path_script, 'lib', 'general_functions.R'))
source(file.path(main_path_script, 'lib', 'dif_expression_packages.R'))
source(file.path(main_path_script, 'lib', 'qc_and_benchmarking_functions.R'))
source(file.path(main_path_script, 'lib', 'correlation_packages.R'))
source(file.path(main_path_script, 'lib', 'correlation_packages_PCIT.R'))
source(file.path(main_path_script, 'lib', 'plotting_functions.R'))



# Prepare command line input 
option_list <- list(
  make_option(c("-i", "--input_file"), type="character", default=NULL,
    help="Input file with read counts"),
  make_option(c("-C", "--Control_columns"), type="character", default=NULL,
    help="Control columns. Please indicate column names of control samples separated by commas"),
  make_option(c("-T", "--Treatment_columns"), type="character", default=NULL,
    help="Treatment columns. Please indicate column names of treatment samples separated by commas"),  
  make_option(c("-r", "--reads"), type="numeric", default=2,
    help="Used in filtering. Reads (counts per million mapped reads) per gene per sample must be higher than this value to count towards the --minlibraries value. Default=%default.
    0 = No filtering"),
  make_option(c("-l", "--minlibraries"), type="integer", default=2,
    help="For each gene, the minimum number of samples that must have more than --reads. Default=%default"),
  make_option(c("-F", "--filter_type"), type="character", default="separate",
    help="Should filtering be done taking into account treatment and control samples separately (value=separate), or accross all samples (global). Default=%default"),
  make_option(c("-o", "--output_files"), type="character", default="hunter_DE_results",
    help="Output path. Default=%default"),
  make_option(c("-p", "--p_val_cutoff"), type="double", default=0.05,
    help="Adjusted p-value cutoff for the differential expression analysis. Default=%default"),
  make_option(c("-f", "--lfc"), type="double", default=1,
    help="Minimum log2 fold change in expression. Note this is on a log2 scale, so a value of 1 would mean a 2 fold change. Default=%default"),
  make_option(c("-m", "--modules"), type="character", default=c("DELNW"), #D = DESeq2, E = edgeR, L = limma, N = NOISeq W = WGCNA.
    help="Differential expression packages to able/disable (D = DESeq2, E = edgeR, L = limma, N = NOISeq, W = WGCNA, P = PCIT.).
    By default the following modules Default=%default are performed"),
  make_option(c("-c", "--minpack_common"), type="integer", default=4,
    help="Number of minimum package to consider a gene as a 'PREVALENT' DEG"),
  make_option(c("-t", "--target_file"), type="character", default=NULL,
    help="Sample descriptions: one column must be named treat and contain values of Treat or Ctrl. This file will take precedent over the -T and -C sample flags."),
  make_option(c("-v", "--model_variables"), type="character", default="",
    help="Variables to include in the model. Must be comma separated and each variable must be a column in the target_file, or the model can be specified precisely if the custom_model flag is TRUE"),
  make_option(c("-n", "--numerics_as_factors"), type="logical", default=TRUE,
    help="If TRUE, numeric variables in the target file specified in the DEG model be treated as distinct factors, i.e. categories. If FALSE, they will be considered as numeric and their values will be taken into account in the DEG model."),
  make_option(c("-M", "--custom_model"), type="logical", default=FALSE,
    help="If true, text in the model_variables variable will be passed directly to the model construction"),
  make_option(c("-S", "--string_factors"), type="character", default="",
    help="Columns in the target file to be used as categorical factors for the correlation analysis. If more than one to be used, should be comma separated"),
  make_option(c("-N", "--numeric_factors"), type="character", default="",
    help="Columns in the target file to be used as numeric (continuous) factors for the correlation analysis. If more than one to be used, should be comma separated"),
  make_option(c("-b", "--WGCNA_memory"), type="integer", default=5000,
    help="Maximum block size value, to be passed to the blockwiseModules function of WGCNA as the maxBlockSize argument"),
  make_option("--WGCNA_deepsplit", type="integer", default=2,
    help="This option control the module building process and is defined as 1,2,3 and 4 values. 1 for rough clustering and 4 for accurate clustering "),
  make_option("--WGCNA_min_genes_cluster", type="integer", default=20,
    help="Mnimum number of genes to keep a cluster"),  
  make_option("--WGCNA_detectcutHeight", type="double", default=0.995,
    help="Cut height to split modules"),
  make_option("--WGCNA_mergecutHeight", type="double", default=0.25,
    help="Value to merge two similar modules: Maximum dissimilarity (i.e., 1-correlation) "),
  make_option(c("-w", "--WGCNA_all"), type="logical", default=TRUE,
    help="Run WGCNA for treated only, control only, and both as 3 separate runs. Needed if using PCIT. If false, WGCNA runs once, on the table including treament and control"),
  make_option(c("--WGCNA_blockwiseNetworkType"), type="character", default="signed",
    help="NetworkType option to be passed to blockwiseModules function"),
  make_option(c("--WGCNA_blockwiseTOMType"), type="character", default="signed",
    help="TOMType option to be passed to blockwiseModules function"),
   make_option(c("--debug"), type="logical", default=FALSE, action = "store_true",
    help="Activate debug mode, which stores RData sessions at different points of the pipeline")
 )
opt <- parse_args(OptionParser(option_list=option_list))
opt_orig <- opt

############################################################
##                       CHECK INPUTS                     ##
############################################################

# Check inputs not allowed values
if (is.null(opt$input_file)){
  stop(cat("No file with RNA-seq counts is provided.\nPlease use -i to submit file"))
}
if (opt$minlibraries < 1){
  stop(cat("Minimum library number to check minimum read counts cannot be less than 1.\nIf you want to avoid filtering, set --reads to 0."))
}
if (opt$reads == 0){
  cat("Minimum reads option is set to zero. Raw count table will not be filtered.")
}

if (opt$model_variables != "" & grepl("N", opt$modules) & nchar(opt$modules) == 1) {
  stop(cat("You cannot run an experimental design that uses the model_variables option with NOISeq only."))
}
if (opt$model_variables != "" & grepl("N", opt$modules) & nchar(opt$modules) > 1) {
  warning("NOISeq will not be run as you have an experimental design that uses the model_variables option.")
  opt$modules <- gsub("N", "", opt$modules)
}

active_modules <- nchar(opt$modules)
if(grepl("W", opt$modules)) {
  active_modules <- active_modules - 1
}
if(grepl("P", opt$modules)) {
  active_modules <- active_modules - 1
}

if(opt$minpack_common > active_modules){
  opt$minpack_common <- active_modules
  warning("The number of active modules is lower than the thresold for tag PREVALENT DEG. The thresold is set to the number of active modules.")
}

# Check either C and T columns or target file.
if( (is.null(opt$Treatment_columns) | is.null(opt$Control_columns)) & is.null(opt$target_file)) {
  stop(cat("You must include either the names of the control and treatment columns or a target file with a treat column."))
}
# In the case of -C/-T AND -t target - give a warning
if( (!is.null(opt$Treatment_columns) | !is.null(opt$Control_columns)) & !is.null(opt$target_file)) {
  warning("You have included at least one -C/-T option as well as a -t option for a target file. The target file will take precedence for assigning samples labels as treatment or control.")
}
# If no -t check there is no -v value
if( is.null(opt$target_file) & opt$model_variables != "") {
  stop(cat("You should not include a -v value if you do not have a target table file."))
}
if(opt$custom_model == TRUE & opt$model_variables == "") {
  stop(cat("If you wish to use a custom model you must provide a value for the model_variables option."))
}
# If factors are specified but WGCNA not selected, throw a warning.
if((opt$string_factors != "" | opt$numeric_factors != "") & (!grepl("W", opt$modules) | is.null(opt$target_file))) {
  warning("If you wish to use factors for the correlation analysis you must also run WGCNA and include a target file. The -S and -N options will be ignored")
}

if(opt$debug){

  # Define only once
  debug_file <- file.path(normalizePath(opt$output_files), paste(c("DHunter_Debug_Session_",format(Sys.Date(),format = "%Y%m%d"),".RData"),collapse = ""))
  # Store session
  debug_point <- function(file, message = "Debug point"){
    debug_message <<- message
    save.image(file)
  }
}

############################################################
##                         I/O DATA                       ##
############################################################

# Create output folder & write original opt
dir.create(opt$output_files)
write.table(cbind(opt_orig), file=file.path(opt$output_files, "opt_input_values.txt"), sep="\t", col.names =FALSE, quote = FALSE)

# Load target file if it exists, otherwise use the -C and -T flags. Note target takes precedence over target.
if(! is.null(opt$target_file)) {
  target <- read.table(opt$target_file, header=TRUE, sep="\t")#, colClasses = "factor", stringsAsFactors=FALSE)
  # Check there is a column named treat
  if(! "treat" %in% colnames(target)) {
    stop(cat("No column named treat in the target file.\nPlease resubmit"))
  }
  index_control_cols <- as.character(subset(target, treat == "Ctrl", select = sample, drop=TRUE))
  index_treatmn_cols <- as.character(subset(target, treat == "Treat", select = sample, drop=TRUE))
  replicatesC <- length(index_control_cols)
  replicatesT <- length(index_treatmn_cols)
} else {
  index_control_cols <- unlist(strsplit(opt$Control_columns, ",")) 
  index_treatmn_cols <- unlist(strsplit(opt$Treatment_columns, ","))
  replicatesC <- length(index_control_cols)
  replicatesT <- length(index_treatmn_cols)
  # Create the target data frame needed in the calls to the DE detection methods
  target <- data.frame(sample=c(index_control_cols, index_treatmn_cols), 
    treat=c(rep("Ctrl", length(index_control_cols)), rep("Treat", length(index_treatmn_cols))))
}

# FOR WGCNA: Check that the appropriate factor columns can be found in the target file and makes a data frame with the specified factor
if(exists("target") & grepl("W", opt$modules)) {
  if(opt$string_factors != "") {
    string_factors <- unlist(strsplit(opt$string_factors, ","))
    string_factors_index <- colnames(target) %in% string_factors 

    if(TRUE %in% string_factors_index) {
      target_string_factors <- target[string_factors_index]
      target_string_factors <-  data.frame(sapply(target_string_factors, as.factor))
    } else {
      stop(cat("Factors specified with the --string_factors option cannot be found in the target file.\nPlease resubmit."))
    }
  } else {
    target_string_factors <- target["treat"] # We checkd this already in load target file code
  }
  if(opt$numeric_factors != "") {
    numeric_factors <- unlist(strsplit(opt$numeric_factors, ","))
    numeric_factors_index <- colnames(target) %in% numeric_factors
    if(TRUE %in% numeric_factors_index) {
      target_numeric_factors <- target[numeric_factors_index]
      # Ensure the factors are numeric
      #invisible(lapply(seq(ncol(target_numeric_factors)), function(i){target_numeric_factors[,i] <<- as.numeric(target_numeric_factors[,i])}))
    } else {
      stop(cat("Factors specified with the --numeric_factors option cannot be found in the target file.\nPlease resubmit."))
    }
  } else {
    target_numeric_factors <- ""
  }
}

# Now coerce the targets to factors for the multifactorial analysis
if(opt$numerics_as_factors == TRUE) {
  target <-  data.frame(sapply(target, as.factor))
}

replicatesC <- length(index_control_cols)
replicatesT <- length(index_treatmn_cols)

# Check if there are enough replicates for specified method
if((replicatesC < 2) | (replicatesT < 2)) stop('At least two replicates per class (i.e. treatment and control) are required\n')

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
  if (opt$filter_type == "separate") {
    # genes with cpm greater than --reads value for at least --minlibrariess samples for either case or control samples
    to_keep_control <- rowSums(edgeR::cpm(raw[index_control_cols]) > opt$reads) >= opt$minlibraries
    to_keep_treatment <- rowSums(edgeR::cpm(raw[index_treatmn_cols]) > opt$reads) >= opt$minlibraries
    keep_cpm <- to_keep_control | to_keep_treatment # Keep if at least minlibraries in either of them
    raw_filter <- raw[keep_cpm,] # Filter out count data frame
  } else if (opt$filter_type == "global") {
    keep_cpm <- rowSums(edgeR::cpm(raw) > opt$reads) >= opt$minlibraries # genes with cpm greater than --reads value for at least --minlibrariess samples
    raw_filter <- raw[keep_cpm,] # Filter out count data frame
  } else {
    warning("Unrecognized minimum read filter type. No filter will be used")
    raw_filter <- raw
  }
} else { # NO FILTER
  raw_filter <- raw
}
# Defining contrast - Experimental design
design_vector <- c(rep("C", replicatesC), rep("T", replicatesT))


############################################################
##                      MODEL TEXT                        ##
############################################################

# Prepare model text
if(opt$model_variables != "") {
  if(opt$custom_model == TRUE) {
    model_formula_text <- opt$model_variables
  } else {
    model_variables_unlist <- unlist(strsplit(opt$model_variables, ","))
    model_formula_text <- paste("~", paste(model_variables_unlist, "+", collapse=" "), "treat")
  }
} else {
  model_formula_text <- "~ treat"
}
cat("Model for gene expression analysis is:", model_formula_text, "\n")

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
                               p_val_cutoff = opt$p_val_cutoff,
                               target = target,
                               model_formula_text = model_formula_text)
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
                              target = target,
                              model_formula_text = model_formula_text)
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
    results <- analysis_limma(data   = raw_filter,
                              target = target,
                              model_formula_text = model_formula_text)
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
                               target = target)
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

##################################################################
##                       CORRELATION ANALYSIS                   ##
##################################################################
#####
################## CASE W: WGCNA
#####

    DESeq2_counts <- counts(package_objects[['DESeq2']][['DESeq2_dataset']], normalize=TRUE)
    DESeq2_counts_treatment <- DESeq2_counts[, index_treatmn_cols]
    DESeq2_counts_control <- DESeq2_counts[, index_control_cols]

if(grepl("W", opt$modules)) {
  if(grepl("D", opt$modules)) { 
    cat('Correlation analysis is performed with WGCNA\n')
    path <- file.path(opt$output_files, "Results_WGCNA")
    dir.create(path)

    if(opt$WGCNA_all == TRUE) {

      WGCNA_treatment_path <- file.path(path, "Treatment_only_data")
      dir.create(WGCNA_treatment_path)

      cat('Performing WGCNA correlation analysis for treated samples\n')
      results_WGCNA_treatment <- analysis_WGCNA(data=DESeq2_counts_treatment,
                     path=WGCNA_treatment_path,
                     target_numeric_factors=target_numeric_factors,
                     target_string_factors=target_string_factors,
                     WGCNA_memory=opt$WGCNA_memory,
                     WGCNA_deepsplit=opt$WGCNA_deepsplit,
                     WGCNA_detectcutHeight=opt$WGCNA_detectcutHeight,
                     WGCNA_mergecutHeight=opt$WGCNA_mergecutHeight,
                     WGCNA_min_genes_cluster=opt$WGCNA_min_genes_cluster,
                     cor_only=TRUE, 
                     blockwiseNetworkType = opt$WGCNA_blockwiseNetworkType, 
                     blockwiseTOMType = opt$WGCNA_blockwiseTOMType
      )

      WGCNA_control_path <- file.path(path, "Control_only_data")
      dir.create(WGCNA_control_path)

      
      cat('Performing WGCNA correlation analysis for control samples\n')
      results_WGCNA_control <- analysis_WGCNA(data=DESeq2_counts_control,
                     path=WGCNA_control_path,
                     target_numeric_factors=target_numeric_factors,
                     target_string_factors=target_string_factors,
                     WGCNA_memory=opt$WGCNA_memory,
                     WGCNA_deepsplit=opt$WGCNA_deepsplit,
                     WGCNA_detectcutHeight=opt$WGCNA_detectcutHeight,
                     WGCNA_mergecutHeight=opt$WGCNA_mergecutHeight,
                     WGCNA_min_genes_cluster=opt$WGCNA_min_genes_cluster,                    
                     cor_only=TRUE, 
                     blockwiseNetworkType = opt$WGCNA_blockwiseNetworkType, 
                     blockwiseTOMType = opt$WGCNA_blockwiseTOMType
      )
    }
    # Need to improve the control, probably by removing PCIT
    # if(results_WGCNA_treatment == "NO_POWER_VALUE" | results_WGCNA_control == "NO_POWER_VALUE") {
    #   warning("WGCNA was unable to generate a suitable power value for at least one of the partial datasets")
    # }


    cat('Performing WGCNA correlation analysis for all samples\n')
    results_WGCNA <- analysis_WGCNA(data=DESeq2_counts,
                                   path=path,
                                   target_numeric_factors=target_numeric_factors,
                                   target_string_factors=target_string_factors,
                                   WGCNA_memory=opt$WGCNA_memory,
                                   WGCNA_deepsplit=opt$WGCNA_deepsplit,
                                   WGCNA_detectcutHeight=opt$WGCNA_detectcutHeight,
                                   WGCNA_mergecutHeight=opt$WGCNA_mergecutHeight,
                                   WGCNA_min_genes_cluster=opt$WGCNA_min_genes_cluster,
                                   cor_only=FALSE, 
                                   blockwiseNetworkType = opt$WGCNA_blockwiseNetworkType, 
                                   blockwiseTOMType = opt$WGCNA_blockwiseTOMType
    )
    if(length(results_WGCNA) == 1) {
      warning("Something went wrong with WGCNA on the full dataset")
      opt$modules <- gsub("W", "", opt$modules)
    }
  } else {
    warning("WGCNA will not be performed as it requires a DESeq2 object")
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
  final_logFC_names, final_FDR_names, opt$p_val_cutoff, opt$lfc, opt$minpack_common)

#####
################## CASE P: PCIT
#####

if(grepl("P", opt$modules)) {
  metrics_pcit <- analysis_diff_correlation(DE_all_genes, DESeq2_counts, DESeq2_counts_control, DESeq2_counts_treatment, PCIT_filter=FALSE)
  DE_all_genes <- transform(merge(DE_all_genes, metrics_pcit, by.x=0, by.y=0), row.names=Row.names, Row.names=NULL)
}

#################################################################################

# Check WGCNA was run and it returned proper results
if(grepl("W", opt$modules) & grepl("D", opt$modules)) {
  DE_all_genes <- transform(merge(DE_all_genes, results_WGCNA[['gene_cluster_info']], by.x=0, by.y="ENSEMBL_ID"), row.names=Row.names, Row.names=NULL)
}

# Add the filtered genes back
DE_all_genes <- add_filtered_genes(DE_all_genes, raw)

# New structure - row names are now actually row names
dir.create(file.path(opt$output_files, "Common_results"))
write.table(DE_all_genes, file=file.path(opt$output_files, "Common_results", "hunter_results_table.txt"), quote=FALSE, row.names=TRUE, sep="\t")

####################### DEBUG POINT #############################
if(opt$debug) debug_point(debug_file,"Full analysis performed")
#################################################################

############################################################
##                    GENERATE REPORT                     ##
############################################################
outf <- file.path(normalizePath(opt$output_files),"DEG_report.html")
rmarkdown::render(file.path(main_path_script, 'templates', 'main_report.Rmd'), 
                  output_file = outf, intermediates_dir = opt$output_files)
