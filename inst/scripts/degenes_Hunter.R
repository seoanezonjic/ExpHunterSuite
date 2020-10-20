#! /usr/bin/env Rscript

#' @author refactor by Fernando Moreno Jabato. Original authors Isabel Gonzalez Gayte

######################################################################################################
############################################# DEgenes Hunter #########################################
######################################################################################################

###############################################
## FOR FUNCTION DEBUGGING
## save(list = ls(all.names = TRUE), file = "~/test.RData", envir = environment())
###############################################

############################################################
##                      SETUP PROGRAM                     ##
############################################################
options(warn=1)
if( Sys.getenv('DEGHUNTER_MODE') == 'DEVELOPMENT' ){
  # Loading libraries
  suppressPackageStartupMessages(library(ggplot2)) 
  suppressPackageStartupMessages(library(limma))
  suppressPackageStartupMessages(library(edgeR))
  suppressPackageStartupMessages(library(DESeq2)) 
  suppressPackageStartupMessages(library(NOISeq))
  suppressPackageStartupMessages(library(VennDiagram))
  suppressPackageStartupMessages(library(gplots))
  suppressPackageStartupMessages(library(stringr))
  suppressPackageStartupMessages(library(plyr))
  suppressPackageStartupMessages(library(knitr))
  suppressPackageStartupMessages(library(FSA))
  suppressPackageStartupMessages(require(rmarkdown))
  suppressPackageStartupMessages(require(reshape2))
  suppressPackageStartupMessages(require(PerformanceAnalytics))
  suppressPackageStartupMessages(require(WGCNA))
  suppressPackageStartupMessages(require(diffcoexp))
  #suppressPackageStartupMessages(require(PCIT))

  # Obtain this script directory
  full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile),  # works when using source
                 error=function(e) # works when using R CMD
                normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', commandArgs())], '='))[2]))
  main_path_script <- dirname(full.fpath)
  root_path <- file.path(main_path_script, '..', '..')
  # Load custom libraries
  custom_libraries <- c('main_degenes_Hunter.R', 'general_functions.R', 'dif_expression_packages.R', 'qc_and_benchmarking_functions.R', 'correlation_packages.R', 'correlation_packages_PCIT.R', 'plotting_functions.R', 'general_functions.R')
  for (lib in custom_libraries){
    source(file.path(root_path, 'R', lib))
  }
  template_folder <- file.path(root_path, 'inst', 'templates')
}else{
  require('DEgenesHunter')
  package_folder <- find.package('DEgenesHunter')
  template_folder <- file.path(package_folder, 'templates')
}

# Prepare command line input 
option_list <- list(
  optparse::make_option(c("-i", "--input_file"), type="character", default=NULL,
    help="Input file with read counts"),
  optparse::make_option(c("-C", "--Control_columns"), type="character", default=NULL,
    help="Control columns. Please indicate column names of control samples separated by commas"),
  optparse::make_option(c("-T", "--Treatment_columns"), type="character", default=NULL,
    help="Treatment columns. Please indicate column names of treatment samples separated by commas"),  
  optparse::make_option(c("-r", "--reads"), type="numeric", default=2,
    help="Used in filtering. Reads (counts per million mapped reads) per gene per sample must be higher than this value to count towards the --minlibraries value. Default=%default.
    0 = No filtering"),
  optparse::make_option(c("-l", "--minlibraries"), type="integer", default=2,
    help="For each gene, the minimum number of samples that must have more than --reads. Default=%default"),
  optparse::make_option(c("-F", "--filter_type"), type="character", default="separate",
    help="Should filtering be done taking into account treatment and control samples separately (value=separate), or accross all samples (global). Default=%default"),
  optparse::make_option(c("-o", "--output_files"), type="character", default="hunter_DE_results",
    help="Output path. Default=%default"),
  optparse::make_option(c("-p", "--p_val_cutoff"), type="double", default=0.05,
    help="Adjusted p-value cutoff for the differential expression analysis. Default=%default"),
  optparse::make_option(c("-f", "--lfc"), type="double", default=1,
    help="Minimum log2 fold change in expression. Note this is on a log2 scale, so a value of 1 would mean a 2 fold change. Default=%default"),
  optparse::make_option(c("-m", "--modules"), type="character", default=c("DELNW"), #D = DESeq2, E = edgeR, L = limma, N = NOISeq W = WGCNA.
    help="Differential expression packages to able/disable (D = DESeq2, E = edgeR, L = limma, N = NOISeq, W = WGCNA, P = PCIT, X = diffcoexp, F = external_DEA_file.).
    By default the following modules Default=%default are performed"),
  optparse::make_option(c("-c", "--minpack_common"), type="integer", default=4,
    help="Number of minimum package to consider a gene as a 'PREVALENT' DEG"),
  optparse::make_option(c("-t", "--target_file"), type="character", default=NULL,
    help="Sample descriptions: one column must be named treat and contain values of Treat or Ctrl. This file will take precedent over the -T and -C sample flags."),
  optparse::make_option(c("-e", "--external_DEA_file"), type="character", default=NULL,
    help="External data file containing preanalysed DE data. Must consist of three columns containing p-value, logFC and FDR/padjust IN THAT ORDER"),
  optparse::make_option(c("-v", "--model_variables"), type="character", default="",
    help="Variables to include in the model. Must be comma separated and each variable must be a column in the target_file, or the model can be specified precisely if the custom_model flag is TRUE"),
  optparse::make_option(c("-n", "--numerics_as_factors"), type="logical", default=TRUE,
    help="If TRUE, numeric variables in the target file specified in the DEG model be treated as distinct factors, i.e. categories. If FALSE, they will be considered as numeric and their values will be taken into account in the DEG model."),
  optparse::make_option(c("-M", "--custom_model"), type="logical", default=FALSE,
    help="If true, text in the model_variables variable will be passed directly to the model construction"),
  optparse::make_option(c("-S", "--string_factors"), type="character", default="",
    help="Columns in the target file to be used as categorical factors for the correlation analysis. If more than one to be used, should be comma separated"),
  optparse::make_option(c("-N", "--numeric_factors"), type="character", default="",
    help="Columns in the target file to be used as numeric (continuous) factors for the correlation analysis. If more than one to be used, should be comma separated"),
  optparse::make_option(c("-b", "--WGCNA_memory"), type="integer", default=5000,
    help="Maximum block size value, to be passed to the blockwiseModules function of WGCNA as the maxBlockSize argument"),
  optparse::make_option(c("--WGCNA_norm_method"), type="character", default="DESeq2",
    help="Method used to normalized the table of counts for WGCNA. Must also run this method in the --modules argument. Default=%default"),
  optparse::make_option("--WGCNA_deepsplit", type="integer", default=2,
    help="This option control the module building process and is defined as 1,2,3 and 4 values. 1 for rough clustering and 4 for accurate clustering "),
  optparse::make_option("--WGCNA_min_genes_cluster", type="integer", default=20,
    help="Mnimum number of genes to keep a cluster"),  
  optparse::make_option("--WGCNA_detectcutHeight", type="double", default=0.995,
    help="Cut height to split modules"),
  optparse::make_option("--WGCNA_mergecutHeight", type="double", default=0.25,
    help="Value to merge two similar modules: Maximum dissimilarity (i.e., 1-correlation) "),
  optparse::make_option(c("-w", "--WGCNA_all"), type="logical", default=TRUE,
    help="Run WGCNA for treated only, control only, and both as 3 separate runs. Needed if using PCIT. If false, WGCNA runs once, on the table including treament and control"),
  optparse::make_option(c("--WGCNA_blockwiseNetworkType"), type="character", default="signed",
    help="NetworkType option to be passed to blockwiseModules function"),
  optparse::make_option(c("--WGCNA_blockwiseTOMType"), type="character", default="signed",
    help="TOMType option to be passed to blockwiseModules function"),
  optparse::make_option(c("--debug"), type="logical", default=FALSE, action = "store_true",
    help="Activate debug mode, which stores RData sessions at different points of the pipeline"),
  optparse::make_option(c("--Debug"), type="character", default=NULL,
    help="Activate debug mode and uses given filename. File must have '.RData' extension")
 )
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

#DRAFT INPUT/OUTPUT
dir.create(opt$output_files)
write.table(cbind(opt), file=file.path(opt$output_files, "opt_input_values.txt"), sep="\t", col.names =FALSE, quote = FALSE)

main_degenes_Hunter(
  input_file=opt$input_file,
  Control_columns=opt$Control_columns,
  Treatment_columns=opt$Treatment_columns,
  reads=opt$reads,
  minlibraries=opt$minlibraries,
  filter_type=opt$filter_type,
  output_files=opt$output_files,
  p_val_cutoff=opt$p_val_cutoff,
  lfc=opt$lfc,
  modules=opt$modules,
  minpack_common=opt$minpack_common,
  target_file=opt$target_file,
  external_DEA_file=opt$external_DEA_file,
  model_variables=opt$model_variables,
  numerics_as_factors=opt$numerics_as_factors,
  custom_model=opt$custom_model,
  string_factors=opt$string_factors,
  numeric_factors=opt$numeric_factors,
  WGCNA_memory=opt$WGCNA_memory,
  WGCNA_norm_method=opt$WGCNA_norm_method,
  WGCNA_deepsplit=opt$WGCNA_deepsplit,
  WGCNA_min_genes_cluster=opt$WGCNA_min_genes_cluster,
  WGCNA_detectcutHeight=opt$WGCNA_detectcutHeight,
  WGCNA_mergecutHeight=opt$WGCNA_mergecutHeight,
  WGCNA_all=opt$WGCNA_all,
  WGCNA_blockwiseNetworkType=opt$WGCNA_blockwiseNetworkType,
  WGCNA_blockwiseTOMType=opt$WGCNA_blockwiseTOMType,
  debug=opt$debug,
  Debug=opt$Debug,
  opt=opt, # To allow markdown report output input data. Think how fit this in the new package
  template_folder=template_folder
)

#DRAFT OUTPUT
#TODO