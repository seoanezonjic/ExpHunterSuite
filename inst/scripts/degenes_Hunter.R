#! /usr/bin/env Rscript

#' @author refactor by Fernando Moreno Jabato. 
#' Original authors Isabel Gonzalez Gayte

################################################################################
############################### DEgenes Hunter #################################
################################################################################


############################################################
##                      SETUP PROGRAM                     ##
############################################################
options(warn=1)
if( Sys.getenv('DEGHUNTER_MODE') == 'DEVELOPMENT' ){
  # Obtain this script directory
  full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile), 
                 error=function(e) # works when using R CMD
                normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', 
                  commandArgs())], '='))[2]))
  main_path_script <- dirname(full.fpath)
  root_path <- file.path(main_path_script, '..', '..')
  # Load custom libraries
  custom_libraries <- c('main_degenes_Hunter.R', 'io_handling.R', 
    'general_functions.R', 'dif_expression_packages.R', 
    'qc_and_benchmarking_functions.R', 'correlation_packages.R', 
    'plotting_functions.R', 'write_report.R', "statistics_functions.R", "factor_mining.R")
  for (lib in custom_libraries){
    source(file.path(root_path, 'R', lib))
  }
  template_folder <- file.path(root_path, 'inst', 'templates')
}else{
  require('ExpHunterSuite')
  root_path <- find.package('ExpHunterSuite')
  template_folder <- file.path(root_path, 'templates')
}

# Prepare command line input 
option_list <- list(
  optparse::make_option(c("-i", "--input_file"), type="character", default=NULL,
    help="Input file with read counts"),
  optparse::make_option(c("-C", "--Control_columns"), type="character", 
    default=NULL,
    help=paste0("Control columns. Please indicate column names of control",
      " samples separated by commas")),
  optparse::make_option(c("-T", "--Treatment_columns"), type="character", 
    default=NULL,
    help=paste0("Treatment columns. Please indicate column names of treatment",
      " samples separated by commas")),  
  optparse::make_option(c("-r", "--reads"), type="numeric", default=2,
    help=paste0("Used in filtering. Reads (counts per million mapped reads)",
      " per gene per sample must be higher than this value to count towards",
      " the --minlibraries value. Default=%default. 0 = No filtering")),
  optparse::make_option(c("--count_var_quantile"), type="double", default=0, 
    help=paste0("Proportion (0-1 value) of low variance genes to filter out. ",
      "Default=%default")),
  optparse::make_option(c("-l", "--minlibraries"), type="integer", default=2,
    help=paste0("For each gene, the minimum number of samples that must have",
      " more than --reads. Default=%default")),
  optparse::make_option(c("-F", "--filter_type"), type="character", 
    default="separate",
    help=paste0("Should filtering be done taking into account treatment and",
      " control samples separately (value=separate), possible combinations of factors ",
      "(combined), or across all factors (global). In the case of combined factors, ",
      "the factors should be specified using \"combined:factor1,factor2,...\". Default=%default")),
  optparse::make_option(c("-o", "--output_files"), type="character", 
    default="hunter_DE_results",
    help="Output path. Default=%default"),
  optparse::make_option(c("-p", "--p_val_cutoff"), type="double", default=0.05,
    help=paste0("Adjusted p-value cutoff for the differential expression",
      " analysis. Default=%default")),
  optparse::make_option(c("-f", "--lfc"), type="double", default=1,
    help=paste0("Minimum log2 fold change in expression. Note this is on a",
      " log2 scale, so a value of 1 would mean a 2 fold change.",
      " Default=%default")),
  optparse::make_option(c("-m", "--modules"), type="character", 
    default=c("DELNW"), #D = DESeq2, E = edgeR, L = limma, N = NOISeq W = WGCNA.
    help=paste0("Differential expression packages to able/disable ",
      "(D = DESeq2, E = edgeR, L = limma, N = NOISeq, W = WGCNA, P = PCIT,",
      " X = diffcoexp, F = external_DEA_file.). By default the following ",
      "modules Default=%default are performed")),
  optparse::make_option(c("-c", "--minpack_common"), type="integer", default=4,
    help="Number of minimum package to consider a gene as a 'PREVALENT' DEG"),
  optparse::make_option(c("-t", "--target_file"), type="character",
    default=NULL,
    help=paste0("Sample descriptions: one column must be named treat and",
      " contain values of Treat or Ctrl. This file will take precedent over",
      " the -T and -C sample flags.")),
  optparse::make_option(c("-e", "--external_DEA_file"), type="character", 
    default=NULL,
    help=paste0("External data file containing preanalysed DE data.",
      " Must consist of four columns with the following names and",
      " corresponding information: Entrez (or other gene id supported by",
      " functional hunter), P.Value, logFC and adj.P.Val IN THAT ORDER")),
  optparse::make_option(c("-v", "--model_variables"), type="character", 
    default="",
    help=paste0("Variables to include in the model. Must be comma separated",
      " and each variable must be a column in the target_file")),
  optparse::make_option(c("-n", "--numerics_as_factors"), type="logical", 
    default=FALSE,
    help=paste0("If TRUE, numeric variables in the target file specified",
      " in the DEG model be treated as distinct factors, i.e. categories.",
      " If FALSE, they will be considered as numeric and their values will be",
      " taken into account in the DEG model.")),
  optparse::make_option(c("-S", "--string_factors"), type="character", 
    default="",
    help=paste0("Columns in the target file to be used as categorical",
      " factors for the correlation analysis. If more than one to be used,",
      " should be comma separated")),
  optparse::make_option(c("-N", "--numeric_factors"), type="character", 
    default="",
    help=paste0("Columns in the target file to be used as numeric (continuous)",
      " factors for the correlation analysis. If more than one to be used,",
      " should be comma separated")),
  optparse::make_option(c("-b", "--WGCNA_memory"), type="integer", default=5000,
    help=paste0("Maximum block size value, to be passed to the",
      " blockwiseModules function of WGCNA as the maxBlockSize argument.",
      " Default=%default")),
  optparse::make_option(c("--WGCNA_norm_method"), type="character", 
    default="DESeq2",
    help=paste0("Method used to normalized the table of counts for WGCNA.",
      " Must also run this method in the --modules argument.",
      " Default=%default")),
  optparse::make_option("--WGCNA_deepsplit", type="integer", default=2,
    help=paste0("This option control the module building process and is",
      " defined as 1,2,3 and 4 values. 1 for rough clustering and 4 for",
      " accurate clustering. Default=%default")),
  optparse::make_option("--WGCNA_min_genes_cluster", type="integer", default=20,
    help="Mnimum number of genes to keep a cluster. Default=%default"),  
  optparse::make_option("--WGCNA_detectcutHeight", type="double", default=0.995,
    help="Cut height to split modules. Default=%default"),
  optparse::make_option("--WGCNA_mergecutHeight", type="double", default=0.25,
    help=paste0("Value to merge two similar modules: Maximum dissimilarity",
      " (i.e., 1-correlation). Default=%default")),
  optparse::make_option(c("-w", "--WGCNA_all"), type="logical", default=FALSE, 
    action = "store_true",
    help=paste0("Run WGCNA for treated only, control only, and both as 3",
      " separate runs. Needed if using PCIT. If false, WGCNA runs once, on the",
      " table including treament and control")),
  optparse::make_option(c("--WGCNA_blockwiseNetworkType"), type="character", 
    default="signed",
    help=paste0("NetworkType option to be passed to blockwiseModules function",
      " (unsigned, signed, signed hybrid).")),
  optparse::make_option(c("--WGCNA_blockwiseTOMType"), type="character", 
    default="signed",
    help=paste0("TOMType option to be passed to blockwiseModules function ",
      "(none, unsigned, signed, signed Nowick, unsigned 2, signed 2 and",
      " signed Nowick 2). If none, adjacency will be used for clustering.",
      " Default=%default")),
  optparse::make_option(c("--WGCNA_minCoreKME"), type="double", default=0.7, 
    help=paste0("Minimum module membership threshold to define module core.",
      " Modules under the threshold will not be considered. Default=%default")),
  optparse::make_option(c("--WGCNA_minCoreKMESize"), type="integer", 
    default=NULL, 
    help=paste0("Minimun genes that belong to module core. Modules under",
      " the threshold will not be considered. Default min(20,",
      " WGCNA_min_genes_cluster/3).")),
  optparse::make_option(c("--WGCNA_minKMEtoStay"), type="double", default=0.5, 
    help=paste0("Minimun module membership of a gene to be kept in module.",
      " Default=%default")),
  optparse::make_option("--WGCNA_corType", type="character", default="pearson",
    help="Set the correlation method to perform WGCNA_algorithm. 'pearson' and 'bicor' are the allowed options"),
  optparse::make_option(c("--multifactorial"), type="character", default="", 
        help=paste0("Currently only a 2x2 or factorial design is possible for interactions, and 2xn for group effects. Nested designs can",
          "also be specified (i.e. Group-specific condition effects, individuals nested within groups) The required contrast must be specified in the following manner: ",
          "FactorA,FactorB:contrast. Contrast can be either",
      "\"interaction,baseA,baseB\" if we are interested in the interaction between the two factors, where baseA and baseB should be the base levels for each factor. ",
      "FC would represent [numA_numB - baseA_numB] - [numA_baseB - baseA_baseB] with numA/B representing the non-base levels for the factorA.",
      "Alternatively, Contrast can be specificed in the form \"effect,baseA,groupB\", where the baseA should be the level in FactorA that should be used as the base for FC calculation, ",
      "and groupB represents the level in Factor B that is the group we are looking for the change in. For effect, FactorB can have more than 2 groups, allowing 2xn designs. ",
      "Finally, if nested is selected, we can perform the same comparisons using a nested design with: \"nested_int,Ctrl,groupA\" for interaction, ",
      "\"nested_effect,Ctrl,groupA\" or \"nested_effect,Ctrl,groupB\" for a group. ")),
  optparse::make_option(c("-s", "--library_sizes"), type="character",
    default=NULL, help="Path to file containing library sizes. If missing,
      certain DROP QC plots will not be drawn, and sample ranks will be defined
      by total counts.")
 )
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
#############################################################################
#DRAFT INPUT/OUTPUT
############################################################################
dir.create(opt$output_files)
write.table(cbind(opt), file=file.path(opt$output_files, 
  "opt_input_values.txt"), sep="\t", col.names =FALSE, quote = FALSE)

#############################################################################
#CHECKINGS
############################################################################

# Don't check if we have given an exernal DEA file
if (is.null(opt$external_DEA_file)) {
  # Check inputs not allowed values
  if (is.null(opt$input_file)){
    stop(cat(paste0("No file with RNA-seq counts is provided.",
      "\nPlease use -i to submit file")))
  }
# Check either C and T columns or target file.
  if( (is.null(opt$Treatment_columns) | is.null(opt$Control_columns)) & 
       is.null(opt$target_file)) {
    stop(cat(paste0("You must include either the names of the control and",
      " treatment columns or a target file with a treat column.")))
  }
  # In the case of -C/-T AND -t target - give a warning
  if( (!is.null(opt$Treatment_columns) | !is.null(opt$Control_columns)) &
       !is.null(opt$target_file)) {
    warning(paste0("You have included at least one -C/-T option as well as a",
      " -t option for a target file. The target file will take precedence for",
      " assigning samples labels as treatment or control."))
  }
}

#############################################################################
#EXPRESSION ANALYSIS
############################################################################
#
if(! is.null(c(opt$target_file, opt$Control_columns, opt$Treatment_columns))) {
  target <- target_generation(from_file=opt$target_file, 
      ctrl_samples=opt$Control_columns, treat_samples=opt$Treatment_columns)
} else {
  target <- NULL
}
if(! is.null(opt$input_file)) {
  raw_count_table <- read.table(opt$input_file, 
    header=TRUE, row.names=1, sep="\t")
} else {
  raw_count_table <- NULL
}
external_DEA_data <- NULL
if (grepl("F", opt$modules)) { 
  # Open the external data file for pre-analysed deg data
  external_DEA_data <- read.table(opt$external_DEA_file,
   header=TRUE, sep="\t", row.names=1)
}
library_sizes <- opt$library_sizes
if(! is.null(library_sizes)) {
  library_sizes <- read.table(library_sizes, header=TRUE)
}


final_results <- main_degenes_Hunter(
  target=target,
  raw=raw_count_table,
  count_var_quantile=opt$count_var_quantile,
  external_DEA_data=external_DEA_data,
  output_files=opt$output_files,
  reads=opt$reads,
  minlibraries=opt$minlibraries,
  filter_type=opt$filter_type,
  p_val_cutoff=opt$p_val_cutoff,
  lfc=opt$lfc,
  modules=opt$modules,
  minpack_common=opt$minpack_common,
  model_variables=opt$model_variables,
  numerics_as_factors=opt$numerics_as_factors,
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
  WGCNA_minCoreKME=opt$WGCNA_minCoreKME,
  WGCNA_minCoreKMESize=opt$WGCNA_minCoreKMESize,
  WGCNA_minKMEtoStay = opt$WGCNA_minKMEtoStay,
  WGCNA_corType = opt$WGCNA_corType,
  multifactorial = opt$multifactorial,
  library_sizes=library_sizes
)

#############################################################################
#WRITE OUTPUT
############################################################################


  write_expression_data(final_results, opt$output_files, template_folder, opt)

  write_expression_report(final_results, opt$output_files, template_folder, opt)


