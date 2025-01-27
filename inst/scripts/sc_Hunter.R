#! /usr/bin/env Rscript


##########################################
## LOAD LIBRARIES
##########################################

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
    'plotting_functions.R', 'write_report.R', "statistics_functions.R", "factor_mining.R",
    "sc_library.R", "main_sc_Hunter.R")
  for (lib in custom_libraries){
    source(file.path(root_path, 'R', lib))
  }
  template_folder <- file.path(root_path, 'inst', 'templates')
}else{
  require('ExpHunterSuite')
  root_path <- find.package('ExpHunterSuite')
  template_folder <- file.path(root_path, 'templates')
}

# Parallelisation options

library(future)
options(future.globals.maxSize = 18 * 1024^6)

##########################################
## OPTPARSE
##########################################

option_list <- list(
  optparse::make_option(c("-n", "--name"), type = "character", default = NULL,
              help = "Name of analysis."),
  optparse::make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output folder."),
  optparse::make_option(c("--doublet_file"), type = "character", default = NULL,
              help = "File containing vector of barcodes to be tagged as doublet and removed from analysis."),
  optparse::make_option("--filter", type = "character", default = NULL,
              help = "TRUE for using only detected cell-associated barcodes,
                      FALSE for using all detected barcodes."),
  optparse::make_option("--mincells", type = "integer", default = NULL,
              help = "Min number of cells for which a feature was recorded."),
  optparse::make_option("--minfeats", type = "integer", default = NULL,
              help = "Min number of features for which a cell was recorded."),
  optparse::make_option("--minqcfeats", type = "integer", default = NULL,
              help = "Min number of features for which a cell was selected in QC."),
  optparse::make_option("--percentmt", type = "integer", default = NULL,
              help = "Max percentage of reads mapped to mitochondrial genes for which a cell is recorded."),
  optparse::make_option("--normalmethod", type = "character", default = NULL,
              help = "Method for normalization. LogNormalize, CLR or RC."),
  optparse::make_option("--scalefactor", type = "integer", default = NULL,
              help = "Scale factor for cell-level normalization."),
  optparse::make_option("--hvgs", type = "integer", default = NULL,
              help = "Number of HVG to be selected."),
  optparse::make_option("--ndims", type = "integer", default = NULL,
              help = "Number of PC to be used for clustering / UMAP / tSNE"),
  optparse::make_option("--dimheatmapcells", type = "integer", default = NULL,
              help = "Heatmap plots the 'extreme' cells on both ends of the spectrum."),
  optparse::make_option("--resolution", type = "double", default = NULL,
              help = "Granularity of the clustering."),
  optparse::make_option(c("-d", "--exp_design"), type = "character", default = NULL,
              help = "Input file with the experiment design."),
  optparse::make_option("--input", type = "character", default = NULL,
            help = "Count results folder. Input will be read from here."),
  optparse::make_option("--suffix", type = "character", help = "Suffix to specific file"),
  optparse::make_option("--samples_to_integrate", type = "character", default = "",
            help = "Path to file containing samples to be processed."),
  optparse::make_option("--int_columns", type = "character", default = "",
            help = "Comma-separated list of conditions by which to perform integration
                    analysis. If empty, all conditions will be analysed."),
  optparse::make_option("--cluster_annotation", type = "character", default = "",
            help = "Clusters annotation file."),
  optparse::make_option("--target_genes", type = "character", default = "",
            help = "Path to target genes table, or comma-separated list of target genes."),
  optparse::make_option("--cpu", type = "double", default = 1,
            help = "Provided CPUs."),
  optparse::make_option("--imported_counts", type = "character", default = "",
            help = "Imported counts directory."),
  optparse::make_option("--DEG_columns", type = "character", default = "",
            help = "Columns for DEG analysis."),
  optparse::make_option("--cell_annotation", type = "character", default = "",
            help = "Cell types annotation file. Will be used to dynamically
                    annotate clusters."),
  optparse::make_option("--p_adj_cutoff", type = "numeric", default = "5e-3",
            help = "Adjusted p-value cutoff."),
  optparse::make_option("--verbose", type = "logical", default = FALSE, action = "store_true",
            help = "Verbosity of base Seurat and harmony function calls."),
  optparse::make_option("--reduce", type = "logical", default = FALSE, action = "store_true",
            help = "Randomly subset seurat object to 3000 cells, for quick testing."),
  optparse::make_option("--SingleR_ref", type = "character", default = "",
            help = "Path to reference to use in SingleR annotation."),
  optparse::make_option("--ref_version", type = "character", default = "",
            help = "SingleR reference version."),
  optparse::make_option("--ref_label", type = "character", default = "main",
            help = "Column of reference metadata to use for annotation."),
  optparse::make_option("--ref_de_method", type = "character", default = "",
            help = "Method to use for markercalculation in single-cell reference."),
  optparse::make_option("--ref_n", type = "integer", default = 25,
            help = "Top N reference markers to consider in annotation. Higher values provide a more
                    accurate annotation, but increase noise and computational time. Will not be used
                    if ref_de_method is empty."),
  optparse::make_option("--integrate", type = "logical", default = FALSE, action = "store_true",
            help = "Activate integrative analysis. If FALSE (the default), script will assume
                    only one sample."),
  optparse::make_option("--saveRDS", type = "logical", default = FALSE, action = "store_true",
            help = "Save final RDS object."),
  optparse::make_option("--loadRDS", type = "logical", default = FALSE, action = "store_true",
            help = "Load RDS object instead of re-processing the entire experiment.
            Loads it from default pipeline saving location.")
)  

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

plan("multicore", workers = opt$cpu)
if(opt$cpu > 1) {
  BiocParallel::register(BiocParallel::MulticoreParam(opt$cpu))
  BPPARAM <- BiocParallel::MulticoreParam(opt$cpu)
} else {
  BPPARAM <- BiocParallel::SerialParam()
}

##########################################
## MAIN
##########################################
if(opt$cluster_annotation != "") {
  cluster_annotation <- read.table(opt$cluster_annotation, sep = "\t", header = TRUE)
} else {
  cluster_annotation <- NULL
}
if(opt$cell_annotation != "") {
  cell_annotation <- read.table(opt$cell_annotation, sep = "\t", header = TRUE)
} else {
  cell_annotation <- NULL
}
if(!is.null(opt$doublet_file)) {
  if (file.exists(opt$doublet_file)) {
    doublet_list <- read.table(opt$doublet_file)[[1]]
  } else {
    stop(paste0("Doublet file does not exist. File supplied was \"",
                opt$doublet_file, "\""))
  }
} else {
  doublet_list <- NULL
}

if(opt$integrate) {
  out_suffix <- "integration_report.html"
} else {
  out_suffix <- "sample_report.html"
}

if(opt$target_genes == ""){
  warning("No target genes provided")
  target_genes <- NULL
} else if(file.exists(opt$target_genes)) {
  target_genes <- read_and_format_targets(opt$target_genes)
} else {
  target_genes <- list(Custom = strsplit(opt$target_genes, split = ";")[[1]])
}

exp_design <- read.table(opt$exp_design, sep = "\t", header = TRUE)

if(opt$integrate) {
  if(opt$int_columns == "") {
  warning("No conditions specified for analysis. Analysing every condition")
  int_columns <- tolower(colnames(exp_design)[!colnames(exp_design)=="sample"])
  } else {
    int_columns <- tolower(unlist(strsplit(opt$int_columns, ",")))
  }
  message(paste0("Selected ", length(int_columns), " condition(s) for analysis: ", paste0(int_columns, collapse = ", ")))
} else {
  message("Starting non-integrative analysis")
  int_columns <- NULL
}

if(opt$DEG_columns == "") {
  DEG_columns <- int_columns
} else {
  DEG_columns <- opt$DEG_columns
}

# Input reading and integration variables setup
if(opt$samples_to_integrate == "") {
  samples <- exp_design$sample
} else {
  samples <- read.table(opt$samples_to_integrate, sep = "\t", header = FALSE)[[1]]
}

exp_design <- exp_design[exp_design$sample %in% samples, ]

if(!opt$loadRDS) {
  split_path <- strsplit(opt$SingleR_ref, "/")[[1]]
  if(split_path[length(split_path)] != "") {
    path_to_ref <- opt$SingleR_ref
    if(opt$ref_version != "") {
      path_to_ref <- paste(path_to_ref, opt$ref_version, sep = "_")
    }
    message("Loading provided SingleR reference")
    SingleR_ref <- HDF5Array::loadHDF5SummarizedExperiment(dir = path_to_ref, prefix = "")
  } else {
    SingleR_ref <- NULL
  }
  if(opt$ref_de_method == "") {
  ref_de_method <- NULL
  ref_n <- NULL
  } else {
    ref_de_method <- opt$ref_de_method
    ref_n <- opt$ref_n
  }
  if(opt$integrate) {
    if(opt$imported_counts == "") {
      seu <- merge_seurat(project_name = opt$name, exp_design = exp_design,
                          suffix = opt$suffix, count_path = opt$input)  
    } else {
      seu <- Seurat::CreateSeuratObject(counts = Seurat::Read10X(opt$imported_counts, gene.column = 1),
                                        project = opt$name, min.cells = 1, min.features = 1)
      seu_meta <- read.table(file.path(opt$imported_counts, "meta.tsv"), sep = "\t", header = TRUE)
      rownames(merged_seu_meta) <- colnames(merged_seu)
      seu <- Seurat::AddMetaData(merged_seu, merged_seu_meta, row.names("Cell_ID"))
    }
  } else {
    input <- file.path(opt$input, ifelse(opt$filter, "filtered_feature_bc_matrix",
                                                   "raw_feature_bc_matrix"))
    seu <- read_sc_counts(name = opt$name, input = input, mincells = opt$mincells,
                      minfeats = opt$minfeats, exp_design = exp_design)
  }
  if(opt$reduce) {
    message('Downsampling seurat object')
    seu <- downsample_seurat(seu, cells = 500, features = 5000, keep = unlist(target_genes))
  }
}

if(opt$loadRDS) {
  file <- file.path(opt$output, paste0(opt$name, ".final_results.rds"))
  message(paste0("Loading processed object from ", file))
  final_results <- readRDS(file)
} else {
  message("Analyzing seurat object")
  final_results <- main_sc_Hunter(seu = seu, cluster_annotation = cluster_annotation, name = opt$name,
                                  ndims = opt$ndims, resolution = opt$resolution, int_columns = int_columns,
                                  cell_annotation = cell_annotation, DEG_columns = DEG_columns,
                                  minqcfeats = opt$minqcfeats, percentmt = opt$percentmt, hvgs = opt$hvgs,
                                  scalefactor = opt$scalefactor, normalmethod = opt$normalmethod,
                                  p_adj_cutoff = opt$p_adj_cutoff, verbose = opt$verbose, sigfig = 2,
                                  output = opt$output, integrate = opt$integrate, query = unlist(target_genes),
                                  reduce = opt$reduce, save_RDS = opt$saveRDS, SingleR_ref = SingleR_ref,
                                  ref_label = opt$ref_label, ref_de_method = ref_de_method, ref_n = ref_n,
                                  BPPARAM = BPPARAM, doublet_list = doublet_list)
}

message("--------------------------------------------")
message("-------------Writing QC report--------------")
message("--------------------------------------------")

write_seurat_report(final_results = final_results, template_folder = template_folder,
                    template = "sc_quality_control.txt", output = file.path(opt$output, "report"),
                    target_genes = target_genes, name = opt$name, out_name = "qc_report.html",
                    use_canvas = TRUE)

message("--------------------------------------------")
message("----------Writing analysis report-----------")
message("--------------------------------------------")

write_seurat_report(final_results = final_results, template_folder = template_folder,
                    output = file.path(opt$output, "report"),
                    target_genes = target_genes, name = opt$name,
                    int_columns = int_columns, cell_annotation = cell_annotation,
                    template = "sc_analysis.txt", out_name = out_suffix)
