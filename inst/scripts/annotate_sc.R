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
  if(Sys.getenv("sketch") == "TRUE") {
    devtools::load_all(Sys.getenv("HTMLREPORT_PATH"))
    custom_libraries <- c('sc_library.R', 'main_sc_functions.R', 'general_functions.R')
    source_folder <- file.path(find.package("htmlreportR"), "inst")
    for (lib in custom_libraries){
      source(file.path(root_path, 'R', lib))
    }
  } else {
    source_folder <- NULL
    devtools::load_all(root_path)
  }
  template_folder <- file.path(root_path, 'inst', 'templates')
}else{
  require('ExpHunterSuite')
  root_path <- find.package('ExpHunterSuite')
  template_folder <- file.path(root_path, 'templates')
}

##########################################
## OPTPARSE
##########################################

option_list <- list(
  optparse::make_option(c("-n", "--name"), type = "character", default = NULL,
              help = "Name of analysis."),
  optparse::make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output folder."),
  optparse::make_option(c("--doublet_file"), type = "character", default = "",
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
  optparse::make_option("--resolution", type = "double", default = NULL,
              help = "Granularity of the clustering."),
  optparse::make_option(c("-d", "--exp_design"), type = "character", default = NULL,
              help = "Input file with the experiment design."),
  optparse::make_option("--input", type = "character", default = NULL,
            help = "Count results folder. Input will be read from here."),
  optparse::make_option("--suffix", type = "character", help = "Suffix to specific file"),
  optparse::make_option("--samples_to_integrate", type = "character", default = "",
            help = "Path to file containing samples to be processed."),
  optparse::make_option("--subset_by", type = "character", default = "",
            help = "Comma-separated list of conditions by which to perform integration
                    analysis. If empty, all conditions will be analysed."),
  optparse::make_option("--extra_columns", type = "character", default = "",
            help = "Comma-separated list of extra conditions to represent in certain plots."),
  optparse::make_option("--int_method", type = "character", default = "RPCA",
            help = "Integration method. Valid methods: \"CCA\", \"RPCA\", \"Harmony\", \"FastMNN\", \"scVI\"."),
  optparse::make_option("--k_weight", type = "integer", default = 100,
            help = "Number of neighbors to consider when weighting anchors. Used in integration."),
  optparse::make_option("--cluster_annotation", type = "character", default = "",
            help = "Clusters annotation file."),
  optparse::make_option("--cpu", type = "integer", default = 1,
            help = "Provided CPUs."),
  optparse::make_option("--imported_counts", type = "character", default = "",
            help = "Imported counts directory."),
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
            help = paste0("Activate integrative analysis. If FALSE (the default), script will assume ",
                    "only one sample.")),
  optparse::make_option("--filter_dataset", type = "character", default = "",
            help = "Filter imported counts according to string."),
  optparse::make_option("--ref_filter", type = "character", default = "",
            help = "Filter SingleR reference according to string."),
  optparse::make_option("--sketch", type = "logical", default = FALSE, action = "store_true",
            help = "Sketch experiment."),
  optparse::make_option("--force_ncells", type = "integer", default = NA_integer_,
            help = "An integer. Skip all sketching tests and force to sketch this many cells from each sample."),
  optparse::make_option("--sketch_pct", type = "numeric", default = 12,
            help = paste0("A numeric. Percentage of total cells to consider representative of the",
  " experiment. Default 12, as suggested by sketching tutorial.")),
  optparse::make_option("--sketch_method", type = "character", default = 12,
            help = "Score calculation method to select cells in sketch."),
  optparse::make_option("--genome", type = "character", default = "Unspecified",
            help = paste0("Genome version. Optional, included as extra ",
              "information in output and reports."))
)

params <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

options(future.globals.maxSize = 15e+09)
options(Seurat.object.assay.version = 'v5')

##########################################
## LOADING INPUT FILES
##########################################

updated_params <- process_sc_params(params, mode = "annotation")
opt <- updated_params$opt
doublet_list <- updated_params$doublet_list
out_suffix <- updated_params$out_suffix

if(opt$cpu > 1) {
  BiocParallel::register(BiocParallel::MulticoreParam(opt$cpu))
  BPPARAM <- BiocParallel::MulticoreParam(opt$cpu)
} else {
  BPPARAM <- BiocParallel::SerialParam()
}
message("CPU provided to BiocParallel: ", opt$cpu)

##########################################
## MAIN
##########################################
final_counts_path <- file.path(opt$output, "counts/matrix.mtx.gz")
if(!file.exists(final_counts_path)) {
  # Input parser function
  # Reference loader function
  SingleR_ref <- NULL
  if(opt$SingleR_ref != "/" & file.exists(opt$SingleR_ref)) {
    SingleR_ref <- load_SingleR_ref(path = opt$SingleR_ref, version = opt$ref_version, filter = opt$ref_filter)
  }
  if(opt$integrate) {
    if(opt$imported_counts == "") {
      seu <- merge_seurat(project_name = opt$name, exp_design = opt$exp_design,
                          suffix = opt$suffix, count_path = opt$input)
    } else {
      seu <- Seurat::CreateSeuratObject(counts = Seurat::Read10X(opt$imported_counts, gene.column = 1),
                                        project = opt$name, min.cells = 1, min.features = 1)
      seu_meta <- read.table(file.path(opt$imported_counts, "meta.tsv"), sep = "\t", header = TRUE)
      rownames(seu_meta) <- colnames(seu)
      seu <- Seurat::AddMetaData(seu, seu_meta, row.names("Cell_ID"))
    }
  } else {
    input <- file.path(opt$input, ifelse(opt$filter, "filtered_feature_bc_matrix",
                                                   "raw_feature_bc_matrix"))
    seu <- read_sc_counts(name = opt$name, input = input, mincells = opt$mincells,
                      minfeats = opt$minfeats, exp_design = opt$exp_design)
  }
  message(paste0("Total cells in dataset: ", ncol(seu), "."))
  if(!is.null(opt$filter_dataset)) {
    message("Filtering input data")
    expressions <- strsplit(opt$filter_dataset, "&|\\|")[[1]]
    expressions <- gsub(" $", "", expressions)
    expressions <- gsub("^ ", "", expressions)
    if(grepl("&", opt$filter_dataset)) {
      operator <- "&"
    } else {
      operator <- "|"
    }
    filter <- vector(mode = "list", length = length(expressions))
    for(i in seq(expressions)) {
      filter[[i]] <- parse_filter(object = "seu@meta.data",
                                  expression = expressions[i]) 
    }
    filter <- Reduce(operator, filter)
    message(paste0("Selected ", sum(filter), " cells for input data subsetting."))
    seu <- seu[, filter]
  }
  if(opt$reduce) {
    message('Downsampling seurat object')
    seu <- downsample_seurat(seu, cells = 500, features = 5000)
  }
  # Function end
}

if(file.exists(final_counts_path)) {
  message("Reconstructing Seurat object from directory ", opt$output)
  seu <- Seurat::CreateSeuratObject(counts = Seurat::Read10X(file.path(opt$output, "counts"), gene.column = 1),
                                    project = opt$name, min.cells = 1, min.features = 1)
  seu_meta <- read.table(file.path(opt$output, "counts/meta.tsv"), sep = "\t", header = TRUE)
  rownames(seu_meta) <- colnames(seu)
  seu <- Seurat::AddMetaData(seu, seu_meta, row.names("Cell_ID"))
  seu$RNA$data <- seu$RNA$counts
  markers <- read.table(file.path(opt$output, "markers.tsv"), sep = "\t", header = TRUE)
  embeddings <- read.table(file.path(opt$output, "embeddings", "cell_embeddings.tsv"), header = TRUE)
  seu$umap <- Seurat::CreateDimReducObject(embeddings = as.matrix(embeddings), key = 'umap_', assay = 'RNA')
  expr_metrics <- get_expression_metrics(seu = seu, sigfig = 2)
  final_results <- list(seu = seu, markers = markers, integrate = TRUE, clusters_pct = expr_metrics$clusters_pct,
                        sample_qc_pct = expr_metrics$sample_qc_pct)
} else {
  message("Analyzing seurat object")
  final_results <- main_annotate_sc(seu = seu, cluster_annotation = opt$cluster_annotation, name = opt$name,
                    ndims = opt$ndims, resolution = opt$resolution, subset_by = opt$subset_by,
                    cell_annotation = opt$cell_annotation, minqcfeats = opt$minqcfeats, percentmt = opt$percentmt,
                    hvgs = opt$hvgs, scalefactor = opt$scalefactor, normalmethod = opt$normalmethod,
                    p_adj_cutoff = opt$p_adj_cutoff, verbose = opt$verbose, sigfig = 2,
                    output = opt$output, integrate = opt$integrate,
                    reduce = opt$reduce, SingleR_ref = SingleR_ref, ref_label = opt$ref_label,
                    ref_de_method = opt$ref_de_method, ref_n = opt$ref_n,
                    BPPARAM = BPPARAM, doublet_list = opt$doublet_list, integration_method = opt$int_method,
                    sketch = opt$sketch, sketch_pct = opt$sketch_pct,
                    sketch_method = opt$sketch_method, force_ncells = opt$force_ncells,
                    k_weight = opt$k_weight)
}

if(!file.exists(final_counts_path) & opt$integrate) {
  message("--------------------------------------------")
  message("----------SAVING RESULTS TO DISK------------")
  message("--------------------------------------------")
  write_annot_output(final_results = final_results, opt = opt)
}

if(file.exists(final_counts_path)) {
  message("Processing reconstruction of seurat object. Not launching QC report.")
} else {
  message("--------------------------------------------")
  message("------------WRITING QC REPORT---------------")
  message("--------------------------------------------")
  write_sc_report(final_results = final_results, template_folder = template_folder, source_folder = source_folder,
                  template = "sc_quality_control.txt", output = file.path(opt$output, "report"), out_suffix = "qc_report.html",
                  use_canvas = TRUE, opt = opt, params = params)
}

message("--------------------------------------------")
message("---------WRITING ANNOTATION REPORT----------")
message("--------------------------------------------")

write_sc_report(final_results = final_results, template_folder = template_folder,
                output = file.path(opt$output, "report"), source_folder = source_folder,
                template = "sc_annotation.txt", out_suffix = out_suffix, opt = opt, params = params)
