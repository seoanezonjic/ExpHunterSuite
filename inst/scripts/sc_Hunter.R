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

option_list <- list(
  optparse::make_option("--DEG_target", type = "character", default = "",
            help = "Columns for DEG analysis."),
  optparse::make_option("--mincells", type = "integer", default = NULL,
            help = "Min number of cells for which a feature was recorded."),
  optparse::make_option("--minfeats", type = "integer", default = NULL,
            help = "Min number of features for which a cell was recorded."),
  optparse::make_option("--p_adj_cutoff", type = "numeric", default = "5e-3",
            help = "Adjusted p-value cutoff."),
  optparse::make_option("--min_avg_log2FC", type = "numeric", default = 1e-3,
            help = "Avg log2fc cutoff for significant DEGs."),
  optparse::make_option("--min_cell_pct", type = "numeric", default = 0.1,
            help = "Min percentage of cells expressing DEG in each group.")
  optparse::make_option(c("-o", "--output"), type = "character", default = NULL,
            help = "Output folder."),
  optparse::make_option("--input", type = "character", default = NULL,
            help = "Count results folder. Input will be read from here."),
  optparse::make_option("--subset_by", type = "character", default = "",
            help = "Comma-separated list of conditions by which to perform integration
                    analysis. If empty, all conditions will be analysed."),
  optparse::make_option("--target_genes", type = "character", default = "",
            help = "Path to target genes table, or comma-separated list of target genes."),
  optparse::make_option("--cpu", type = "double", default = 1,
            help = "Provided CPUs."),
  optparse::make_option("--verbose", type = "logical", default = FALSE, action = "store_true",
            help = "Verbosity of base Seurat and harmony function calls."),

)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

if(opt$target_genes == ""){
  warning("No target genes provided")
  target_genes <- NULL
} else if(file.exists(opt$target_genes)) {
  target_genes <- read_and_format_targets(opt$target_genes)
} else {
  target_genes <- list(Custom = strsplit(opt$target_genes, split = ";")[[1]])
}

if(opt$subset_by == "") {
  message("No conditions specified for analysis. Data will not be subset.")
} else {
  subset_by <- tolower(unlist(strsplit(opt$subset_by, ",")))
  message(paste0("Selected ", length(subset_by), " subset condition(s): ", paste0(subset_by, collapse = ", ")))    
}

message("Reconstructing Seurat object from directory ", opt$input)

seu <- Seurat::CreateSeuratObject(counts = Seurat::Read10X(opt$input, gene.column = 1),
                                  project = opt$name, min.cells = 1, min.features = 1)
seu_meta <- read.table(file.path(opt$imported_counts, "meta.tsv"), sep = "\t", header = TRUE)
rownames(seu_meta) <- colnames(seu)
seu <- Seurat::AddMetaData(seu, seu_meta, row.names("Cell_ID"))
