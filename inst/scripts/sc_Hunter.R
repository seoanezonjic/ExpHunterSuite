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
  devtools::load_all(root_path)
  template_folder <- file.path(root_path, 'inst', 'templates')
}else{
  require('ExpHunterSuite')
  root_path <- find.package('ExpHunterSuite')
  template_folder <- file.path(root_path, 'templates')
}

option_list <- list(
  optparse::make_option(c("-n", "--name"), type = "character", default = NULL,
            help = "Experiment name."),
  optparse::make_option("--DEG_target", type = "character", default = "",
            help = "A string describing DEG analyses to perform."),
  optparse::make_option("--p_val_cutoff", type = "numeric", default = "5e-3",
            help = "Adjusted p-value cutoff."),
  optparse::make_option("--min_avg_log2FC", type = "numeric", default = 1e-3,
            help = "Avg log2fc cutoff for significant DEGs."),
  optparse::make_option("--min_cell_proportion", type = "numeric", default = 0.1,
            help = "Min percentage of cells expressing DEG in each group."),
  optparse::make_option("--min_counts", type = "numeric", default = 10,
            help = "Min counts to consider a gene is expressed in a cell."),
  optparse::make_option(c("-o", "--output"), type = "character", default = NULL,
            help = "Output folder."),
  optparse::make_option("--input", type = "character", default = NULL,
            help = "Directory containing counts"),
  optparse::make_option("--subset_by", type = "character", default = "",
            help = "Comma-separated list of conditions by which to perform integration
                    analysis. If empty, all conditions will be analysed."),
  optparse::make_option("--target_genes", type = "character", default = "",
            help = "Path to target genes table, or comma-separated list of target genes."),
  optparse::make_option("--cpu", type = "double", default = 1,
            help = "Provided CPUs."),
  optparse::make_option("--verbose", type = "logical", default = FALSE, action = "store_true",
            help = "Verbosity of base Seurat and harmony function calls.")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

if(opt$DEG_target == "") {
  message("DEG analysis disabled.")
  DEG_target <- NULL
} else {
  DEG_target <- unlist(strsplit(opt$DEG_target, ";"))
  DEG_target <- unlist(lapply(DEG_target, strsplit, split = ":"), recursive = FALSE)
  DEG_target <- lapply(DEG_target, function(vector) {
    if(length(vector) < 2) {
      res <- c(vector, "all")
    } else {
      res <- vector
    }
    return(res)
  })
  DEG_names <- unlist(lapply(DEG_target, `[[`, 1))
  DEG_names <- strsplit(DEG_names, ">")
  DEG_columns <- unlist(lapply(DEG_names, `[[`, 2))
  DEG_names <- unlist(lapply(DEG_names, `[[`, 1))
  DEG_values <- unlist(lapply(DEG_target, `[[`, 2))
  DEG_target <- data.frame(column = tolower(DEG_columns), values = DEG_values)
  rownames(DEG_target) <- DEG_names
}

if(opt$target_genes == ""){
  warning("No target genes provided")
  target_genes <- NULL
} else if(file.exists(opt$target_genes)) {
  target_genes <- read_and_format_targets(opt$target_genes)
} else {
  target_genes <- list(Custom = strsplit(opt$target_genes, split = ";")[[1]])
}

subset_by <- NULL
if(opt$subset_by == "") {
  message("No conditions specified for analysis. Data will not be subset.")
} else {
  subset_by <- tolower(unlist(strsplit(opt$subset_by, ",")))
  message(paste0("Selected ", length(subset_by), " subset condition(s): ", paste0(subset_by, collapse = ", ")))    
}

message("Reconstructing Seurat object from directory ", opt$input)

seu <- Seurat::CreateSeuratObject(counts = Seurat::Read10X(opt$input, gene.column = 1),
                                  project = opt$name, min.cells = 1, min.features = 1)
seu_meta <- read.table(file.path(opt$input, "meta.tsv"), sep = "\t", header = TRUE)
rownames(seu_meta) <- colnames(seu)
seu <- Seurat::AddMetaData(seu, seu_meta, row.names("Cell_ID"))
seu$RNA$data <- seu$RNA$counts
DEG_results <- main_sc_Hunter(seu = seu, DEG_target = DEG_target, p_val_cutoff = opt$p_val_cutoff,
                              min_avg_log2FC = opt$min_avg_log2FC, min_cell_proportion = opt$min_cell_proportion,
                              query = unlist(target_genes), subset_by = subset_by, output_path = opt$output,
                              min_counts = opt$min_counts, verbose = opt$verbose)
