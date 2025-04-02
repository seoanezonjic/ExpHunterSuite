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
  optparse::make_option("--targets_folder", type = "character", default = "",
            help = "Directory containing target files"),
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
            help = "Directory containing processed single-cell counts and metadata."),
  optparse::make_option("--target_genes", type = "character", default = "",
            help = "Path to target genes table, or comma-separated list of target genes."),
  optparse::make_option("--verbose", type = "logical", default = FALSE, action = "store_true",
            help = "Verbosity of base Seurat and harmony function calls."),
  optparse::make_option("--cpu", type = "integer", default = 1,
            help = "Provided CPUs.")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

if(opt$targets_folder == "" | !file.exists(opt$targets_folder)) {
  stop(paste0("Provided targets directory does not exist. Was ", opt$targets_folder))
} else {
  target_files <- dir(opt$targets_folder, full.names = TRUE)
  DEG_targets <- vector(mode = "list", length = length(target_files))
  names(DEG_targets) <- target_files
  for(target_file in target_files) {
    DEG_targets[[target_file]] <- read.table(target_file, header = TRUE)
  }
  names(DEG_targets) <- tools::file_path_sans_ext(basename(target_files))
}

if(opt$cpu > 1) {
  BiocParallel::register(BiocParallel::MulticoreParam(opt$cpu))
  BPPARAM <- BiocParallel::MulticoreParam(opt$cpu)
} else {
  BPPARAM <- BiocParallel::SerialParam()
}

if(opt$target_genes == ""){
  message("No target genes provided")
  target_genes <- NULL
} else if(file.exists(opt$target_genes)) {
  target_genes <- read_and_format_targets(opt$target_genes)
} else {
  target_genes <- strsplit(opt$target_genes, split = ";")[[1]]
}

message("Reconstructing Seurat object from directory ", opt$input)

seu <- Seurat::CreateSeuratObject(counts = Seurat::Read10X(opt$input, gene.column = 1),
                                  project = opt$name, min.cells = 1, min.features = 1)
seu_meta <- read.table(file.path(opt$input, "meta.tsv"), sep = "\t", header = TRUE)
rownames(seu_meta) <- colnames(seu)
seu <- Seurat::AddMetaData(seu, seu_meta, row.names("Cell_ID"))
seu$RNA$data <- seu$RNA$counts
save.image('testing.RData')
stop('test')
DEG_results <- BiocParallel::bplapply(X = DEG_targets, FUN = main_sc_Hunter, BPPARAM = BPPARAM, seu = seu,
                                      p_val_cutoff = opt$p_val_cutoff, min_avg_log2FC = opt$min_avg_log2FC,
                                      min_cell_proportion = opt$min_cell_proportion, query = target_genes,
                                      output_path = opt$output, min_counts = opt$min_counts, verbose = opt$verbose)
names(DEG_results) <- unlist(strsplit(names(DEG_targets), "_target"))
