#! /usr/bin/env Rscript


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
  optparse::make_option("--simple_DEG_pct", type = "logical", default = FALSE, action = "store_true",
            help = "Do not override Seurat's default handling of min.pct argument (passed as min_cell_proportion in this script)."), 
  optparse::make_option("--top_N", type = "numeric", default = 10,
            help = "Top N DEGs to represent in certain plots"),
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

params <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

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

opt <- process_sc_params(params, mode = "DEG")$opt
if(!file.exists(opt$input)) {
  stop("No counts file found, please ensure the annotation module has been run. Path searched: ", opt$input)
}

if(opt$targets_folder == "" | !file.exists(opt$targets_folder)) {
  stop(paste0("Provided targets directory does not exist. Was ", opt$targets_folder))
} else {
  target_files <- dir(opt$targets_folder, full.names = TRUE)
  DEG_targets <- vector(mode = "list", length = length(target_files))
  names(DEG_targets) <- target_files
  for(target_file in target_files) {
    DEG_targets[[target_file]] <- read.table(target_file, sep = "\t", header = TRUE)
  }
  names(DEG_targets) <- tools::file_path_sans_ext(basename(target_files))
}

load_targets <- list()
for(target in names(DEG_targets)) {
  path <- file.path(opt$output, "DEG", gsub("_target", "", target))
  if(file.exists(path)) {
    load_targets <- c(load_targets, DEG_targets[names(DEG_targets) == target])
    DEG_targets <- DEG_targets[names(DEG_targets) != target]
  }
}

DEG_list <- NULL
if(length(DEG_targets) > 0 | isTRUE(opt$recalc_query)) {
  message("Reconstructing Seurat object from directory ", opt$input)
  seu <- SeuratObject::CreateSeuratObject(counts = Seurat::Read10X(opt$input, gene.column = 1),
                                          project = opt$name, min.cells = 1, min.features = 1)
  seu_meta <- read.table(file.path(opt$input, "meta.tsv"), sep = "\t", header = TRUE)
  rownames(seu_meta) <- colnames(seu)
  seu <- Seurat::AddMetaData(seu, seu_meta, row.names("Cell_ID"))
  seu$RNA$data <- seu$RNA$counts
} else {
  seu <- NULL
}

if(length(DEG_targets) > 0) {
  DEG_list <- parallel_list(X = DEG_targets, FUN = main_sc_Hunter, workers = opt$cpu, seu = seu,
                            p_val_cutoff = opt$p_val_cutoff, min_avg_log2FC = opt$min_avg_log2FC,
                            min_cell_proportion = opt$min_cell_proportion, top = opt$top_N,
                            query = opt$target_genes, output_path = opt$output,
                            min_counts = opt$min_counts, verbose = opt$verbose)
  names(DEG_list) <- unlist(strsplit(names(DEG_targets), "_target"))
  parallel_list(X = names(DEG_list), FUN = write_DEG_output, workers = opt$cpu, DEG_list = DEG_list, opt = opt)
  message("Execution successful. Writing parameters to ", params$output, "/execution_parameters.txt")
  write_opt(opt = params, output = file.path(params$output, "execution_parameters.txt"))
}

load_DEG_list <- NULL
if(length(load_targets) > 0) {
  message("Reading DEG results from disk")
  target_names <- unlist(strsplit(names(load_targets), "_target"))
  load_DEG_list <- load_DEG_output(targets = target_names, load_path = opt$output, recalc_query = opt$recalc_query,
                                   min_avg_log2FC = opt$min_avg_log2FC, min_cell_proportion = opt$min_cell_proportion,
                                   p_val_cutoff = opt$p_val_cutoff, seu = seu)
}

DEG_list <- c(DEG_list, load_DEG_list)
DEG_targets <- c(DEG_targets, load_targets)
write_opt(opt = params, output = file.path(params$output, "execution_parameters.txt"))

message("--------------------------------------------")
message("---------WRITING SC HUNTER REPORTS----------")
message("--------------------------------------------")
for(target_name in names(DEG_targets)) {
  DEG_name <- unlist(strsplit(target_name, "_target"))
  message(paste0("Writing ", DEG_name, " report"))
  DEG_results <- list(DEG_list = DEG_list[[DEG_name]], opt = opt, target = DEG_targets[[target_name]], target_name = DEG_name)
  write_sc_report(final_results = DEG_results, opt = opt, template_folder = template_folder, analysis = "Single-Cell Differential Expression",
                  output = file.path(opt$output, "report"), template = "sc_DEGs.txt", out_suffix = paste0(DEG_name, "_DEG_report.html"), params = params)
}
