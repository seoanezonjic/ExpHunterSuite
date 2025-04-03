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
  optparse::make_option("--input", type = "character", default = NULL,
            help = "Directory containing processed single-cell counts and metadata."),
  optparse::make_option("--embeddings", type = "character", default = NULL,
            help = "Directory containing embedding information."),
  optparse::make_option("--target_genes", type = "character", default = "",
            help = "Path to target genes table, or comma-separated list of target genes."),
  optparse::make_option("--sigfig", type = "integer", default = 3,
            help = "Significant figures by which to round off results."),
  optparse::make_option(c("-o", "--output"), type = "character", default = NULL,
            help = "Output folder."),
  optparse::make_option("--cpu", type = "double", default = 1,
            help = "Provided CPUs."),
  optparse::make_option("--verbose", type = "logical", default = FALSE, action = "store_true",
            help = "Verbosity of base Seurat and harmony function calls."),
  optparse::make_option("--extra_columns", type = "character", default = "",
            help = "Comma-separated list of extra conditions to represent in certain plots.")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

if(opt$target_genes == ""){
  warning("No target genes provided")
  target_genes <- NULL
} else if(file.exists(opt$target_genes)) {
  target_genes <- read_and_format_targets(opt$target_genes)
} else {
  target_genes <- strsplit(opt$target_genes, split = ";")[[1]]
}

if(opt$extra_columns == "") {
  extra_columns <- NULL
} else {
  extra_columns <- tolower(unlist(strsplit(opt$extra_columns, ",")))
}

seu <- Seurat::CreateSeuratObject(counts = Seurat::Read10X(opt$input, gene.column = 1),
                                  project = opt$name, min.cells = 1, min.features = 1)
seu_meta <- read.table(file.path(opt$input, "meta.tsv"), sep = "\t", header = TRUE)
rownames(seu_meta) <- colnames(seu)
seu <- Seurat::AddMetaData(seu, seu_meta, row.names("Cell_ID"))
seu$RNA$data <- seu$RNA$counts
embeddings <- read.table(file.path(opt$embeddings, "cell_embeddings.txt"), header = TRUE)
seu$UMAP_full <- Seurat::CreateDimReducObject(embeddings = as.matrix(embeddings), key = 'UMAPfull_', assay = 'RNA')
query_data <- main_analyze_sc_query(seu = seu, query = target_genes, sigfig = opt$sigfig)

query_results <- list(seu = seu, query_data = query_data)
message("--------------------------------------------")
message("-------Writing query analysis report--------")
message("--------------------------------------------")

write_sc_report(final_results = query_results, template_folder = template_folder, output = file.path(opt$output, "report"),
                query = unlist(target_genes), extra_columns = extra_columns,
                template = "sc_analysis.txt", out_name = "query_report.html", opt = opt)
