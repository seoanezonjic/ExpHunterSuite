#! /usr/bin/env Rscript


##########################################
## OPTION PARSER
##########################################

option_list <- list(
  optparse::make_option(c("-n", "--name"), type = "character", default = NULL,
            help = "Experiment name."),
  optparse::make_option("--input", type = "character", default = NULL,
            help = "Directory containing processed single-cell counts, metadata, embeddings and markers."),
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

opt <- process_sc_params(params, mode = "query")$opt

message("Reconstructing Seurat object from directory ", opt$input)
seu <- Seurat::CreateSeuratObject(counts = Seurat::Read10X(file.path(opt$input, "counts"), gene.column = 1),
                                  project = opt$name, min.cells = 1, min.features = 1)
seu_meta <- read.table(file.path(opt$input, "counts", "meta.tsv"), sep = "\t", header = TRUE)
rownames(seu_meta) <- colnames(seu)
seu <- Seurat::AddMetaData(seu, seu_meta, row.names("Cell_ID"))
seu$RNA$data <- seu$RNA$counts
embeddings <- read.table(file.path(opt$input, "embeddings", "cell_embeddings.tsv"), header = TRUE)
seu$UMAP_full <- Seurat::CreateDimReducObject(embeddings = as.matrix(embeddings), key = 'UMAPfull_', assay = 'RNA')
markers <- read.table(file.path(opt$input, "markers.tsv"), sep = "\t", header = TRUE)
query_data <- main_analyze_sc_query(seu = seu, query = opt$target_genes, sigfig = opt$sigfig)

query_results <- list(seu = seu, query_data = query_data)

message("--------------------------------------------")
message("------WRITING QUERY EXPRESSION REPORT-------")
message("--------------------------------------------")

write_sc_report(final_results = query_results, template_folder = template_folder,
                output = file.path(opt$output, "report"), template = "sc_query.txt",
                out_suffix = "query_report.html", opt = opt, params = params)
