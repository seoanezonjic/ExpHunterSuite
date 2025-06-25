#! /usr/bin/env Rscript


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

# Auxiliary function

col_to_table <- function(column, col_names) {
  res <- as.data.frame(table(column))
  colnames(res) <- col_names
  return(res)
}

option_list <- list(
  optparse::make_option(c("-r", "--reference"), type = "character",
    help = "SingleR reference to use."),
  optparse::make_option(c("-v", "--version"), type = "character",
    help = "Celldex version of reference."),
  optparse::make_option(c("-o", "--output"), type = "character",
    help = "Output directory."),
  optparse::make_option("--replace", type = "logical", default = FALSE, action = "store_true",
  	help = "Replace existing directory"),
  optparse::make_option("--verbose", type = "logical", default = TRUE, action = "store_true",
  	help = "Display progress"),
  optparse::make_option("--quiet", type = "logical", default = FALSE, action = "store_false",
  	dest = "verbose", help = "Display progress"),
  optparse::make_option("--database", type = "character", help = "Database to consult
    (\"celldex\" or \"scRNAseq\") or a path to local reference"),
  optparse::make_option("--ref_label", type = "character", default = NULL,
    help = "Column of reference metadata to use for annotation. Only used in scRNAseq mode"),
  optparse::make_option("--cpu", type = "integer", default = NULL,
    help = "CPUs to use when calling data.table::fread"),
  optparse::make_option("--only_showcase", type = "logical", default = FALSE, action = "store_true",
    help = "Simply load reference and produce showcase report")
)  

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
opt$name <- basename(opt$reference)
opt$output <- file.path(opt$output, basename(opt$reference))
if(opt$version != "") {
  opt$output <- paste(output, opt$version, sep = "_")
}

if(!opt$only_showcase) {
  if(is.null(opt$reference)) {
    stop('No reference provided. Please see get_SingleR_ref.R --help')
  }
  if(is.null(opt$output)) {
    stop('No output path provided. Please see get_SingleR_ref.R --help')
  }
  if(file.exists(opt$database)) {
    message("Generating reference from specified database")
    meta_file <- Sys.glob(file.path(opt$database, "metadata/meta*.txt"))
    if(length(meta_file) > 1) {
      stop('More than one match for metadata file. Please ensure only one metadata
        file matches the expression "meta*.txt"')
    }
    metadata <- read.table(meta_file, sep = '\t', header = TRUE)
    expr_files <- Sys.glob(paste0(opt$database, "/expression/*txt*"))
    if(length(expr_files) > 0) {
      data.table::setDTthreads(threads = opt$CPU)
      expr_matrices <- vector(mode = "list", length = length(expr_files))
      for(file in seq(expr_files)) {
        message("Loading expression file ", file, " of ", length(expr_files))
        expr_matrix <- as.matrix(data.table::fread(expr_files[file], verbose = TRUE))
        message("File loaded. Processing and converting to sparse matrix.")
        row_names <- expr_matrix[, 1]
        expr_matrix <- expr_matrix[, -1]
        expr_matrix <- apply(expr_matrix, 2, as.numeric)
        rownames(expr_matrix) <- row_names
        sparse_matrix <- Matrix::Matrix(expr_matrix, sparse = TRUE)
        sparse_matrix <- sparse_matrix[, Matrix::colSums(sparse_matrix) > 0]
        expr_matrices[[file]] <- sparse_matrix
        message(paste0("File ", file, " converted"))
      }
      read_sparse_matrix <- do.call(cbind, expr_matrices)
    } else {
      read_sparse_matrix <- NULL
    }
    tenX_dirs <- dirname(Sys.glob(paste0(opt$database, "/expression/*/*mtx*")))
    if(length(tenX_dirs) > 0) {
      tenX_matrices <- vector(mode = "list", length = length(tenX_dirs))
      for(tenX_dir in seq(tenX_dirs)) {
        message(paste0("Loading 10X directory ", tenX_dir, " of ", length(tenX_dirs)))
        tenX_matrix <- Seurat::Read10X(tenX_dirs[tenX_dir], gene.column = 1)
        tenX_matrices[[tenX_dir]] <- tenX_matrix
      }
      message("10X matrices loaded. Merging (this may take a while)")
      merged_tenX_matrix <- SeuratObject::RowMergeSparseMatrices(tenX_matrices[[1]], tenX_matrices[-1])
    } else {
      read_sparse_matrix <- NULL
    }
    if(!is.null(read_sparse_matrix) & !is.null(merged_tenX_matrix)) {
     final_sparse_matrix <- SeuratObject::RowMergeSparseMatrices(read_sparse_matrix, merged_tenX_matrix)  
    } else {
      if(is.null(read_sparse_matrix)) {
        final_sparse_matrix <- merged_tenX_matrix
      } else {
        final_sparse_matrix <- read_sparse_matrix
      }
    }
    final_sparse_matrix <- final_sparse_matrix[, order(colnames(final_sparse_matrix))]
    metadata <- metadata[metadata$NAME %in% colnames(final_sparse_matrix), ]
    metadata <- metadata[order(metadata$NAME), ]
    rownames(metadata) <- metadata$NAME
    metadata <- metadata[, -1]
    ref <- SummarizedExperiment::SummarizedExperiment(assays=list(counts = final_sparse_matrix), colData = metadata)
    ref <- scater::logNormCounts(ref)
  }
  if(opt$database == "scRNAseq") {
    ref <- scRNAseq::fetchDataset(opt$reference, opt$version)
    # Removing unlabelled cells or cells without a clear label.
    ref <- ref[, !is.na(ref[[opt$ref_label]]) & ref[[opt$ref_label]]!="unclear"] 
    ref <- scater::logNormCounts(ref)
  }
  if(opt$database == "celldex") {
    ref <- celldex::fetchReference(opt$reference, opt$version)
  }
  HDF5Array::saveHDF5SummarizedExperiment(x = ref, dir = opt$output, verbose = opt$verbose, replace = opt$replace)
  message(paste0("Reference saved successfully in ", opt$output))
  message("Converting reference to seurat object to calculate UMAP")
  sce <- as(ref, "SingleCellExperiment")
  sce@assays@data$counts <- as(sce@assays@data$counts, "dgCMatrix")
  sce@assays@data$logcounts <- as(sce@assays@data$logcounts, "dgCMatrix")
  seu <- Seurat::as.Seurat(sce)
  seu <- Seurat::NormalizeData(object = seu, verbose = FALSE, normalization.method = "LogNormalize", scale.factor = 1e6)
  seu <- Seurat::FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE, selection.method = "vst", assay = "originalexp")
  seu <- Seurat::ScaleData(object = seu, verbose = FALSE, features = rownames(seu))
  seu <- Seurat::FindNeighbors(object = seu, dims = seq(10), assay = "originalexp", reduction = "PCA", verbose = FALSE)
  seu <- Seurat::RunUMAP(object = seu, dims = seq(10), verbose = FALSE, reduction = "PCA", return.model = TRUE)
  message("Saving results to disk")
  dir.create(file.path(opt$output, "counts"))
  dir.create(file.path(opt$output, "embeddings"))
  write_annot_output(final_results = list(seu = seu), opt = opt, assay = "originalexp")
} else {
  message("Loading SingleR ref for showcase report")
  seu <- Seurat::CreateSeuratObject(counts = Seurat::Read10X(file.path(opt$output, "counts"), gene.column = 1),
                                  project = opt$name, min.cells = 1, min.features = 1)
  seu_meta <- read.table(file.path(opt$output, "counts/meta.tsv"), sep = "\t", header = TRUE)
  rownames(seu_meta) <- colnames(seu)
  seu <- Seurat::AddMetaData(seu, seu_meta, row.names("Cell_ID"))
  seu$RNA$data <- seu$RNA$counts
  markers <- read.table(file.path(opt$output, "markers.tsv"), sep = "\t", header = TRUE)
  embeddings <- read.table(file.path(opt$output, "embeddings", "cell_embeddings.tsv"), header = TRUE)
  seu$umap <- Seurat::CreateDimReducObject(embeddings = as.matrix(embeddings), key = 'umap_', assay = 'RNA')
}

message("Preparing report showcasing selected reference")

meta <- seu@meta.data
fields_to_remove <- grepl("orig.ident|nCount|nFeature|sizeFactor", names(meta))
meta <- meta[, !fields_to_remove]
tables <- lapply(meta, col_to_table, col_names = c("Type", "Frequency"))
names(tables) <- colnames(meta)

write_sc_report(final_results = list(tables = tables, seu = seu), template_folder = template_folder, output = getwd(),
                template = "sc_ref_showcase.txt", out_suffix = "ref_showcase.html", opt = opt)
