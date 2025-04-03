
#' main_annotate_sc
#' is the main seurat analysis function, which can be used for integrative or
#' non-integrative analysis.
#'
#' @param seu A seurat object.
#' @param name Project name. Default NULL (no project name)
#' @param minqcfeats An integer. Minimum features to consider a cell valid
#' @param percentmt A float. Maximun MT percentage to consider a cell valid
#' @param query A string vector. List of genes to explore in dataset
#' @param sigfig An integer. Significant figures to output
#' @param resolution An integer. Controls clustering granularity. Default 0.5
#' @param p_adj_cutoff A float. Adjusted p-value cutoff by which to consider a
#' marker valid for cell type annotation.
#' @param integrate A boolean.
#'   * `TRUE`: Integrate seurat object using the `harmony` package.
#'   * `FALSE` (the default): Do not perform integration.
#' @param cluster_annotation A data frame. Table to use to rename clusters.
#' @param cell_annotation A data frame. Table of markers to use for cell type
#' annotation.
#' @param scalefactor An integer. Factor by which to scale data in normalisation
#' @param hvgs An integer. Number of highly-variable features to select.
#' default 2000.
#' @param subset_by A string vector. Categories to consider in integrative
#' analysis.
#' @param normalmethod A string. Method to use in normalisation. Default is
#' "LogNormalize", the Seurat default.
#' @param ndims An integer. Target dimensions in dimensionality reduction.
#' Default 10.
#' @param verbose A boolean.
#'   * `TRUE`: Prints progress bars and messages to closely monitor progress.
#'   * `FALSE` (the default): Print minimal progress messages.
#' @param output A string. Path to output directory.
#' @param save_RDS A boolean.
#'   * `TRUE`: Save output as an rds file.
#'   * `FALSE` (the default): Do not save output as an rds file.
#' @param reduce A boolean.
#'   * `TRUE`: Skip QC filtering. Intended for development and testing.
#'   * `FALSE` (the default): QC filtering will be performed.
#' @param SingleR_ref SummarizedExperiment object to use as reference
#' for SingleR cell type annotation. If NULL (default value),
#' SingleR will not be used.
#' @param ref_label Reference column to use in annotation.
#' @param ref_de_method Method to use for marker calculation in single-cell
#' reference.
#' @param ref_n Top N reference markers to consider in annotation. Higher values
#' provide a more accurate annotation, but increase noise and computational 
#' time. Will not be used if ref_de_method is empty.
#' @param doublet_list A vector containing barcodes to be marked as doublets.
#' NULL by default. Per-sample analysis finds this vector
#' for every sample, integrative mode requires this vector.
#' @param BPPARAM Parameters to pass to BiocParallel framework.
#' @param integration_method A string. Method to use in integration. "Harmony"
#' (the default), "RPCA", "CCA", "FastMNN" or "scVI".
#' @param sketch A boolean. If TRUE, data will be sketched. If FALSE
#' (the default), data will be analysed as-is.
#' @param sketch_pct A numeric. Percentage of total cells to consider
#' representative of the experiment. Default 12, as suggested by sketching
#' tutorial.
#' @param sketch_method A string. Method to use in score calculation for
#' sketching.
#' @param k_weight Number of neighbors to consider when weighting anchors. Used
#' in integration.
#' @export
#' @examples
#'  \dontrun{
#'    main_annotate_sc(seu = seurat_object, minqcfeats = 500, percenmt = 5,
#'                   query = "TREM2", sigfig = 2, resolution = 0.5,
#'                   p_adj_cutoff = 5e-3, name = "project_name",
#'                   integrate = TRUE, cluster_annotation = NULL,
#'                   cell_annotation = cell_types, DEG_target = "genotype",
#'                   scalefactor = 10000, hvgs = 2000, subset_by = "genotype",
#'                   normalmethod = "LogNormalize", ndims = 10, verbose = FALSE,
#'                   output = getwd(), save_RDS = FALSE, reduce = FALSE,
#'                   ref_label = NULL, SingleR_ref = NULL, ref_de_method = NULL,
#'                   ref_n = NULL, BPPARAM = NULL, doublet_list = NULL,
#'                   integration_method = "Harmony", sketch_pct = 12,
#'                   sketch_ncells = 5000, DEG_p_val_cutoff = 5e-3,
#'                   min_avg_log2FC = 0.5, k_weight = 90)
#'  }
#' @returns invisible(final_results list. Contains multiple items:)
#' * qc: seurat object prior to filtering and analysis.
#' * seu: processed seurat object.
#' * sample_qc_pct: Gene expression matrix of union of top N expressed genes in
#'                  all samples.
#' * clusters_pct: A data frame with expression levels for query genes in each
#'                 cluster (or cell type, if annotated).
#' * query_exp: A data frame with expression levels for query genes in each
#'              sample.
#' * query_pct: A data frame with expression levels for query genes in each 
#'              sample.
#' * query_cluster_pct: Same as query_pct, but subset by one or even two
#'                      conditions.
#' * markers: Data frame of marker genes.
#' * SingleR_annotation: Trained SingleR annotation object.
#' * DEG_list: List of differentially expressed genes across speficied
#'             conditions. Performed globally and in each cluster or cell type
#'             subset.
#' * subset_seu: Analyzed seurat objects, subset by integration condition.
#' * subset_DEGs: Same as DEG_list, but subset by integration condition.
#' * integrate: Whether or not to trigger integrative analysis. Workflow
#'              contains a few key differences depending on this argument.
#' @details subset_seu and subset_DEGs only trigger if two conditions have
#' been provided for integrative analysis, and they will contain two objects,
#' each only containing samples corresponding to one of the two possible values
#' in the second experimental condition. That allows isolating the effect
#' of each experimental condition to be analyzed and compared separately.

main_annotate_sc <- function(seu, minqcfeats, percentmt, query, sigfig = 2,
  resolution = 0.5, p_adj_cutoff = 5e-3, name = NULL, integrate = FALSE,
  cluster_annotation = NULL, cell_annotation = NULL, scalefactor = 10000,
  hvgs = 2000, subset_by = NULL, ndims = 10, normalmethod = "LogNormalize",
  verbose = FALSE, output = getwd(), save_RDS = FALSE, reduce = FALSE,
  ref_label = NULL, SingleR_ref = NULL, ref_de_method = NULL, ref_n = NULL,
  BPPARAM = NULL, doublet_list = NULL, integration_method = "Harmony",
  sketch = FALSE, sketch_pct = 25, k_weight = 100, sketch_ncells = 5000,
  force_ncells = NA_integer_, sketch_method = "LeverageScore",
  doublet_path = getwd()){
  #check_sc_input(metadata = seu@meta.data, DEG_target = DEG_target)
  qc <- tag_qc(seu = seu, minqcfeats = minqcfeats, percentmt = percentmt,
               doublet_list = doublet_list)
  if(!reduce) {
    seu <- subset(qc, subset = qc == 'Pass')
    aggr.ref <- FALSE
    fine.tune <- TRUE
  } else {
    message(paste0("Reduce argument is set to TRUE. Skipping QC subsetting. ",
                   "Updating SingleR configuration"))
    seu <- qc
    aggr.ref <- TRUE
    fine.tune <- FALSE
    k_weight <- 30
  }
  if(integrate) {
    message(paste0("Splitting seurat object by sample."))
    seu[["RNA"]] <- split(seu[["RNA"]], f = seu$sample)
  }
  message('Normalizing data')
  norm_start <- Sys.time()
  seu <- Seurat::NormalizeData(object = seu, verbose = verbose,
  									           normalization.method = normalmethod,
                               scale.factor = scalefactor)
  norm_end <- Sys.time()
  if(verbose) {
    message(paste0("Normalization time: ", norm_end-norm_start))
  }
  message('Finding variable features')
  seu <- Seurat::FindVariableFeatures(seu, nfeatures = hvgs, verbose = verbose,
                                      selection.method = "vst", assay = "RNA")
  if(sketch) {
    seu <- process_sketch(seu = seu, sketch_method = sketch_method,
                          sketch_pct = sketch_pct, force_ncells = force_ncells,
                          hvgs = hvgs, verbose = verbose)
  }
  assay <- Seurat::DefaultAssay(seu)
  message('Scaling data')
  scale_start <- Sys.time()
  seu <- Seurat::ScaleData(object = seu, verbose = verbose,
                           features = rownames(seu))
  scale_end <- Sys.time()
  if(verbose) {
    message(paste0("Scaling time: ", scale_end - scale_start))
  }
  message('Reducing dimensionality')  
  seu <- Seurat::RunPCA(seu, assay = assay, npcs = ndims, verbose = verbose)
  reduction <- "pca"
  if(integrate) {
  	message('Integrating seurat object')
    seu <- Seurat::IntegrateLayers(object = seu, orig.reduction = "pca",
      new.reduction = integration_method, verbose = FALSE, dims = seq(1, ndims),
      method = paste0(integration_method, "Integration"), assay = assay,
      scale.layer = "scale.data", k.weight = k_weight)
    reduction <- integration_method
  }
  seu <- Seurat::FindNeighbors(object = seu, dims = seq(1, ndims),
                               assay = assay, reduction = reduction,
                               verbose = verbose)
  if(is.null(SingleR_ref) | !integrate) {
    seu <- Seurat::FindClusters(seu, resolution = resolution, verbose = verbose)
    # Seurat starts counting clusters from 0, which is the source of many
    # headaches when working in R, which starts counting from 1. Therefore,
    # we introduce this correction. Weirdly enough, coercing it to numeric
    # already adds 1.
    seu@meta.data$seurat_clusters <- as.numeric(seu@meta.data$seurat_clusters)
    Seurat::Idents(seu) <- seu@meta.data$seurat_clusters
  } else {
    message(paste0("Annotation by clusters not active and multiple samples",
                   " detected. Skipping clustering"))
  }
  seu <- Seurat::RunUMAP(object = seu, dims = seq(ndims), reduction = reduction,
                         return.model = T, verbose = verbose)
  if(!integrate) {
    processed_doublets <- process_doublets(seu = seu, qc = qc,
                                           doublet_path = doublet_path)
  }
  seu <- SeuratObject::JoinLayers(seu)
  annot_start <- Sys.time()
  annotation <- annotate_seurat(seu = seu, SingleR_ref = SingleR_ref,
    cell_annotation = cell_annotation, ref_n = ref_n, subset_by = subset_by,
    cluster_annotation = cluster_annotation, p_adj_cutoff = p_adj_cutoff,
    BPPARAM = BPPARAM, ref_de_method = ref_de_method, aggr.ref = aggr.ref,
    fine.tune = fine.tune, assay = assay, integrate = integrate,
    verbose = verbose)
  annot_end <- Sys.time()
  if(verbose) {
    message(paste0("Time to annotate: ", annot_end - annot_start))
  }
  seu <- annotation$seu
  markers <- annotation$markers
  if("sketch" %in% names(seu@assays)) {
    seu <- project_sketch(seu = seu, reduction = reduction, ndims = ndims)
  }
  assay <- "RNA"
  message("Extracting expression quality metrics")
  sample_qc_pct <- get_qc_pct(seu, by = "sample")
  message("Extracting query expression metrics. This might take a while.")
  clusters_pct <- get_clusters_distribution(seu = seu, sigfig = sigfig)
  if(!is.null(query)) {
    query_data <- analyze_query(seu = seu, query = query, sigfig = sigfig)
  } else {
    query_data <- NULL
  }
  final_results <- list()
  final_results$qc <- qc
  final_results$seu <- seu
  final_results$sample_qc_pct <- sample_qc_pct
  final_results$clusters_pct <- clusters_pct
  final_results$markers <- markers
  final_results$SingleR_annotation <- annotation$SingleR_annotation
  final_results$integrate <- integrate
  if(save_RDS){
    message('Writing results to disk.')
    saveRDS(final_results, file.path(output,
                           paste0(seu@project.name, ".final_results.rds")))
  }
  return(final_results)
}

#' main_sc_Hunter
#' `main_sc_Hunter` is the main seurat analysis function. Can be used
#' for integrative or non-integrative analysis.
#'
#' @param seu Seurat object to analyze.
#' @param DEG_target A string vector. Categories by which DEG analysis will be
#' performed
#' @param p_val_cutoff Adjusted p-val cutoff for significant DEGs.
#' Default 5e-3.
#' @param min_avg_log2FC Average log2fc cutoff for significant DEGs.
#' Default 0.5.
#' @param min_cell_proportion Minimum threshold of percentage of cells
#' expressing DEG. Despite Seurat::FindMarkers being called "min.pct", it is
#' not a percentage but a value from 0 to 1.
#' @param min_counts Min counts to consider a gene is expressed in a cell.
#' @param verbose Print extra execution information.
#' @param output_path Path where output will be written.
#' @export

main_sc_Hunter <- function(DEG_target = stop("No DEG target supplied"),
                           seu = stop("No seurat object supplied."),
                           p_val_cutoff = 1e-3, min_avg_log2FC = 0.5,
                           min_counts = 10, verbose = FALSE, query = NULL,
                           min_cell_proportion = 1, output_path = getwd()) {
  message('Starting DEG analysis.')
  seu <- seu[, which(seu$sample %in% DEG_target$sample)]
  seu$deg_group <- DEG_target$treat[match(seu$sample, DEG_target$sample)]
  subset_target <- ifelse("cell_type" %in% colnames(seu@meta.data),
                           yes = "cell_type", no = "seurat_clusters")
  DEGs <- get_sc_markers(seu = seu, cond = "deg_group", DEG = TRUE,
                         subset_by = subset_target, verbose = verbose,
                         values = "Ctrl,Treat", min.pct = min_cell_proportion)
  message("Extracting DEG cell metrics")
  ## Here we need to control a list of FALSE data frames being returned
  ## from get_sc_markers, it's a possible case if cluster-conditon
  ## combinations do not allow a differential analysis for specified
  ## condition.
  DEG_metrics <- get_fc_vs_ncells(seu = seu, DEG_list = DEGs$markers,
                  min_avg_log2FC = min_avg_log2FC, p_val_cutoff = p_val_cutoff,
                  min_counts = min_counts)
  if(!is.null(query)) {
    DEG_query <- get_fc_vs_ncells(seu = seu, DEG_list = DEGs$markers,
                                  min_counts = min_counts, query = query)
  }
  return(list(DEGs = DEGs, DEG_metrics = DEG_metrics, DEG_query = DEG_query))
}

main_analyze_sc_query <- function(seu, query, sigfig, output_path) {
  query_data <- analyze_sc_query(seu = seu, query = query, sigfig = sigfig)
  query_results <- list()
  query_results$query_exp <- query_data$query_exp
  query_results$query_pct <- query_data$query_pct
  query_results$query_cluster_pct <- query_data$query_cluster_pct
  return(query_results)
}

#' write_annot_output
#' writes final counts matrix and annotated metadata of a seurat experiment.
#'
#' @param final_results final_results object output from main_annotate_sc
#' @param assay,layer Assay and layer from where counts will be retrieved.
#' @param opt Options used in script call.
#' @export
#' @returns invisible(NULL)

write_annot_output <- function(final_results = stop("Missing results object"),
                               opt = NULL, assay = "RNA", layer = "data") {
  if(is.null(opt$output)) {
    opt$output <- getwd()
  }
  seu <- final_results$seu
  reduction <- ifelse("umap.full" %in% names(seu), yes = "umap.full",
                       no = "umap")
  metadata <- seu@meta.data
  reduction <- seu[[reduction]]
  counts <- Seurat::GetAssayData(seu, assay = assay, layer = layer)
  DropletUtils::write10xCounts(path = file.path(opt$output, "counts"),
                              x = counts, overwrite = TRUE, genome = opt$genome,
                              gene.type = "Gene Expression")
  write.table(metadata, sep = "\t", quote = FALSE, row.names = TRUE,
              file = file.path(opt$output, "counts", "meta.tsv"))
  write.table(reduction@cell.embeddings, sep = "\t", quote = FALSE,
              row.names = TRUE,
              file = file.path(opt$output, "embeddings", "cell_embeddings.txt"))
  message(paste0("Counts matrix saved to ", file.path(opt$output, "counts")))
  return(invisible(NULL))
}

#' write_sc_report
#' Write integration HTML report
#'
#' @inheritParams main_annotate_sc
#' @param final_results Output from main_sc_Hunter function
#' @param output directory where report will be saved
#' @param name experiment name, will be used to build output file name
#' @param out_name suffix to add to output file name. Useful when rendering
#' different templates with this function (as is the case in our workflow)
#' @param template_folder directory where template is located
#' @param template Template to render
#' @param source_folder htmlreportR source folder
#' @param subset_by factors present in experiment design
#' @param use_canvas Parameter to select whether or not CanvasXpress plots will
#' be triggered in templates where this control parameter has been implemented.
#' Setting it to FALSE can be useful for big datasets, as CanvasXpress might
#' have trouble in certain plots
#'
#' @keywords preprocessing, write, report
#' 
#' @returns invisible(NULL)

write_sc_report <- function(final_results, output = getwd(),
                            template_folder, source_folder = NULL,
                            query = NULL, subset_by = NULL, opt,
                            cell_annotation = NULL, template = NULL,
                            out_name = NULL, use_canvas = TRUE,
                            extra_columns = NULL, DEG_p_val_cutoff = NA_real_,
                            min_avg_log2FC = NA_real_){
  if(is.null(template_folder)) {
    stop("No template folder was provided.")
  }
  if(is.null(source_folder)) {
    source_folder <- find.package("htmlreportR")
  }
  if(!file.exists(source_folder)) {
    stop(paste0("Source folder not found. Was ", source_folder, " ."))
  }
  if(any(is.null(final_results))) {
    stop("ERROR: final results object contains NULL fields. Analysis
       is not complete.")
  }
  if(is.null(template)) {
    stop("Please specify a template to render")
  }
  template <- file.path(template_folder, template)
  if(!file.exists(template)) {
    stop("Specified template does not exist in template folder.")
  }
  out_file <- file.path(output, paste0(opt$name, "_", out_name))
  tmp_folder <- file.path(output, paste0(opt$name, "_tmp_lib"))
  dir.create(tmp_folder)
  container <- list(opt = opt, seu = final_results$seu, qc = final_results$qc,
                    subset_by = subset_by, use_canvas = use_canvas,
                    DEG_list = final_results$DEG_list, query = query,
                    DEG_metrics_list = final_results$DEG_metrics_list,
                    DEG_query_list = final_results$DEG_query_list,
                    DEG_p_val_cutoff = opt$DEG_p_val_cutoff,
                    min_avg_log2FC = opt$min_avg_log2FC,
                    marker_meta = final_results$marker_meta,
                    subset_seu = final_results$subset_seu,
                    subset_DEGs = final_results$subset_DEGs,
                    sample_qc_pct = final_results$sample_qc_pct,
                    clusters_pct = final_results$clusters_pct,
                    query_exp = final_results$query_exp,
                    query_pct = final_results$query_pct,
                    query_cluster_pct = final_results$query_cluster_pct,
                    markers = final_results$markers, use_canvas = use_canvas,
                    cell_annotation = cell_annotation,
                    extra_columns = extra_columns,
                    integrate = final_results$integrate)
  plotter <- htmlreportR:::htmlReport$new(title_doc = paste0("Single-Cell ",
                            opt$name, " report"), container = container,
                            tmp_folder = tmp_folder, src = source_folder,
                            compress_obj = FALSE)
  plotter$build(template)
  plotter$write_report(out_file)
  message(paste0("Report written in ", out_file))
  return(invisible(NULL))
}

#' check_sc_input
#' check input of main SC analysis function

check_sc_input <- function(metadata, DEG_target){
  colnames(metadata) <- tolower(colnames(metadata))
  PASS <- rep(TRUE, length(DEG_target))
  names(PASS) <- DEG_target
  for(DEG_column in DEG_target) {
    uniques <- length(unique(metadata[[DEG_column]]))
    if(uniques != 2) {
      PASS[DEG_column] <- FALSE
    }
  }
  if(any(!PASS)) {
    stop(paste0("ERROR. Please check DEG_target have exactly two unique ",
      "values. DEG_column(s) \"", paste(names(PASS[PASS==FALSE]),
        collapse = "\", \""), "\" contain an invalid number."))
  }
}
