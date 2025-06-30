
#' main_annotate_sc
#' is the main seurat analysis function, which can be used for integrative or
#' non-integrative analysis.
#'
#' @inheritParams process_sketch
#' @inheritParams process_doublets
#' @inheritParams tag_qc
#' @importFrom BiocParallel SerialParam
#' @param seu A seurat object.
#' @param name Project name. Default NULL (no project name)
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
#' NULL by default. Per-sample analysis finds this vector for every sample,
#' integrative mode requires this vector.
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
#'                   DEG_p_val_cutoff = 5e-3,
#'                   min_avg_log2FC = 0.5, k_weight = 90)
#'  }
#' @returns invisible(final_results list. Contains multiple items:
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
#' of each experimental condition to be analyzed and compared separately.)

main_annotate_sc <- function(seu, minqcfeats, percentmt, query, sigfig = 2,
    resolution = 0.5, p_adj_cutoff = 5e-3, name = NULL, integrate = FALSE,
    cluster_annotation = NULL, cell_annotation = NULL, scalefactor = 10000,
    hvgs = 2000, subset_by = NULL, ndims = 10, normalmethod = "LogNormalize",
    verbose = FALSE, output = getwd(), reduce = FALSE,
    ref_label = NULL, SingleR_ref = NULL, ref_de_method = NULL, ref_n = NULL,
    BPPARAM = SerialParam(), doublet_list = NULL, k_weight = 100,
    integration_method = "Harmony", sketch = FALSE, sketch_pct = 25, 
    force_ncells = NA_integer_, sketch_method = "LeverageScore",
    doublet_path = getwd()){
    main_start <- Sys.time()
    new_opt <- check_sc_input(integrate = integrate, sketch = sketch,
                              SingleR_ref = SingleR_ref, reduce = reduce)
    annotate <- TRUE
    qc <- tag_qc(seu = seu, minqcfeats = minqcfeats, percentmt = percentmt,
                 doublet_list = doublet_list)
    if(length(unique(qc$sample)) == 1) {
      qc <- process_doublets(seu = qc, name = name, doublet_path = doublet_path,
                             assay = "RNA", nfeatures = hvgs, BPPARAM = BPPARAM,
                             includePCs = seq(1, ndims))
    }
    if(!reduce) {
      seu <- subset(qc, subset = qc == 'Pass')
    } else {
      seu <- qc
    }
    if(new_opt$integrate) {
      message("Splitting seurat object by sample.")
      seu[["RNA"]] <- split(seu[["RNA"]], f = seu$sample)
    }
    message('Normalizing data')
    norm_start <- Sys.time()
    seu <- Seurat::NormalizeData(object = seu, verbose = verbose,
      normalization.method = normalmethod, scale.factor = scalefactor)
    message("Normalization time: ", Sys.time() - norm_start)
    message('Finding variable features')
    seu <- Seurat::FindVariableFeatures(seu, nfeatures = hvgs,
                                        verbose = verbose,
                                        selection.method = "vst", assay = "RNA")
    if(new_opt$sketch) {
      seu <- process_sketch(seu = seu, sketch_method = sketch_method,
                            sketch_pct = sketch_pct, hvgs = hvgs,
                            force_ncells = force_ncells, verbose = verbose)
    }
    assay <- Seurat::DefaultAssay(seu)
    if(!is.null(SingleR_ref)) {
      message("SingleR reference provided. Annotating cells.")
      annotation <- annotate_SingleR(seu = seu, SingleR_ref = SingleR_ref,
                    BPPARAM = BPPARAM, ref_n = ref_n, ref_label = ref_label,
                    verbose = verbose, ref_de_method = ref_de_method,
                    aggr.ref = new_opt$aggr.ref, fine.tune = new_opt$fine.tune)
      annotate <- FALSE
      seu <- annotation$seu
      markers <- annotation$markers
    }
    message('Scaling data')
    scale_start <- Sys.time()
    seu <- Seurat::ScaleData(object = seu, verbose = verbose,
                             features = rownames(seu))
    message("Scaling time: ", Sys.time() - scale_start)
    message('Reducing dimensionality')
    seu <- Seurat::RunPCA(seu, assay = assay, npcs = ndims, verbose = verbose)
    reduction <- "pca"
    if(new_opt$integrate) {
      message('Integrating seurat object')
      seu <- Seurat::IntegrateLayers(object = seu, orig.reduction = "pca",
        new.reduction = integration_method, verbose = FALSE, assay = assay,
        dims = seq(1, ndims), method = paste0(integration_method,
        "Integration"), scale.layer = "scale.data", k.weight = k_weight)
      reduction <- integration_method
    }
    seu <- Seurat::FindNeighbors(object = seu, dims = seq(1, ndims),
                                 assay = assay, reduction = reduction,
                                 verbose = verbose)
    if(is.null(SingleR_ref)) {
      seu <- Seurat::FindClusters(seu, resolution = resolution,
                                  verbose = verbose)
      # Seurat starts counting clusters from 0, which is the source of many
      # headaches when working in R, which starts counting from 1. Therefore,
      # we introduce this correction. Weirdly enough, coercing it to numeric
      # already adds 1.
      seu@meta.data$seurat_clusters <- as.numeric(seu@meta.data$seurat_clusters)
      Seurat::Idents(seu) <- seu@meta.data$seurat_clusters
    } else {
      message("Annotation by clusters not active, skipping clustering.")
    }
    seu <- Seurat::RunUMAP(object = seu, dims = seq(ndims), verbose = verbose,
                           reduction = reduction, return.model = TRUE)
    seu <- SeuratObject::JoinLayers(seu)
    if(annotate) {
      annot_start <- Sys.time()
      annotation <- annotate_seurat(seu = seu, cell_annotation = cell_annotation,
        subset_by = subset_by, cluster_annotation = cluster_annotation,
        p_adj_cutoff = p_adj_cutoff, assay = assay,
        integrate = new_opt$integrate, verbose = verbose)
        message("Time to annotate: ", Sys.time() - annot_start)
      seu <- annotation$seu
      markers <- annotation$markers
    }
    if("sketch" %in% names(seu@assays)) {
      seu <- project_sketch(seu = seu, reduction = reduction, ndims = ndims)
    }
    assay <- "RNA"
    expr_metrics <- get_expression_metrics(seu = seu, sigfig = sigfig)
    message("Total processing time: ", Sys.time() - main_start)
    final_results <- list()
    final_results$qc <- qc
    final_results$seu <- seu
    final_results$sample_qc_pct <- expr_metrics$sample_qc_pct
    final_results$clusters_pct <- expr_metrics$clusters_pct
    final_results$markers <- markers
    final_results$SingleR_annotation <- annotation$SingleR_annotation
    final_results$integrate <- new_opt$integrate
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
#' @param query A string vector. List of genes to focus on DEG analysis in
#' addition to regular DEG analysis.
#' @param output_path Path where output will be written.
#' @returns A list containing DEG results, metrics and a query subset of
#' DEG results.
#' @export
#' @examples
#' data(pbmc_tiny)
#' seu <- pbmc_tiny
#' seu$sample <- "a"
#' seu$sample[8:15] <- "b"
#' seu$seurat_clusters <- c("0", "1")
#' DEG_target <- data.frame(sample = seu$sample,
#'                          treat = c(rep("Ctrl", 7), rep("Treat", 8)))
#' DEG_results <- main_sc_Hunter(DEG_target = DEG_target, seu = seu,
#'                         p_val_cutoff = 1, min_avg_log2FC = 0, min_counts = 0,
#'                         min_cell_proportion = 0)
#' print(DEG_results)

main_sc_Hunter <- function(DEG_target, seu, p_val_cutoff = 1e-3,
                           min_avg_log2FC = 0.5, min_counts = 10, query = NULL,
                           verbose = FALSE, min_cell_proportion = 0.01,
                           output_path = getwd()) {
  message('Starting DEG analysis.')
  DEG_query <- NULL
  seu <- .add_target_info(seu = seu, DEG_target = DEG_target)
  subset_target <- ifelse("cell_type" %in% colnames(seu@meta.data),
                           yes = "cell_type", no = "seurat_clusters")
  DEGs <- get_sc_markers(seu = seu, cond = "deg_group", DEG = TRUE,
                         subset_by = subset_target, verbose = verbose,
                         values = "Ctrl,Treat", min.pct = min_cell_proportion)
  message("Extracting DEG cell metrics")
  DEG_metrics <- get_fc_vs_ncells(seu = seu, DEG_list = DEGs$markers,
                  min_avg_log2FC = min_avg_log2FC, p_val_cutoff = p_val_cutoff,
                  min_counts = min_counts)
  if(!is.null(query)) {
    DEG_query <- get_fc_vs_ncells(seu = seu, DEG_list = DEGs$markers,
                                  min_counts = min_counts, query = query)
  }
  return(list(DEGs = DEGs, DEG_metrics = DEG_metrics, DEG_query = DEG_query))
}

#' analyze_sc_query
#' is a wrapper for the three query analysis steps:
#' `get_query_distribution`, `get_query_pct` by samples and `get_query_pct` by
#' samples and cell types.
#'
#' @importFrom Seurat GetAssayData
#' @inheritParams analyze_sc_query
#' @param layer Seurat object layer to analyze.
#' @returns A list with the three query analysis objects.
#' @examples
#' \dontrun{
#'    data(pbmc_tiny)
#'    pbmc_tiny$seurat_clusters <- c(rep(1, 7), rep(2, 8))
#'    main_analyze_sc_query(pbmc_tiny, c("PPBP", "CA2"), 2,
#'    sample_col = "orig.ident", layer = "counts")
#' }
#' @export

main_analyze_sc_query <- function(seu, query, sigfig = 2, layer = "counts",
                                  sample_col = "sample") {
  query_data <- analyze_sc_query(seu = seu, query = query, sigfig = sigfig,
                                 layer = layer, sample_col = sample_col)
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
            row.names = TRUE, file = file.path(opt$output, "embeddings",
            "cell_embeddings.tsv"))
    if(!is.null(final_results$markers)) {
      write.table(final_results$markers, sep = "\t", quote = FALSE,
                  row.names = TRUE, file = file.path(opt$output, "markers.tsv"))
    }
    if(!is.null(final_results$SingleR_annotation)) {
      data.table::fwrite(final_results$SingleR_annotation, quote = FALSE,
        file = file.path(opt$output, "SingleR_annotation.tsv"), sep = "\t")
    }
    message("Counts matrix saved to ", file.path(opt$output, "counts"))
    return(invisible(NULL))
}

#' write_sc_report
#' Write integration HTML report
#'
#' @param final_results Output from main_sc_Hunter function
#' @param output directory where report will be saved
#' @param out_suffix suffix to add to output file name. Useful when rendering
#' different templates with this function (as is the case in our workflow)
#' @param template_folder directory where template is located
#' @param template Template to render
#' @param source_folder htmlreportR source folder
#' @param use_canvas Parameter to select whether or not CanvasXpress plots will
#' be triggered in templates where this control parameter has been implemented.
#' Setting it to FALSE can be useful for big datasets, as CanvasXpress might
#' have trouble in certain plots
#' @param opt Processed ptions list, will be consulted in report.
#' @param params Input options list, will be included as execution parameters.
#' @param analysis Analysis type. Will only be used to build report name inside
#' html document.
#' @param files_css Path to css file to inject. By default it injects styles.css
#' found in ExpHunterSuite's templates directory.
#'
#' @keywords preprocessing, write, report
#' 
#' @returns invisible(NULL)
#' @export

write_sc_report <- function(final_results, analysis = "Single-Cell",
                output = getwd(), out_suffix = NULL, template_folder,
                files_css = file.path(template_folder, "styles.css"),
                source_folder = NULL, opt, params = NULL, template = NULL,
                use_canvas = TRUE){
    if(is.null(template_folder)) {
      stop("No template folder was provided.")
    }
    if(is.null(source_folder)) {
      source_folder <- find.package("htmlreportR")
    }
    if(!file.exists(source_folder)) {
      stop("Source folder not found. Was ", source_folder, " .")
    }
    if(any(is.null(final_results))) {
      # This check is stupid @me. If a list element is set to NULL, the element
      # will simply not exist. It would be best to check for an element
      # that should ALWAYS be there. Other specific checks should go in
      # reporting templates.
      stop("Final results object contains NULL fields. Analysis
         is not complete.")
    }
    if(is.null(template)) {
      stop("Please specify a template to render")
    }
    template <- file.path(template_folder, template)
    if(!file.exists(template)) {
      stop("Specified template does not exist in template folder.")
    }
    out_file <- file.path(output, paste0(opt$name, "_", out_suffix))
    tmp_folder <- file.path(output, paste0(opt$name, "_tmp_lib"))
    dir.create(tmp_folder)
    container <- list(params = params, seu = final_results$seu,
    qc = final_results$qc, subset_by = opt$subset_by, use_canvas = use_canvas,
    DEG_list = final_results$DEG_list, query = opt$target_genes, opt = opt,
    p_val_cutoff = opt$DEG_p_val_cutoff, min_avg_log2FC = opt$min_avg_log2FC,
    tables = final_results$tables, marker_meta = final_results$marker_meta,
    sample_qc_pct = final_results$sample_qc_pct, target = final_results$target,
    clusters_pct = final_results$clusters_pct, markers = final_results$markers,
    query_exp = final_results$query_data$query_exp,
    SingleR_annotation = final_results$SingleR_annotation,
    query_pct = final_results$query_data$query_pct,
    query_cluster_pct = final_results$query_data$query_cluster_pct,
    cell_annotation = opt$cell_annotation, extra_columns = opt$extra_columns,
    target_name = final_results$target_name,
    integrate = final_results$integrate)
    plotter <- htmlreportR::htmlReport$new(title_doc = paste0(opt$name,
                            analysis, " report"), container = container,
                            tmp_folder = tmp_folder, src = source_folder,
                            compress_obj = FALSE, files_css = files_css)
    plotter$build(template)
    plotter$write_report(out_file)
    message("Report written in ", out_file)
    return(invisible(NULL))
}

#' check_sc_input
#' checks and pre-processes input of main SC analysis functions. Outdated, needs
#' overhaul.
#'
#' @param metadata Seurat object metadata. 
#' @param DEG_target A data frame describing DEG analysis to perform.
#' @param integrate A boolean, or NULL. Whether or not to perform integration.
#' @param sketch A boolean, or NULL. Controls sketching of seurat object.
#' @param SingleR_ref SingleR_reference to use. Default NULL.
#' @param reduce A boolean, or NULL. Whether or not to reduce input for testing.
#' @returns A markers data frame with a new column for cell type assigned to
#' cluster.
#' @examples
#'  \dontrun{
#'    match_cell_types(markers_df = markers_df, p_adj_cutoff = 1e-5,
#'                     cell_annotation = markers_celltypes_df)
#'  }
#' @export

check_sc_input <- function(metadata = NULL, DEG_target = NULL, integrate = NULL,
                           sketch = NULL, SingleR_ref = NULL, reduce = NULL){
    if(reduce) {
      aggr.ref <- TRUE
      fine.tune <- FALSE
    } else {
      aggr.ref <- FALSE
      fine.tune <- TRUE
    }
    if(!is.null(DEG_target)) {
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
        not_pass <- paste(names(PASS[PASS==FALSE]), collapse = "\", \"")
        stop("Please check DEG_target has exactly two unique ",
          "values. DEG_column(s) \"", not_pass, "\" contain an invalid number.")
      }
    }
    if(!is.null(SingleR_ref)) {
      message("SingleR reference provided. Disabling integrative analysis and ",
              "sketching.")
      integrate <- FALSE
      sketch <- FALSE
    }
    return(list(integrate = integrate, sketch = sketch, aggr.ref = aggr.ref,
                fine.tune = fine.tune))
}
