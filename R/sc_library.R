
#' read_sc_counts
#' creates a seurat object from cellranger counts
#'
#' @importFrom Seurat Read10X CreateSeuratObject
#' @param name sample name
#' @param input path to cellranger counts
#' @param mincells min number of cells for which a feature is recorded
#' @param minfeats min number of features for which a cell is recorded
#' @param exp_design experiment design table
#' 
#' @return Seurat object
#' @examples
#'  \dontrun{
#'    read_sc_counts(name = "experiment_name", mincells = 1, minfeats = 1,
#'                   input = "path\to\counts\matrices",
#'                   exp_design = exp_design)
#'  }
#' @export

read_sc_counts <- function(name, input, mincells = 1, minfeats = 1, exp_design){
  mtx <- Seurat::Read10X(input)
  seu <- Seurat::CreateSeuratObject(counts = mtx, project = name,
                                    min.cells = mincells,
                                    min.features = minfeats)
  seu <- add_sc_design(seu = seu, name = name, exp_design = exp_design)
  return(seu)
}

#' tag_qc
#' tags cell IDs not passing QC filters
#'
#' @importFrom Seurat PercentageFeatureSet
#' @param seu Seurat object to tag
#' @param minqcfeats Min number of features for which a cell is selected.
#' Default 500
#' @param percentmt Max percentage of reads mapped to mitochondrial genes for
#' which a cell is selected. Default 5
#' @param doublet_list Vector of barcodes to be marked as doublet. Default NULL
#'
#' @keywords preprocessing, qc
#' 
#' @return Seurat object with tagged metadata
#' @examples
#' data(pbmc_tiny)
#' # Integrative experiment, previously loaded doublet list:
#' tag_qc(seu = pbmc_tiny, minqcfeats = 500, percentmt = 5,
#'        doublet_list = c("ATGCCAGAACGACT", "CATGGCCTGTGCAT"))
#' # Per-sample experiment, no prior doublet knowledge
#' tag_qc(seu = pbmc_tiny, minqcfeats = 500, percentmt = 5,
#'        doublet_list = NULL)
#' @export

tag_qc <- function(seu, minqcfeats = 500, percentmt = 5, doublet_list = NULL){
  seu@meta.data$percent.mt <- Seurat::PercentageFeatureSet(seu,
                                                        pattern = "(?i)^MT-")
  seu@meta.data$percent.rb <- Seurat::PercentageFeatureSet(seu,
                                                        pattern = "(?i)^RP[SL]")
  seu@meta.data$qc <- vector(mode = "character", length = nrow(seu@meta.data))
  seu@meta.data$qc[seu@meta.data$nFeature_RNA < minqcfeats] <- "Low_nFeature"
  high_mt <- seu@meta.data$percent.mt > percentmt
  seu@meta.data$qc[high_mt] <- paste(seu@meta.data$qc[high_mt], "High_MT",
                                     sep = ",")
  if(!is.null(doublet_list)) {
     message("Doublet list provided. Marking barcodes")
     seu <- tag_doublets(seu, doublet_list)
  }
  seu@meta.data$qc[seu@meta.data$qc == ""] = "Pass"
  commas <- grep("^,", seu@meta.data$qc)
  seu@meta.data$qc[commas] <- sub(",", "", seu@meta.data$qc[commas])
  colnames(seu@meta.data) <- tolower(colnames(seu@meta.data))
  return(seu)
}

# tag_doublets
#' tags barcodes in a seurat object that appear in a vector of doublets
#'
#' @inheritParams tag_qc
#'
#' @return Seurat object with tagged doublets
#' @examples
#' data(pbmc_tiny)
#' tag_doublets(seu = pbmc_tiny,
#'              doublet_list = c("ATGCCAGAACGACT", "CATGGCCTGTGCAT"))
#' @export

tag_doublets <- function(seu, doublet_list) {
  doublet <- rownames(seu@meta.data) %in% doublet_list
  seu@meta.data$qc[doublet] <- paste(seu@meta.data$qc[doublet], "Doublet",
                                        sep = ",")
  seu@meta.data$qc <- gsub("Pass,Doublet", "Doublet", seu@meta.data$qc)
  return(seu)
}

#' add_sc_design
#' adds experimental condition to single-sample Seurat metadata
#'
#' @param seu Seurat object
#' @param name Sample name
#' @param exp_design Data frame containing experiment design
#' 
#' @keywords preprocessing, subsetting, integration
#' 
#' @return Seurat object with the experimental conditions added as metadata
#' @examples
#' data(pbmc_tiny)
#' pbmc_tiny$sample <- "sampleA"
#' exp_design <- data.frame(sample = c("sampleA", "sampleB", "sampleC"),
#'                          new_field = c("info_A", "info_B", "info_C"))
#' add_sc_design(seu = pbmc_tiny, name = "sampleA", exp_design = exp_design)
#' @export

add_sc_design <- function(seu, name, exp_design){
  if(length(unique(seu@meta.data$orig.ident)) > 1) {
    stop(paste0("Seurat object contains more than one sample. ",
                "Please use Seurat::AddMetaData"))
  }
  exp_design <- as.list(exp_design[exp_design$sample == name,])
  seu <- Seurat::AddMetaData(object = seu, metadata = exp_design)
  return(seu)
}

##########################################################################

#' merge_seurat
#' loads single-cell count matrices and creates a merged
#' seurat object.
#'
#' @importFrom Seurat Read10X CreateSeuratObject
#' @importFrom SeuratObject JoinLayers
#' @param project_name Will appear in output project.name slot.
#' @param exp_design Experiment design table in TSV format.
#' @param count_path Directory that will be searched for count matrices.
#' @param suffix Suffix to paste with main count_path. Avoids unnecessary
#' globbing.
#' 
#' @keywords preprocessing, merging, integration
#' 
#' @return Merged Seurat object
#' @examples
#'  \dontrun{
#'    merge_seurat(project_name = "experiment_name", exp_design = exp_design,
#'                 count_path = "path/to/counts/folder", suffix = "path/suffix")
#'  }
#' @export

merge_seurat <- function(project_name, exp_design, count_path,
                         suffix=''){
  full_paths <- Sys.glob(paste(count_path, suffix, sep = "/"))
  seu.list <- sapply(exp_design$sample, function(sample) {
    sample_path <- grep(sample, full_paths, value = TRUE)
    d10x <- Seurat::Read10X(sample_path)
    seu <- Seurat::CreateSeuratObject(counts = d10x, project = sample, min.cells = 1,
                              min.features = 1)
    seu <- add_sc_design(seu = seu, name = sample, exp_design = exp_design)
    })
  merged_seu <- merge(seu.list[[1]], y = seu.list[-1],
                      add.cell.ids = exp_design$sample, project = project_name)
  merged_seu <- SeuratObject::JoinLayers(merged_seu)
  return(merged_seu)
}

#' rename_clusters
#' renames seurat clusters according to dictionary, and stores result in new
#' cell_type metadata column.
#' @importFrom Seurat RenameIdents Idents
#' @param seu Non-annotated seu with markers
#' @param new_clusters Vector of names to assign to clusters.
#'
#' @return Annotated seu object
#' @examples
#' data(pbmc_tiny)
#' pbmc_tiny$seurat_clusters <- 1
#' rename_clusters(seu = pbmc_tiny, new_clusters = "typeA")
#' @export

rename_clusters <- function(seu, new_clusters = NULL ) {
  Seurat::Idents(seu) <- seu$seurat_clusters
  names(new_clusters) <- levels(seu)
  seu <- Seurat::RenameIdents(seu, new_clusters)
  seu@meta.data$cell_type <- Seurat::Idents(seu)
  return(seu)
}

#' collapse_markers
#' takes list of marker gene data frames and collapses it
#' into a cluster-markers data frame.
#'
#' @importFrom plyr rbind.fill
#' @param markers_list A list containing marker gene data frames.
#' @return A data frame. Column `cluster` contains element names of original
#' list (seurat clusters) and column `genes` contains, for each row, the top
#' markers of that element from the original list separated by commas.
#' @examples
#'  \dontrun{
#'    collapse_markers(markers_list = list_of_marker_dfs)
#'  }
#' @export

collapse_markers <- function(markers_list) {
  df_list <- vector(mode = "list", length = length(markers_list))
  for(i in seq(length(markers_list))) {
    df_list[[i]] <- as.data.frame(markers_list[[i]])
    df_list[[i]]$gene <- rownames(df_list[[i]])
    rownames(df_list[[i]]) <- NULL
  }
  merged_df <- do.call(plyr::rbind.fill, df_list)
  merged_df <- cbind(merged_df$gene, merged_df[, colnames(merged_df) != "gene"])
  fcols <- grep("log2FC", colnames(merged_df))
  if(length(fcols) > 1) {
    merged_df$avg_log2FC <- rowMeans(merged_df[, fcols])
  } else {
    merged_df$avg_log2FC <- merged_df$fcols
  }
  colnames(merged_df)[1] <- "gene"
  return(merged_df)
}

#' match_cell_types
#' takes a cluster-marker gene data frame and a cell type
#' marker file. It then looks for matches between the two and assigns a cell
#' type to each cluster of the data frame.
#'
#' @param markers_df Data frame of markers, clusters and p-values
#' @param cell_annotation Table of cell types and their associated markers
#' @param p_adj_cutoff Minimum adjusted p-value of markers to consider.
#' @return A markers data frame with a new column for cell type assigned to
#' cluster.
#' @examples
#'  \dontrun{
#'    match_cell_types(markers_df = markers_df, p_adj_cutoff = 1e-5,
#'                     cell_annotation = markers_celltypes_df)
#'  }
#' @export

match_cell_types <- function(markers_df, cell_annotation, p_adj_cutoff = 1e-5) {
  canon_types <- unique(cell_annotation$type)
  if(any(markers_df$seurat_clusters == 0)) {
    markers_df$seurat_clusters <- markers_df$seurat_clusters + 1
  }
  clusters <- unique(markers_df$seurat_clusters)
  subset_list <- vector(mode = "list", length = length(clusters))
  pcols <- grep("p_val_adj", colnames(markers_df))
  if(any(is.na(markers_df[, pcols]))) {
    warning("WARNING: NAs detected in marker p-values. Coercing to 1.",
            immediate. = TRUE)
    markers_df[, pcols][is.na(markers_df[, pcols])] <- 1
  }
  if(length(pcols) > 1) {
    pvals <- unlist(markers_df[, pcols])
    pvals <- pvals[pvals != 0]
    min_pval <- min(pvals) # Small correction for machine-zeroes
    int_p_val <- vector(mode = "double", length = nrow(markers_df))
    p_vals_1 <- markers_df[[pcols[1]]]
    p_vals_2 <- markers_df[[pcols[2]]]
    for(i in seq(1, length(int_p_val))) {
      int_p_val[i] <- corto::fisherp(c(max(min_pval, p_vals_1[i]),
                                       max(min_pval, p_vals_2[i])))
    }
    markers_df$p_val_adj <- int_p_val
  } else {
    colnames(markers_df)[pcols] <- "p_val_adj"
  }
  fcols <- grep("log2FC", colnames(markers_df))
  if(any(is.na(markers_df[, fcols]))) {
    warning("WARNING: NAs detected in marker log2FC. Coercing to 0.",
            immediate. = TRUE)
    markers_df[, fcols][is.na(markers_df[, fcols])] <- 0
  }
  if(length(fcols) > 1) {
    markers_df$avg_log2FC <- (markers_df[[fcols[1]]] +
                              markers_df[[fcols[2]]]) / 2
  } else {
    colnames(markers_df)[fcols] <- "avg_log2FC"
  }
  for(cluster in unique(markers_df$seurat_clusters)) {
    subset <- markers_df[markers_df$seurat_clusters == cluster &
                         markers_df$p_val_adj <= p_adj_cutoff, ]
    if(nrow(subset) < 1) {
      warning(paste("WARNING: cluster", cluster, "contains no significant",
                     "markers", sep = " "), immediate. = TRUE)
      subset <- markers_df[markers_df$seurat_clusters == cluster, ][1,]
      subset$gene <- "None"
      subset$avg_log2FC <- 1
      subset$p_val_adj <- 1
      subset$seurat_clusters <- cluster
      subset$cell_type <- paste0(cluster, ". Unknown")
    } else {
      subset <- subset[order(subset$p_val_adj), ]
      max_log2FC <- max(subset$avg_log2FC)
      scores <- vector(mode = "numeric", length = length(canon_types))
      names(scores) <- canon_types
      for(type in canon_types) {
        type_markers <- cell_annotation[cell_annotation$type == type, ]$marker
        found_markers <- which(subset$gene %in% type_markers)
        scores[[type]] <- sum(subset$avg_log2FC[found_markers] / max_log2FC)
      }
      if(max(unlist(scores)) == 0 || is.na(max(unlist(scores)))) {
        cluster_match <- "Unknown"
      } else {
        cluster_match <- names(scores[which(scores == max(scores))])
        cluster_match <- paste0(cluster_match, collapse = " / ")
      }
      subset$cell_type <- paste0(subset$seurat_clusters, ". ", cluster_match)
    }
    subset_list[[as.numeric(cluster)]] <- subset
  }
  stats_table <- do.call(rbind, subset_list)
  stats_table <- stats_table[order(stats_table$seurat_clusters), ]
  columns <- colnames(stats_table)
  anno_types <- strsplit(stats_table$cell_type, "\\. ")
  anno_types <- sapply(anno_types, `[`, 2)
  types <- unique(anno_types)
  for(type in types) {
    matches <- which(anno_types == type)
    type_clusters <- stats_table$seurat_clusters[matches]
    uniques <- unique(type_clusters)
    if(length(uniques) > 1) {
      indices <- vector(mode = "double", length = length(uniques))
        for (i in seq(type_clusters)) {
          indices[i] <- which(type_clusters[i] == uniques)
        }
      dupe_cluster <- paste0(stats_table$cell_type[matches],
                             " (", letters[indices], ")")
      stats_table[matches, ]$cell_type <- dupe_cluster
    }
  }
  sum_columns <- c("gene", "p_val_adj", "avg_log2FC", "seurat_clusters", "cell_type")
  res <- list(stats_table = stats_table,
              cell_types = unique(stats_table$cell_type),
              summary = stats_table[, sum_columns])
  return(res)
}

#' get_sc_markers
#' performs differential expression analysis on OR selects
#' conserved markers from a seurat object.
#'
#' @importFrom Seurat Idents FindMarkers FindConservedMarkers
#' @param seu Seurat object to analyse.
#' @param DEG A boolean.
#'   * `TRUE`: Function will calculate differentally expressed genes.
#'   * `FALSE` (the default): Function will calculate cluster markers.
#' @param cond A string. Condition by which to perform DEG analysis, or by which
#' to group data to find conserved markers.
#' @param subset_by Metadata column by which seurat object will be subset for
#' marker calculation.
#' @param verbose A boolean. Will be passed to Seurat function calls.
#' @param assay A string. Assay whose markers will be calculated. Default "RNA",
#' as per usual workflow.
#' @param values A comma-separated string delimiting values of cond column by
#' which to perform the comparison. Must specify two conditions. Default NULL,
#' will use all values. If there are more than two, function will return an
#' error.
#' @return A list containing one marker or DEG data frame per cluster, plus
#' an additional one for global DEGs if performing differential analysis.
#' @examples
#' data(pbmc_tiny)
#' pbmc_tiny$seurat_clusters <- 1
#' pbmc_tiny$seurat_clusters[8:15] <- 2
#' get_sc_markers(seu = pbmc_tiny, cond = "groups", DEG = TRUE, verbose = TRUE,
#'                subset_by = "seurat_clusters")
#' @export

get_sc_markers <- function(seu, cond = NULL, subset_by, DEG = FALSE,
                           verbose = FALSE, assay = "RNA", values = NULL,
                           min.pct = 0.1) {
  if(is.null(values)) {
    conds <- unique(seu@meta.data[[cond]])
  } else if (values == "all") {
    conds <- unique(seu@meta.data[[cond]])
  } else {
    conds <- c(unlist(strsplit(values, ",")))
  }
  if(length(conds) != 2) {
    stop(paste0("ERROR: get_sc_markers only works with condition with two ",
                 "values. Provided condition has ", length(conds)))
  }
  marker_meta <- list(high = paste0(cond, ": ", conds[1]),
                      low = paste0(cond, ": ", conds[2]))
  sub_values <- as.character(sort(unique(seu@meta.data[[subset_by]])))
  sub_markers <- vector(mode = "list", length = length(sub_values))
  names(sub_markers) <- as.character(sub_values)
  for (i in seq(length(sub_values))) {
    message(paste0("Analysing cluster ", i, "/", length(sub_values)))
    if(DEG) {
      subset_seu <- subset_seurat(seu, subset_by, sub_values[i])
      meta <- as.character(subset_seu@meta.data[[cond]])
      ncells <- c(sum(meta==conds[1]), sum(meta==conds[2]))
      if(any(ncells < 3)) {
        warning(paste0('Cluster ', i, ' contains fewer than three cells for',
                        ' condition \'', conds[which(ncells < 3)],
                        '\'. Skipping DEG analysis\n', collapse = ""),
                immediate. = TRUE)
        markers <- data.frame(FALSE)
      } else {
        Seurat::Idents(subset_seu) <- cond  
        markers <- Seurat::FindMarkers(subset_seu, ident.1 = conds[1],
                                       ident.2 = conds[2], verbose = verbose,
                                       min.pct = min.pct)
        markers$gene <- rownames(markers)
      }
    } else {
      markers <- Seurat::FindConservedMarkers(seu, ident.1 = sub_values[i],
                                              grouping.var = cond,
                                              verbose = verbose, assay = assay)
      markers[[subset_by]] <- sub_values[i]
    }
    nums <- sapply(markers, is.numeric)
    markers[nums] <- lapply(markers[nums], signif, 2)
    sub_markers[[as.character(sub_values[i])]] <- markers
  }
  if(DEG) {
    message("Calculating global DEGs")
    Seurat::Idents(seu) <- seu@meta.data[, tolower(cond)]
    global_markers <- Seurat::FindMarkers(seu, ident.1 = conds[1],
                                ident.2 = conds[2], verbose = verbose)
    global_markers$gene <- rownames(global_markers)
    nums <- sapply(global_markers, is.numeric)
    global_markers[nums] <- lapply(global_markers[nums], signif, 2)
    sub_markers[["global"]] <- global_markers
  }
  if(all(c("cell_type", "seurat_clusters") %in% colnames(seu@meta.data))) {
    metadata <- unique(seu@meta.data[, c("seurat_clusters", "cell_type")])
    cell_type <- metadata[order(metadata$seurat_clusters), ]$cell_type
    if(length(sub_markers) == length(cell_type) + 1){
      names(sub_markers) <- c(as.character(cell_type), "global")
    } else {
      message("Cell types have not been annotated by cluster. Clusters cannot
        be renamed.")
      } 
  }
  res <- list(meta = marker_meta, markers = sub_markers)
  return(res)
}

#' calculate_markers
#'
#' is a wrapper that dynamically decides the best strategy to find marker genes
#' in a seurat dataset
#' @importFrom Seurat Idents FindAllMarkers
#' @inheritParams get_sc_markers
#' @param verbose A boolean. Will be passed to Seurat function calls.
#' @param idents Identity class to which to set seurat object before calculating
#' markers in case conserved mode cannot be triggered.
#' @param min.pct,logfc.threshold See ?Seurat::FindAllMarkers
#' @param subset_by Metadata column by which seurat object will be subset for
#' marker calculation. Default NULL.
#' @param integrate Whether or not integrative analysis is active. Does not
#' bother with subsetting if FALSE.
#' @param assay A string. Assay whose markers will be calculated. Default "RNA",
#' as per usual workflow.
#' @return Marker or DEG data frame.
#' @examples
#'  \dontrun{
#'    data(pbmc_tiny)
#'    calculate_markers(seu = pbmc_tiny, subset_by = "groups",
#'                      verbose = TRUE, idents = "letter.idents")
#'  }
#' @export

calculate_markers <- function(seu, subset_by = NULL, verbose = FALSE, idents = NULL,
                              integrate = FALSE, min.pct = 0.25,
                              logfc.threshold = 0.25, assay = "RNA") {
  test <- FALSE
  if(!is.null(subset_by)) {
    test <- length(subset_by) == 1 & integrate
    test <- test & length(unique(seu@meta.data[[subset_by]])) == 2
  }
  run_conserved <- ifelse(test = test, no = FALSE,
                          yes = !.has_exclusive_idents(seu = seu,
                                  idents = idents, cond = tolower(subset_by)))
  if(run_conserved) {
    markers <- get_sc_markers(seu = seu, cond = subset_by, DEG = FALSE,
                              subset_by = idents, verbose = verbose,
                              assay = assay)
    markers <- collapse_markers(markers$markers)
  }else{
    Seurat::Idents(seu) <- seu@meta.data[[idents]]
    markers <- Seurat::FindAllMarkers(seu, only.pos = TRUE, min.pct = min.pct,
                                      logfc.threshold = logfc.threshold,
                                      verbose = verbose, assay = assay)
    colnames(markers)[colnames(markers) == "cluster"] <- "seurat_clusters"
    rownames(markers) <- NULL
  }
  markers <- cbind(markers$gene, markers[, -grep("gene", colnames(markers))])
  colnames(markers)[1] <- "gene"
  return(markers)
}

#' analyze_sc_query
#' is a wrapper for the three query analysis steps:
#' `get_query_distribution`, `get_query_pct` by samples and `get_query_pct` by
#' samples and cell types.
#'
#' @importFrom Seurat GetAssayData
#' @inheritParams get_query_distribution
#' @return A list with the three query analysis objects.
#' @examples
#' data(pbmc_tiny)
#' pbmc_tiny$seurat_clusters <- c(rep(1, 7), rep(2, 8))
#' analyze_sc_query(pbmc_tiny, c("PPBP", "CA2"), 2, sample_col = "orig.ident")
#' @export

analyze_sc_query <- function(seu, query, sigfig, sample_col = "sample",
  layer = "data") {
  if(all(!query %in% rownames(seu))) {
    warning("None of the query genes are expressed in the dataset",
             immediate. = TRUE)
    res <- NULL
  } else {
    query_exp <- get_query_distribution(seu = seu, query = query, layer = layer,
                                     sigfig = sigfig, sample_col = sample_col)
    query_pct <- get_query_pct(seu = seu, query = query, by = sample_col,
                               sigfig = sigfig, layer = layer)
    if("cell_type" %in% colnames(seu@meta.data)) {
      get_by <- c(sample_col, "cell_type")
    } else {
      get_by <- c(sample_col, "seurat_clusters")
    }
    query_cluster_pct <- get_query_pct(seu = seu, query = query, by = get_by,
                                       sigfig = sigfig, layer = layer)
    res <- list(query_exp = query_exp, query_pct = query_pct,
                query_cluster_pct = query_cluster_pct)
  }
  return(res)
}



#' get_clusters_distribution
#' calculates the percentage of cells that make up each cluster for each
#' different sample in a seurat object. If clusters are annotated, it will
#' show cell types instead of cluster number.
#'
#' @inheritParams get_query_distribution
#' @param seu Clustered seurat object.
#' @param sigfig Significant figure cutoff
#' @return A data frame with cell type distribution in each sample.
#' @examples
#' data(pbmc_tiny)
#' pbmc_tiny$seurat_clusters <- c(rep(1, 7), rep(2, 8))
#' get_clusters_distribution(pbmc_tiny, 2, "orig.ident")
#' @export

get_clusters_distribution <- function(seu, sigfig = 3, sample_col = "sample") {
  clusters_column <- ifelse("cell_type" %in% colnames(seu@meta.data), 
                            "cell_type", "seurat_clusters")
  clusters_table <- table(seu@meta.data[, c(sample_col, clusters_column)])
  percent_table <- signif(clusters_table/rowSums(clusters_table)*100, sigfig)
  percent_table <- as.data.frame.matrix(percent_table)
  return(percent_table)
}

#' get_query_distribution
#' builds a table of query genes expression levels
#' across all samples of a Seurat object
#'
#' @importFrom stats aggregate
#' @param seu Seurat object
#' @param query Vector of query genes whose expression to analyse.
#' @param sigfig Significant figure cutoff
#' @param sample_col Name of column that specifies sample. Default "sample",
#' as per our workflow.
#' @return A data frame with expression levels for query genes in each sample.
#' @examples
#' data(pbmc_tiny)
#' get_query_distribution(pbmc_tiny, c("PPBP", "CA2"), 2, "orig.ident")
#' @export

get_query_distribution <- function(seu, query, sigfig = 3, layer = "data",
  sample_col = "sample") {
  genes <- SeuratObject::FetchData(seu, query, layer = layer)
  genes <- cbind(seu@meta.data[sample_col], genes)
  colnames(genes)[1] <- "sample"
  gene_distribution <- stats::aggregate(genes[, -1], list(genes$sample), FUN = sum)
  gene_distribution[, -1] <- signif(gene_distribution[, -1], sigfig)
  rownames(gene_distribution) <- gene_distribution[, 1]
  gene_distribution <- gene_distribution[, -1, drop = FALSE]
  # In the case where query vector is of length one, its names are dropped.
  # Therefore, we need to set them forcefully.
  colnames(gene_distribution) <- query[query %in% colnames(genes)]
  return(gene_distribution)
}

#' get_query_pct
#' gets the percentage of cells in each sample of a seurat
#' object which expresses genes specified in a list of queries.
#'
#' @param seu Seurat object to analyze.
#' @param query Vector of query genes whose expression to analyse.
#' @param sigfig Significant figure cutoff
#' @param by Factor (metadata column name) by which query will be broken down.
#' @param assay Seurat object assay to break down
#' @param layer Layer of seurat object assay to break down
#' @return A data frame with percentage of cells expressing query genes in each
#'         sample.
#' @examples
#' data(pbmc_tiny)
#' pbmc_tiny$seurat_clusters <- c(rep(1, 7), rep(2, 8))
#' get_query_pct(pbmc_tiny, c("PPBP", "CA2"), "seurat_clusters", 2)
#' @export

get_query_pct <- function(seu, query, by, sigfig = 2, assay = "RNA",
                          layer = "data") {
  if(length(by) < 1 || 2 < length(by)) {
    stop("Invalid 'by' length. Must be 1 or 2")
  }
  items <- unique(seu@meta.data[[by[1]]])
  subset_list <- vector(mode = "list", length = length(items))
  names(subset_list) <- items
  for(i in seq(length(items))) {
    message(paste("Subsetting", by[1],  paste0(i, "/", length(items)), sep = " "))
    subset <- subset_seurat(seu, by[1], as.character(items[i]))
    subset_list[[as.character(items[i])]] <- subset
  }
  if(length(by) == 2){
    pct_list <- vector(mode = "list", length = length(subset_list))
    for(element in seq(subset_list)) {
      message(paste("Subdividing", by[1], paste0(element, "/", length(items)),
              sep = " "))
      sec_items <- unique(subset_list[[element]]@meta.data[[by[2]]])
      new_subset <- vector(mode = "list", length = length(sec_items))
      names(new_subset) <- sec_items
      for(j in seq(length(sec_items))) {
        sublist_name <- as.character(sec_items[j])
        message(paste("Re-subsetting by", by[2],
                       paste0(j, "/", length(sec_items)), sep = " "))
        # Expr set to TRUE because double subset can result in only selecting
        # cells that have not been selected by sketch, which breaks seurat
        # object subsetting.
        new_subset[[sublist_name]] <- subset_seurat(column = by[2], expr = TRUE,
               seu = subset_list[[element]], value = as.character(sec_items[j]),
               layer = layer)
      }
      subset_list[[element]] <- new_subset
      pct_list[[element]] <- lapply(X = subset_list[[element]],
                                    FUN = breakdown_query, query = query,
                                    assay = assay, layer = layer)
      pct_list[[element]] <- do.call(rbind, pct_list[[element]]) * 100
      pct_list[[element]] <- signif(pct_list[[element]], sigfig)
      if("cell_type" %in% by) {
        rows <- strsplit(rownames(pct_list[[element]]), "\\.")
        rows <- unlist(lapply(rows, `[[`, 1))
        pct_list[[element]] <- pct_list[[element]][order(unlist(
                                                          as.numeric(rows))), ]
      }
    }
    names(pct_list) <- names(subset_list)
  } else {
    pct_list <- lapply(X = subset_list, FUN = breakdown_query, query = query,
                       assay = assay, layer = layer)
    pct_list <- do.call(rbind, pct_list) * 100
    pct_list <- signif(pct_list, sigfig)
  }
  return(pct_list)
}

#' get_top_genes
#' extracts the top N genes expressed in the highest percentage
#' of cells for each sample in a seurat object, and returns a vector with the
#' union of these genes.
#'
#' @importFrom Seurat GetAssayData
#' @inheritParams get_query_pct
#' @param seu Seurat object to analyze
#' @param top Top N genes to retrieve.
#' @param sample_col Seurat object metadata column containing sample
#' information. Default "sample", as per our workflow.
#'
#' @return A vector containing the union of the top N genes of each sample of
#' input Seurat object.
#' @examples
#' data(pbmc_tiny)
#' get_top_genes(seu = pbmc_tiny, top = 5, assay = "RNA", layer = "counts",
#'               sample_col = "orig.ident")
#' @export

get_top_genes <- function(seu, top = 20, assay = "RNA", layer = "counts",
                          sample_col = "sample") {
  if(top < 1) {
    stop(paste0("Invalid \"top\" argument. Must be greater than 1, was ",
                 top))
  }
  samples <- as.character(unique(seu[[sample_col, drop = TRUE]]))
  top_samples <- vector(mode = "list", length = length(samples))
  names(top_samples) <- samples
  for(sample in samples) {
    subset <- subset_seurat(seu, sample_col, sample)
    genes <- Seurat::GetAssayData(subset, assay, layer)
    expressed_genes <- vector(mode = "integer", length = nrow(genes))
    names(expressed_genes) <- rownames(genes)
    expressed_genes <- Matrix::rowSums(genes!=0) / ncol(genes)
    if(length(expressed_genes) > top) {
      expressed_genes <- sort(expressed_genes, decreasing = TRUE)[1:top]
    }
    top_samples[[sample]] <- names(expressed_genes)
  }
  top_genes <- unique(unlist(top_samples))
  return(top_genes)
}

#' get_qc_pct
#' creates a gene expression matrix of the union of the top N genes
#' expressed in every sample in a seurat object.
#' @inheritParams get_query_pct
#' @inheritParams get_top_genes
#' @return Gene expression matrix of union of top N expressed genes in all
#' samples.
#' @examples
#' data(pbmc_tiny)
#' pbmc_tiny$seurat_clusters <- c(rep(1, 7), rep(2, 8))
#' get_qc_pct(seu = pbmc_tiny, top = 5, assay = "RNA", layer = "counts",
#'            sample_col = "orig.ident", by = "seurat_clusters", sigfig = 2)
#' @export

get_qc_pct <- function(seu, top = 20, assay = "RNA", layer = "counts", by,
                       sigfig = 2, sample_col = "sample") {
  top_genes <- get_top_genes(seu = seu, top = top, assay = assay, layer = layer,
                             sample_col = sample_col)
  res <- get_query_pct(seu = seu, query = top_genes, by = by, sigfig = sigfig)
  return(res)
}

#' get_fc_vs_ncells
#' takes a seurat object and a list of differentially expressed gene data
#' frames, and converts them to two heatmaps with the same dimensions. The first
#' contains the amount of cells of each type that express each gene, and the
#' second the log2FC of said gene for every given cell type.
#'
#' @param seu Input seurat object.
#' @param DEG_list A data frame. Table of differentially expressed genes per
#' cell type, as calculated with Seurat.
#' @param min_avg_log2FC A numeric. Minimum absolute log2FC cutoff.
#' @param p_val_cutoff A numeric. Max adjusted p value cutoff.
#' @param min_counts An integer. Minimum gene counts to consider a gene as expressed.
#' @param query A character vector. List of genes to analyze. Default
#' NULL.

get_fc_vs_ncells<- function(seu, DEG_list, min_avg_log2FC = 0.2,
                             p_val_cutoff = 0.01, min_counts = 1,
                             query = NULL) {
  res <- NULL
  return_output <- TRUE
  if(length(DEG_list) < 1) {
    warning("Empty DEG_list provided.")
    return_output <- FALSE
  } else {
    if(!is.null(query)) {
      message("Starting analysis for target gene list. DEG cutoffs disabled.")
      if(all(!query %in% rownames(seu))) {
        warning("None of the target genes are expressed in seurat object.")
        return_output <- FALSE
      }
      min_avg_log2FC <- 0
      p_val_cutoff <- 1
    }
    if("cell_type" %in% colnames(seu@meta.data)) {
      meta <- seu$cell_type
    } else {
      meta <- seu$seurat_clusters
    }
    DEG_list$global <- NULL
    empty_DEGs <- unlist(lapply(DEG_list, function(x) return(isFALSE(unlist(x)))))
    DEG_list[empty_DEGs] <- NULL
    if(length(DEG_list) < 1) {
      warning("No per-identity DEG analysis present in DEG list.")
      return_output <- FALSE
    } else {
      if(is.null(query)) {
        gene_list <- .get_union_DEGenes(DEG_list = DEG_list)
      } else {
        gene_list <- query
        if(any(!query %in% rownames(seu))) {
          warning(paste0("Target gene(s) \"",
          paste(query[!query %in% rownames(seu)], collapse = "\", \""),
          "\" are not expressed in dataset."))
        }
        query <- query[query %in% rownames(seu)]
      }
      if(return_output) {
        matrix_list <- .get_matrices(seu = seu, meta = meta, DEG_list = DEG_list,
                        genes = gene_list, min_avg_log2FC = min_avg_log2FC,
                        p_val_cutoff = p_val_cutoff)
        DEG_df <- matrix_list$DEG_df
        if(any(!query %in% colnames(DEG_df))) {
          warning(paste0("Target gene(s) \"",
                  paste(query[!query %in% colnames(DEG_df)],
                  collapse = "\", \""), "\" are not differentially expressed in ",
                  "dataset."))
        }
        if(any(dim(DEG_df)) < 1) {
          warning("None of specified target genes are differentially expressed.")
          return_output <- FALSE
        }
      }
    }
  }
  if(return_output) {
    res <- .add_ncell_df(DEG_df = DEG_df, matrices = matrix_list$matrices,
                         min_counts = min_counts)
  }
  return(res)
}

.get_matrices <- function(seu = stop("No seurat object provided"),
                          meta = stop("No identities vector provided"),
                          genes = stop("No gene list provided."),
                          DEG_list = stop("No DEG_list provided"),
                          min_avg_log2FC = 0.2, p_val_cutoff = 0.01) {
  idents <- as.character(unique(meta))
  matrices <- vector(mode = "list", length = length(idents))
  names(matrices) <- idents
  DEG_matrices <- matrices
  for(ident in idents) {
    ident_seu <- seu[, meta == ident]
    matrix <- SeuratObject::GetAssayData(ident_seu, layer = "data")
    matrix <- matrix[rownames(matrix) %in% genes, , drop = FALSE]
    matrices[[ident]] <- matrix
    DEG_matrices[[ident]] <- .process_DEG_matrix(DEG_matrix = DEG_list[[ident]],
                   min_avg_log2FC = min_avg_log2FC, p_val_cutoff = p_val_cutoff,
                   genes_list = genes, ident = ident)
  }
  DEG_df <- do.call(rbind, DEG_matrices)
  ## Have to remove zeroes again, not in process_DEG_matrix, because we want to
  ## remove genes set to zero in ALL idents, not ident by ident
  DEG_df <- DEG_df[, colSums(abs(DEG_df)) > 0, drop = FALSE]
  DEG_df <- DEG_df[rowSums(abs(DEG_df)) > 0, , drop = FALSE ]
  return(list(DEG_df = DEG_df, matrices = matrices))
}

.get_union_DEGenes <- function(DEG_list) {
  res <- unique(do.call(c, lapply(DEG_list,
                function(matrix) return(rownames(matrix)))))
  return(res)
}

.add_ncell_df <- function(matrices = stop("No counts matrices supplied"),
                          DEG_df = stop("No DEG data frame supplied"),
                          min_counts = 1) {
  ncell_df <- .process_matrix_list(matrices, .is_expressed_matrix,
              min_counts = min_counts)
  # Different genes or cell types/clusters might be discarded due to DEG
  # thresholds or count thresholds. We must update both data frames accordingly.
  # This could be a function instead of this repetitive code.
  col_intersection <- intersect(colnames(ncell_df), colnames(DEG_df))
  row_intersection <- intersect(rownames(ncell_df), rownames(DEG_df))
  DEG_df <- DEG_df[row_intersection, col_intersection, drop = FALSE]
  ncell_df <- ncell_df[row_intersection, col_intersection, drop = FALSE]
  return(list(DEG_df = DEG_df, ncell_df = ncell_df))
}
    

.process_DEG_matrix <- function(DEG_matrix, min_avg_log2FC = 0.2, ident = NULL,
                                p_val_cutoff = 0.01, genes_list = NULL) {
  if(is.null(DEG_matrix) | FALSE %in% DEG_matrix) {
    res <- NULL
  } else {
    DEG_matrix <- DEG_matrix[rownames(DEG_matrix) %in% genes_list, ,
                             drop = FALSE]
    DEG_matrix <- DEG_matrix[, c("avg_log2FC", "p_val_adj", "gene")]
    DEG_matrix$ident <- ident
    rownames(DEG_matrix) <- NULL
    DEG_matrix <- DEG_matrix[abs(DEG_matrix$avg_log2FC) >= min_avg_log2FC, ]
    DEG_matrix <- DEG_matrix[DEG_matrix$p_val_adj <= p_val_cutoff, ]
    if(!is.null(genes_list)) {
      missing <- which(!genes_list %in% DEG_matrix$gene)
      if(length(missing) > 0) {
      missing_df <- data.frame(matrix(data = 0, ncol = ncol(DEG_matrix),
                                      nrow = length(missing)))
      colnames(missing_df) <- colnames(DEG_matrix)
      missing_df$gene <- genes_list[missing]
      missing_df$ident <- ident
      DEG_matrix <- rbind(DEG_matrix, missing_df)
      }
    }
    DEG_matrix$p_val_adj <- NULL
    DEG_matrix <- DEG_matrix[order(DEG_matrix$gene), ]
    res <- data.frame(matrix(ncol = nrow(DEG_matrix), nrow = 1))
    colnames(res) <- DEG_matrix$gene
    rownames(res) <- ident
    res[1, ] <- DEG_matrix$avg_log2FC
  }
  return(res)
}

.is_expressed_matrix <- function(matrix, min_counts = 1) {
  matrix@x[matrix@x < min_counts] <- 0
  matrix@x[matrix@x > 0] <- 1
  return(matrix)
}

.process_matrix_list <- function(matrix_list, processing_function = NULL,
                                 min_counts = 1) {
  if(!is.null(processing_function)) {
    matrix_list <- lapply(matrix_list, processing_function, min_counts = min_counts)
  }
  res <- lapply(matrix_list, Matrix::rowSums)
  res <- do.call(rbind, res)
  res <- res[, Matrix::colSums(res) > 0, drop = FALSE]
  return(as.data.frame(res))
}

#' breakdown_query
#' breaks down the expression of a list of query genes by
#' the specified parameter.
#' @inheritParams get_query_pct
#' @return A data frame of the proportion of cells (between 0 and 1) that
#' express each query gene.
#' @examples
#' data(pbmc_tiny)
#' breakdown_query(seu = pbmc_tiny, query = c("PPBP", "CA2"))
#' @export

breakdown_query <- function(input, query, assay = "RNA", layer = "data") {
  # This pseudo-method dispatch comes from an error when subsetting a seurat
  # object only by cells NOT selected by sketch. Once it is solved, this will
  # always assume a seurat object.
  if("Seurat" %in% class(input)) {
    if(nrow(input) > 1) {
      genes <- SeuratObject::GetAssayData(input, assay = assay, layer = layer)
    } else {
      genes <- input[[assay]][[, layer]]
    }
  } else {
    genes <- input
  }
  if(is.null(names(genes))) {
    features <- rownames(genes)
  } else {
    features <- names(genes)
  }
  missing <- !query %in% rownames(genes)
  if(any(missing)) {
    warning(paste0("Query gene(s) ", paste0(query[missing], collapse = ", "),
                   " not present in seurat object."), immediate. = TRUE)
    query <- query[!missing]
  }
  if(!is.vector(genes)) {
    queries <- genes[query, , drop = FALSE]
    pct <- Matrix::rowSums(!(queries == 0))/ncol(queries)
    pct[is.na(pct)] <- 0
  } else {
    queries <- genes[query]
    pct <- as.numeric(!(queries == 0))
    names(pct) <- query
  }
  return(pct)   
}

#' .has_exclusive_idents
#' checks whether any condition-identity pairs in
#' seurat object has less than three occurrences, which makes certain analyses
#' impossible. Not exported
#'
#' @param seu Seurat object
#' @param cond Condition to check
#' @param idents Identity class to which to set seurat object before calculating
#' markers in case conserved mode cannot be triggered.
#'
#' @return A boolean. `TRUE` if it contains exclusive pairs, `FALSE` otherwise.

.has_exclusive_idents <- function(seu, cond, idents) {
  meta <- seu@meta.data[, c(cond, idents)]
  meta <- meta[complete.cases(meta), ]
  groups <- unique(meta[[cond]])
  clusters <- unique(meta[[idents]])
  pairs <- expand.grid(groups, clusters)
  sum_matches <- vector(mode = "integer", length = nrow(pairs))
  for(pair in seq(nrow(pairs))) {
    matches <- apply(meta, 1, function(x) x == pairs[pair, ])
    sum_matches[pair] <- sum(colSums(matches) == 2)
  }
  if(any(sum_matches < 3)) {
    mismatch <- pairs[which(sum_matches < 3), ]
    warning('One or more identities contain fewer than three cells for one or ',
            'more categories. Affected pair(s): ',
            paste(apply(mismatch, 1, paste, collapse = "-"), collapse = ", "),
            ". \nDefaulting to general marker analysis.")
    res <- TRUE
  }else{
    res <- FALSE
  }
  return(res)
}

#' subset_seurat
#' subsets a seurat object by a specified value of provided
#' column.
#'
#' @importFrom Seurat FetchData GetAssayData
#' @param seu Seurat object.
#' @param column Column with value by which to subset input.
#' @param value Value to search within column.
#' @param expr Directly subset expression matrix instead of seurat object.
#' @return A subset of the seurat object, which itself is a seurat object or an
#' expression matrix if expr is set to TRUE.
#' @examples
#' data(pbmc_tiny)
#' head(subset_seurat(seu = pbmc_tiny, column = "groups", value = "g2",
#' expr = FALSE))
#' @export

subset_seurat <- function(seu, column, value, expr = FALSE, layer = "data") {
  # Argument "expr" comes from an error when subsetting a seurat
  # object only by cells NOT selected by sketch. Once it is solved, this will
  # always return a seurat object.
  expr_vec <- Seurat::FetchData(seu, vars = column)
  if(expr) {
    expr_vec <- expr_vec[expr_vec == value, , drop = FALSE]
    cells <- rownames(expr_vec)
    expr_data <- Seurat::GetAssayData(seu, layer = layer)
    subset <- expr_data[, colnames(expr_data) %in% cells, drop = FALSE]
  } else {
    subset <- seu[, which(expr_vec == value)]  
  }
  return(subset)
}

#' downsample_seurat
#' `downsample_seurat` takes a seurat object as input, and downsamples it to
#' specified number of cells and features. You can also input specific lists
#' if you know which cells and/or features you want to retrieve.
#'
#' @inheritParams get_query_pct
#' @param cells,features Lists or integers. An integer will trigger random
#' downsampling by its respective variable.
#' @param keep A vector of genes to keep when downsampling features.
#' @return A downsampled seurat object.
#' @examples
#' data(pbmc_tiny)
#' print(pbmc_tiny)
#' print(downsample_seurat(seu = pbmc_tiny, cells = 2, features = 5))
#' @export

downsample_seurat <- function(seu, cells = 500, features = 5000, keep = "",
                              assay = "RNA", layer = "counts") {
  if(length(unique(seu$orig.ident)) > 1) {
    seu <- SeuratObject::JoinLayers(seu)
  }
  counts <- SeuratObject::GetAssayData(seu, assay = assay, layer = layer)
  if(is.numeric(features)) {
    gene_list <- sample(rownames(counts), size = features, replace = F)
  } else {
    gene_list <- features
  }
  if(paste0(keep, collapse = "") != "") {
    gene_list <- unique(c(gene_list), keep)
  }
  counts <- counts[gene_list, ]
  seu <- subset(seu, features = rownames(counts))
  if(is.numeric(cells)) {
    cell_list <- sample(colnames(seu), size = cells, replace = F)
  } else {
    cell_list <- cells
  }
  seu <- seu[, cell_list]
  return(seu)
}

#' read_and_format_targets
#' `read_and_format_targets` formats a marker-celltype table into a list
#'
#' @param file Path to target gene file
#' @return A list with one element per cell type
#' @examples
#'  \dontrun{
#'    read_and_format_markers("path/to/markers/table")
#'  }
#' @export

read_and_format_targets <- function(file) {
  markers_df <- read.table(file, sep = "\t", header = FALSE,
                           stringsAsFactors = FALSE)
  cell_types <- markers_df[, 1]
  markers <- strsplit(markers_df[, 2], ",")
  names(markers) <- cell_types
  return(markers)
}

#' extract_metadata
#' Extract metadata dataframe from Seurat objects
#' 
#' @param seu Seurat object / list of Seurat objects 
#' @keywords preprocessing, report, metadata
#' @return Dataframe with metadata
#' @examples 
#' data(pbmc_tiny)
#' tiny_list <- list(pbmc_tiny1 = pbmc_tiny, pbmc_tiny2 = pbmc_tiny)
#' extract_metadata(seu = tiny_list)
#' @export

extract_metadata <- function(seu){
  if (!is.list(seu)){
    seu <- seu[[]]
  } else {
    seu <- lapply(seu, "[[")
    seu <- do.call(rbind, seu)
    }
  return(seu)
}

#' find_doublets
#' is a wrapper for the recommended steps for doublet
#' calculation in package DoubletFinder
#'
#' @param seu Seurat object
#' @return A list. Item "seu" contains seurat object with tagged doublets
#' in metadata slot. Item "barcodes" contains a vector of cell barcodes
#' corresponding to barcodes.
#' @examples
#'  \dontrun{
#'    find_doublets(seu = pbmc_tiny)
#'  }
#' @export

find_doublets <- function(seu) {
  nExp <- round(ncol(seu) * 0.04)  # Expect 4% doublets BUT WHY???
  seu <- DoubletFinder::doubletFinder(seu, pN = 0.25, pK = 0.09, nExp = nExp,
                                      PCs = 1:10)
  doublet_col <- grep('DF.classifications*', colnames(seu@meta.data))
  doublets <- seu@meta.data[seu@meta.data[, doublet_col] == "Doublet", ]
  doublet_barcodes <- rownames(doublets)
  seu <- subset(seu, cells = doublet_barcodes, invert = TRUE)
  res <- list(seu = seu, barcodes = doublet_barcodes)
  return(res)
}

#' sketch_sc_experiment
#' determines optimal cell number to sketch seurat assay. If it is larger than
#' specified minimum, it proceeds with sketching, else it is skipped.
#' @param seu Seurat object to sketch.
#' @param min.ncells An integer. If estimated cell number optimal value for
#' sketching is smaller than this number, sketching will not be performed.
#' Default value of 5000, recommended by Seurat tutorials.
#' @param force.ncells An integer. Number of cells to forcibly use in sketch.
#' NA by default.
#' @param assay A string. Name of assay to sketch. Default "RNA".
#' @param method A string. Method to use for sketching. Default "LeverageScore".
#' @param sketched.assay A string. Name of assay where sketch will be saved.
#' @param cell.pct A numeric. Percentage of total cells to consider
#' representative of the experiment. Default 25, as suggested by sketching
#' tutorial.
#' @return A seurat object with a new sketched assay.
#' @examples
#'  \dontrun{
#'    sketched_experiment <- sketch_sc_experiment(seu = nonsketched_experiment,
#'           assay = "RNA", method = "LeverageScore", sketched.assay = "sketch",
#'           cell.pct = 50)
#'  }

sketch_sc_experiment <- function(seu, assay = "RNA", method = "LeverageScore",
                        sketched.assay = "sketch", cell.pct = 25,
                        force.ncells = NA_integer_) {
  perform_sketch <- FALSE
  if(!is.na(force.ncells)) {
    message("force.ncells has non-empty value. Forcing sketch.")
    perform_sketch <- TRUE
    ncells <- force.ncells
    min.ncells <- 0
  } else {
    input_cells <- ncol(seu)
    if(input_cells < 60000) {
      message(paste0("Fewer than thirty thousand cells detected in experiment.",
                     " Skipping sketch."))
    } else {
      ncells <- ceiling(mean(table(seu$sample)) * cell.pct / 100)
      if(ncells < 5000) {
        warning("WARNING: Fewer than 5000 cells selected for sketching.")
      }
      if(ncells * length(unique(seu$sample)) < 30000) {
        message("Fewer than thirty thousand cells selected. Disabling sketch.")
      } else {
        perform_sketch <- TRUE
      }
    }
  } 
  if(perform_sketch) {
    seu <- Seurat::SketchData(object = seu, ncells = ncells, method = method,
                              assay = assay, sketched.assay = sketched.assay)
    Seurat::DefaultAssay(seu) <- "sketch" 
  }
  return(seu)
}

annotate_seurat <- function(seu = stop(paste0("Please provide a seurat object",
                            " to annotate")), cell_annotation = NULL,
                            subset_by = NULL, cluster_annotation = NULL,
                            p_adj_cutoff = 1e-5, verbose = FALSE, assay = "RNA",
                            integrate = FALSE) {
  res <- NULL
  if(!is.null(cell_annotation)) {
    message("Dynamically annotating clusters.")
    res <- annotate_clusters(seu = seu, subset_by = subset_by, assay = assay,
      integrate = integrate, idents = "seurat_clusters", verbose = verbose,
      p_adj_cutoff = p_adj_cutoff)
  } else {
  if(!is.null(cluster_annotation)){
    message("Clusters annotation file provided. Renaming clusters.")
    seu <- rename_clusters(seu = seu,
                           new_clusters = cluster_annotation$name)
    idents <- "cell_type"
  } else {
    warning("No data provided for cluster annotation.", immediate. = TRUE)
    idents <- "seurat_clusters"
  }
  markers <- calculate_markers(seu = seu, verbose = verbose, assay = assay,
                integrate = integrate, idents = idents, subset_by = subset_by)
  res <- list(seu = seu, markers = markers)
  }
  return(res)
}

annotate_SingleR <- function(seu = stop(paste0("Please provide a seurat object",
                             " to annotate")), SingleR_ref = NULL, ref_n = 25,
                             BPPARAM = NULL, ref_de_method = "wilcox",
                             ref_label = ref_label, aggr.ref = FALSE,
                             fine.tune = TRUE, assay = "RNA", verbose = FALSE,
                             save_pdf = getwd(), subset_by = NULL){
    counts_matrix <- Seurat::GetAssayData(seu, assay = assay)
    SingleR_annotation <- SingleR::SingleR(test = counts_matrix,
      ref = SingleR_ref, labels = SingleR_ref[[ref_label]],
      assay.type.test = "scale.data", de.method = ref_de_method,
      de.n = ref_n, BPPARAM = BPPARAM, aggr.ref = aggr.ref,
      fine.tune = fine.tune)
    seu@meta.data$cell_type <- SingleR_annotation$labels
    pdf(file.path(save_pdf, "DeltaDistribution.pdf"), width = 20, height = 10)
    print(SingleR::plotScoreHeatmap(SingleR_annotation))
    print(SingleR::plotDeltaDistribution(SingleR_annotation))
    message("Calculating cell type markers")
    markers <- calculate_markers(seu = seu, subset_by = subset_by,
                                 integrate = integrate, verbose = verbose,
                                 idents = "cell_type")
    return(list(seu = seu, markers = markers,
            SingleR_annotation = SingleR_annotation))
}

annotate_clusters <- function(seu = stop(paste0("Please provide a seurat ",
                              "object to annotate")), subset_by = NULL,
                              assay = "RNA", integrate = FALSE, verbose = FALSE,
                              idents = "seurat_clusters", p_adj_cutoff = 1e-5){
  message("Calculating cluster markers")
      markers <- calculate_markers(seu = seu, subset_by = subset_by,
                                   integrate = integrate, verbose = verbose,
                                   idents = "seurat_clusters", assay = assay)
      message("Annotating clusters")
      annotated_clusters <- match_cell_types(markers_df = markers,
                                             cell_annotation = cell_annotation,
                                             p_adj_cutoff = p_adj_cutoff)
      markers <- annotated_clusters$summary
      seu <- rename_clusters(seu = seu,
                             new_clusters = annotated_clusters$cell_types)
      return(list(seu = seu, markers = markers))
}

generate_sc_target <- function(DEG_target = stop("Missing target string")){
  processed_target <- unlist(strsplit(DEG_target, ";"))
  processed_target <- unlist(lapply(processed_target, strsplit, split = ":"),
                       recursive = FALSE)
  processed_target <- lapply(processed_target, function(vector) {
    if(length(vector) < 2) {
      processed_target <- c(vector, "all")
    } else {
      processed_target <- vector
    }
    return(processed_target)
  })
  DEG_names <- unlist(lapply(processed_target, `[[`, 1))
  DEG_names <- strsplit(DEG_names, ">")
  DEG_columns <- unlist(lapply(DEG_names, `[[`, 2))
  DEG_names <- unlist(lapply(DEG_names, `[[`, 1))
  DEG_values <- unlist(lapply(processed_target, `[[`, 2))
  processed_target <- data.frame(column = tolower(DEG_columns),
                                 values = DEG_values)
  rownames(processed_target) <- DEG_names
  return(processed_target)
}

project_sketch <- function(seu, reduction, ndims){
  refdata <- list(seurat_clusters = "seurat_clusters", cell_type = "cell_type")
  refdata <- refdata[names(refdata) %in% colnames(seu@meta.data)]
  message("Projecting sketched data")
  seu[["sketch"]] <- split(seu[["sketch"]], f = seu$sample)
  seu <- Seurat::ProjectIntegration(object = seu, sketched.assay = "sketch",
    assay = "RNA", reduction = reduction)
  seu <- Seurat::ProjectData(object = seu, sketched.assay = "sketch",
    assay = "RNA", dims = seq(1, ndims), refdata = refdata,
    full.reduction = paste0(reduction, ".full"),
    sketched.reduction = paste0(reduction, ".full")) ## WHY 30 dims???
  seu <- Seurat::RunUMAP(seu, reduction = paste0(reduction, ".full"),
    dims = seq(1, ndims), reduction.name = "umap.full",
    reduction.key = "UMAP_full_")
  Seurat::DefaultAssay(seu) <- "RNA"
  seu <- SeuratObject::JoinLayers(seu, assay = "RNA")
  seu[["sketch"]] <- NULL
  return(seu)
}

process_doublets <- function(seu = stop("Missing seurat object"),
                             qc = stop("Missing QC object"),
                             doublet_path = getwd()){
  doublets <- find_doublets(seu)
  seu <- doublets$seu
  doublet_list <- doublets$barcodes
  message("Removing doublets from Seurat object")
  seu <- subset(seu, cells = doublet_list, invert = TRUE)
  file_conn <- file(file.path(doublet_path, paste0(unique(seu$sample),
                    "_doublet_list.txt")))
  # Sample name is included in barcode when merging, we include it here
  # so it can be recognised in merged object.
  sample_doublet_list <- paste(name, doublet_list, sep = "_")
  writeLines(sample_doublet_list, file_conn)
  close(file_conn)
  qc <- tag_doublets(seu = qc, doublet_list = doublet_list)
  return(list(seu = seu, qc = qc))
}

process_sketch <- function(seu, sketch_method, sketch_pct, force_ncells, hvgs,
                           verbose = FALSE){
  message("Sketching sample data")
  sketch_start <- Sys.time()
  seu <- sketch_sc_experiment(seu = seu, assay = "RNA",
    method = sketch_method, sketched.assay = "sketch", cell.pct = sketch_pct,
    force.ncells = force_ncells)
  sketch_end <- Sys.time()
  if(verbose) {
    message(paste0("Sketching time: ", sketch_end-sketch_start))
  }
  if("sketch" %in% names(seu@assays)) {
    Seurat::DefaultAssay(seu) <- "sketch"
    seu <- Seurat::FindVariableFeatures(seu, nfeatures = hvgs,
                   verbose = verbose, selection.method = "vst", assay = "sketch")
  } else {
    Seurat::DefaultAssay(seu) <- "RNA"
  }
  return(seu)
}