

#' read_sc_counts
#' Create seurat object from cellranger counts
#'
#' @param name sample name
#' @param input path to cellranger counts
#' @param mincells min number of cells for which a feature is recorded
#' @param minfeats min number of features for which a cell is recorded
#' @param exp_design experiment design table
#' 
#' @return Seurat object
read_input <- function(name, input, mincells, minfeats, exp_design){
  mtx <- Seurat::Read10X(input)
  seu <- Seurat::CreateSeuratObject(counts = mtx, project = name,
                                    min.cells = mincells,
                                    min.features = minfeats)
  seu <- add_exp_design(seu = seu, name = name, exp_design = exp_design)
  return(seu)
}

#' tag_gc
#' Perform Quality Control
#'
#' @param seu Seurat object
#' @param minqcfeats Min number of features for which a cell is selected
#' @param percentmt Max percentage of reads mapped to mitochondrial genes for which a cell is selected
#'
#' @keywords preprocessing, qc
#' 
#' @return Seurat object
tag_qc <- function(seu, minqcfeats, percentmt){
  seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
  seu[["percent.rb"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^RP[SL]")
  seu[['qc']] <- ifelse(seu@meta.data$nFeature_RNA < minqcfeats,
                        'Low_nFeature', 'Pass')
  seu[['qc']] <- ifelse(seu@meta.data$nFeature_RNA < minqcfeats &
                        seu@meta.data$qc != 'Pass' &
                        seu@meta.data$qc != 'Low_nFeature',
                        paste('Low_nFeature', seu@meta.data$qc, sep = ','),
                              seu@meta.data$qc)
  seu[['qc']] <- ifelse(seu@meta.data$percent.mt > percentmt &
                        seu@meta.data$qc == 'Pass','High_MT', seu@meta.data$qc)
  seu[['qc']] <- ifelse(seu@meta.data$nFeature_RNA < minqcfeats &
                        seu@meta.data$qc != 'Pass' &
                        seu@meta.data$qc != 'High_MT',
                        paste('High_MT', seu@meta.data$qc, sep = ','),
                              seu@meta.data$qc)
  return(seu)
}


##########################################################################


#' do_dimred
#' Perform linear (PCA) and non-linear (UMAP/tSNE) dimensionality reduction
#'
#' @param seu Seurat object
#' @param ndims Number of PC to be used for UMAP / tSNE
#' @param dimreds character vector with the dimensional reductions to perform. E.g. c("pca", "tsne", "umap")
#' @param reduction Dimensional reduction to use for UMAP /tSNE. "pca" if no integration, or "harmony" if integration
#'
#' @keywords preprocessing, dimensionality, reduction, PCA, UMAP, tSNE
#' 
#' @return Seurat object
do_dimred <- function(seu, ndims, dimreds, reduction = "pca"){
  if ("pca" %in% dimreds){
    seu <- RunPCA(seu, features = VariableFeatures(object = seu))
  }
  if ("umap" %in% dimreds){
    seu <- RunUMAP(seu, dims = 1:ndims, reduction = reduction)
  }
  if ("tsne" %in% dimreds){
    seu <- RunTSNE(seu, dims = 1:ndims, reduction = reduction)
  }
  return(seu)
}


##########################################################################


#' do_clustering
#' Perform clustering of cells
#'
#' @param seu Seurat object
#' @param ndims Number of PC to be used for clustering
#' @param resolution Granularity of the downstream clustering (higher values -> greater number of clusters)
#' @param reduction Dimensional reduction to use for clustering. "pca" if no integration, or "harmony" if integration
#' 
#' @keywords preprocessing, clustering
#' 
#' @return Seurat object
do_clustering <- function(seu, ndims, resolution, reduction){
  seu <- FindNeighbors(seu, dims = 1:ndims, reduction = reduction)
  seu <- FindClusters(seu, resolution = resolution)
  return(seu)
}

#' add_exp_design
#' Add experimental condition to Seurat metadata
#'
#' @param seu Seurat object
#' @param name Sample name
#' @param exp_design Experiment design table in CSV format
#' 
#' @keywords preprocessing, subsetting, integration
#' 
#' @return Seurat object with the experimental conditions added as metadata
add_exp_design <- function(seu, name, exp_design){
  exp_design <- as.list(exp_design[exp_design$sample == name,])
  for (i in names(exp_design)){
    seu@meta.data[[i]] <- c(rep(exp_design[[i]], nrow(seu@meta.data)))
  }
  return(seu)
}


##########################################################################

#' merge_seurat
#' `merge_seurat` loads single-cell count matrices and creates a merged
#' seurat object.
#'
#' @param project_name Will appear in output project.name slot.
#' @param exp_design Experiment design table in TSV format
#' @param count_path Directory that will be searched for count matrices.
#' @param suffix Suffix to paste with main count_path. Avoids unnecessary
#' globbing.
#' 
#' @keywords preprocessing, merging, integration
#' 
#' @return Merged Seurat object
merge_seurat <- function(project_name, exp_design, count_path,
                         suffix=''){
  full_paths <- Sys.glob(paste(count_path, suffix, sep = "/"))
  seu.list <- sapply(exp_design$sample, function(sample) {
    sample_path <- grep(sample, full_paths, value = TRUE)
    d10x <- Seurat::Read10X(sample_path)
    seu <- Seurat::CreateSeuratObject(counts = d10x, project = sample, min.cells = 1,
                              min.features = 1)
    seu <- add_exp_design(seu = seu, name = sample, exp_design = exp_design)
    seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
    seu <- subset(seu, subset = nFeature_RNA > 500 & nCount_RNA > 1000
                                & percent.mt < 20)
    })
  merged_seu <- merge(seu.list[[1]], y = seu.list[-1], add.cell.ids = samples,
                      project = project_name)
  return(merged_seu)
}

#' annotate_clusters
#' `annotate_clusters` renames seurat clusters according to dictionary
#'
#' @param seu Non-annotated seu with markers
#' @param new_clusters Vector of names to assign to clusters.
#'
#' @return Annotated seu object
#'

annotate_clusters <- function(seu, new_clusters = NULL ) {
  names(new_clusters) <- levels(seu)
  seu <- Seurat::RenameIdents(seu, new_clusters)
  seu@meta.data$cell_type <- Seurat::Idents(seu)
  return(seu)
}

#' collapse_markers
#' `collapse_markers` takes list of marker gene data frames and collapses it
#' into a cluster-markers data frame.
#'
#' @param marker_list A list containing marker gene data frames.
#' @returns A data frame. Column `cluster` contains element names of original
#' list (seurat clusters) and column `genes` contains, for each row, the top
#' markers of that element from the original list separated by commas.

collapse_markers <- function(markers_list) {
  df_list <- vector(mode = "list", length = length(markers_list))
  for(i in seq(length(markers_list))) {
    df_list[[i]] <- as.data.frame(markers_list[[i]])
    df_list[[i]]$cluster <- i - 1
    df_list[[i]]$gene <- rownames(df_list[[i]])
    rownames(df_list[[i]]) <- NULL
  }
  merged_df <- do.call(plyr::rbind.fill, df_list)
  merged_df <- cbind(merged_df$gene, merged_df[, colnames(merged_df) != "gene"])
  fcols <- grep("log2FC", colnames(merged_df))
  merged_df$avg_log2FC <- rowMeans(merged_df[, fcols])
  colnames(merged_df)[1] <- "gene"
  return(merged_df)
}

#' match_cell_types
#' `match_cell_types` takes a cluster-marker gene data frame and a cell type
#' marker file. It then looks for matches between the two and assigns a cell
#' type to each cluster of the data frame.
#'
#' @param markers_df Data frame of markers, clusters and p-values
#' @param cell_annotation Table of cell types and their associated markers
#' @param top Top markers by p-value to use in cell type assignment
#' @returns A markers data frame with a new column for cell type assigned to
#' cluster.

match_cell_types <- function(markers_df, cell_annotation, p_adj_cutoff = 1e-5) {
  canon_types <- unique(cell_annotation$type)
  clusters <- unique(markers_df$cluster)
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
  markers_df <- markers_df[markers_df$p_val_adj <= p_adj_cutoff, ]
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
  for(cluster in unique(markers_df$cluster)) {
    subset <- markers_df[markers_df$cluster == cluster, ]
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
    subset$cell_type <- paste0(subset$cluster, ". ", cluster_match)
    subset_list[[as.numeric(cluster) + 1]] <- subset
  }
  stats_table <- do.call(rbind, subset_list)
  stats_table <- stats_table[order(stats_table$cluster), ]
  columns <- colnames(stats_table)
  anno_types <- strsplit(stats_table$cell_type, "\\. ")
  anno_types <- sapply(anno_types, `[`, 2)
  types <- unique(anno_types)
  for(type in types) {
    matches <- which(anno_types == type)
    type_clusters <- stats_table$cluster[matches]
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
  sum_columns <- c("gene", "p_val_adj", "avg_log2FC", "cluster", "cell_type")
  res <- list(stats_table = stats_table,
              cell_types = unique(stats_table$cell_type),
              summary = stats_table[, sum_columns])
  return(res)
}

#' get_sc_markers
#' `get_sc_markers` performs differential expression analysis on OR selects
#' conserver markers from a seurat object.
#'
#' @param seu Seurat object to analyse.
#' @param DEG A boolean.
#'   * `TRUE`: Function will calculate differentally expressed genes.
#'   * `FALSE` (the default): Function will calculate cluster markers.
#' @param cond A string. Condition by which to perform DEG analysis, or by which
#' to group data to find conserved markers.
#' @param subset_by Metadata column by which seurat object will be subset for
#' marker calculation.
#' @param verbose A boolean. Will be passed to Seurat function calls.
#' @returns A list containing one marker or DEG data frame per cluster, plus
#' an additional one for global DEGs if performing differential analysis.

get_sc_markers <- function(seu, cond = NULL, subset_by, DEG = FALSE,
                           verbose = FALSE) {
  conds <- unique(seu@meta.data[[cond]])
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
        warning(paste0('Cluster ', i, ' contains less than three cells for',
                        ' condition \'', conds[which(ncells < 3)],
                        '\'. Skipping DEG analysis', collapse = ""),
                immediate. = TRUE)
        markers <- data.frame(FALSE)
      } else {
        Seurat::Idents(subset_seu) <- cond  
        markers <- Seurat::FindMarkers(subset_seu, ident.1 = conds[1],
                                       ident.2 = conds[2], verbose = verbose)
        markers$gene <- rownames(markers)
      }
    } else {
      markers <- Seurat::FindConservedMarkers(seu, ident.1 = sub_values[i],
                                              grouping.var = cond,
                                              verbose = verbose)
    }
    nums <- sapply(markers, is.numeric)
    markers[nums] <- lapply(markers[nums], signif, 2)
    sub_markers[[as.character(sub_values[i])]] <- markers
  }
  if(DEG) {
    message("Calculating global DEGs")
    Seurat::Idents(seu) <- cond
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
#' @inheritParams get_sc_markers
#' @param verbose A boolean. Will be passed to Seurat function calls.
#' @param idents Identity class to which to set seurat object before calculating
#' markers in case conserved mode cannot be triggered.

calculate_markers <- function(seu, int_columns, verbose = FALSE, idents = NULL,
                              DEG = FALSE, integrate = FALSE) {
  run_conserved <- ifelse(test = length(int_columns) == 1 & integrate,
                          no = FALSE,
                          yes = !has_exclusive_idents(seu = seu,
                                  idents = idents, cond = tolower(int_columns)))
  if(run_conserved) {
    markers <- get_sc_markers(seu = seu, cond = int_columns, DEG = FALSE,
                              subset_by = idents, verbose = verbose)
    markers <- collapse_markers(markers$markers)
  }else{
    Seurat::Idents(seu) <- idents
    markers <- Seurat::FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25,
                                      logfc.threshold = 0.25, verbose = verbose)
    rownames(markers) <- NULL
  }
  markers <- cbind(markers$gene, markers[, -grep("gene", colnames(markers))])
  colnames(markers)[1] <- "gene"
  return(markers)
}

#' analyze_query
#' `analyze_query` is a wrapper for the three query analysis steps:
#' `get_query_distribution`, `get_query_pct` by samples and `get_query_pct` by
#' samples and cell types.
#'
#' @inheritParams get_query_distribution
#' @returns A list with the three query analysis objects.

analyze_query <- function(seu, query, sigfig) {
  if(all(!query %in% rownames(Seurat::GetAssayData(seu)))) {
    warning("None of the query genes are expressed in the dataset",
             immediate. = TRUE)
    res <- NULL
  } else {
    query_exp <- get_query_distribution(seu = seu, query = query, sigfig = sigfig)
    query_pct <- get_query_pct(seu = seu, query = query, by = "sample",
                           sigfig = sigfig)
    if("cell_type" %in% colnames(seu@meta.data)) {
      get_by <- c("sample", "cell_type")
    } else {
      get_by <- c("sample", "seurat_clusters")
    }
    query_cluster_pct <- get_query_pct(seu = seu, query = query, by = get_by,
                                       sigfig = sigfig)
    res <- list(query_exp = query_exp, query_pct = query_pct,
                query_cluster_pct = query_cluster_pct)
  }
  return(res)
}



#' get_clusters_distribution
#' `get_clusters_distribution` calculates the percentage of cells that make up
#' each cluster for each different sample in a seurat object. If clusters are
#' annotated, it will show cell types instead of cluster number.
#'
#' @param seu Clustered seurat object.
#' @param sigfig Significant figure cutoff
#' @returns A data frame with cell type distribution in each sample.

get_clusters_distribution <- function(seu, sigfig = 3) {
  clusters_column <- ifelse("cell_type" %in% colnames(seu@meta.data), 
                            "cell_type", "seurat_clusters")
  clusters_table <- table(seu@meta.data[, c("sample", clusters_column)])
  percent_table <- signif(clusters_table/rowSums(clusters_table)*100, sigfig)
  percent_table <- as.data.frame.matrix(percent_table)
  return(percent_table)
}

#' get_query_distribution
#' `get_query_distribution` builds a table of query genes expression levels
#' across all samples of a Seurat object
#'
#' @param seu Seurat object
#' @param query Vector of query genes whose expression to analyse.
#' @param sigfig Significant figure cutoff
#' @return A data frame with expression levels for query genes in each sample.

get_query_distribution <- function(seu, query, sigfig = 3) {
  genes <- SeuratObject::FetchData(seu, query)
  genes <- cbind(seu@meta.data$sample, genes)
  colnames(genes)[1] <- "sample"
  gene_distribution <- aggregate(genes[, -1], list(genes$sample), FUN = sum)
  gene_distribution[, -1] <- signif(gene_distribution[, -1], sigfig)
  rownames(gene_distribution) <- gene_distribution[, 1]
  gene_distribution <- gene_distribution[, -1, drop = FALSE]
  # In the case where query vector is of length one, its names are dropped.
  # Therefore, we need to set them forcefully.
  colnames(gene_distribution) <- query[query %in% colnames(genes)]
  return(gene_distribution)
}

#' get_query_pct
#' `get_query_pct` gets the percentage of cells in each sample of a seurat
#' object which expresses genes specified in a list of queries.
#'
#' @inheritParams breakdown_query
#' @param query Vector of query genes whose expression to analyse.
#' @param sigfig Significant figure cutoff
#' @return A data frame with expression levels for query genes in each sample.

get_query_pct <- function(seu, query, by, sigfig = 2, assay = "RNA",
                          layer = "counts") {
  if(length(by) < 1 || 2 < length(by)) {
    stop("Invalid 'by' length. Must be 1 or 2")
  }
  items <- unique(seu@meta.data[[by[1]]])
  subset_list <- vector(mode = "list", length = length(items))
  names(subset_list) <- items
  for(i in seq(length(items))) {
    message(paste("Subsetting", by[1],  paste0(i, "/", length(items)), sep = " "))
    subset <- subset_seurat(seu, by[1], items[i])
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
        new_subset[[sublist_name]] <- subset_seurat(subset_list[[element]],
                                                    by[2],
                                                    as.character(sec_items[j]))
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
#' `get_top_genes` extracts the top N genes expressed in the highest percentage
#' of cells for each sample in a seurat object, and returns a vector with the
#' union of these genes.
#'
#' @inheritParams qc_pct
#' @return A vector containing the union of the top N genes of each sample of
#' input Seurat object.

get_top_genes <- function(seu, top = 20, assay = "RNA", layer = "counts") {
  if(top < 1) {
    stop(paste0("Invalid \"top\" argument. Must be greater than 1, was ",
                 top))
  }
  samples <- unique(seu@meta.data$sample)
  top_samples <- vector(mode = "list", length = length(samples))
  names(top_samples) <- samples
  for(sample in samples) {
    subset <- subset_seurat(seu, "sample", sample)
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
#' `get_qc_pct` creates a gene expressio matrix of the union of the top N genes
#' expressed in every sample in a seurat object.
#' @param top Top N genes to take from each sample.
#' @inheritParams breakdown_query get_query_pct

get_qc_pct <- function(seu, top = 20, assay = "RNA", layer = "counts", by,
                   sigfig = 2) {
  top_genes <- get_top_genes(seu = seu, top = top, assay = assay, layer = layer)
  res <- get_query_pct(seu = seu, query = top_genes, by = by, sigfig = sigfig)
  return(res)
}

#' breakdown_query
#' `breakdown_query` breaks down the expression of a list of query genes by
#' the specified parameter.
#'
#' @param seu Seurat object
#' @param query Vector of query genes whose expression to analyse.
#' @param sigfig Significant figure cutoff, default 2
#' @param assay Seurat assay from which to extract data. Default is "RNA",
#' the default assay.
#' @param layer Layer of Seurat object from which to extract data. Default is
#' "counts", normalised assay data.
#'
#' @returns A data frame of the proportion of cells (between 0 and 1) that
#' express each query gene.

breakdown_query <- function(seu, query, assay = "RNA", layer = "counts") {
  if(nrow(seu@meta.data) > 1) {
    genes <- SeuratObject::GetAssayData(seu, assay = assay, layer = layer)
    missing <- !(query %in% rownames(genes))
  } else {
    genes <- seu[[assay]][[, layer]]
    missing <- !(query %in% names(genes))
  }
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

#' has_exclusive_idents
#' `has_exclusive_idents` checks whether any condition-identity pairs in
#' seurat object has less than three occurrences, which makes certain analyses
#' impossible.
#'
#' @param seu Seurat object
#' @param cond Condition to check
#' @param idents Identity class to which to set seurat object before calculating
#' markers in case conserved mode cannot be triggered.
#'
#' @returns A boolean. `TRUE` if it contains exclusive pairs, `FALSE` otherwise.

has_exclusive_idents <- function(seu, cond, idents) {
  meta <- seu@meta.data[, c(cond, idents)]
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
    warning('One or more identities contain less than three cells for one or ',
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
#' `subset_seurat` subsets a seurat object by a specified value of provided
#' column.
#'
#' @param seu Seurat object
#' @param column Column with value by which to subset input
#' @param value Value to search within column
#' @returns A subset of the seurat object, which itself is a seurat object.

subset_seurat <- function(seu, column, value) {
  expr <- Seurat::FetchData(seu, vars = column)
  subset <- seu[, which(expr == value)]
  return(subset)
}

#' downsample_seurat
#' `downsample_seurat` takes a seurat object as input, and downsamples it to
#' specified number of cells and features. You can also input specific lists
#' if you know which cells and/or features you want to retrieve.
#'
#' @param seu Seurat object.
#' @param cells,features Lists or integers. An integer will trigger random
#' downsampling by its respective variable.
#' @param keep A vector of genes to keep when downsampling features.
#' @returns A downsampled seurat object.

downsample_seurat <- function(seu, cells = NULL, features = NULL, keep = "",
                              assay = "RNA", layer = "counts") {
  if(!is.null(features)) {
    seu <- SeuratObject::JoinLayers(seu)
    counts <- SeuratObject::GetAssayData(seu, assay = assay, layer = layer)
    if(is.numeric(features)) {
      gene_list <- sample(rownames(counts), size = features, replace = F)
      if(!"" %in% keep) {
        gene_list <- unique(c(gene_list), keep)
      }
    } else {
      gene_list <- features
    }
    counts <- counts[gene_list, ]
    seu <- subset(seu, features = rownames(counts))
  }
  if(!is.null(cells)) {
    if(is.numeric(cells)) {
      cell_list <- sample(colnames(seu), size = cells, replace = F)
    } else {
      cell_list <- cells
    }
    seu <- seu[, cell_list]
  }
  return(seu)
}

#' read_and_format_targets
#' `read_and_format_targets` formats a marker-celltype table into a list
#' 
#' @param file Path to target gene file
#' 
#' @return A list with one element per cell type

read_and_format_targets <- function(file) {
  markers_df <- read.table(file, sep = "\t", header = FALSE,
                           stringsAsFactors = FALSE)
  cell_types <- markers_df[, 1]
  markers <- strsplit(markers_df[, 2], ",")
  names(markers) <- cell_types
  return(markers)
}


##########################################################################


#' extract_metadata
#' Extract metadata dataframe from Seurat objects
#' 
#' @param seu Seurat object / list of Seurat objects
#' 
#' @keywords preprocessing, report, metadata
#' 
#' @return Dataframe with metadata
extract_metadata <- function(seu){
  if (!is.list(seu)){
    seu <- seu[[]]
  } else {
    seu <- lapply(seu, "[[")
    seu <- do.call(rbind, seu)
    }
return(seu)
}

##########################################################################

#' make_vln
#' Make Violin plot
#' 
#' @param seu Seurat object / list of Seurat objects
#' @param feature metadata feature to plot
#' 
#' @keywords preprocessing, report, plot, violin
#' 
#' @return nothing
make_vln <- function(seu, feature){
  seu <- extract_metadata(seu)
  seu <- seu[, c("orig.ident", feature)]
  colnames(seu)[2] <- "values"
  ggplot() + 
    geom_point(seu,
               mapping = aes(orig.ident, values, fill = orig.ident),
               shape = 21,
               colour = "white",
               size = 2,
               stroke = 0.5,
               position = "jitter",
               alpha = 0.3) +
    geom_violin(seu,
                mapping = aes(orig.ident, values, fill = orig.ident),
                width = 0.5,
                color = "black") +
    xlab(NULL) +
    ylab(feature) +
    labs(fill = NULL) +
    theme_bw()
}

##########################################################################


#' ensure_list
#' Makes sure you have list of Seurat objects (even if you have only one)
#' 
#' @param seu Seurat object / list of Seurat objects
#' 
#' @keywords preprocessing, report, list
#' 
#' @return List of Seurat objects
ensure_list <- function(seu){
  if (!is.list(seu)){
    seu <- list(seu)
  }
  return(seu)
}

##########################################################################

