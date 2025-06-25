###############################################################################
############### Plot functions
###############################################################################

#' get_plot_df
#'
#' `get_plot_df` Gets geneSets as df from enrichment object
#'
#' @param enrich_obj WGCNA results object.
#' @param showCategory Number of categories to show
#' @returns geneSets data.frame for printing in reports
#' @examples
#'  \dontrun{
#'    get_plot_df(enrich_MF)
#'  }
#' @export
get_plot_df <- function(enrich_obj, showCategory = 30) {
  extract_geneSets <- get_unexported_function("enrichplot", "extract_geneSets")
  list2df <- get_unexported_function("enrichplot", "list2df")
  update_n <- get_unexported_function("enrichplot", "update_n")
  geneSets <- extract_geneSets(enrich_obj, update_n(enrich_obj, showCategory))
  geneSets <- list2df(geneSets)
  return(geneSets)
}


#' get_clusters_count
#'
#' `get_clusters_count` count how many clusters are found in a WGCNA object
#'
#' @param results_WGCNA WGCNA results object.
#' @returns number of clusters in a WGCNA object
#' @examples
#'  \dontrun{
#'    get_clusters_count(results_WGCNA)
#'  }
#' @export
get_clusters_count <- function(results_WGCNA){
  cl_count <- length(unique(
                    results_WGCNA[['gene_cluster_info']][['Cluster_ID']]))
  return(cl_count)
}

#' get_features_count
#'
#' `get_features_count` Extracts the number of traits presents in a WGCNA
#' results object.
#'
#' @param results_WGCNA WGCNA results object.
#' @returns An integer. Number of traits in object.
#' @examples
#'  \dontrun{
#'    get_features_count(results_WGCNA)
#'  }
#' @export
get_features_count <- function(results_WGCNA){
  trait_count <- length(colnames(
          results_WGCNA[['package_objects']][['module_trait_cor']]))
  return(trait_count)
}

get_clusters <- function(compareCluster_obj){ # TODO: what's happening here?
    all_clusters <- compareCluster_obj@compareClusterResult$Cluster
    return(unique(all_clusters))
}

get_categories <- function(enrich_obj, showCategory = 30){
  if(is.null(enrich_obj)) return(NULL)
  enrich_obj_class <- class(enrich_obj)[1]
  if(enrich_obj_class == "enrichResult"){
      categories <- get_plot_df(enrich_obj, showCategory)$categoryID
    }else if(enrich_obj_class == "compareClusterResult"){
      categories <- enrich_obj@compareClusterResult$Description
    }else{
      warning("Not valid enrich object")
      return(NULL)
    }
    return(unique(categories))
}

#' get_genes
#'
#' `get_genes` Extracts the vector of unique genes from an enrichment object.
#'
#' @param enrich_obj Enrichment object.
#' @param showCategory showCategory.
#' @returns An integer. Number of traits in object.
#' @examples
#'  \dontrun{
#'    get_genes(enrich_obj, showCategory = 30)
#'  }
#' @export

get_genes <- function(enrich_obj, showCategory = 30){
  genes <- get_plot_df(enrich_obj, showCategory)$Gene
  return(unique(genes))
}


#' set_default_width
#'
#' `set_default_width` Calculates optimal default with for plot object.
#'
#' @param enrich_obj Enrichment object.
#' @param default Default width, will be modified by new calculations (and NOT
#' simply replaced). Default 12.
#' @param showCategory showCategory. Default 30.
#' @param threshold threshold. Default 30.
#' @param character_size An integer. Size of characters in plot. Default 0.04.
#' @returns Optimal width.
#' @examples
#'  \dontrun{
#'    set_default_width(enrich_obj, default = 12, showCategory = 30,
#'                      threshold = 30, character_size = 0.04)
#'  }
#' @export

set_default_width <- function(enrich_obj, 
  default = 12, 
  showCategory = 30, 
  threshold = 30, 
  character_size = 0.04){
    longer_category <- max(nchar(as.character(get_categories(enrich_obj, 
                                                           showCategory))))
    if(longer_category > threshold){
        default_width <- default + character_size * longer_category
    }else{
        default_width <- default
    }
    return(default_width)
}




#' gg_heatmap
#'
#' `gg_heatmap` Builds a ggplot2 heatmap from a data table.
#'
#' @importFrom ggplot2 ggplot aes_string geom_tile theme_minimal theme 
#' element_blank element_text scale_fill_gradient2 geom_text
#' @param data_table Input data table.
#' @param input Input string. Only used if it is equal to "matrix", in which
#' case it will convert input data to a data frame.
#' @param x_axis,y_axis Names of X and Y axes.
#' @param transpose A boolean.
#'   * `TRUE`: Data will be transposed.
#'   * `FALSE` (the default): Data will NOT be transposed.
#' @param text_plot Label to add to plot.
#' @param labs A boolean.
#'   * `TRUE` (the default): Data will be transposed.
#'   * `FALSE`: Data will NOT be transposed.
#' @param x_angle Angle by which x axis labels will be rotated.
#' @param text_size,text_colour Font size and colour for text in plot
#' @param col Color pallete to use when plotting.
#' @param na_col Color for NA values.
#' @inheritParams ggplot2::aes_
#' @inheritParams ggplot2::element_text
#' @inheritParams ggplot2::scale_fill_gradient2
#' @inheritParams ggplot2::geom_text
#' @returns A ggplot2 heatmap.
#' @examples
#'  \dontrun{
#'    gg_heatmap(matrix, input = "matrix")
#'  }
#' @export
gg_heatmap <- function(data_table, input = "", transpose = FALSE,
                y_axis= "y_axis", x_axis = "x_axis", fill = "Freq",
                text_plot = NULL, labs = TRUE, x_angle = 25, text_size = 2.5,
                text_colour = "black", col = c("#0000D5","#FFFFFF","#D50000"),
                na_col = "grey50"){
  if(transpose){
    data_table <- t(data_table)
  }
  if(input == "matrix"){
    data_table <- as.data.frame(as.table(as.matrix(data_table)))
    colnames(data_table) <- c(y_axis, x_axis, fill)
  }

  pp <- ggplot2::ggplot(data_table, ggplot2::aes_(x = as.name(x_axis), 
        y = as.name(y_axis), fill = as.name(fill))) +
  ggplot2::geom_tile(show.legend = TRUE) +
  ggplot2::theme_minimal() +
  ggplot2::theme(panel.grid.major = ggplot2::element_blank())+
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = x_angle, 
      face = "bold", hjust = 1),
    axis.text.y = ggplot2::element_text(face = "bold")) +
  ggplot2::scale_fill_gradient2(
    low = col[1],
    high = col[3],
    mid= col[2],
    na.value = na_col,
    guide = "colourbar",
    aesthetics = "fill"
  )
  if(!labs){
    pp <- pp + ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank())
  }
  if(!is.null(text_plot)){
    pp <- pp + ggplot2::geom_text(ggplot2::aes_(label=as.name(text_plot)), 
      colour = text_colour, size = text_size) 
  }
  return(pp)
}


#' Function to generate scatter plot for each gene in a hunter table
#' and show logFC per package
#' taking into account the variablity between them
#' @param ht : hunter table dataframe
#' @param var_filter : variability threshold to show gene into this
#' graph (low variability will be removed)
#' @param title : plot title
#' @param top : plots only the top N items with more variance.
#' NULL will not filter
#' @param alpha : transparency of dots
#' @return plot ready to be rendered
#' @author Fernando Moreno Jabato <jabato(at)uma(dot)com>
#' @importFrom ggplot2 ggplot aes_string geom_point scale_shape_manual ylim 
#' xlab ggtitle theme theme_classic element_blank
#' @importFrom stats var
#' @export
ht2logFCPlot <- function(ht,
  var_filter = 0.001, 
  title = "Filtered logFC", 
  top = 50, 
  alpha = 0.5){

  logFC_columns <- grep("logFC_", colnames(ht), value=TRUE)
  logFCs_tab <- data.frame(gene_name=row.names(ht), ht[,logFC_columns])
  logFCs_tab <- logFCs_tab[! apply(logFCs_tab, 1, anyNA), ]
  logFCs_tab$var <- apply(logFCs_tab[,-1], 1, var)
  variable_genes <- logFCs_tab[logFCs_tab$var > var_filter, "gene_name"]
  logFCs_tab <- logFCs_tab[variable_genes, ]
  logFCs_tab <- logFCs_tab[order(logFCs_tab$var, decreasing=TRUE), ]
  if(nrow(logFCs_tab) > top) logFCs_tab <- logFCs_tab[seq(top), ]  
  logFCs_tab$gene_name <- factor(logFCs_tab$gene_name, 
    levels = logFCs_tab$gene_name[order(logFCs_tab$var, decreasing=TRUE)])

  if( length(variable_genes) == 0) return(paste("No genes have variance greater than", 
                                                var_filter))

  plot_tab <- reshape2::melt(logFCs_tab, id="gene_name", measure.vars=logFC_columns)
  names(plot_tab) <- c("Gene", "Package", "logFC")
  plot_tab$Package <- gsub("logFC_", "", plot_tab$Package)

  pp <- ggplot2::ggplot(plot_tab, ggplot2::aes_string(x = "Gene", 
      y = "logFC", colour = "Package", shape = "Package")) + 
        ggplot2::geom_point(alpha = alpha) 
  pp <- pp + 
    ggplot2::xlab(paste0("Gene (",100 - round(length(variable_genes)/length(
      unique(rownames(ht))),4)*100,
    "% of genes have lower variability than the threshold)")) + 
    ggplot2::ggtitle(title) + 
    ggplot2::theme_classic() + 
    ggplot2::theme(text = ggplot2::element_text(size = 14),axis.text.x=ggplot2::element_blank(),
        axis.ticks.x=ggplot2::element_blank()) +
    ggplot2::geom_hline(yintercept = 0,linetype="dashed", color = "#636363")
  return(pp)
}



