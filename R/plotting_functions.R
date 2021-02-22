###############################################################################
############### Plot functions
###############################################################################
#' @importFrom knitr knit knit_expand
#' @importFrom stats runif
resize <- function(g, fig_height=5, fig_width=12) {
  g_deparsed <- paste0(deparse(function() {g}), collapse = '')
  sub_chunk <- paste0("\n```{r sub_chunk_", floor(stats::runif(1) * 10000), 
    ", fig.height=", fig_height, ", fig.width=", fig_width, 
    ", echo=FALSE}", "\n(", g_deparsed, ")()\n```\n\n\n")
  cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = TRUE))
}

#' @importFrom ggplot2 fortify
#' @importFrom DT datatable
plot_enrResult_DT <- function(ER){
  toshow <- ggplot2::fortify(ER)[,!colnames(ER@compareClusterResult) %in% 
                                c("pvalue","qvalue")]
  toshow$Cluster <- gsub("[\n,\t,\r]{0,1}\\(.*\\)","",toshow$Cluster)
  DT::datatable(toshow, filter = 'top', rownames = FALSE, 
                      extensions = c('Buttons','ColReorder'),
                      options = list(
                        colReorder = TRUE,
                        dom = 'lftBip',
                          buttons = c('copy', 'csv', 'excel')
  ))
}

#' @importFrom knitr knit knit_expand
#' @importFrom stats runif
plot_in_div <- function(g, fig_height=7, fig_width=12, 
## height and width in inches
  cex = 1, #size multiplier
  max_size = 50, # max size for plots in inches
  min_size = 10, 
  counter = NULL) {
  cat('\n<div class="plot_real_size">\n')
  fig_height <- fig_height * cex
  fig_width <- fig_width * cex

  if(fig_height > max_size){
    fig_height <- max_size
    }else if(fig_height < min_size){
    fig_height <- min_size
  }
  if(fig_width > max_size){
    fig_width <- max_size
    }else if(fig_width < min_size){
    fig_width <- min_size  
  }
  g_deparsed <- paste0(deparse(function() {g}), collapse = '')
  # set.seed(Sys.time())
  if (!is.null(counter)){
    chunk_name <- paste0("sub_chunk_", counter)
    counter <- counter + 1
  } else {
    chunk_name <- paste0("sub_chunk_", floor(stats::runif(1) * 10000))
  }
  sub_chunk <- paste0("\n```{r ", chunk_name, ", fig.height=", 
    fig_height, ", fig.width=", fig_width, ", echo=FALSE}", 
    "\n(", g_deparsed, ")()\n```\n\n\n") 
    # sub_chunk_", floor(runif(1) * 1000000), "
  cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = TRUE))
  cat('\n</div>\n')
  return(counter)
}
 
get_plot_df <- function(enrich_obj, showCategory = 30) {
  extract_geneSets <- get_unexported_function("enrichplot", "extract_geneSets")
  list2df <- get_unexported_function("enrichplot", "list2df")
  update_n <- get_unexported_function("enrichplot", "update_n")
  geneSets <- extract_geneSets(enrich_obj, update_n(enrich_obj, showCategory))
  geneSets <- list2df(geneSets)
  return(geneSets)
}

check_categories <- function(enrichplot_obj, min_categories = 1) {
  aux <- enrichplot_obj@result
  categories_count <- length(aux$Description[aux$p.adjust <= 
                                            enrichplot_obj@pvalueCutoff])  
  check <- categories_count > min_categories 
  return(check)
}

get_clusters_count <- function(results_WGCNA){
  cl_count <- length(unique(
                    results_WGCNA[['gene_cluster_info']][['Cluster_ID']]))
  return(cl_count)
}

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

get_genes <- function(enrich_obj, showCategory = 30){
  genes <- get_plot_df(enrich_obj, showCategory)$Gene
  return(unique(genes))
}

calc_width_clusters <- function(elements, multiplier = 0.3){
 width <- elements * multiplier
 return(width)
}

calc_height_clusters <- function(elements, multiplier = 0.3){
  height <- elements * multiplier
  return(height)
}

calc_width <- function(enrich_obj, 
  showCategory = 30, 
  category_character_size = 0.035, 
  character_width = 3, 
  legend_size = 1){
  enrich_obj_class <- class(enrich_obj)[1]
  if(enrich_obj_class == "enrichResult"){
    longer_x_element <- length(get_genes(enrich_obj, showCategory))
  }else if(enrich_obj_class == "compareClusterResult"){
      longer_x_element <- length(get_clusters(enrich_obj))
    }else{
      warning("Not valid enrich object")
  }

  width_size <- (legend_size + (category_character_size * 
    max(nchar(as.character(get_categories(enrich_obj, showCategory))))) + 
    (character_width * longer_x_element))

  if(width_size > 100){
    width_size <- 50
  }
  return(width_size)
}

calc_height <- function(enrich_obj, 
  showCategory = NA, 
  min_size = 7, 
  gene_character_size = 3, 
  category_name_size = 0.1){
  enrich_obj_class <- class(enrich_obj)[1]
  if(enrich_obj_class == "enrichResult"){
    plot_area_elements <- max(nchar(as.character(get_genes(enrich_obj, 
      showCategory = showCategory))))
  }else if(enrich_obj_class == "compareClusterResult"){
      plot_area_elements <- 2
    }else{
      warning("Not valid enrich object")
  }
  if(!is.na(showCategory)){
    categories <- showCategory
  } else {
    categories <- length(get_categories(enrich_obj))

  }
    height_size <- (gene_character_size* plot_area_elements) + 
        (category_name_size * categories)
  if(height_size < min_size){
    height_size <- min_size
  } else if(height_size > 100){
    height_size <- 50
  }

  return(height_size)
}

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
# Prepare resize function

#' @importFrom ggplot2 ggplot aes_string geom_tile
clusters_heatplot <- function(compare_cluster_obj){
  pp <- ggplot2::ggplot(compare_cluster_obj, ggplot2::aes_string(x = "Cluster", 
        y = "Description", fill = "p.adjust")) + 
  ggplot2::geom_tile() 
}


#' @importFrom ggplot2 theme_light theme element_text element_blank
set_standard_size <- function(pp){
  pp <- pp + ggplot2::theme_light() + 
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, 
          face = 'bold'), axis.title = ggplot2::element_blank()) 
  return(pp)
}

#' @importFrom ggplot2 ggplot aes_string geom_tile theme_minimal theme 
#' element_blank element_text scale_fill_gradient2 geom_text
gg_heatmap <- function(data_table, 
  input = "", 
  traspose = FALSE, 
  y_axis= "y_axis",
  x_axis = "x_axis", 
  fill = "Freq",
  text_plot = NULL,
  labs = TRUE,
  x_angle = 25,
  text_size = 2.5,
  text_colour = "black",
#  dendro = NULL,
  col = c("#0000D5","#FFFFFF","#D50000"),
  na_col = "grey50"){

    if(traspose){
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

#' @importFrom ggplot2 ggplot_gtable ggplot_build
extract_legend <- function(a.gplot){ 
# https://stackoverflow.com/questions/13649473/(...)
# (...)add-a-common-legend-for-combined-ggplots  
# tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(a.gplot))
  leg <- which(unlist(lapply(tmp$grobs, function(x) x$name) == "guide-box"))
  legend <- tmp$grobs[[leg]]
  return(legend)
}


#' Function to generate scatter plot for each gene in a hunter table
#' and show logFC per package
#' taking into account the variablity between them
#' @param ht : hunter table dataframe
#' @param var_filter : variability threshold to show gene into this
#' graph (low variability will be removed)
#' @param title : plot title
#' @param y_range : y limit to be applied. Only a number must be provided.
#' NULL will not fix the axis
#' @param top : plots only the top N items with more variance.
#' NULL will not filter
#' @param alpha : transparency of dots
#' @return plot ready to be rendered
#' @author Fernando Moreno Jabato <jabato(at)uma(dot)com>
#' @importFrom ggplot2 ggplot aes_string geom_point scale_shape_manual ylim 
#' xlab ggtitle theme theme_classic element_blank
#' @importFrom stats var
ht2logFCPlot <- function(ht,
  var_filter = 0.001, 
  title = "Filtered logFC", 
  y_range = NULL, 
  top = 50, 
  alpha = 0.5){
  gene_names <- rownames(ht)
  target_cols <- which(grepl("logFC_",colnames(ht)))
  aux_pack <- gsub("logFC_","",colnames(ht))
  df_logfc <- data.frame(Gene = character(),
              Package = character(),
              logFC = numeric(),
              stringsAsFactors = FALSE)
  invisible(lapply(target_cols,function(j){
    df_logfc <<- rbind(df_logfc,data.frame(Gene = gene_names,
                        Package = rep(aux_pack[j],length(gene_names)),
                        logFC = ht[,j],
                        stringsAsFactors = FALSE))
  }))
  nas <- which(is.na(df_logfc$logFC))
  if(length(nas) > 0) df_logfc <- df_logfc[-nas,]
  # Calculate var
  gene_names <- unique(df_logfc$Gene)
  vars <- unlist(lapply(seq_along(gene_names),function(i){
    stats::var(df_logfc$logFC[which(df_logfc$Gene == gene_names[i])])
  }))
  names(vars) <- gene_names
  vars <- sort(vars, decreasing = TRUE)
  vars <- vars[vars > var_filter]
  # Check
  if(length(vars) <= 0){
    return(NULL)
  }
  aux_vars <- length(vars)
  if(!is.null(top)){
    if(top <= length(vars)) vars <- vars[seq(top)]
  }
  df_logfc$Gene <- factor(df_logfc$Gene, levels = names(vars))
  df_logfc <- df_logfc[-which(is.na(df_logfc$Gene)),]
  # Check special cases
  if(!is.null(y_range)){
    df_logfc$Type <- unlist(lapply(df_logfc$logFC,function(lgfc){
      ifelse(abs(lgfc) <= abs(y_range),"Regular","Outlier")
    }))
    invisible(lapply(which(df_logfc$Type == "Outlier"),function(i){
      if(df_logfc$logFC[i] < 0){
        df_logfc$logFC[i] <<- -abs(y_range) 
      }else{
        df_logfc$logFC[i] <<- abs(y_range)
      }
    }))
    df_logfc$Type <- factor(df_logfc$Type, levels = c("Regular","Outlier"))
  }

  # Plot
  # Check special cases
  if(!is.null(y_range)){
    pp <- ggplot2::ggplot(df_logfc, ggplot2::aes_string(x = "Gene", 
      y = "logFC", colour = "Package")) + 
        ggplot2::geom_point(alpha = alpha, 
          ggplot2::aes_string(shape = "Type")) +
        ggplot2::scale_shape_manual(values=c(16, 17)) +
        ggplot2::ylim(c(-abs(y_range),abs(y_range))) 
  }else{
    pp <- ggplot2::ggplot(df_logfc, ggplot2::aes_string(x = "Gene", 
      y = "logFC", colour = "Package")) + 
        ggplot2::geom_point(alpha = alpha) 
  }
  pp <- pp + 
    ggplot2::xlab(paste0("Gene (",100 - round(aux_vars/length(
      unique(rownames(ht))),4)*100,
    "% of genes has not significant variability)")) + 
    ggplot2::ggtitle(title) + 
    ggplot2::theme_classic() + 
    ggplot2::theme(axis.text.x=ggplot2::element_blank(),
        axis.ticks.x=ggplot2::element_blank())
  # Return
  return(pp)
}


