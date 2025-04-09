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


#' @importFrom stats runif
rechunk <- function(code, counter = NULL,chunk_options = "") {

   if (chunk_options != ""){
    chunk_options <- paste0(", ", chunk_options)
   } 
   if (!is.null(counter)){
    counter <- counter + 1
   } else {
    counter <- floor(stats::runif(1) * 10000)
   }
  if (!is.character(code)){
   code_deparsed <- paste0(deparse(function() {code}), collapse = '')
   code_deparsed <- paste0("(", code_deparsed, ")()")
  } else {
    code_deparsed <- code
  }
   sub_chunk <- paste0("\n```{r sub_chunk_", counter, chunk_options, 
   "}", "\n", code_deparsed, "\n```\n\n\n")

   cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = TRUE))
   return(counter)
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

#' @export
get_clusters_count <- function(results_WGCNA){
  cl_count <- length(unique(
                    results_WGCNA[['gene_cluster_info']][['Cluster_ID']]))
  return(cl_count)
}

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

#' @export
get_genes <- function(enrich_obj, showCategory = 30){
  genes <- get_plot_df(enrich_obj, showCategory)$Gene
  return(unique(genes))
}

force_max_genes_enrobj <- function(enrich_obj, maxGenes = 200, showCategory = 30) {
  plot_df <- get_plot_df(enrich_obj, showCategory)
  for(i in seq_len(nrow(plot_df))) {
    sub_df <- plot_df[1:i,]
    unique_genes <- unique(sub_df$Gene)
    if(length(unique_genes) == maxGenes) {
      cat_ids <- plot_df[1:i, "categoryID"]
      return(enrich_obj[enrich_obj$Description %in% cat_ids, asis=TRUE])
    }
  }
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
  if(enrich_obj_class %in% c("enrichResult", "gseaResult")){
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
  if(enrich_obj_class  %in% c("enrichResult", "gseaResult")){
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
#' @export
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
  tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(a.gplot))
  leg <- which(unlist(lapply(tmp$grobs, function(x) x$name) == "guide-box"))
  legend <- tmp$grobs[[leg]]
  return(legend)
}


make_top_n_expression_table <- function(count_data, n=5) {
  top_n_index <- order(rowSums(count_data), decreasing=TRUE)[1:n]
  sample_totals <- colSums(count_data)
  top_n_count <- count_data[top_n_index, ]
  top_n_perc <- apply(top_n_count, 1, function(x) { 
    round(x / sample_totals * 100, 3)
  })
  knitr::kable(t(top_n_perc))
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


#' @importFrom ggplot2 aes_string aes ggplot geom_point scale_colour_gradientn guide_colorbar guides
#' @importFrom ggrepel geom_text_repel 
plot_odds_ratio <- function(cont_tables, OR_col = "Odds_ratio", pval_col = "Pvalue", y_col = "strategy", text_col = "TP") {
  cont_tables <- as.data.frame(cont_tables)
  idx <- order(cont_tables[,OR_col], decreasing = TRUE)
  cont_tables[,y_col] <- factor(cont_tables[,y_col],
                            levels=rev(unique(cont_tables[,y_col][idx])))
 

  OR <- ggplot2::ggplot(cont_tables, ggplot2::aes_string(x = OR_col, y = y_col)) +
  ggplot2::geom_point(ggplot2::aes_string(color = pval_col, size = "TP")) + 
  ggrepel::geom_text_repel(ggplot2::aes(label = as.character(cont_tables[,text_col]))) +
  ggplot2::scale_colour_gradientn(
    colours = c("red", "white", "blue"),
    values = c(0, 0.05, 1),
    name = pval_col, 
    guide=ggplot2::guide_colorbar(reverse=TRUE))+
  ggplot2::guides(size = FALSE)
  return(OR)
}

get_pie_clusters <- function(enrich_DF, enrichplot){
  showed_categories <- enrichplot$data$name
  all_categories <- enrich_DF@compareClusterResult$Description
  plotted_clusters <- enrich_DF@compareClusterResult[all_categories %in% showed_categories,"Cluster"]
  return(unique(plotted_clusters))
}

#' @importFrom enrichplot cnetplot
clnetplot <- function(compareCluster, ...) {
  mod_enrich_obj <- compareCluster
  mod_enrich_obj@compareClusterResult$geneID <- as.character(mod_enrich_obj@compareClusterResult$Cluster)
  clnet_plot <- enrichplot::cnetplot(mod_enrich_obj, ...) 
  return(clnet_plot)
}



parse_dim_table <- function(dimension, dim_data, tag = "##"){
  pca_dim_data <- dim_data[[dimension]]

  cat(paste0(tag, " Valiable and factors inpection for ", gsub("Dim.", "dimension ", dimension), "\n"))
  quanti_data <- pca_dim_data$quanti
  cat(paste0(tag,"# **Quantitative** Variables and Dimension Relationship\n"))
  cat("This table shows the relationship between the dimension and the numerical variables. 
    It includes the Correlation coefficient as a association metric and P value as significance mentric. \n")
  quanti_table <- DT::datatable(quanti_data, filter = 'top', rownames = TRUE, extensions = c('Buttons','ColReorder'),
    options = list(paging = TRUE,
                   colReorder = TRUE,
                   dom = 'lftBip',
                   buttons = c('copy', 'csv', 'excel')))
  invisible(rechunk(quanti_table, chunk_options = "echo=FALSE, results = 'asis', warning = FALSE, message = FALSE"))

  cat(paste0(tag,"# **Qualitative** Variables and Dimension Relationship\n"))
  cat("This table shows the relationship between the dimension and the categorical variables. 
       It includes the P value as significance mentric and R2 as a association metric.
       R2 represents the proportion (in ranges from 0 to 1) of the total variability in the categorical variable that is accounted for by the PCA dimension.
       A higher R2 value indicates a better fit or a stronger association between the variable and the dimension. \n")
   quali_table <- DT::datatable(pca_dim_data$quali, filter = 'top', rownames = TRUE, extensions = c('Buttons','ColReorder'),
    options = list(paging = TRUE,
                   colReorder = TRUE,
                   dom = 'lftBip',
                   buttons = c('copy', 'csv', 'excel')))
  invisible(rechunk(quali_table, chunk_options = "echo=FALSE, results = 'asis', warning = FALSE, message = FALSE"))

  cat(paste0(tag,"# **Qualitative** Factors and Dimension Relationship\n"))
   quali_factors_table <- DT::datatable(pca_dim_data$category, filter = 'top', rownames = TRUE, extensions = c('Buttons','ColReorder'),
    options = list(paging = TRUE,
                   colReorder = TRUE,
                   dom = 'lftBip',
                   buttons = c('copy', 'csv', 'excel')))
  invisible(rechunk(quali_factors_table, chunk_options = "echo=FALSE, results = 'asis', warning = FALSE, message = FALSE"))


}