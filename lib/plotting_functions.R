#################################################################################
############### Plot functions
###############################################################################
resize <- function(g, fig_height=5, fig_width=12) {
  g_deparsed <- paste0(deparse(function() {g}), collapse = '')
  sub_chunk <- paste0("\n```{r sub_chunk_", floor(runif(1) * 10000), ", fig.height=", fig_height, ", fig.width=", fig_width, ", echo=FALSE}", "\n(", g_deparsed, ")()\n```\n\n\n")
  cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = TRUE))
}

plot_in_div <- function(g, fig_height=7, fig_width=12) {
  cat('\n<div class="plot_real_size">\n')
  if(fig_height > 50){
	fig_height <- 50
	}
  if(fig_width > 50){
	fig_wigth <- 50
	}
  g_deparsed <- paste0(deparse(function() {g}), collapse = '')
  sub_chunk <- paste0("\n```{r sub_chunk_", floor(runif(1) * 10000), ", fig.height=", fig_height, ", fig.width=", fig_width, ", echo=FALSE}", "\n(", g_deparsed, ")()\n```\n\n\n") 
  cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = TRUE))
  cat('\n</div>\n')
}

get_plot_df <- function(enrich_obj, showCategory = 30) {
  geneSets <- enrichplot:::extract_geneSets(enrich_obj, enrichplot:::update_n(enrich_obj, showCategory))
  geneSets <- enrichplot:::list2df(geneSets)
  return(geneSets)
}

get_clusters_count <- function(results_WGCNA){
  cl_count <- length(unique(results_WGCNA[['gene_cluster_info']][['Cluster_ID']]))
  return(cl_count)
}

get_features_count <- function(results_WGCNA){
  trait_count <- length(colnames(results_WGCNA[['package_objects']][['module_trait_cor']]))
  return(trait_count)
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

calc_width <- function(enrich_obj, showCategory = 30, category_character_size = 0.035, character_width = 3, legend_size = 1){
  enrich_obj_class <- class(enrich_obj)[1]
  if(enrich_obj_class == "enrichResult"){
    longer_x_element <- length(get_genes(enrich_obj, showCategory))
  }else if(enrich_obj_class == "compareClusterResult"){
      longer_x_element <- length(get_clusters(enrich_obj))
    }else{
      warning("Not valid enrich object")
  }

  width_size <- (legend_size + (category_character_size * max(nchar(as.character(get_categories(enrich_obj, showCategory))))) + (character_width * longer_x_element))

  if(width_size > 100){
    width_size <- 50
  }
  return(width_size)
}

calc_height <- function(enrich_obj, showCategory = NA, min_size = 7, gene_character_size = 3, category_name_size = 0.1){
  enrich_obj_class <- class(enrich_obj)[1]
  if(enrich_obj_class == "enrichResult"){
    plot_area_elements <- max(nchar(as.character(get_genes(enrich_obj, showCategory = showCategory))))
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
    height_size <- (gene_character_size* plot_area_elements) + (category_name_size * categories)
  if(height_size < min_size){
    height_size <- min_size
  } else if(height_size > 100){
    height_size <- 50
  }

  return(height_size)
}

set_default_width <- function(enrich_obj, default = 12, showCategory = 30, threshold = 30, character_size = 0.04){
	longer_category <- max(nchar(as.character(get_categories(enrich_obj, showCategory))))
	if(longer_category > threshold){
		default_width <- default + character_size * longer_category
	}else{
		default_width <- default
	}
	return(default_width)
}
# Prepare resize function

clusters_heatplot <- function(compare_cluster_obj){
  pp <- ggplot(compare_cluster_obj, aes(x = Cluster, y = Description, fill = p.adjust)) + 
  geom_tile() 
}
set_standard_size <- function(pp){
  pp <- pp + theme_light() + 
        theme(axis.text.y = element_text(size = 10, face = 'bold'), axis.title = element_blank()) 
  return(pp)
}







