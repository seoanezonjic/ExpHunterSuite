#############################################
### QC AND BENCHMARKING FUNCTIONS 
#############################################


generate_heatmap <- function(norm_data){  ## A) Hierarchical clustering routine
  h_pearson <- hclust(as.dist(1-cor(t(norm_data), method="pearson")), method="complete") #Generates row and column dendrograms
  h_spearman <- hclust(as.dist(1-cor(norm_data, method="spearman")), method="complete")
  
  mycl <- cutree(h_pearson, h=max(h_pearson$height)/1.5)
  mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9)
  mycolhc <- mycolhc[as.vector(mycl)] #creates color vector for clusters
 
  write_data(mycl, file.path(paths[["Common_results"]]),"cluster_jerarquico_DE.txt")
    
  pdf(file.path(paths[['Common_results']], "Heatmap_DE.pdf"), height=12, width=10) #Creates heatmap for entire data set where the obtained clusters are indicated in the color bar.
    heatmap.2(norm_data, Rowv=as.dendrogram(h_pearson), Colv=as.dendrogram(h_spearman), col=redgreen(75), scale="row", density.info="none", trace="none", RowSideColors=mycolhc)
  dev.off()  
}


clustering_with_zscore <- function(norm_data){
  data_z <- genescale(norm_data, axis=1, method="Z")
  cl_z <- kmeans(data_z,3) ## B) Clustering kmeans from normalized data by z-score

  write_data(data_z, file.path(paths[["Common_results"]]),"rows-zscore.txt")
  write_data(cl_z$cluster, file.path(paths[["Common_results"]]),"cluster_kmeans_DE.txt")
  
  pdf(file.path(paths[['Common_results']], "plot_groups.pdf"), height=12, width=10) # plot de los agrupamientos
    plot(data_z, col=cl_z$cluster)
    points(cl_z$centers, pch="x", cex=1, col="magenta")
  dev.off()

  pdf(file.path(paths[['Common_results']], "expression_cluster.pdf"), height=12, width=10)
    par(mfrow=c(2,3)) # 2 x 3 pictures on one plot
    for(i in 1:3) {
      matplot(t(data_z[cl_z$cluster==i,]), type = "l" , main=paste("cluster:", i), ylab="Z.score.Exp", xlab="sample")
    }
  dev.off()
}


heatmapping_and_clustering <- function(norm_data, genes){     
  norm_data <- get_specific_dataframe_names(norm_data, rownames(norm_data), genes)
  norm_data <- as.matrix(na.omit(norm_data)) # Assigns row indices and converts the data into a matrix object.

  generate_heatmap(norm_data)

  clustering_with_zscore(norm_data)    
}

calculate_percentage_DEGs_in_intersection <- function(raw, x_all){
  genes_raw <- nrow(raw)
  intersection_genes <- length(x_all)
  percentage_DEGs <- intersection_genes/genes_raw*100
  return(percentage_DEGs)
}

calculate_percentage_DEGs_in_package_results <- function(raw, all_results_in_package){
  genes_raw <- nrow(raw)
  percentage_DEGs_in_package <- all_results_in_package/genes_raw*100
  return(percentage_DEGs_in_package)
}


calculate_lfc_trend <- function(extraction_lfcs, all_results_in_package){     
    positive_lfc_index <- which(extraction_lfcs>0)
    positive_lfc_genes <- length(positive_lfc_index)
    up_regulated_percentage <- positive_lfc_genes/all_results_in_package*100
    return(up_regulated_percentage)
}

calculate_mean_logFC <- function(extraction_lfcs){
    abs_lfcs <- abs(extraction_lfcs)
    abs_FCs <- 2^abs_lfcs
    mean_FCs <- mean(abs_FCs)
    return(mean_FCs)
}


generate_report <- function(all_data, all_LFC_names, genes){
  vector_names <- names(all_data)
  print(vector_names)
  statistics_report <- NULL
  union_names <- unite_result_names(all_data)
  print(head(union_names))   
  for (i in c(1:length(all_data))){
    all_results_in_package <- nrow(all_data[[i]])
    percentage_DEGs_in_package <- calculate_percentage_DEGs_in_package_results(raw, all_results_in_package)
    percentage_intersection <- calculate_percentage_DEGs_in_intersection(raw, x_all)
    intersection_number <- length(x_all)
    union_number <- length(union_names)

    percentage_union <- calculate_percentage_DEGs_in_intersection(raw, union_names)
    extraction_lfcs <-c(all_data[[i]][[all_LFC_names[[i]]]])
    extraction_fdrs <-c(all_data[[i]][[all_FDR_names[[i]]]])
    package_name <- vector_names[[i]]
    percentage_DEGs <- calculate_percentage_DEGs_in_intersection(raw, genes)
    mean_FCs <- calculate_mean_logFC(extraction_lfcs)
    up_regulated_percentage <- calculate_lfc_trend(extraction_lfcs, all_results_in_package)
    median_fdr <- median(extraction_fdrs) 
    statistics_report <- rbind(statistics_report, data.frame(name = package_name, 
      pDEGs = percentage_DEGs_in_package, FC =mean_FCs, UPreg_DEGs = up_regulated_percentage, 
      FDR = median_fdr, p_common_DEGs = percentage_intersection, n_common_DEGs = intersection_number,
      n_union_DEGs = union_number, p_union_DEGs = percentage_union))  
  }
 
  write.table(statistics_report, file=file.path(paths$root, "statistics_report.txt"), quote=F, sep="\t", row.names = FALSE)  
}

 

