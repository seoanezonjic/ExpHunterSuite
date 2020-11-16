#' Function to generate target
#'
#' Load target file if it exists, otherwise use the -C and -T flags. Note target takes precedence over target.
#' @param from_file 
#' @param ctrl_samples
#' @param treat_samples
#' @keywords design
#' @export
#' @examples
#' target_generation()
target_generation <- function(from_file=NULL, ctrl_samples=NULL, treat_samples=NULL){
	target <- NULL
	if(! is.null(from_file)) {
	  target <- read.table(from_file, header=TRUE, sep="\t", check.names = FALSE)
	  # Check there is a column named treat
	  if(! "treat" %in% colnames(target)) {
	    stop(cat("No column named treat in the target file.\nPlease resubmit"))
	  }
	} else {
	  index_control_cols <- unlist(strsplit(ctrl_samples, ",")) 
	  index_treatmn_cols <- unlist(strsplit(treat_samples, ","))
	  # Create the target data frame needed in the calls to the DE detection methods
	  target <- data.frame(sample=c(index_control_cols, index_treatmn_cols), 
	    treat=c(rep("Ctrl", length(index_control_cols)), rep("Treat", length(index_treatmn_cols))))
	}	
	return(target)
}	


#' Function to write expression package results
#'
#' Write to disk the results of each diff expression package
#' @param df_list 
#' @param prefix
#' @param root
#' @keywords output
#' @export
#' @examples
#' write_df_list_as_tables()

write_df_list_as_tables <- function(df_list, prefix, root) {
  invisible(
    lapply(1:length(df_list), function(i) {
      pack <- names(df_list)[i]
      folder <- file.path(root, paste0("Results_", pack))
      if(!dir.exists(folder)){dir.create(folder)}
      write.table(df_list[[pack]],
      file=file.path(folder, paste0(prefix, pack, '.txt')), quote=FALSE, col.names = NA, sep="\t")
    })
  )
}


#' Loads stored results in DEgenes Hunter expression analysis folder and returns in correct object format
#' @param path of DGH folder
#' @export
load_hunter_folder <- function(path){
	dgh_exp_results <- list()
	dgh_exp_results[["raw_filter"]] <- read.table(file.path(path,"filtered_count_data.txt"), stringsAsFactors = FALSE) 
	dgh_exp_results[["sample_groups"]] <- read.table(file.path(path,"control_treatment.txt"), stringsAsFactors = FALSE, header = TRUE)
	dgh_exp_results[["index_control_cols"]] <- dgh_exp_results[["sample_groups"]][dgh_exp_results[["sample_groups"]][1] == "C",2]
	dgh_exp_results[["index_treatmn_cols"]] <- dgh_exp_results[["sample_groups"]][dgh_exp_results[["sample_groups"]][1] == "T",2]
	dgh_exp_results[["design_vector"]] <- dgh_exp_results[["sample_groups"]][[1]]
	aux <- table(dgh_exp_results[["sample_groups"]][[1]])
	dgh_exp_results[["replicatesC"]] <- aux[["C"]]
	dgh_exp_results[["replicatesT"]] <- aux[["T"]]
	
	packages_results <- list.dirs(path,recursive=FALSE)
	packages_results <- packages_results[grepl("Results_",packages_results)]
	if(any(grepl("WGCNA$",packages_results))){
		packages_results <- packages_results[!grepl("WGCNA$",packages_results)]
	}
	
	dgh_exp_results[["DE_all_genes"]] <- read.table(file.path(path,"Common_results","hunter_results_table.txt"), stringsAsFactors = FALSE, header = TRUE)
	dgh_exp_results[["all_FDR_names"]] <- c("padj", "FDR")
	dgh_exp_results[["all_LFC_names"]] <- c("log2FoldChange", "logFC")
	dgh_exp_results[["all_pvalue_names"]] <- c("pvalue", "PValue")
	aux <- colnames(dgh_exp_results[["DE_all_genes"]])	

	dgh_exp_results[["all_data_normalized"]] <- list()
	dgh_exp_results[["all_counts_for_plotting"]] <- list()
	dgh_exp_results[["final_logFC_names"]] <- c()
	dgh_exp_results[["final_FDR_names"]] <- c()
	dgh_exp_results[["final_pvalue_names"]] <- c()
	dgh_exp_results[["DEG_pack_columns"]] <- c() 
	invisible(lapply(packages_results,function(pack_path){
		pack_res <- load_package_result(pack_path)
		pack_name <- tail(unlist(strsplit(pack_path,"_")),1)
		if(length(pack_res) > 0){
			for(n in names(pack_res)){
				dgh_exp_results[[n]][[pack_name]] <<- pack_res[[n]]
			}
			dgh_exp_results[["final_logFC_names"]] <<- c(dgh_exp_results[["final_logFC_names"]],c(unlist(lapply(dgh_exp_results[["all_LFC_names"]],function(regex){aux[grepl(paste0("^",regex,"*_",pack_name,"$"),aux)]}))))
			dgh_exp_results[["final_FDR_names"]] <<- c(dgh_exp_results[["final_FDR_names"]],c(unlist(lapply(dgh_exp_results[["all_FDR_names"]],function(regex){aux[grepl(paste0("^",regex,"*_",pack_name,"$"),aux)]}))))
			dgh_exp_results[["final_pvalue_names"]] <<- c(dgh_exp_results[["final_pvalue_names"]],c(unlist(lapply(dgh_exp_results[["all_pvalue_names"]],function(regex){aux[grepl(paste0("^",regex,"*_",pack_name,"$"),aux)]}))))
			dgh_exp_results[["DEG_pack_columns"]] <<- c(dgh_exp_results[["DEG_pack_columns"]], aux[grepl(paste0(pack_name,"_DEG"),aux)])
		}
	}))

	aux <- load_WGCNA_results(path,dgh_exp_results[["DE_all_genes"]])
	if(length(aux) > 0) dgh_exp_results[["WGCNA_all"]] <- aux
	
	# dgh_exp_results[["all_package_results"]] <- 
	# dgh_exp_results[["package_objects"]] <- 

	return(dgh_exp_results)
}

#' Loads stored results in DEgenes Hunter expression analysis package folder and returns in compact object
#' @param pack_path of expression package folder
load_package_result <- function(pack_path){
	info <- list()
	pack_name <- tail(unlist(strsplit(pack_path,"_")),1)
	if(file.exists(file.path(pack_path,paste0("Normalized_counts_",pack_name,".txt")))){
		info[["all_data_normalized"]] <- read.table(file.path(pack_path,paste0("Normalized_counts_",pack_name,".txt")), stringsAsFactors = FALSE)
	}
	if(file.exists(file.path(pack_path,paste0("allgenes_",pack_name,".txt")))){
		info[["all_counts_for_plotting"]] <- read.table(file.path(pack_path,paste0("allgenes_",pack_name,".txt")), stringsAsFactors = FALSE)
	}
	return(info)
}

#' Loads stored WGCNA results in DEgenes Hunter expression analysis package folder and returns in compact object
#' @param path of expression analysis results (Hunter folder)
#' @param main_deg_table dataframe with main DEG analysis results
load_WGCNA_results <- function(path, main_deg_table){
	info <- list()
	if(dir.exists(file.path(path,"Results_WGCNA"))){
		# Package objects
		package_objects <- list()
		package_objects[["gene_module_cor"]] <- read.table(file.path(path,"Results_WGCNA","gene_MM.txt"))
		package_objects[["gene_module_cor_p"]] <- read.table(file.path(path,"Results_WGCNA","gene_MM_p_val.txt"))
		package_objects[["gene_trait_cor"]] <- read.table(file.path(path,"Results_WGCNA","gene_trait.txt"))
		package_objects[["gene_trait_cor_p"]] <- read.table(file.path(path,"Results_WGCNA","gene_trait_p_val.txt"))
		package_objects[["module_trait_cor"]] <- read.table(file.path(path,"Results_WGCNA","module_trait.txt"))
		package_objects[["module_trait_cor_p"]] <- read.table(file.path(path,"Results_WGCNA","module_trait_p_val.txt"))
		# gene_cluster_info
		gene_cluster_info <- cbind(rownames(main_deg_table), main_deg_table[,c("Cluster_ID", "Cluster_MM", "Cluster_MM_pval")])
		colnames(gene_cluster_info) <- c("ENSEMBL_ID", "Cluster_ID", "Cluster_MM", "Cluster_MM_pval")
		gene_cluster_info <- gene_cluster_info[rownames(package_objects[["gene_module_cor"]]),]

		# plot_objects
		plot_objects <- list()

		sample_module <- read.table(file.path(path,"Results_WGCNA","eigen_values_per_samples.txt"))
		sample_trait <- read.table(file.path(path,"Results_WGCNA","sample_trait.txt"))
		plot_objects[["trait_and_module"]] <- cbind(sample_module,sample_trait)


		# Store into info
		info$package_objects <- package_objects
		info$gene_cluster_info <- gene_cluster_info
		info$plot_objects <- plot_objects
	}

	return(info)
}
