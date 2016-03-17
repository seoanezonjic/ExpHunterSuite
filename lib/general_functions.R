 #############################################
### GENERAL FUNCTIONS 
#############################################


calculate_lfc <- function(opt){
  lfc <- log2(opt$fc)
  return(lfc)
}

#creates main folder(root)
create_folder <- function(folder, arg_paths) {
	local_path = file.path(arg_paths$root, folder)
	dir.create(local_path)
	e <- globalenv() 
  e $ paths[[folder]] <- local_path
}

#creates subfolders
create_subfolders <- function(subfolders, paths){
	for (i in 1:length(subfolders)){
		local_folder = subfolders[[i]]
		create_folder(local_folder, paths)
	}
}

write_data <- function(data, path, name){
  write.table(data, file=file.path(path, name), quote=F, col.names = NA, sep="\t")
}

   
filter_gene_expression <- function(dataframe, p_val_cutoff, lfc, pval_col_name, log2_val_col_name){
  dataframe[which(dataframe[[pval_col_name]] <= p_val_cutoff & abs(dataframe[[log2_val_col_name]]) >= lfc), ]
}


get_vector_names <- function(all_data){
  for (i in c(1:length(all_data))){
    all_data[[i]] <- row.names(all_data[[i]])
  }
  return(all_data)
}


calculate_intersection <- function(all_package_results){
  x_all <- all_package_results[[1]]
  for (i in c(2:length(all_package_results))){
    x_all <- intersect(x_all, all_package_results[[i]])
  }
  return(x_all)
}


get_subset_for_fdr_df <- function(all_data, x_all, all_FDR_names){
  all_fdrs <- c()
  all_names <- c()
  vector_names <- names(all_data)
  for (i in c(1:length(all_data))){
    reduced_dataframe <- get_specific_dataframe_names(all_data[[i]], rownames(all_data[[i]]), x_all)
    package_name <- vector_names[[i]]
    fdr_column <-reduced_dataframe[[all_FDR_names[[i]]]]
    all_fdrs <- c(all_fdrs, fdr_column)
    fdr_names <- rep(package_name, length(fdr_column))
    all_names <- c(all_names, fdr_names)
    }
    intersection_FDR_data <- data.frame(package_name = all_names, fdr= all_fdrs)
  return(intersection_FDR_data)
}

get_all_fdr_df <- function(all_data, x_all, all_FDR_names){
  all_fdrs <- c()
  all_names <- c()
  vector_names <- names(all_data)
  for (i in c(1:length(all_data))){
    dataframe <- all_data[[i]]
    package_name <- vector_names[[i]]
    fdr_column <-dataframe[[all_FDR_names[[i]]]]
    all_fdrs <- c(all_fdrs, fdr_column)
    fdr_names <- rep(package_name, length(fdr_column))
    all_names <- c(all_names, fdr_names)
    }
    all_FDR_data <- data.frame(package_name = all_names, fdr= all_fdrs)
  return(all_FDR_data)
}


write_data_frames_list <- function(dataframe_list, prefix, paths){
  dataframe_names <- names(dataframe_list)
  for (i in c(1:length(dataframe_list))){
    dataframe <- dataframe_list[[i]]
    name <- dataframe_names[[i]]
    write_data(as.data.frame(dataframe), paths[[paste('Results_', name, sep='')]], paste(prefix, name, '.txt', sep=''))
  }
}

mapping_results <- function(raw_filter, x_all){
  read_mapped <- subset(raw_filter, rownames(raw_filter) %in% x_all)
  write_data(read_mapped, file.path(paths[["Common_results"]]),"Mapped_genes.txt")
}


annot_DEgenes <- function(raw_filter, all_data, annot, x_all, replicates){ #faltan los ECs y pathname
  raw_filter_reduced <- get_specific_dataframe_names(raw_filter, raw_filter[,1], x_all)
  rm_columns <- replicates*2+1 #'1' is the column with the row names
  raw_filter_reduced <- raw_filter[ -c(2:rm_columns) ]

  DESeq2_res_df <- get_specific_dataframe_names(all_data[[1]], rownames(all_data[[1]]), x_all)
  lfc_column_name <- all_LFC_names[[1]]  
  DESeq2_logFCs <- DESeq2_res_df[lfc_column_name]
  print(head(DESeq2_logFCs))
  annotated_results <- merge(raw_filter_reduced, DESeq2_logFCs, by.x=1, by.y="row.names")
  print(head(annotated_results))
  write_data(annotated_results, file.path(paths[["Common_results"]]),"Annotated_genes.txt")
}

 
separate_intersection_logFCs_by_sign <- function(all_data, raw_filter, x_all, all_LFC_names){
  reduced_dataframe <- get_specific_dataframe_names(all_data[[1]], rownames(all_data[[1]]), x_all)
  reduced_dataframe <- reduced_dataframe[order(rownames(reduced_dataframe)), , drop = FALSE]  
  lfc_column_name <- all_LFC_names[[1]]
  positive_common_lfcs <- reduced_dataframe[reduced_dataframe[[lfc_column_name]] > 0, ]
  pos_common_lfcs <- rownames(positive_common_lfcs)     
  negative_common_lfcs <- reduced_dataframe[reduced_dataframe[[lfc_column_name]] < 0, ]
  neg_common_lfcs <- rownames(negative_common_lfcs)
  return(list(pos_common_lfcs, neg_common_lfcs))
}


unite_result_names <- function(all_data){
  dataframe_names <- as.character(rownames(all_data[[1]]))
  for (i in c(2:length(all_data))){
      names <- as.character(rownames(all_data[[i]]))
      united_result_names <- union(dataframe_names, names)
  }
  return(united_result_names)
}


get_specific_dataframe_names <- function(dataframe, dataframe_names, selected_names){
  dataframe_selection<- subset(dataframe, dataframe_names %in% selected_names)
  return(dataframe_selection)
}


handling_errors <- function(a){
  normalized_counts <- NULL
  expres_diff <- NULL
  all_genes_df <- NULL
  return(list(expres_diff, normalized_counts, all_genes_df))
}


get_all_data_lfcs_results <- function(all_data, all_LFC_names){ #REVISAR!!!
  all_data_first_element <- as.data.frame(all_data[[3]]) #aqui estaba puesto un 1, lo que solo recogia los DEGs y no todos...por lo que al final los resultados de la union serian los mismos...revisar!!!
  
  first_df_positive_lfcs <- all_data_first_element[all_data_first_element[[all_LFC_names[1]]] > 0, ]
  first_positive_lfcs <- rownames(first_df_positive_lfcs)
  
  first_df_negative_lfcs <- all_data_first_element[all_data_first_element[[all_LFC_names[1]]] < 0, ]
  first_negative_lfcs <- rownames(first_df_negative_lfcs)

  for (i in c(2:length(all_data))){
    union_df <- as.data.frame(all_data[[i]])
    
    union_df_positive_lfcs <- union_df[union_df[[all_LFC_names[i]]] > 0, ]
    union_names_pos <- rownames(union_df_positive_lfcs)
    pos_union_lfcs <- union(first_positive_lfcs, union_names_pos)

    union_df_negative_lfcs <- union_df[union_df[[all_LFC_names[i]]] < 0, ]
    union_names_neg <- rownames(union_df_negative_lfcs)
    neg_union_lfcs <- union(first_negative_lfcs, union_names_neg)
  }
  return(list(pos_union_lfcs, neg_union_lfcs))
}


get_specific_condition_dataframe_names <- function(dataframe, selected_names_condition){
  met_condition_df <- subset(dataframe, selected_names_condition)
  condition_names <- as.character(met_condition_df[,1])
  unique_condition_names <- unique(condition_names)
  return(unique_condition_names)
}
