################################# GENERAL FUNCTIONS ##############################################
checking_input <- function(opt){
  if (is.null(opt$input_file)){
    stop(cat("No file with RNA-seq counts is provided. \n Please use -i to submit file"))
  }
  if (is.null(opt$Control_columns)){
    stop(cat("No control samples are indicated. \n Please select column names of the count table with -C"))
  }
  if (is.null(opt$Treatment_columns)){
    stop(cat("No treatment samples are indicated. \n Please select column names of the count table with -T"))
  }
}


defining_subfolders <- function(replicatesC, replicatesT, opt){
  subfolders <- c()
  if (sum(replicatesC, replicatesT)<=3){
    subfolders <- c('Results_DESeq2')
    } else if (sum(replicatesC, replicatesT)<=5){
       subfolders <- c(subfolders, 'Common_results')
        if (grepl("D", opt$modules)){
          subfolders <- c(subfolders, 'Results_DESeq2')
        }
        if (grepl("E", opt$modules)){
          subfolders <- c(subfolders, 'Results_edgeR')
        }
    } else if ((replicatesC > 2) & (replicatesT >2)){
        subfolders <- c(subfolders, 'Common_results')
        if (grepl("D", opt$modules)){
          subfolders <- c(subfolders, 'Results_DESeq2')
        }
        if (grepl("E", opt$modules)){
          subfolders <- c(subfolders, 'Results_edgeR')
        }
        if (grepl("L", opt$modules)){    
          subfolders <- c(subfolders, 'Results_limma')  
        }
        if (grepl("N", opt$modules)){
          subfolders <- c(subfolders, 'Results_NOISeq')
        }
    }
  return(subfolders)  
}

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

calculate_replicates_number <- function(opt){
  ccolumns <- unlist(strsplit(opt, ","))
  replicatesC <- length(ccolumns)
  return(replicatesC)
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


annot_DEgenes <- function(raw_filter, all_data, annot, x_all, replicates){
  raw_filter_reduced <- get_specific_dataframe_names(raw_filter, raw_filter[,1], x_all)
  rm_columns <- replicates*2+1 #'1' is the column with the row names
  raw_filter_reduced <- raw_filter[ -c(2:rm_columns) ]

  DESeq2_res_df <- get_specific_dataframe_names(all_data[[1]], rownames(all_data[[1]]), x_all)
  lfc_column_name <- all_LFC_names[[1]]  
  DESeq2_logFCs <- DESeq2_res_df[lfc_column_name]
  annotated_results <- merge(raw_filter_reduced, DESeq2_logFCs, by.x=1, by.y="row.names")
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


get_all_data_lfcs_results <- function(all_data, all_LFC_names){ 
  all_data_first_element <- as.data.frame(all_data[[3]])
  
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




####################### preparing complete_genes_statistics #################################

unite_all_list_dataframes <- function(all_counts_for_plotting, all_FDR_names, all_LFC_names, all_pvalue_names, final_pvalue_names, final_logFC_names, final_FDR_names){
    all_genes_df <- as.data.frame(all_counts_for_plotting[[1]])
    all_genes_df <- all_genes_df[c(all_LFC_names[1], all_FDR_names[1], all_pvalue_names[1])]
    colnames(all_genes_df)[colnames(all_genes_df)==all_LFC_names[1]] <- final_logFC_names[1]
    colnames(all_genes_df)[colnames(all_genes_df)==all_FDR_names[1]] <- final_FDR_names[1]
    colnames(all_genes_df)[colnames(all_genes_df)==all_pvalue_names[1]] <- final_pvalue_names[1]

    for (i in c(2:length(all_counts_for_plotting))){  
         next_all_genes_df <- as.data.frame(all_counts_for_plotting[[i]])
         next_all_genes_df <- next_all_genes_df[c(all_LFC_names[i], all_FDR_names[i], all_pvalue_names[i])]
         colnames(next_all_genes_df)[colnames(next_all_genes_df)==all_LFC_names[i]] <- final_logFC_names[i] 
         colnames(next_all_genes_df)[colnames(next_all_genes_df)==all_FDR_names[i]] <- final_FDR_names[i]
         colnames(next_all_genes_df)[colnames(next_all_genes_df)==all_pvalue_names[i]] <- final_pvalue_names[i]
         all_genes_df <- merge(all_genes_df, next_all_genes_df, by="row.names") #o bin by=0, all=TRUE
         rownames_all <- all_genes_df[,1]
         all_genes_df[1] <- NULL
         rownames(all_genes_df) <- rownames_all 
    }
    return(all_genes_df) 
}



unite_all_rownames_from_dataframes_list <- function(all_data){
   alldata_df_names <- rownames(all_data[[1]])
    for (i in c(2:length(all_data))){
      next_alldata_df_names <- rownames(all_data[[i]])
      alldata_df_names <- c(alldata_df_names, next_alldata_df_names)
      complete_alldata_df <- unique(alldata_df_names)
    }
    return(complete_alldata_df) 
}



check_deg_in_pck <- function(all_counts_for_plotting, all_data, all_genes_df, DEG_pack_columns){
  DESeq2_DEG <- rownames(all_counts_for_plotting[[1]]) %in% rownames(all_data[[1]])  
  names(DESeq2_DEG) <- rownames(all_counts_for_plotting[[1]])
  DESeq2_DEG <- as.data.frame(DESeq2_DEG)
  all_genes_df <- merge(all_genes_df, DESeq2_DEG, by="row.names")
  for (i in c(2:length(all_counts_for_plotting))){
    pck_DEG <- rownames(all_counts_for_plotting[[i]]) %in% rownames(all_data[[i]])
    names(pck_DEG) <- rownames(all_counts_for_plotting[[i]])
    pck_DEG <- as.data.frame(pck_DEG)
    colnames(pck_DEG) <- DEG_pack_columns[i]
    all_genes_df <- merge(all_genes_df, pck_DEG, by.x="Row.names", by.y="row.names")
  }
  return(all_genes_df)
}


plotting_FDR_values <- function(FDRs, graphname, yaxis){
  pdf(file.path(paths[["Common_results"]], graphname), w=11, h=8.5)
    p_seguros_Int <- ggplot(FDRs, aes(x = package_name, y = fdr, color = package_name))
    plot(p_seguros_Int + geom_boxplot(outlier.colour = rgb(0, 0, 0, 0)) + theme_bw(base_size = 30) + geom_point(position = position_jitter(w = 0.1), color = "grey50", size = 1) + geom_hline(aes(yintercept = opt$p_val_cutoff)) + ylab("1 - precision (FDR)") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("") + scale_colour_discrete(guide = "none") + coord_cartesian(ylim = c(0, yaxis)))
  dev.off()    
}


count_coincidences <- function(is_a_DEG, coincid_counter){  
  sum_DEGs <- sum(is_a_DEG)
  if (sum_DEGs <= 1){
    coincid_counter[i] <- length(is_a_DEG)-sum_DEGs
  } else {
    coincid_counter[i] <- sum_DEGs
  }
  return(coincid_counter)
}


adding_filtered_transcripts <- function(raw, raw_filter, final_BIG_table){
  filtered_rows <- subset(raw, !( rownames(raw) %in% rownames(raw_filter) ))
  df_for_filtered_rows <- as.data.frame(matrix(NA, ncol = ncol(final_BIG_table), nrow = nrow(filtered_rows)))
  cols_final_BIG_table <- colnames(final_BIG_table)
  colnames(df_for_filtered_rows) <- cols_final_BIG_table
  df_for_filtered_rows["Row.names"] <- rownames(filtered_rows)
  df_for_filtered_rows["is_union_genenames"] <- c(rep("FILTERED_OUT", nrow(df_for_filtered_rows)))
  final_BIG_table <- rbind(final_BIG_table, df_for_filtered_rows)

  final_BIG_table <- final_BIG_table[order(final_BIG_table["combined_pvalues"]),]
  return(final_BIG_table)
}



generate_DE_report <- function(){
  template_path_DEreport <- file.path(main_path_script, 'templates', 'generate_report.tex')
  current_template_path_DE <- file.path(paths$root, 'generate_report.tex')
  latex_file_path_DE <- file.path(paths$root, 'DE_report.tex')
  file.copy(template_path_DEreport, paths$root, overwrite = TRUE)
  opts_chunk$set(echo=FALSE)
  knit(current_template_path_DE , output = latex_file_path_DE)
  cmd <- paste('pdflatex', latex_file_path_DE, sep=' ')
  system(cmd)
}


creating_final_BIG_table <- function(all_counts_for_plotting, complete_alldata_df, all_FDR_names, all_LFC_names, all_pvalue_names, final_pvalue_names, final_logFC_names, final_FDR_names, all_data, DEG_pack_columns){
  
  is_union_genenames <- rownames(all_genes_df) %in% complete_alldata_df  #checking which genes are union results
  names(is_union_genenames) <- rownames(all_genes_df)
  all_genes_df <- check_deg_in_pck(all_counts_for_plotting, all_data, all_genes_df, DEG_pack_columns)


  ############### Establishing gene label requirements and assigning (POSSIBLE, PREVALENT, NOTDEG and FILTERED_OUT) ##########
  labeling_results <- labeling_genes(all_genes_df, DEG_pack_columns, is_union_genenames) 
  coincid_counter <- labeling_results[[1]]
  is_union_genenames <- labeling_results[[2]]

  un_common_rej_df <- as.data.frame(is_union_genenames)
  final_BIG_table <- merge(un_common_rej_df, all_genes_df, by.x="row.names", by.y="Row.names")

  coincidences <- as.data.frame(coincid_counter)
  final_BIG_table <- cbind(final_BIG_table, coincidences)

  combined_pvalues <- calculating_combined_pvalue_per_geneID(final_BIG_table, final_FDR_names)
  combined_pvalues_column <- as.data.frame(combined_pvalues)
  final_BIG_table <- merge(final_BIG_table, combined_pvalues_column, by.x="Row.names", by.y="row.names")
  print(head(final_BIG_table))

  combined_nominal_pvalues <- calculating_combined_nominal_pvalue_per_geneID(final_BIG_table, final_pvalue_names)
  combined_nominal_pvalues_column <- as.data.frame(combined_nominal_pvalues)
  final_BIG_table <- merge(final_BIG_table, combined_nominal_pvalues_column, by.x="Row.names", by.y="row.names")
  print(head(final_BIG_table))

  final_BIG_table$BH = p.adjust(final_BIG_table$combined_nominal_pvalues, method = "BH")
  print(head(final_BIG_table))

  final_BIG_table$Bonferroni = p.adjust(final_BIG_table$combined_nominal_pvalues, method = "bonferroni")
  print(head(final_BIG_table))

  final_BIG_table$Holm = p.adjust(final_BIG_table$combined_nominal_pvalues, method = "holm")
  print(head(final_BIG_table))

  final_BIG_table$Hochberg = p.adjust(final_BIG_table$combined_nominal_pvalues, method = "hochberg")
  print(head(final_BIG_table))

  final_BIG_table$Hommel = p.adjust(final_BIG_table$combined_nominal_pvalues, method = "hommel")
  print(head(final_BIG_table))

  final_BIG_table$BY = p.adjust(final_BIG_table$combined_nominal_pvalues, method = "BY")
  print(head(final_BIG_table))  

  ################### combined pvalue labeling ##############
  pval_labeling <- labeling_comb_pvalue(combined_pvalues) #checks if combined pvalue is < 0.05 (significative)
  pval_labeling_col <- as.data.frame(pval_labeling)
  final_BIG_table <- cbind(final_BIG_table, pval_labeling_col)

  ################### calculating mean logFC per gene ##############
  mean_logFCs <- calculating_logFC_mean(final_BIG_table)
  logFC_means_column <- as.data.frame(mean_logFCs)
  final_BIG_table <- merge(final_BIG_table, logFC_means_column, by.x="Row.names", by.y="row.names")

  ########################################################

  ################## Adding filtered genes ###########
  final_BIG_table <- adding_filtered_transcripts(raw, raw_filter, final_BIG_table)
  return(final_BIG_table)  
}



