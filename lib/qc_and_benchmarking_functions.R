#############################################
### QC AND BENCHMARKING FUNCTIONS 
#############################################
preparing_rlog_PCA <- function(raw_filter,design_vector){
  require(DESeq2) # DESeqDataSetFromMatrix, rlog
  coldata_df <- data.frame(cond = design_vector,
                           each = colnames(raw_filter),
                           row.names = colnames(raw_filter))
  # rownames(coldata_df) <- colnames(raw_filter)
  dds <- DESeqDataSetFromMatrix(countData = raw_filter, colData = coldata_df, design = ~ cond)
  rld <- rlog(dds, blind=FALSE)
  return(rld)
}


counting_trues <- function(all_genes_df, DEG_pack_columns){
  summary_DEGs <- c()
  for (i in c(1:nrow(all_genes_df))){
    is_a_DEG <- as.logical(all_genes_df[i, DEG_pack_columns])
    sum_DEGs <- sum(is_a_DEG)
    summary_DEGs[i] <- sum_DEGs
    }
  return(summary_DEGs)
}

tagging_genes <- function(final_BIG_table, opt, DEG_counts, DEG_pack_columns){
  activated_packages <- length(DEG_pack_columns)
  tag <- c()
  counts <- final_BIG_table$DEG_counts
  for (i in c(1:nrow(final_BIG_table))){
    if (counts[i] == opt$minpack_common){
      tag[i] <- "PREVALENT_DEG"
    }
    else if (counts[i] > 0){
      tag[i] <- "POSSIBLE_DEG"
    }
    else {
      tag[i] <- "NOT_DEG"
    }
  }
  return(tag)
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
  statistics_report <- NULL
  union_names <- unite_result_names(all_data) 
  for (i in c(1:length(all_data))){
    all_results_in_package <- nrow(all_data[[i]])
    if(length(all_results_in_package == 0)) next
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


combine_log_fischer <- function(adjp_values){
  value_range <- (adjp_values > 0) & (adjp_values <= 1)
  log_adjp <- log(adjp_values[value_range])
  xi_squared <- (-2) * sum(log_adjp)
  freedom_degree <- 2 * length(log_adjp)
  combined_pvalue <- pchisq(xi_squared, freedom_degree, lower.tail = FALSE)
  if (length(log_adjp) != length(adjp_values)) {
      warning("Some adjusted p-values probably are set cero in some package")
  }
  return(combined_pvalue)
}


calculating_combined_nominal_pvalue_per_geneID <-function(final_BIG_table, final_pvalue_names){
  genenames <- c()
  combined_nominal_pvalues <- c()
    for (i in c(1:nrow(final_BIG_table))){
      genenames[i] <- final_BIG_table[i, "Row.names"]

      nominalp_values <- as.numeric(final_BIG_table[i, final_pvalue_names])
      if (any(is.na(nominalp_values))){
        combined_nominal_pvalues[i] <- "NA"
      } else {
      combined_nominal_pvalues[i] <- combine_log_fischer(nominalp_values)
    }
  }
  
  combined_nominal_pvalues <- as.numeric(combined_nominal_pvalues)
  names(combined_nominal_pvalues) <- genenames
  return(combined_nominal_pvalues)
}


labeling_comb_pvalue <- function(combined_pvalues){
  pval_labeling <- c()
  for (i in c(1:length(combined_pvalues))){
    if (is.na(combined_pvalues[i])){
      pval_labeling[i] <- "NOTSIGN"
    } else{
    if (combined_pvalues[i] < 0.05){
      pval_labeling[i] <- "SIGN"
    }
    if (combined_pvalues[i] > 0.05){    
      pval_labeling[i] <- "NOTSIGN"
    }
  }}
  return(pval_labeling)
}


calculate_sensitivity <- function(True_Pos, False_Neg){
  sensitivity <- True_Pos/(True_Pos + False_Neg)
  return(sensitivity)
}

calculate_specificity <- function(True_Neg, False_Pos){
  specificity <- True_Neg/(False_Pos + True_Neg)
  return(specificity)
} 

calculate_positive_predictive_value_PPV <- function(True_Pos, False_Pos){
  positive_predictive_value_PPV <- True_Pos/(True_Pos + False_Pos)
  return(positive_predictive_value_PPV)
}

calculate_negative_predictive_value_NPV <- function(True_Neg, False_Neg){
  negative_predictive_value_NPV <- True_Neg/(False_Neg + True_Neg)
  return(negative_predictive_value_NPV)
}

calculate_accuracy <- function(True_Pos, True_Neg, False_Pos, False_Neg){
  negative_predictive_value_NPV <- True_Neg/(False_Neg + True_Neg)
  return(negative_predictive_value_NPV)
}




calculating_logFC_mean <- function(final_BIG_table){
  geneids <- c()
  mean_logFCs <- c()
  for (i in c(1:nrow(final_BIG_table))){
      geneids[i] <- final_BIG_table[i, "Row.names"]
      logFC_values <- as.numeric(final_BIG_table[i, final_logFC_names])
      mean_logFCs[i] <- mean(logFC_values)
  }
  names(mean_logFCs) <- geneids
  return(mean_logFCs)
} 



creating_genenumbers_barplot <- function(raw, raw_filter, all_data, x_all){
  complete_alldata_df <- unite_all_rownames_from_dataframes_list(all_data)
  gene_numbers <- c(nrow(raw), nrow(raw_filter), length(complete_alldata_df), length(x_all))
  names <- c("Raw counts", "Filtered raw counts","All possible DEGs","Prevalent DEGs")
  barplot_df <- data.frame(numbers=gene_numbers, cat=names)
  barplot_df$cat <- factor(barplot_df$cat, levels = barplot_df$cat[order(barplot_df$numbers)])
  return(barplot_df)
}


creating_top20_table <- function(final_BIG_table){
  final_BIG_table <- final_BIG_table[order(final_BIG_table["combined_pvalues"]),]
  best_20_genes <- final_BIG_table[1:20,]
  write.table(best_20_genes, file=file.path(paths$root, "top20_genes.txt"), quote=F, col.names=NA, sep="\t")
}
