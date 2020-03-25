################################# GENERAL FUNCTIONS ##############################################

handling_errors <- function(a){
  normalized_counts <- NULL
  expres_diff <- NULL
  all_genes_df <- NULL
  return(list(expres_diff, normalized_counts, all_genes_df))
}

########################################################
# Functions to generate output files
unite_DEG_pack_results <- function(all_DE, all_FDR_names, all_LFC_names, all_pvalue_names, final_pvalue_names, final_logFC_names, final_FDR_names, p_val_cutoff, lfc, minpack_common) {
  # Reorder all output by gene id so all equal
  all_DE <- lapply(all_DE, function(x) x[order(row.names(x)),])
  gene_ids <- row.names(all_DE[[1]])
  
  # Initialise output dataframe and add p/FDR/logFC values from the different packages
  all_DE_df <- data.frame(row.names=gene_ids)
  for(i in 1:length(all_DE)) {
    all_DE_df[final_logFC_names[i]] <- all_DE[[i]][, all_LFC_names[i]]
    all_DE_df[final_FDR_names[i]] <- all_DE[[i]][, all_FDR_names[i]]
    all_DE_df[final_pvalue_names[i]] <- all_DE[[i]][, all_pvalue_names[i]]
  }

  # Add TRUE/FALSE for each DEG package
  for(i in 1:length(all_DE)) {
    all_DE_df[DEG_pack_columns[i]] <- all_DE_df[, final_FDR_names[i]] < p_val_cutoff & abs(all_DE_df[, final_logFC_names[i]]) >= lfc
    all_DE_df[DEG_pack_columns[i]][is.na(all_DE_df[DEG_pack_columns[i]])] <- FALSE
    # Check: if no DE genes for package, give warning
    if(sum(all_DE_df[, DEG_pack_columns[i]]) == 0) warning(paste("No significant", DEG_pack_columns[i], "found"))
  }
  # Get DEG_counts
  all_DE_df["DEG_counts"] <- rowSums(all_DE_df[DEG_pack_columns], na.rm = TRUE)
 
  # Calc and add combined FDR-values
  log_FDR <- log(all_DE_df[final_FDR_names]) # Log all final p-values
  log_FDR[is.na(log_FDR)] <- 0 # any NAs made to 0 so as not to contribute to the combined score
  if("FDR_NOISeq" %in% final_FDR_names){ # NOISeq can give FDR values of 0 - these become 0 when -logged:
    log_FDR[,"FDR_NOISeq"][log_FDR[,"FDR_NOISeq"] == -Inf] <- sort(unique(log_FDR[, "FDR_NOISeq"]))[2]  # Give NOISeq -Inf values smallest possible value
  }

  xi_squared <-  -2 * rowSums(log_FDR)
  degrees_freedom <- 2 * length(final_FDR_names)
  combined_FDR <- pchisq(xi_squared, degrees_freedom, lower.tail = FALSE)
  all_DE_df[,"combined_FDR"] <- combined_FDR

  # Reorder by combined FDR value
  all_DE_df <- all_DE_df[order(all_DE_df[,"combined_FDR"]), ]
  
  # Label as significant or not using combined FDR values
  all_DE_df[, "FDR_labeling"] <- ifelse(all_DE_df[, "combined_FDR"] < p_val_cutoff, "SIGN", "NOTSIGN")

  # Calculate average fold changes
  all_DE_df[, "mean_logFCs"] <- rowMeans(all_DE_df[final_logFC_names])

  # Add PREVALENT_DEG tag if as many as minpack_common; POSSIBLE_DEG if less but > 0; NOT_DEG if == 0
  genes_tag <- ifelse(all_DE_df[, "DEG_counts"] >= minpack_common, yes = "PREVALENT_DEG",
         no = ifelse(all_DE_df[, "DEG_counts"] > 0, yes = "POSSIBLE_DEG", no = "NOT_DEG")
  )
  all_DE_df[, "genes_tag"] <- genes_tag

  return(all_DE_df)
}

add_filtered_genes <- function(all_DE_df, raw) {
  filtered_genes <- row.names(raw)[! row.names(raw) %in% row.names(all_DE_df)]
  filtered_df <- as.data.frame(matrix(NA, ncol = ncol(all_DE_df), nrow = length(filtered_genes)))
  colnames(filtered_df) <- colnames(all_DE_df)
  row.names(filtered_df) <- filtered_genes
  filtered_df[, "genes_tag"] <- rep("FILTERED_OUT", dim(filtered_df)[1])
  all_DE_df <- rbind(all_DE_df, filtered_df)

  return(all_DE_df)
}

write_df_list_as_tables <- function(df_list, prefix, root) {
  invisible(
    lapply(1:length(df_list), function(i) {
      pack <- names(df_list)[i]
      write.table(df_list[[pack]],
      file=file.path(root, paste0("Results_", pack), paste0(prefix, pack, '.txt')), quote=FALSE, col.names = NA, sep="\t")
    })
  )
}

debug_point <- function(file, message = "Debug point"){
    debug_message <<- message
    save.image(file)
}


save_times <- function(time_control, output="times_control.txt", plot_name = "time_control.pdf"){
 spent_times <- list()
    invisible(lapply(seq(2,length(time_control)), function(time_control_i){
      spent_times[[names(time_control)[time_control_i]]] <<- as.numeric(unlist(time_control[time_control_i]) - unlist(time_control[time_control_i - 1]))
    }))
    spent_times_df <- do.call(rbind.data.frame, spent_times)
    colnames(spent_times_df) <- c("time")
    spent_times_df$control <- names(spent_times)
    pp <- ggplot(spent_times_df, aes(x= control, y = time)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x=element_text(angle = 25, hjust=1))
    pdf(file.path(dirname(output), plot_name))    
      plot(pp)
    dev.off()
    write.table(spent_times_df, file=output, quote=FALSE, row.names=FALSE, sep="\t")
}

