rtable2measures <- function(
  htfile,
  realprediction,
  aucfile = NULL,
  experiment = NULL,
  logFC_threshold = 1,
  pval_threshold = 0.01
  ){

  # Load hunter table
  dgh_res_raw <- read.table(file = htfile, sep = "\t", header = T, stringsAsFactors = FALSE, row.names = 1)

  # Load real predictions
  true_pred <- read.table(file = realprediction, sep = "\t", header = T, stringsAsFactors = FALSE)

  # Sort
  dgh_res_raw <- dgh_res_raw[true_pred[,1],]

  # Prepare target columns
  pck_names <- c("DESeq2_DEG","edgeR_DEG","limma_DEG","NOISeq_DEG")
  pck_names <- pck_names[pck_names %in% colnames(dgh_res_raw)]
  dgh_res <- dgh_res_raw[,pck_names]
  dgh_res$Gene <- rownames(dgh_res)

  # CALCULATE AUCS
  if(!is.null(aucfile)){
    pck_fdrnames <- c("FDR_DESeq2","FDR_edgeR","FDR_limma","FDR_NOISeq", "combined_FDR")
    dgh_fdr <- dgh_res_raw[, pck_fdrnames]
    dgh_fdr[is.na(dgh_fdr)] <- 1
    dgh_fdr <- 1 - dgh_fdr
    real_df <- as.data.frame(replicate(ncol(dgh_fdr), true_pred[,2]))
    pred <- prediction(dgh_fdr, real_df)
    perf_auc <- performance(pred, "auc")
    auc_vals <- data.frame(Method = gsub("_FDR","",gsub("FDR_","",colnames(dgh_fdr))), 
                 Set = rep("All",ncol(dgh_fdr)),
                 Measure = rep("auc", ncol(dgh_fdr)),
                 Value = unlist(perf_auc@y.values), stringsAsFactors = FALSE) 
    # Write
    write.table(auc_vals, file = aucfile, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)  
  }



  # Calculate values for combined values
  dgh_res$Combined <- (abs(dgh_res_raw$mean_logFCs) >= logFC_threshold) & (dgh_res_raw$combined_FDR <= pval_threshold)


  # Combine into unique matrix
  preds_res <- merge(x = dgh_res, y = true_pred, by = "Gene")
  genes <- preds_res$Gene
  preds_res <- as.matrix(preds_res[,-1])
  preds_res[is.na(preds_res)] <- FALSE# Substitute NAs

  # Prepare matrix for vote cuts
  votes <- unlist(lapply(seq(nrow(preds_res)),function(i){
    sum(head(preds_res[i,],-1)) 
  }))

  # Add cuts
  cut_names <- unlist(lapply(seq(length(pck_names)),function(pckcut){
    preds_res <<- cbind(preds_res,votes >= pckcut)
    return(paste0("Cut_",pckcut))
  }))


  colnames(preds_res) <- c(pck_names,"Combined","Prediction",cut_names)
  rownames(preds_res) <- genes

  cut_names <- c(cut_names,"Combined")

  #############################################
  ### CALC MEASURES
  #############################################

  # Calculate TP,TN,FP,FN
  df_cuts <- as.data.frame(do.call(rbind,lapply(cut_names,function(cut){
    return(data.frame(Cut = cut,
              TP = sum(preds_res[,cut] & preds_res[,"Prediction"]),
              FP = sum(preds_res[,cut] & !preds_res[,"Prediction"]),
              TN = sum(!preds_res[,cut] & !preds_res[,"Prediction"]),
              FN = sum(!preds_res[,cut] & preds_res[,"Prediction"])))
  })))

  # Calculate values for

  # Calculate measures
  df_cuts$ACC <- acc(df_cuts)           # Accuracy
  df_cuts$Precision <- ppv(df_cuts)     # Precision
  df_cuts$Recall <- recall(df_cuts)     # Recall / sensitivity
  df_cuts$Specificity <- spc(df_cuts)   # Specificity
  df_cuts$FMeasure <- fmeasure(df_cuts) # F-measure

  #############################################
  ### EXPORT
  #############################################
  # Concat experiment info (if proceed)
  if(!is.null(experiment)){
    exp_info <- read.table(file = experiment, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    exp_info <- as.data.frame(do.call(cbind,lapply(seq(ncol(exp_info)),function(i){
      aux <- data.frame(Value = rep(exp_info[1,i],nrow(df_cuts)))
      colnames(aux) <- colnames(exp_info)[i]
      return(aux)
    })))
    df_cuts <- cbind(exp_info,df_cuts)  
  }

  return(df_cuts)

}
