#! /usr/bin/env Rscript

#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>

#############################################
### CONFIGURE 
#############################################

#Loading libraries
suppressPackageStartupMessages(library(optparse)) 

option_list <- list(
  make_option(c("-i", "--inputfile"), type="character", 
    help="DEGenesHunter expression resulta table file"),
  make_option(c("-r", "--realprediction"), type="character",
    help="Real prediction values file"),
  make_option(c("-e", "--experiment"), type="character", default = NULL,
    help="[OPTIONAL] Experiment metrics file to be included to output."),
  make_option(c("-f", "--logFC_threshold"), type="double", default=1,
    help="Log Fold Change threshold. Default : %default"),
  make_option(c("-p", "--pval_threshold"), type="double", default=0.01,
    help="P-value threshold. Default : %default"),
  make_option(c("-o", "--outfile"), type="character",
    help="Output file")
)

opt <- parse_args(OptionParser(option_list=option_list))

#############################################
### FUNCTIONS 
#############################################

# ACC
acc <- function(df){
	(df$TP + df$TN) / (df$TP + df$TN + df$FP + df$FN)
}

# PPV
ppv <- function(df){
	df$TP / (df$TP + df$FP)
}

# Recall
recall <- function(df){
	df$TP / (df$TP + df$FN)
}

# SPC
spc <- function(df){
	df$TN / (df$TN + df$FP)
}

# F-measure
fmeasure <- function(df){
	(df$Precision * df$Recall) / (df$Precision + df$Recall)
}


#############################################
### LOAD & PREPARE
#############################################

# Load hunter table
dgh_res_raw <- read.table(file = opt$inputfile, sep = "\t", header = T, stringsAsFactors = FALSE, row.names = 1)

# Prepare target columns
pck_names <- c("DESeq2_DEG","edgeR_DEG","limma_DEG","NOISeq_DEG")
pck_names <- pck_names[pck_names %in% colnames(dgh_res_raw)]
dgh_res <- dgh_res_raw[,pck_names]
dgh_res$Gene <- rownames(dgh_res)

# Calculate values for combined values
dgh_res$Combined <- (abs(dgh_res_raw$mean_logFCs) >= opt$logFC_threshold) & (dgh_res_raw$combined_FDR <= opt$pval_threshold)

# Load real predictions
true_pred <- read.table(file = opt$realprediction, sep = "\t", header = T, stringsAsFactors = FALSE)

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
if(!is.null(opt$experiment)){
	exp_info <- read.table(file = opt$experiment, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
	exp_info <- as.data.frame(do.call(cbind,lapply(seq(ncol(exp_info)),function(i){
		aux <- data.frame(Value = rep(exp_info[1,i],nrow(df_cuts)))
		colnames(aux) <- colnames(exp_info)[i]
		return(aux)
	})))
	df_cuts <- cbind(exp_info,df_cuts)	
}

# Store results table
write.table(df_cuts, file = opt$outfile, quote = FALSE, col.names = TRUE, sep = "\t", row.names = FALSE)