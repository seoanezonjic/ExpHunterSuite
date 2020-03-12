#! /usr/bin/env Rscript

#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>

#############################################
### CONFIGURE 
#############################################

#Loading libraries
suppressPackageStartupMessages(library(optparse)) 

option_list <- list(
  make_option(c("-i", "--inputfile"), type="character", default = NULL,
    help="DEGenesHunter expression resulta table file"),
  make_option(c("-f", "--logFC_threshold"), type="double", default=1,
    help="Log Fold Change threshold. Default : %default"),
  make_option(c("-F", "--logFC_max"), type="double", default=NULL,
    help="[OPTIONAL] Maximum Log Fold Change. Default : maximum observed"),
  make_option(c("-p", "--pval_threshold"), type="double", default=0.01,
    help="P-value threshold. Default : %default"),
  make_option(c("-o", "--outfile"), type="character",
    help="Output file")
)

opt <- parse_args(OptionParser(option_list=option_list))

#############################################
### FUNCTIONS 
#############################################

#' Funtion used to scale a vector of numeric values from [x_min,x_max] range to [y_min,y_max] range.
#' 	f : R -> R
#'      [x_min, x_max] -> [y_min, y_max]
#'                         vect - x_min
#'  f(vect,y_min,y_max) =  ------------- x (y_max - y_min) + y_min
#'                         x_max - x_min
#'
#' @param vect vector to be transformated
#' @param nmin new minimum
#' @param nmax new maximum
#' @return transformated vector to new range
scale_range <- function(vect,nmin,nmax){
	((vect - min(vect))/(max(vect)-min(vect)))*(nmax-nmin)+nmin
}


#'
#'
accordion_range <- function(vect,nmin,nmax,ncenter = NULL,vthreshold){
	# Divide vector in over/under center items
	indx_under <- which(vect < vthreshold)
	indx_over <- which(vect >= vthreshold)
	# Check center
	if(is.null(ncenter)) ncenter <- nmin + (nmax - nmin)/2
	# Re-range vectors
	vect_under <- scale_range(vect[indx_under],nmin,ncenter)
	vect_over <- scale_range(vect[indx_over],ncenter,nmax)
	# Reorganize
	vect[indx_under] <- vect_under
	vect[indx_over] <- vect_over
	# END
	return(vect)
}


#############################################
### LOAD & PREPARE 
#############################################

# Load DEGenesHunter result table
dgh_res <- read.table(file = opt$inputfile, sep = "\t", header = T, stringsAsFactors = FALSE, row.names = 1)

# Prepare target columns
pck_names <- c("DESeq2","edgeR","limma","NOISeq","Combined")
fc_cols <- c("logFC_DESeq2","logFC_edgeR","logFC_limma","logFC_NOISeq","mean_logFCs")
pv_cols <- c("FDR_DESeq2","FDR_edgeR","FDR_limma","FDR_NOISeq","combined_FDR")
# Check
aux <- fc_cols %in% colnames(dgh_res)
pck_names <- pck_names[aux]
fc_cols <- fc_cols[aux]
pv_cols <- pv_cols[aux]

# Config new sclae
score_min <- 0
score_max <- 1
score_center <- 0.5

#############################################
### CALC SCORES 
#############################################

# Calculate scores for each set
scores <- as.data.frame(do.call(cbind,lapply(seq_along(pck_names),function(i){
	# Check
	if(any(dgh_res[,fc_cols[i]] > opt$logFC_max, na.rm = TRUE)) dgh_res[(dgh_res[,fc_cols[i]] > opt$logFC_max & !is.na(dgh_res[,fc_cols[i]])),] <- opt$logFC_max
	# Apply rescale to score
	scrs <- accordion_range(dgh_res[,fc_cols[i]],score_min,score_max,score_center,opt$logFC_threshold)
	# Adjust
	scrs <- scrs - (dgh_res[,pv_cols[i]] > opt$pval_threshold & dgh_res[,fc_cols[i]] < opt$logFC_threshold)*score_center
	if(any(scrs < score_min)) scrs[scrs < score_min] <- score_min
	# Return
	aux <- data.frame(Scores = scrs)
	colnames(aux) <- pck_names[i]
	return(aux)
})))

# Add gene names
scores <- cbind(list(Gene = rownames(dgh_res)), scores)

#############################################
### EXPORT
#############################################
write.table(scores, file = opt$outfile, quote = FALSE, col.names = TRUE, sep = "\t", row.names = FALSE)