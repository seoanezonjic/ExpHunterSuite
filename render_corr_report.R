#! /usr/bin/env Rscript
#############################################
############## FUNCTIONAL HUNTER ###########
#############################################

# this is wrapped in a tryCatch. The first expression works when source executes, the
# second expression works when R CMD does it.
full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile),  # works when using source
               error=function(e) # works when using R CMD
              normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', commandArgs())], '='))[2]))
main_path_script <- dirname(full.fpath)


#Loading libraries  
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(knitr))


#############################################
### MAIN 
#############################################

# Parse command line
#------------------------------------------------

option_list <- list(
  make_option(c("-i", "--input_hunter_folder"), type="character",
    help="DEgenes Hunter's differential expression analysis output folder"), 
  make_option(c("-o", "--output_files"), type="character", default="results",
    help="Output path. Default=%default")
)
opt <- parse_args(OptionParser(option_list=option_list))


############ CREATE FOLDERS #########3
paths <- list()
dir.create(opt$output_files)
paths$root <-opt$output_files



#############################################
### LOAD AND PARSE 
#############################################

DEGH_results <- read.table(file.path(opt$input_hunter_folder, "Common_results", "hunter_results_table.txt"), header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE)
aux <- which(DEGH_results$genes_tag == "FILTERED_OUT") 
if(length(aux) > 0){
	DEGH_results <- DEGH_results[-aux,]
}


#############################################
### PREPARE AND TRANSFORM DATA
#############################################

####
# LOAD NORMALIZED COUNTS
norm_counts_raw <- as.matrix(read.table(file.path(opt$input_hunter_folder, "Results_DESeq2", "Normalized_counts_DESeq2.txt"), header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE))

# Normalize
norm_counts_raw_gnorm <- norm_counts_raw
invisible(lapply(seq_along(norm_counts_raw[,1]),function(i){
	m <- min(norm_counts_raw[i,])
	M <- max(norm_counts_raw[i,])
	dff <- M - m
	norm_counts_raw_gnorm[i,] <<- (norm_counts_raw[i,] - m) / dff
}))
# Modify to plot better later
norm_counts <- as.data.frame(as.table(norm_counts_raw))
colnames(norm_counts) <- c("Gene","Sample","Count")
norm_counts_gnorm <- as.data.frame(as.table(norm_counts_raw_gnorm))
colnames(norm_counts_gnorm) <- c("Gene","Sample","Count")
# norm_counts_gnorm <- cbind(norm_counts_gnorm,list(Type = rep("Regular",nrow(norm_counts_gnorm))))

####
# LOAD WGCNA clusters representative profiles with samples
cl_eigvalues <- as.matrix(read.table(file.path(opt$input_hunter_folder, "Results_WGCNA", "eigen_values_per_samples.txt"), header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE))
cl_eigvalues <- as.data.frame(as.table(cl_eigvalues),stringsAsFactors = FALSE)
colnames(cl_eigvalues) <- c("Sample","Cluster_ID","Count") 
cl_eigvalues_gnorm <- cl_eigvalues
cl_eigvalues_gnorm$Count <- (cl_eigvalues_gnorm$Count + 1) / 2 

####
# LOAD WGCNA - PVal (Cluster - Trait)
wgcna_pval_cl_trait <- as.matrix(read.table(file.path(opt$input_hunter_folder, "Results_WGCNA", "module_trait_p_val.txt"), header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE))
wgcna_corr_cl_trait <- as.matrix(read.table(file.path(opt$input_hunter_folder, "Results_WGCNA", "module_trait.txt"), header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE))
####
# LOAD WGCNA - Correlation (Sample - Trait)
wgcna_count_sample_trait <- as.matrix(read.table(file.path(opt$input_hunter_folder, "Results_WGCNA", "sample_trait.txt"), header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE))
wgcna_count_sample_trait_gnorm <- as.data.frame(do.call(cbind,lapply(seq(ncol(wgcna_count_sample_trait)),function(j){
	(wgcna_count_sample_trait[,j] - min(wgcna_count_sample_trait[,j],na.rm = TRUE)) / (max(wgcna_count_sample_trait[,j],na.rm = TRUE) - min(wgcna_count_sample_trait[,j],na.rm = TRUE))
})))
colnames(wgcna_count_sample_trait_gnorm) <- colnames(wgcna_count_sample_trait)


# Obtain clusters
cls <- unique(DEGH_results$Cluster_ID)
if(any(c(0,"grey") %in% cls)){
	cls <- cls[!cls %in% c(0,"grey")]
}else{
	warning("Cluster Zero/Grey not found")
}
clgenes <- lapply(cls,function(cl){unique(rownames(DEGH_results[which(DEGH_results$Cluster_ID == cl),]))}) # Find
names(clgenes) <- cls


############################################################
##                    GENERATE REPORT                     ##
############################################################
results_path <- normalizePath(paths$root)

invisible(lapply(cls,function(cl){
	# Take output name
	aux <- paste0("cl_func_",cl,".html")
	outf_cls_i <- file.path(results_path, aux)
	# Generate report
	rmarkdown::render(file.path(main_path_script, 'templates', 'corrprofiles_report.Rmd'), output_file = outf_cls_i, intermediates_dir = results_path)
}))
