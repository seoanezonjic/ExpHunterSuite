#! /usr/bin/env Rscript


#' @author refactor by Fernando Moreno Jabato. Original authors Isabel Gonzalez Gayte


######################################################################################################
############################################# DEgenes Hunter #########################################
######################################################################################################


############################################################
##                      SETUP PROGRAM                     ##
############################################################

# Loading libraries
suppressPackageStartupMessages(library(ggplot2)) 
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(DESeq2)) 
suppressPackageStartupMessages(library(NOISeq))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(FSA))

suppressPackageStartupMessages(require(rmarkdown))


# Obtain this script directory
full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile),  # works when using source
               error=function(e) # works when using R CMD
              normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', commandArgs())], '='))[2]))
main_path_script <- dirname(full.fpath)

# Load custom libraries
source(file.path(main_path_script, 'lib', 'general_functions.R'))
source(file.path(main_path_script, 'lib', 'dif_expression_packages.R'))
source(file.path(main_path_script, 'lib', 'qc_and_benchmarking_functions.R'))



# Prepare command line input 

option_list <- list(
  make_option(c("-i", "--input_file"), type="character", default=NULL,
    help="Input file with read counts"),
  make_option(c("-C", "--Control_columns"), type="character", default=NULL,
    help="Control columns. Please indicate column names of control samples separated by commas"),
  make_option(c("-T", "--Treatment_columns"), type="character", default=NULL,
    help="Treatment columns. Please indicate column names of treatment samples separated by commas"),  
  make_option(c("-r", "--reads"), type="integer", default=2,
    help="Number of minimum reads required. Lesser number of reads are discarded. Default=%default.
    0 = No filtering"),
  make_option(c("-l", "--minlibraries"), type="integer", default=2,
    help="Number of minimum reads required. Lesser number of reads are discarded. Default=%default"),
  make_option(c("-o", "--output_files"), type="character", default="hunter_DE_results",
    help="Output path. Default=%default"),
  make_option(c("-p", "--p_val_cutoff"), type="double", default=0.05,
    help="Adjusted p-value for the differential expression analysis. Default=%default"),
  make_option(c("-f", "--fc"), type="double", default=1.5,
    help="Minimum fold expression change . Default=%default"),
  make_option(c("-q", "--q_value"), type="double", default=0.95,
    help="q value for NOISeq package. Default=%default"),
  make_option(c("-a", "--adjust_method"), type="character", default=c("BH"), #D = DESeq2, E = edgeR, L = limma, N = NOISeq
    help="Method selection to adjust the combined nominal p-values. By default method Default=%default is performed"),
  make_option(c("-n", "--name_exp"), type="character", default="experiment1",
    help="Type the name of your experiment."),
  make_option(c("-m", "--modules"), type="character", default=c("DELN"), #D = DESeq2, E = edgeR, L = limma, N = NOISeq
    help="Differential expression packages to able/disable (D = DESeq2, E = edgeR, L = limma, N= NOISeq.).
    By default the following modules Default=%default are performed"),
  make_option(c("-c", "--minpack_common"), type="integer", default=4,
    help="Number of minimum package to consider a gene as a 'PREVALENT' DEG")
)
opt <- parse_args(OptionParser(option_list=option_list))

active_modules <- nchar(opt$modules)
if(opt$minpack_common > active_modules){
	opt$minpack_common <- active_modules
	warning("The number of active modules is lower than the thresold for tag PREVALENT DEG. The thresold is set to the number of active modules.")
}





############################################################
##                       CHECK INPUTS                     ##
############################################################

# Check inputs not allowed values
if (is.null(opt$input_file)){
  stop(cat("No file with RNA-seq counts is provided.\nPlease use -i to submit file"))
}
if (is.null(opt$Control_columns)){
  stop(cat("No control samples are indicated.\nPlease select column names of the count table with -C"))
}
if (is.null(opt$Treatment_columns)){
  stop(cat("No treatment samples are indicated.\nPlease select column names of the count table with -T"))
}
if (opt$minlibraries < 1){
  stop(cat("Minimum library number to check minimum read counts cannot be less than 1."))
}
if (opt$reads == 0){
  cat("Minimum reads option is set to cero. Raw count table will not be filtered.")
}


# Control replicate number VS method selection
index_control_cols <- unlist(strsplit(opt$Control_columns, ",")) 
index_treatmn_cols <- unlist(strsplit(opt$Treatment_columns, ","))
replicatesC <- length(index_control_cols)
replicatesT <- length(index_treatmn_cols)

# Chekc if there are not enough replicates for specified method
if((sum(replicatesC, replicatesT)<3) &
    (((grepl("E", opt$modules)) | 
      (grepl("L", opt$modules)) | 
      (grepl("N", opt$modules))))){
  stop("Not enough replicates to perform an analysis with the selected method. Select 'D' with parameter -m")
}else if((sum(replicatesC, replicatesT)<=5) &
        ((grepl("L", opt$modules)) | 
          (grepl("N", opt$modules)))){
  stop("Not enough replicates to perform an analysis with the selected method")
}








############################################################
##                         I/O DATA                       ##
############################################################

# Create output container folders 
dir.create(opt$output_files)
paths <- list()
paths$root <- opt$output_files

# Creating Subfolders
subfolders <- defining_subfolders(replicatesC, replicatesT, opt$modules)
paths <- create_subfolders(subfolders, paths)

# Load raw count data
raw <- read.table(opt$input_file, header=T, row.names=1, sep="\t")
raw <- raw[c(index_control_cols,index_treatmn_cols)] #Indexing selected columns from input count file

# Substitute NA values
raw[is.na(raw)] <- 0







############################################################
##                          FILTER                        ##
############################################################
# Prepare filtered set
if(opt$reads != 0){
  keep_cmp <- rowSums(cpm(raw) > opt$reads) >= opt$minlibraries # two reads at least in two libraries
  raw_filter <- raw[keep_cmp,] # Filter out count data frame
  write.table(raw_filter, file=file.path(paths$root, "filtered_count_data.txt"),
                                          quote=F, col.names=NA, sep="\t")
}else{ # NO FILTER
  raw_filter <- raw
}




# Defining contrast - Experimental design
design_vector <- c(rep("C", replicatesC), rep("T", replicatesT)) 




############################################################
##                        QC GRAPHS                       ## 
##                 (before normalization)                 ##
############################################################

# # Simple boxplot to see overall counts' distribution
# pdf(file.path(paths$root, "boxplot_rawcounts_distribution.pdf"), w=8.5, h=6) 
#   boxplot(raw_filter)
# dev.off()

# pdf(file.path(paths$root, "boxplot_before_normalization.pdf"), w=8.5, h=6)
#   max_mean <- max(apply(raw_filter, MARGIN = 2, function(x) mean(x, na.rm=TRUE)))
#   boxplot(raw_filter,  ylim=c(0, max_mean*10), cex.lab=0.8, cex.axis=0.8, notch=TRUE, col=(c(rep("gold",replicatesC),rep("darkgreen",replicatesT))))
# dev.off()

# if('Results_DESeq2' %in% names(paths)){
#   rld <- preparing_rlog_PCA(raw_filter, design_vector)
#   pdf(file.path(paths[['Results_DESeq2']], "PCAplot.pdf"))
#     plotPCA(rld, intgroup=c("cond", "each")) # DESeq2 package
#   dev.off()  
# }








############################################################
##                       D.EXP ANALYSIS                   ##
############################################################

# Prepare results containers
all_data                <- list()
all_data_normalized     <- list()
all_counts_for_plotting <- list()
all_package_results     <- list()
all_FDR_names           <- c()
all_LFC_names           <- c()
all_pvalue_names        <- c()
final_logFC_names       <- c()
final_FDR_names         <- c()
final_pvalue_names      <- c()
DEG_pack_columns        <- c()

all_plots               <- list()

# Calculate global parameters (log fold change value)
lfc <- calculate_lfc(opt$fc)



############## Verbose point
if((replicatesC == 1) & (replicatesT == 1)){ 
  warning('There are one replicate available. Only DESeq2 will be performed.\n')
}else{
  cat('There are', replicatesC, 'replicates in the control condition and', replicatesT, 'replicates in the treatment condition.\n')
}



#####
################## CASE D: DESeq2
#####
if(grepl("D",opt$modules)){
  # Verbose
  cat('Gene expression analysis is performed with DESeq2.\n')

  # Calculate results
  results <- analysis_DESeq2(data   = raw_filter, 
                             num_controls  = replicatesC, 
                             num_treatmnts = replicatesT, 
                             p_val_cutoff  = opt$p_val_cutoff, 
                             lfc    = lfc, 
                             groups = design_vector)
  # Store results
  all_data[['DESeq2']] <- results[[1]]
  all_data_normalized[['DESeq2']] <- results[[2]]
  all_counts_for_plotting[['DESeq2']] <- results[[3]]

  # Result Plot Visualization
  if (!is.null(all_counts_for_plotting[['DESeq2']])){
    all_FDR_names      <- c(all_FDR_names, 'padj')
    all_LFC_names      <- c(all_LFC_names, 'log2FoldChange')
    all_pvalue_names   <- c(all_pvalue_names, 'pvalue')
    final_pvalue_names <- c(final_pvalue_names, 'pvalue_DESeq2')
    final_logFC_names  <- c(final_logFC_names, 'logFC_DESeq2')
    final_FDR_names    <- c(final_FDR_names, 'FDR_DESeq2')
    DEG_pack_columns   <- c(DEG_pack_columns, 'DESeq2_DEG')

    # pdf(file.path(paths[['Results_DESeq2']], "MA_plot_DESeq2.pdf"), w=11, h=8.5)
    #   plotMA(all_counts_for_plotting[['DESeq2']], cex.lab=1.6, cex.axis=1.5) 
    # dev.off()
  } 
}




#####
################## CASE E: edgeR   (AT LEAST 2 REPLICLATES)
#####
if((replicatesC >= 2) & 
   (replicatesT >= 2) & 
   grepl("E", opt$modules)){ 

  # Verbose point
  cat('Gene expression analysis is performed with edgeR.\n')

  # Calculate results
  results <- tryCatch(
    # CODE
    analysis_edgeR(data   = raw_filter, 
                   p_val_cutoff = opt$p_val_cutoff, 
                   lfc    = lfc, 
                   paths  = paths,
                   groups = design_vector),
    # CATCH
    error = handling_errors, warning = handling_errors)
  # Store results
  all_data[['edgeR']] <- results[[1]]
  all_data_normalized[['edgeR']] <- results[[2]]
  all_counts_for_plotting[['edgeR']] <- results[[3]]
 
  # Result Plot Visualization
  if (!is.null(all_counts_for_plotting[['edgeR']])){
    all_FDR_names      <- c(all_FDR_names, 'FDR')
    all_LFC_names      <- c(all_LFC_names, 'logFC')
    all_pvalue_names   <- c(all_pvalue_names, 'PValue')
    final_pvalue_names <- c(final_pvalue_names, 'pvalue_edgeR')
    final_logFC_names  <- c(final_logFC_names, 'logFC_edgeR')
    final_FDR_names    <- c(final_FDR_names, 'FDR_edgeR')
    DEG_pack_columns   <- c(DEG_pack_columns, 'edgeR_DEG')

    # pdf(file.path(paths[['Results_edgeR']], "MA_plot_edgeR.pdf"), w=11, h=8.5)
    #   with(all_counts_for_plotting[['edgeR']], plot(logCPM, logFC, pch=20, main="edgeR: Fold change vs abundance", cex.lab=1.5, cex.axis=1.5))
    #   with(subset(all_counts_for_plotting[['edgeR']], FDR < opt$p_val_cutoff), points(logCPM, logFC, pch=20, col="red"))
    #   abline(h=c(-1,1), col="blue")
    # dev.off()
  }
}






#####
################## CASE L: limma   (AT LEAST 3 REPLICLATES)
#####
if((replicatesC >= 3) &
   (replicatesT >= 3) &
   grepl("L", opt$modules)){ 

  # Verbose
  cat('Gene expression analysis is performed with limma.\n')

  # Calculate results
  results <- analysis_limma(data = raw_filter, 
                            num_controls  = replicatesC, 
                            num_treatmnts = replicatesT, 
                            p_val_cutoff  = opt$p_val_cutoff, 
                            lfc  = lfc)

  # Store results
  all_data[['limma']]                <- results[[1]]    
  all_data_normalized[['limma']]     <- results[[2]]
  all_counts_for_plotting[['limma']] <- results[[3]]

  # Result Plot Visualization
  if (!is.null(all_counts_for_plotting[['limma']])){
    all_FDR_names      <- c(all_FDR_names, 'adj.P.Val')
    all_LFC_names      <- c(all_LFC_names, 'logFC')
    all_pvalue_names   <- c(all_pvalue_names, 'P.Value')
    final_pvalue_names <- c(final_pvalue_names, 'pvalue_limma')
    final_logFC_names  <- c(final_logFC_names, 'logFC_limma')
    final_FDR_names    <- c(final_FDR_names, 'FDR_limma')
    DEG_pack_columns   <- c(DEG_pack_columns, 'limma_DEG')

    k_limma <- rownames(all_counts_for_plotting[['limma']]) %in% rownames(all_data[['limma']])
    # pdf(file.path(paths[['Results_limma']], "Volcanoplot_limma.pdf"), w=11, h=8.5)
    #   plot(x=all_counts_for_plotting[['limma']]$logFC, y=-log10(all_counts_for_plotting[['limma']]$P.Value), xlab="logFC", ylab="logOdds", col=c("blue", "red") [k_limma+1], pch=20, main= c("groupsB-groupsA"), cex.lab=1.6, cex.axis=1.5)
    #   abline(v= opt$lfc, col="cyan")
    #   limit.pval_limma <- -log10(max(all_data[['limma']]$P.Value)) 
    #   abline(h=limit.pval_limma, col="green")
    #   abline(h=-log10(opt$p_val_cutoff), col="red", lty="dashed")
    # dev.off()
  }
}




#####
################## CASE N: NOISeq  (AT LEAST 3 REPLICLATES)
#####
if((replicatesC >= 3) &
   (replicatesT >= 3) &
   grepl("N", opt$modules)){ 

  # Verbose
  cat(paste('\n Gene expression analysis is performed with NOISeqBIO function within NOISeq.'))

  # Calculate results
  results <- analysis_NOISeq(data    = raw_filter, 
                             num_controls  = replicatesC, 
                             num_treatmnts = replicatesT, 
                             q_value = opt$q_value, 
                             groups  = design_vector)

  # Store results
  all_data[['NOISeq']]                <- results[[1]]
  all_data_normalized[['NOISeq']]     <- results[[2]]
  all_counts_for_plotting[['NOISeq']] <- results[[3]]

  #Result Plot Visualization
  if (!is.null(all_counts_for_plotting[['NOISeq']])){
    all_FDR_names      <- c(all_FDR_names, 'adj.p')
    all_LFC_names      <- c(all_LFC_names, 'log2FC')
    all_pvalue_names   <- c(all_pvalue_names, 'prob')
    final_pvalue_names <- c(final_pvalue_names, 'pvalue_NOISeq')
    final_logFC_names  <- c(final_logFC_names, 'logFC_NOISeq')
    final_FDR_names    <- c(final_FDR_names, 'FDR_NOISeq')
    DEG_pack_columns   <- c(DEG_pack_columns, 'NOISeq_DEG')
  }
}  








############################################################
##                   FINAL RESULTS TABLE                  ##
############################################################

#### Preparing and creating final table
all_genes_df <- unite_all_list_dataframes(all_counts_for_plotting, all_FDR_names, all_LFC_names, all_pvalue_names, final_pvalue_names, final_logFC_names, final_FDR_names)
all_genes_df <- check_deg_in_pck(all_counts_for_plotting, all_data, all_genes_df, DEG_pack_columns)

DEG_counts   <- counting_trues(all_genes_df, DEG_pack_columns)
DEG_counts   <- as.data.frame(DEG_counts)
all_genes_df <- cbind(all_genes_df, DEG_counts)

final_BIG_table <- creating_final_BIG_table(all_genes_df, all_FDR_names, all_LFC_names, all_pvalue_names, final_pvalue_names, final_logFC_names, final_FDR_names, opt)

tag <- as.data.frame(tagging_genes(final_BIG_table, opt, DEG_counts, DEG_pack_columns))
colnames(tag)   <- "genes_tag"
final_BIG_table <- cbind(final_BIG_table, tag)
final_BIG_table <- adding_filtered_transcripts(raw, raw_filter, final_BIG_table)










############################################################
##                       EXPORT RESULTS                   ##
############################################################

# Export results tables
write_data_frames_list(dataframe_list = all_data, prefix = 'DEgenes_', paths = paths)
write_data_frames_list(dataframe_list = all_data_normalized, prefix = 'Normalized_counts_', paths = paths)
write_data_frames_list(dataframe_list = all_counts_for_plotting, prefix = 'allgenes_', paths = paths)

# # Plot final QC graphs after normalization
# pdf(file.path(paths$root, "boxplot_normcounts_distribution.pdf"), w=8.5, h=6) #simple boxplot to see overall normalized counts' distribution
#   boxplot(all_data_normalized[[1]])
# dev.off()

# pdf(file.path(paths$root, "boxplot_normalized_data.pdf"), w=8.5, h=6)
#   max_mean <- max(apply(all_data_normalized[[1]], MARGIN = 2, function(x) mean(x, na.rm=TRUE)))
#   boxplot(all_data_normalized[[1]],  ylim=c(0, max_mean*10), cex.lab=0.8, cex.axis=0.8, notch=TRUE, col=(c(rep("gold",replicatesC),rep("darkgreen",replicatesT))))
# dev.off()


# Expor "BIG FINAL TABLE" 
write.table(final_BIG_table, file=file.path(paths[["Common_results"]], "hunter_results_table.txt"), quote=F, col.names=T, sep="\t", row.names=F)










############################################################
##                PREVALENT RESULTS GRAPHS                ##
############################################################

if (length(all_data) > 1){
  ########### Venn diagram ##############
  all_package_results <- get_vector_names(all_data)

  # pdf(paste(file.path(paths[['Common_results']], "VennDiagram.pdf"))) 
  #   venn_plot <- venn.diagram(all_package_results, cex = 2, cat.fontface = 1, lty = 2, filename = NULL, cat.cex=1.5)
  #   grid.draw(venn_plot)
  # dev.off()

  ########################
  x_all <- calculate_intersection(all_package_results)

  write_data(x_all, file.path(paths[["Common_results"]]),"Prevalent_geneIDs.txt")
  ##############################

  raw_filter_x_all_separate_lfcs <- separate_intersection_logFCs_by_sign(all_data, raw_filter, x_all, all_LFC_names)
  write_data(raw_filter_x_all_separate_lfcs[[1]], file.path(paths[["Common_results"]]),"pos_prevalentDEGs_logFCs.txt")
  write_data(raw_filter_x_all_separate_lfcs[[2]], file.path(paths[["Common_results"]]),"neg_prevalentDEGs_logFCs.txt")
  intersection_data <- get_subset_for_fdr_df(all_data, x_all, all_FDR_names)
  # plotting_FDR_values(intersection_data, "padj_prevalent_DEGs.pdf" , opt$p_val_cutoff)

  all_fdr_data <- get_all_fdr_df(all_data, x_all, all_FDR_names)  
  # plotting_FDR_values(all_fdr_data, "padj_possible_DEGs.pdf" , opt$p_val_cutoff)

  all_fdr_counts_data <- get_all_fdr_df(all_counts_for_plotting, x_all, all_FDR_names)
  # plotting_FDR_values(all_fdr_counts_data, "padj_all_genes.pdf" , 1.00)
}

# generate_report(all_data, all_LFC_names, x_all)

barplot_df <- creating_genenumbers_barplot(raw, raw_filter, all_data, x_all)
# pdf(file=file.path(paths$root, "genenumbers.pdf"), width=7, height=1.2)
#   p <- ggplot(barplot_df, aes(cat, numbers)) + ylab("Number of genes") + xlab("") +
#             geom_bar(position="dodge", stat="identity", fill=c("#000034", "red", "orange", "blue"), show.legend=FALSE) + coord_flip() + 
#             geom_text(aes(label = numbers, y= numbers + 1500))

#   p + theme(text = element_text(face="bold", size=10))
# dev.off()


creating_top20_table(final_BIG_table)

# generate_DE_report()








############################################################
##                    GENERATE REPORT                     ##
############################################################
# message(dirname(normalizePath(paths$root, "report.html")))
outf <- paste(dirname(normalizePath(paths$root,"report.html")),"report.html",sep=.Platform$file.sep)

rmarkdown::render(file.path(main_path_script, 'templates', 'main_report.Rmd'), 
                  output_file = outf, intermediates_dir = outf)