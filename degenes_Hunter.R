#! /usr/bin/env Rscript
######################################################################################################
### Isabel González Gayte, Rocío Bautista Moreno, Pedro Seoane Zonjic and M. Gonzalo Claros , 2016. ##
######################################################################################################

# this is wrapped in a tryCatch. The first expression works when source executes, the
# second expression works when R CMD does it.
full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile),  # works when using source
               error=function(e) # works when using R CMD
              normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', commandArgs())], '='))[2]))
main_path_script <- dirname(full.fpath)


#Loading libraries
suppressPackageStartupMessages(library(ggplot2)) 
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(DESeq2)) 
suppressPackageStartupMessages(library(NOISeq))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(genefilter)) 
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(plyr))


# Load custom libraries
source(file.path(main_path_script, 'lib', 'general_functions.R'))
source(file.path(main_path_script, 'lib', 'dif_expression_packages.R'))
source(file.path(main_path_script, 'lib', 'qc_and_benchmarking_functions.R'))


#############################################
### MAIN 
#############################################

### Input/Output (I/O)
#############################################

# Parse command line
#------------------------------------------------

option_list <- list(
  make_option(c("-i", "--input_file"), type="character", 
    help="Input file with read counts"),
  make_option(c("-C", "--Control_columns"), type="character", default=NULL,
    help="Control columns. Please indicate column names of control samples separated by commas"),
  make_option(c("-T", "--Treatment_columns"), type="character", default=NULL,
    help="Treatment columns. Please indicate column names of treatment samples separated by commas"),  
  make_option(c("-r", "--reads"), type="integer", default=2,
    help="Number of minimum reads required. Lesser number of reads are discarded. Default=%default"),
  make_option(c("-o", "--output_files"), type="character", default="",
    help="Output path. Default=%default"),
  make_option(c("-p", "--p_val_cutoff"), type="double", default=0.05,
    help="Adjusted p-value for the differential expression analysis. Default=%default"),
  make_option(c("-f", "--fc"), type="double", default=1.5,
    help="Fold Change value. Default=%default"),
  make_option(c("-q", "--q_value"), type="double", default=0.99,
    help="q value for NOISeq package. Default=%default"),
  make_option(c("-n", "--name_exp"), type="character", default="experiment1",
    help="Type the name of your experiment."),
  make_option(c("-m", "--modules"), type="character", default=c("DELN"), #D = DESeq2, E = edgeR, L = limma, N = noiseq
    help="Differential expression packages to able/disable (D = DESeq2, E = edgeR, L = limma, N= NOISeq).
    By default the following modules Default=%default are performed")
)
opt <- parse_args(OptionParser(option_list=option_list))

#Calculate global parameters
lfc <- calculate_lfc(opt)

##########


# Parse raw count data
#--------------------------------------------------

raw <- read.table(opt$input_file, header=T, row.names=1, sep="\t")

#replicates <- get_replic_number(raw)
ccolumns <- unlist(strsplit(opt$Control_columns, ","))
tcolumns <- unlist(strsplit(opt$Treatment_columns, ","))
raw <- raw[c(ccolumns,tcolumns)]

replicatesC <- length(ccolumns)
replicatesT <- length(tcolumns)
print(replicatesT)
print(replicatesC)
##########

# Create tree folder structure
#--------------------------------------------------------
paths <- list() #empty list in which all output paths will be stored
dir.create(opt$output_files)
paths$root <-opt$output_files

#
# Subfolders are created depending on the number of replicates available
subfolders <- c('Results_DESeq2')
if((replicatesC==2) & (replicatesT==2)) {
	subfolders <- c(subfolders, 'Results_edgeR')
  subfolders <- c(subfolders, 'Common_results')
} else if((replicatesC > 2) & (replicatesT >2)){
	subfolders <- c(subfolders, 'Results_NOISeq')
	subfolders <- c(subfolders, 'Results_edgeR')
	subfolders <- c(subfolders, 'Results_limma')
  subfolders <- c(subfolders, 'Common_results')
}

create_subfolders(subfolders, paths)


###################################################
### PREPROCESSING (Filtering and normalization) ### 
###################################################

# Filtering data
#------------------------------------------------------------

raw[is.na(raw)] <- 0  # make sure there are no missing values (NAs) in any column
keep_cmp <- rowSums(cpm(raw) > opt$reads) >=2 # two reads at least in two libraries
raw_filter <- raw[keep_cmp,] #Filter out count data frame
write.table(raw_filter, file=file.path(paths$root, "filtered_count_data.txt"), quote=F, col.names=NA, sep="\t")

# Boxplot BEFORE normalization 
#-----------------------------------------
pdf(file.path(paths$root, "boxplot_before_normalization.pdf"), w=11, h=8.5)
  boxplot(raw_filter, las = 2, ylim = c(-500, 1500))
dev.off()

raw_filter_m <- as.matrix(raw_filter)
raw_filter_m_sample <- raw_filter_m[sample(nrow(raw_filter_m), 1000), ]

pdf(file.path(paths[['Common_results']], "Heatmap_DE_raw.pdf"), height=12, width=10)
  heatmap(raw_filter_m_sample)
dev.off()


### Experimental design
#############################################
 
# Defining contrast
#---------------------------------------------------------------------
design_vector <- c(rep("C", replicatesC), rep("T", replicatesT))
index_lev <- levels(design_vector)


#############################################
### Differential expression analysis 
#############################################
all_data <- list()
all_data_normalized <- list()
all_counts_for_plotting <- list()
all_package_results <- list()
all_FDR_names <- c()
all_LFC_names <- c()
all_pvalue_names <- c()
final_logFC_names <- c()


# ANALYSIS BY DEFAULT (NO REPLICATES) -> DESeq2
#------------------------------------------------------
if ((replicatesC == 1)&(replicatesT == 1)){ 
  cat(paste('There are no replicates available. Gene expression analysis is performed with DESeq2 only.'))
}

module_selected <- grepl("D", opt$modules)

if (module_selected == TRUE){
  results <- analysis_DESeq2(raw_filter, replicatesC, replicatesT, opt, lfc, paths)
  all_data[['DESeq2']] <- results[[1]]
  all_data_normalized[['DESeq2']] <- results[[2]]
  all_counts_for_plotting[['DESeq2']] <- results[[3]]

  if (!is.null(all_counts_for_plotting[['DESeq2']])){
    all_FDR_names <- c(all_FDR_names, 'padj')
    all_LFC_names <- c(all_LFC_names, 'log2FoldChange')
    all_pvalue_names <- c(all_pvalue_names, 'pvalue')
    final_logFC_names <- c(final_logFC_names, 'logFC_DESeq2')

    pdf(file.path(paths[['Results_DESeq2']], "MAplot.pdf"), w=11, h=8.5)
    plotMA(all_counts_for_plotting[['DESeq2']]) # alpha is the adjusted p-value (FDR)
    dev.off()
  }
}


################## CUSTOMISED DIFFERENTIAL EXPRESSION ANALYSES #####################

if ((replicatesC >= 2)&(replicatesT >= 2)){ ############## 2 REPLICATES PER COMPARISON GROUP ###############
cat(paste('There are', replicatesC, 'replicates in the control and ',replicatesT,'replicates in the treatment condition.
  Gene expression analysis is performed with R-packages DESeq2 and edgeR'))

  module_selected <- grepl("E", opt$modules)

  if (module_selected == TRUE){
    results <- tryCatch(analysis_edgeR(raw_filter, replicatesC, replicatesT, opt, lfc, paths), error = handling_errors, warning = handling_errors)
    all_data[['edgeR']] <- results[[1]]
    all_data_normalized[['edgeR']] <- results[[2]]
    all_counts_for_plotting[['edgeR']] <- results[[3]]
   
    # Visualizing edgeR result (MA Plot)
    if (!is.null(all_counts_for_plotting[['edgeR']])){
      all_FDR_names <- c(all_FDR_names, 'FDR')
      all_LFC_names <- c(all_LFC_names, 'logFC')
      all_pvalue_names <- c(all_pvalue_names, 'PValue')
      final_logFC_names <- c(final_logFC_names, 'logFC_edgeR')

      pdf(file.path(paths[['Results_edgeR']], "MA_plots_edgeR.pdf"), w=11, h=8.5)
        with(all_counts_for_plotting[['edgeR']], plot(logCPM, logFC, pch=20, main="edgeR: Fold change vs abundance"))
        with(subset(all_counts_for_plotting[['edgeR']], FDR < opt$p_val_cutoff), points(logCPM, logFC, pch=20, col="red"))
        abline(h=c(-1,1), col="blue")
      dev.off()
    }
  }
}


if((replicatesC > 2)&(replicatesT > 2)){ ############## MORE THAN TWO REPLICATES PER COMPARISON GROUP ###############
  
  module_selected <- grepl("L", opt$modules)

  if (module_selected == TRUE){
    cat(paste('There are', replicatesC, 'replicates in the control condition and', replicatesT, 'replicates in the treatment condition.
    Gene expression analysis is performed with R-packages DESeq2, edgeR, limma and NOISeq'))
    results <- analysis_limma(raw_filter, replicatesC, replicatesT, opt, lfc, paths)

    all_data[['limma']] <- results[[1]]    
    all_data_normalized[['limma']] <- results[[2]]
    all_counts_for_plotting[['limma']] <- results[[3]]

    # Visualizing result (MA Plot)
    if (!is.null(all_counts_for_plotting[['limma']])){
      all_FDR_names <- c(all_FDR_names, 'adj.P.Val')
      all_LFC_names <- c(all_LFC_names, 'logFC')
      all_pvalue_names <- c(all_pvalue_names, 'P.Value')
      final_logFC_names <- c(final_logFC_names, 'logFC_limma')

      k_limma <- rownames(all_counts_for_plotting[['limma']]) %in% rownames(all_data[['limma']])
      pdf(file.path(paths[['Results_limma']], "MA_plots_Voom.pdf"), w=11, h=8.5)
        plot(x=all_counts_for_plotting[['limma']]$logFC, y=-log10(all_counts_for_plotting[['limma']]$P.Value), xlab="logFC", ylab="logOdds", col=c("blue", "red") [k_limma+1], pch=20, main= c("groupsB-groupsA"))
        abline(v= opt$lfc, col="cyan")
        limit.pval_limma <- -log10(max(all_data[['limma']]$P.Value)) 
        abline(h=limit.pval_limma, col="green")
        abline(h=-log10(opt$p_val_cutoff), col="red", lty="dashed")
      dev.off()
    }
  }

  module_selected <- grepl("N", opt$modules)

  if (module_selected == TRUE){
    results <- analysis_NOISeq(raw_filter, replicatesC, replicatesT, opt, paths)
    all_data[['NOISeq']] <- results[[1]]
    all_data_normalized[['NOISeq']] <- results[[2]]
    all_counts_for_plotting[['NOISeq']] <- results[[3]]


    #Result Plot Visualization
    if (!is.null(all_counts_for_plotting[['NOISeq']])){
      all_FDR_names <- c(all_FDR_names, 'adj.p')
      all_LFC_names <- c(all_LFC_names, 'log2FC')
      all_pvalue_names <- c(all_pvalue_names, 'P.Value')
      final_logFC_names <- c(final_logFC_names, 'logFC_NOISeq')
    }
  }
}  

write_data_frames_list(all_data, 'DEgenes_', paths)
write_data_frames_list(all_data_normalized, 'Normalized_counts', paths)


##################################################################################################
######### Quality Control - Boxplot and heatmap AFTER normalization (with edgeR TMM method) ######
##################################################################################################

pdf(file.path(paths$root, "boxplot_normalized_data.pdf"), w=11, h=8.5)
  boxplot(all_data_normalized[["edgeR"]], las = 2, ylim = c(-500, 1500))
dev.off()


## Heatmap AFTER normalization
print(str(all_data_normalized))
normalized_raw_counts_m <- as.matrix(all_data_normalized[["edgeR"]])
normalized_counts_sample <- normalized_raw_counts_m[sample(nrow(normalized_raw_counts_m), 1000), ]

pdf(file.path(paths[['Common_results']], "Heatmap_DE_norm.pdf"), height=12, width=10)
  heatmap(normalized_counts_sample)
dev.off()


################################################################
### Common Results (Venn, FDR, DEheadmap, clustering, DEGlists)
################################################################

if (length(all_data) > 1){
  ########### Venn diagram ##############
  all_package_results <- get_vector_names(all_data)

  pdf(paste(file.path(paths[['Common_results']], "Venn.pdf"))) 
    venn_plot <- venn.diagram(all_package_results, cex = 2, cat.fontface = 4, lty = 2, filename = NULL)
    grid.draw(venn_plot)
  dev.off

  x_all <- calculate_intersection(all_package_results)
  write_data(x_all, file.path(paths[["Common_results"]]),"Intersection_geneIDs.txt")
  
  raw_filter_x_all_separate_lfcs <- separate_intersection_logFCs_by_sign(all_data, raw_filter, x_all, all_LFC_names)
  write_data(raw_filter_x_all_separate_lfcs[[1]], file.path(paths[["Common_results"]]),"pos_common_lfcs.txt")
  write_data(raw_filter_x_all_separate_lfcs[[2]], file.path(paths[["Common_results"]]),"neg_common_lfcs.txt")

  intersection_data <- get_subset_for_fdr_df(all_data, x_all, all_FDR_names)
  
  pdf(file.path(paths[["Common_results"]],"ggplot_intersect.pdf"), w=11, h=8.5)
    p_seguros_Int <- ggplot(intersection_data, aes(x = package_name, y = fdr, color = package_name))
    plot(p_seguros_Int + geom_boxplot(outlier.colour = rgb(0, 0, 0, 0)) + theme_bw(base_size = 30) + geom_point(position = position_jitter(w = 0.1), color = "grey50", size = 1) + geom_hline(aes(yintercept = opt$p_val_cutoff)) + ylab("1 - precision (FDR)") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("") + scale_colour_discrete(guide = "none") + coord_cartesian(ylim = c(0, 0.05)))
  dev.off()
   
  all_fdr_data <- get_all_fdr_df(all_data, x_all, all_FDR_names)
  
  pdf(file.path(paths[["Common_results"]],"ggplot_all.pdf"), w=11, h=8.5)
    p_seguros_Int <- ggplot(all_fdr_data, aes(x = package_name, y = fdr, color = package_name))
    plot(p_seguros_Int + geom_boxplot(outlier.colour = rgb(0, 0, 0, 0)) + theme_bw(base_size = 30) + geom_point(position = position_jitter(w = 0.1), color = "grey50", size = 1) + geom_hline(aes(yintercept = opt$p_val_cutoff)) + ylab("1 - precision (FDR)") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("") + scale_colour_discrete(guide = "none") + coord_cartesian(ylim = c(0, 0.05)))
  dev.off()

  ########### Mapping genes ##############
  mapping_results(raw_filter,x_all)

  ##### Generate statistics report ######
  generate_report(all_data, all_LFC_names, x_all)

  ########### Heatmap and clustering candidates ##############
  heatmapping_and_clustering(unlist(all_data_normalized[['edgeR']]), x_all) #with package edgeR as it works with standard TMM normalization

}else{
  generate_report(all_data, all_LFC_names, x_all)
}


############ PREPARING OUTPUT FOR functional_Hunter.R ###################

complete_alldata_df <- as.data.frame(subset(all_data[[1]], select = all_LFC_names[1]))
colnames(complete_alldata_df)[colnames(complete_alldata_df)==all_LFC_names[1]] <- final_logFC_names[1]

for (i in c(2:length(all_data))){
   next_alldata_df <- as.data.frame(subset(all_data[[i]], select = all_LFC_names[i]))
   colnames(next_alldata_df)[colnames(next_alldata_df)==all_LFC_names[i]] <- final_logFC_names[i]
   complete_alldata_df <- merge(complete_alldata_df, next_alldata_df, by=0, all=TRUE)
   rownames_alldata <- complete_alldata_df[,1]
   complete_alldata_df[1] <- NULL
   rownames(complete_alldata_df) <- rownames_alldata
}

common_DEGs <- rownames(complete_alldata_df) %in% x_all
names(common_DEGs) <- rownames(complete_alldata_df)
common_DEGs_df <- as.data.frame(common_DEGs)
final_table <- merge(common_DEGs_df, complete_alldata_df, by.x="row.names", by.y="row.names")
write_data(final_table, file.path(paths[["Common_results"]]),"hunter_DEGs_file.txt")

