#! /usr/bin/env Rscript

#' @author refactor by Fernando Moreno Jabato. Original authors Isabel Gonzalez Gayte

#############################################
### Synthetic Data Generator 
#############################################

#Loading libraries
suppressPackageStartupMessages(library(optparse)) 
suppressPackageStartupMessages(library(TCC))


# Parse command line
#------------------------------------------------

option_list <- list(
  make_option(c("-r", "--replicates"), type="integer", default=3,
    help="Number of replicates for control/treatment group. Default : %default"),
  make_option(c("-n", "--ngenes"), type="integer", default=20000,
    help="Number of genes. Default : %default"),
  make_option(c("-d", "--DEGs_proportion"), type="double", default=0.2,
    help="Proportion of differentially expressed genes (DEGs). Default : %default"),
  make_option(c("-c", "--FC_control"), type="integer", default=1,
    help="Fold Change value for control group. Note: values different of 1 provokes that this group is not a control. Default : %default"),
  make_option(c("-t", "--FC_treat"), type="integer", default=3,
    help="Fold Change value for treatment group. Default : %default"),
  make_option(c("-T", "--P_treat"), type="double", default=1,
    help="Proportion of DEGs that will be up-regulated in treatment group. Rest will be up-regulated at control, making control not to be a real control if !={1,-1}. Default : %default"),
#   make_option(c("-F", "--Filt_nonDEGs"), action="store_true",type="logical", default=FALSE,
#     help="Flag to execute special non-DEG filtering"), 
#   make_option(c("-e", "--name_exp"), type="character", default="experiment1",
#     help="Type the name of your experiment."),
  make_option(c("-o", "--outfile"), type="character",
    help="Output file basenames. A *_scount_ttc and *_predv files will be generated")
)



opt <- parse_args(OptionParser(option_list=option_list))

###### CREATING SIMULATION READ COUNTS ######################

########### Verbose point
tcc_object <- simulateReadCounts(Ngene      = opt$ngenes, 
                                 PDEG       = opt$DEGs_proportion,
                                 DEG.assign = c(1-opt$P_treat, opt$P_treat),
                                 DEG.foldchange = c(opt$FC_control, opt$FC_treat),
                                 replicates = rep(opt$replicates, 2))

# Obtain necessary objects
colnames(tcc_object$count) <- gsub("G2","Treat",gsub("G1","Ctrl",colnames(tcc_object$count)))
names(tcc_object$simulation$trueDEG) <- rownames(tcc_object$count)
prediction_vector <- tcc_object$simulation$trueDEG

# Prepare final prediction vector format
prediction_vector <- replace(prediction_vector, prediction_vector >= 1, "DEG")
prediction_vector <- replace(prediction_vector, prediction_vector == 0, "NOTDEG")
prediction_vector <- prediction_vector[intersect(names(prediction_vector), rownames(tcc_object$count))] 
prediction_vector <- as.data.frame(prediction_vector)
prediction_vector <- cbind(rownames(prediction_vector),prediction_vector)
colnames(prediction_vector) <- c("Gene","Prediction")

# Write output
write.table(tcc_object$count, file=paste0(opt$outfile,"_scount_ttc"), quote = FALSE, col.names = TRUE, sep = "\t")
write.table(prediction_vector, file=paste0(opt$outfile,"_predv"), quote = FALSE, col.names = TRUE, sep = "\t", row.names = FALSE)



