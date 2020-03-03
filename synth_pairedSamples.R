#! /usr/bin/env Rscript

#' @author Fernando Moreno Jabato <jabato(at)uma(dot)com>

#############################################
### CONFIGURE 
#############################################

#Loading libraries
suppressPackageStartupMessages(library(optparse)) 
suppressPackageStartupMessages(library(ssizeRNA))


option_list <- list(
  make_option(c("-r", "--replicates"), type="integer", default=3,
    help="Number of replicates for control/treatment group. Default : %default"),
  make_option(c("-n", "--ngenes"), type="integer", default=20000,
    help="Number of genes. Default : %default"),
  make_option(c("-d", "--DEGs_proportion"), type="double", default=0.2,
    help="Proportion of differentially expressed genes (DEGs). Default : %default"),
  make_option(c("-t", "--FC_treat"), type="integer", default=3,
    help="Fold Change value for treatment group. Default : %default"),
  make_option(c("-T", "--P_treat"), type="double", default=1,
    help="Proportion of DEGs that will be up-regulated in treatment group. Rest will be down-regulated. Default : %default"),
#   make_option(c("-F", "--Filt_nonDEGs"), action="store_true",type="logical", default=FALSE,
#     help="Flag to execute special non-DEG filtering"), 
#   make_option(c("-e", "--name_exp"), type="character", default="experiment1",
#     help="Type the name of your experiment."),
  make_option(c("-o", "--outfile"), type="character",
    help="Output file basenames. A *_scount_ttc and *_predv files will be generated")
)


opt <- parse_args(OptionParser(option_list=option_list))

#############################################
### SIMULATE
#############################################

simu_set <- sim.counts(m    = opt$replicates, # Number of biological samples
                       pi0  = opt$DEGs_proportion, # Diff expressed genes amount
                       nGenes = opt$ngenes, # Number of genes
                       mu   = 10, # Mean counts in control group 
                       disp = 0.1, # Dispersion for all genes
                       up = opt$P_treat, # Proportion of up-regulated genes 
                       fc   = opt$FC_treat) # Fold-change of treatment set

# Prepare set
aux <- ncol(simu_set$counts)/2
colnames(simu_set$counts) <- unlist(lapply(seq(aux*2),function(i){ifelse(i <= aux,paste0("Ctrl_",i),paste0("Treat_",i-aux))}))
rownames(simu_set$counts) <- paste0("gene_",seq(nrow(simu_set$counts)))

# Prepare final predition vector format
prediction_vector <- simu_set$de
prediction_vector <- replace(prediction_vector, prediction_vector != 0, "DEG")
prediction_vector <- replace(prediction_vector, prediction_vector == 0, "NOTDEG")
prediction_vector <- as.data.frame(prediction_vector)
rownames(prediction_vector) <- rownames(simu_set$counts)
prediction_vector <- cbind(rownames(prediction_vector),prediction_vector)
colnames(prediction_vector) <- c("Gene","Prediction")


#############################################
### EXPORT
#############################################
write.table(simu_set$count, file=paste0(opt$outfile,"_scount_ttc"), quote = FALSE, col.names = TRUE, sep = "\t")
write.table(prediction_vector, file=paste0(opt$outfile,"_predv"), quote = FALSE, col.names = TRUE, sep = "\t", row.names = FALSE)
