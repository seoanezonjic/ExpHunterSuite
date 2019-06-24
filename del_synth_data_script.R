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
  make_option(c("-R", "--replicates"), type="integer", default=3,
    help="Replicates number of treatment group"),
  make_option(c("-n", "--number_DEGs"), type="integer", default=20000,
    help="Number of DEGs. Default=%default"),
  make_option(c("-f", "--fraction_DEGs"), type="double", default=0.2,
    help="Number of DEGs. Default=%default"),
  make_option(c("-c", "--FC_groupC"), type="integer", default=3,
    help="Fold Change value groupC. Default=%default"),
  make_option(c("-t", "--FC_groupT"), type="integer", default=3,
    help="Fold Change value groupT. Default=%default"),
  make_option(c("-a", "--DEG_assign_groupC"), type="double", default=0.5, 
    help="q value for NOISeq package. Default=%default"),
  make_option(c("-F", "--Filt_nonDEGs"), action="store_true",type="logical", default=FALSE,
    help="Flag to execute special non-DEG filtering"), 
  make_option(c("-e", "--name_exp"), type="character", default="experiment1",
    help="Type the name of your experiment."),
  make_option(c("-o", "--outfile"), type="character",
    help="Output file basenames. A *_scount_ttc and *_predv files will be generated")
)
opt <- parse_args(OptionParser(option_list=option_list))



###### CREATING SIMULATION READ COUNTS ######################

# Simulate
tcc_object <- simulateReadCounts(Ngene      = opt$number_DEGs, 
                                 PDEG       = opt$fraction_DEGs,
                                 DEG.assign = c(opt$DEG_assign_groupC, 1-opt$DEG_assign_groupC),
                                 DEG.foldchange = c(opt$FC_groupC, opt$FC_groupT),
                                 replicates = rep(opt$replicates, 2))

# Obtain necessary objects
tcc_counts <- tcc_object$count
names(tcc_object$simulation$trueDEG) <- rownames(tcc_counts)
prediction_vector <- tcc_object$simulation$trueDEG


# Special case
if(opt$Filt_nonDEGs){ 
  # print("iniciando filtrado non-degs")
  message("Stating non-degs filtering")
  tcc_sim <- tcc_object$simulation$trueDEG
  truedegs_vector <- tcc_sim>=1
  degnames <- tcc_sim[c(truedegs_vector)]
  degs <- subset(tcc_counts, rownames(tcc_counts) %in% names(degnames))
  number_degs <- nrow(degs)

  notdegs_vector <- tcc_sim==0
  notdegnames <- tcc_sim[c(notdegs_vector)]
  notdegs_df <- subset(tcc_counts, rownames(tcc_counts) %in% names(notdegnames))
  notdegs_sample <- notdegs_df[sample(nrow(notdegs_df), number_degs), ]

  tcc_counts <- rbind(degs, notdegs_sample)
  #prediction_vector <- prediction_vector[nrow(tcc_final_counts)]
  # rows_counts <- rownames(tcc_counts)
}

# Set row counts
# rows_counts <- rownames(tcc_counts)

# Prepare final prediction vector format
prediction_vector <- replace(prediction_vector, prediction_vector ==1, "DEG")
prediction_vector <- replace(prediction_vector, prediction_vector ==2, "DEG")
prediction_vector <- replace(prediction_vector, prediction_vector ==0, "NOTDEG")
# prediction_vector <- prediction_vector[intersect(names(prediction_vector), rows_counts)] 
prediction_vector <- prediction_vector[intersect(names(prediction_vector), rownames(tcc_counts))] 


# (FMJ) WTF is that shit Isabel?
prediction_vector <- as.data.frame(prediction_vector)
colnames(prediction_vector) <- "prediction"

rownames_vector <- rownames(prediction_vector)
rownames(prediction_vector) <- NULL
prediction_table <- data.frame(Row.names=rownames_vector, prediction=prediction_vector[,"prediction"])



# Store final prediction table 



# Write output
write.table(as.data.frame(tcc_counts), file=paste(opt$outfile,"_scount_ttc",sep=""), quote=F, col.names=NA, sep="\t")
write.table(prediction_table, file=paste(opt$outfile,"_predv",sep=""), quote=F, col.names=T, sep="\t", row.names=F)



