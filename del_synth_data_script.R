#! /usr/bin/env Rscript
#############################################
### Synthetic Data Generator 
#############################################

# this is wrapped in a tryCatch. The first expression works when source executes, the
# second expression works when R CMD does it.
full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile),  # works when using source
               error=function(e) # works when using R CMD
              normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', commandArgs())], '='))[2]))
main_path_script <- dirname(full.fpath)

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
  make_option(c("-a", "--DEG_assign_groupC"), type="double", default=NULL, #0.9
    help="q value for NOISeq package. Default=%default"),
  make_option(c("-F", "--Filt_nonDEGs"), type="logical", default=FALSE,
    help="If TRUE, special non-DEG filtering is performed"), 
  make_option(c("-e", "--name_exp"), type="character", default="experiment1",
    help="Type the name of your experiment.")
)
opt <- parse_args(OptionParser(option_list=option_list))


# Create tree folder structure
#--------------------------------------------------------

print(opt$replicates)


DEG_assign_groupT <- 1-opt$DEG_assign_groupC

###### CREATING SIMULATION READ COUNTS ######################

tcc_object <- simulateReadCounts(Ngene = opt$number_DEGs, PDEG = opt$fraction_DEGs,
DEG.assign = c(opt$DEG_assign_groupC, DEG_assign_groupT),
DEG.foldchange = c(opt$FC_groupC, opt$FC_groupT),
replicates = c(opt$replicates, opt$replicates))

tcc_counts <- tcc_object$count
names(tcc_object$simulation$trueDEG) <- rownames(tcc_counts)
prediction_vector <- tcc_object$simulation$trueDEG




if (opt$Filt_nonDEGs==TRUE){ 
  print("iniciando filtrado non-degs")
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
  rows_counts <- rownames(tcc_counts)
}



rows_counts <- rownames(tcc_counts)


prediction_vector <- replace(prediction_vector, prediction_vector ==1, "DEG")
prediction_vector <- replace(prediction_vector, prediction_vector ==2, "DEG")
prediction_vector <- replace(prediction_vector, prediction_vector ==0, "NOTDEG")

prediction_vector <- prediction_vector[intersect(names(prediction_vector), rows_counts)] 



prediction_vector <- as.data.frame(prediction_vector)
colnames(prediction_vector) <- "prediction"

rownames_vector <- rownames(prediction_vector)
rownames(prediction_vector) <- NULL
prediction_table <- data.frame(Row.names=rownames_vector, prediction=prediction_vector[,"prediction"])


write.table(as.data.frame(tcc_counts), file="simulated_counts_tcc.txt", quote=F, col.names=NA, sep="\t")
write.table(prediction_table, file="prediction_vector.txt", quote=F, col.names=T, sep="\t", row.names=F)



