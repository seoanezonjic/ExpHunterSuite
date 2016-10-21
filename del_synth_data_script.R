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
  make_option(c("-e", "--name_exp"), type="character", default="experiment1",
    help="Type the name of your experiment.")
)
opt <- parse_args(OptionParser(option_list=option_list))


# Create tree folder structure
#--------------------------------------------------------
#paths <- list() #empty list in which all output paths will be stored
#dir.create(opt$output_files)
#paths$root <-opt$output_files

print(opt$replicates)


DEG_assign_groupT <- 1-opt$DEG_assign_groupC

###### CREATING SIMULATION READ COUNTS ######################

tcc_object <- simulateReadCounts(Ngene = opt$number_DEGs, PDEG = opt$fraction_DEGs,
DEG.assign = c(opt$DEG_assign_groupC, DEG_assign_groupT),
DEG.foldchange = c(opt$FC_groupC, opt$FC_groupT),
replicates = c(opt$replicates, opt$replicates))
print(dim(tcc_object$count))
print(head(tcc_object$count))
print(str(tcc_object$simulation))

tcc_counts <- tcc_object$count
names(tcc_object$simulation$trueDEG) <- rownames(tcc_counts)
prediction_vector <- tcc_object$simulation$trueDEG

tcc_sim <- tcc_object$simulation$trueDEG
truedegs_vector <- tcc_sim>=1
degnames <- tcc_sim[c(truedegs_vector)]
degs <- subset(tcc_counts, rownames(tcc_counts) %in% names(degnames))
number_degs <- nrow(degs)

notdegs_vector <- tcc_sim==0
notdegnames <- tcc_sim[c(notdegs_vector)]
notdegs_df <- subset(tcc_counts, rownames(tcc_counts) %in% names(notdegnames))
notdegs_sample <- notdegs_df[sample(nrow(notdegs_df), number_degs), ]

tcc_final_counts <- rbind(degs, notdegs_sample)
#prediction_vector <- prediction_vector[nrow(tcc_final_counts)]
rows_counts <- rownames(tcc_final_counts)
print("tcc_counts")
print(head(rows_counts))
print(class(rows_counts))



# names(tcc_object$simulation$trueDEG) <- rownames(tcc_counts)
# print(head(tcc_object$simulation$trueDEG))
# prediction_vector <- tcc_object$simulation$trueDEG[1:8000]
# print(head(prediction_vector))

print("prediction_vector")
prediction_vector <- replace(prediction_vector, prediction_vector ==1, "DEG")
prediction_vector <- replace(prediction_vector, prediction_vector ==2, "DEG")
prediction_vector <- replace(prediction_vector, prediction_vector ==0, "NOTDEG")
print(head(prediction_vector))
print(length(prediction_vector))
print(class(prediction_vector))

prediction_vector <- prediction_vector[intersect(names(prediction_vector), rows_counts)] 
print("prediction_table")
print(head(prediction_vector))
print(length(prediction_vector))
print(class(prediction_vector))


prediction_vector <- as.data.frame(prediction_vector)
colnames(prediction_vector) <- "prediction"

rownames_vector <- rownames(prediction_vector)
rownames(prediction_vector) <- NULL
prediction_table <- data.frame(Row.names=rownames_vector, prediction=prediction_vector[,"prediction"])
print(head(prediction_table))
print(nrow(prediction_table))
print(class(prediction_table))


#stopifnot(1>500)
#tcc_prediction <- tcc_prediction[1:8000,]
#print(head(tcc_prediction))
# print(tcc_final_counts)
# print(nrow(tcc_final_counts))
# print(prediction_table)
# print(nrow(prediction_table))
#stopifnot(1>500)
###### Saving simulated counts and prediction vector ###
write.table(as.data.frame(tcc_final_counts), file="simulated_counts_tcc.txt", quote=F, col.names=NA, sep="\t")
write.table(prediction_table, file="prediction_vector.txt", quote=F, col.names=T, sep="\t", row.names=F)

#write.table(as.data.frame(tcc_object$count), file=file.path(paths$root, "simulated_counts_tcc.txt"), quote=F, col.names=NA, sep="\t")
#write.table(as.data.frame(tcc_object$simulation$trueDEG), file=file.path(paths$root, "prediction_vector.txt"), quote=F, col.names=NA, sep="\t")

