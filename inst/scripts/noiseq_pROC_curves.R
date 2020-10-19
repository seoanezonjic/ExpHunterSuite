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
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(optparse))

# Load custom libraries
source(file.path(main_path_script, 'lib', 'general_functions.R'))
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
    help="Big hunter table"),
  make_option(c("-c", "--FC_groupC"), type="integer", default=3,
    help="Fold Change value groupC. Default=%default"),
  make_option(c("-o", "--output_files"), type="character", default="",
    help="Output path. Default=%default")
)
opt <- parse_args(OptionParser(option_list=option_list))

#Calculate global parameters


##########
paths <- list() #empty list in which all output paths will be stored
dir.create(opt$output_files)
paths$root <-opt$output_files


# Parse raw count data
#--------------------------------------------------
print("CORRRECTO")
big_final_table <- read.table(opt$input_file, header=T, row.names=1, sep="\t")


# subsetting big_final_table in common DEGs and Union DEGs
#--------------------------------------------------
big_final_table <-  na.omit(big_final_table)
print(head(big_final_table))
print(nrow(big_final_table))


drawing_ROC_curves <- function(big_final_table, graphname){ #prediction vector entra dentro de la tabla que sea
	pdf(graphname, height=12, width=10)
		roc(big_final_table$prediction_099, big_final_table$FDR_NOISeq_099, smooth=TRUE, plot=TRUE, col="blue", main="Statistical comparison", percent=FALSE, add=FALSE)
		roc(big_final_table$prediction_NOISeq_3rpl_sin, big_final_table$FDR_NOISeq_3repl_sin, smooth=TRUE, plot=TRUE, col="green", percent=FALSE, add=TRUE)
		roc(big_final_table$prediction_6rpl_sin, big_final_table$FDR_NOISeq_6rpl_sin, smooth=TRUE, plot=TRUE, col="red", percent=FALSE, add=TRUE)
		roc(big_final_table$prediction_001, big_final_table$FDR_NOISeq_post001, smooth=TRUE, plot=TRUE, col="orange", percent=FALSE, add=TRUE)
 		
		AUC_NOISeq_099 <- auc(big_final_table$prediction_099, big_final_table$FDR_NOISeq_099)
 		AUC_NOISeq_3rpl_sin <- auc(big_final_table$prediction_NOISeq_3rpl_sin, big_final_table$FDR_NOISeq_3repl_sin)
 		AUC_NOISeq_6repl_sin <- auc(big_final_table$prediction_6rpl_sin, big_final_table$FDR_NOISeq_6rpl_sin)
 		AUC_NOISeq_post001 <- auc(big_final_table$prediction_001, big_final_table$FDR_NOISeq_post001)

		legend("bottomright", legend=c(paste("NOISeq_099 = ",AUC_NOISeq_099), paste("NOISeq_3rpl_sin = ",AUC_NOISeq_3rpl_sin), paste("NOISeq_6rpl_sin = ",AUC_NOISeq_6repl_sin), paste("NOISeq_postfilt_0.01 = ",AUC_NOISeq_post001)), col=c("blue", "green", "red", "orange"), lwd=2)
	dev.off()
}

drawing_ROC_curves(big_final_table, "ROC_Curves_NOISeq_allstats_4FC_new.pdf")



stats <- c(TP_counter=0, TN_counter=0, FP_counter=0, FN_counter=0)
package_truefalse_counts <- list(NOISeq_prediction_6rpl_sin=stats, NOISeq_DEG_sin3repl=stats, NOISeq_DEG_099=stats, NOISeq_DEG_001=stats)
print(str(package_truefalse_counts))
#stopifnot(1>500)
 deg_in_packages <- c("NOISeq_prediction_6rpl_sin", "NOISeq_DEG_sin3repl", "NOISeq_DEG_099", "NOISeq_DEG_001")
 #deg_in_packages <- c("NOISeq_DEG_6sin", "NOISeq_DEG_3sin", "NOISeq_DEG_q099", "NOISeq_DEG_001") #cambiarlo a esto

 calculate_truefalse_counts_per_package <- function(big_final_table, deg_in_packages){
 for (i in c(1:nrow(big_final_table))){
 	prediction_label <- as.character(big_final_table[["prediction_099"]][i])
 	print(prediction_label)
 	is_a_DEG <- as.logical(big_final_table[i, deg_in_packages])
 	print(is_a_DEG)
 	#stopifnot(1>500)
 	 for (i in c(1:length(deg_in_packages))){
 	 	pck_name <- deg_in_packages[i] #Deseq2
 	 	print(pck_name)
 	 	if ((is_a_DEG[i] == TRUE) & (prediction_label == "DEG")){ #TRUE POSITIVE
 	 		counter <- package_truefalse_counts[[pck_name]][1]
 	 		package_truefalse_counts[[pck_name]][1] <- counter + 1
 	 		print(package_truefalse_counts)
 	 	}
 	 	if ((is_a_DEG[i] == FALSE) & (prediction_label == "NOTDEG")){ #TRUE NEGATIVE
		counter <- package_truefalse_counts[[pck_name]][2]
		package_truefalse_counts[[pck_name]][2] <- counter + 1
		print(package_truefalse_counts)
 	 	}
 	 	if ((is_a_DEG[i] == TRUE) & (prediction_label == "NOTDEG")){ #FALSE POSITIVE
		counter <- package_truefalse_counts[[pck_name]][3]
		package_truefalse_counts[[pck_name]][3] <- counter + 1
		print(package_truefalse_counts)
 	 	}
 	 	if ((is_a_DEG[i] == FALSE) & (prediction_label == "DEG")){  #FALSE NEGATIVE
		counter <- package_truefalse_counts[[pck_name]][4]
		package_truefalse_counts[[pck_name]][4] <- counter + 1
		print(package_truefalse_counts)
 	 	}

 	 }
 }
 return(package_truefalse_counts)
}

package_truefalse_counts <- calculate_truefalse_counts_per_package(big_final_table, deg_in_packages)
package_truefalse_counts <- as.data.frame(package_truefalse_counts)

print(str(package_truefalse_counts))
print(package_truefalse_counts)
write.table(package_truefalse_counts, file="measure_table_noiseq_TP.txt", quote=F, col.names=NA, sep="\t")



########################################### CALCULATING MEASURES ################################

###### FC in data set
set_FC <- opt$FC_groupC
FC_column <- as.data.frame(c(rep(set_FC, length(deg_in_packages))))
rownames(FC_column) <- deg_in_packages
colnames(FC_column) <- "FC"
print(FC_column)
#stopifnot(1>500)


##### mean FDR in each test
	
FDR_cols <- big_final_table[deg_in_packages]
FDR_means <- c()
pck_name <- c()
for (i in c(1:length(FDR_cols))){
	column_FDR <- as.numeric(FDR_cols[,i])
	print(head(column_FDR))
	mean_pck_FDR <- mean(column_FDR)
	print(mean_pck_FDR)
	FDR_means[i] <- mean_pck_FDR
	pck_name[i] <- deg_in_packages[i]
}
print(FDR_means)
names(FDR_means) <- pck_name
print(FDR_means)
FDR_means_df <- as.data.frame(FDR_means)
#stopifnot(1>500)

############################

recopilate_values_table <- function(package_truefalse_counts, deg_in_packages){
measure_table <- NULL
for (i in c(1:ncol(package_truefalse_counts))){
	package_column <- package_truefalse_counts[,i]
	package_name <- deg_in_packages[i]
	True_Pos <- package_column[1]
	True_Neg <- package_column[2]
	False_Pos <- package_column[3]
	False_Neg <- package_column[4]

	
	sensitivity_TPR <- calculate_sensitivity(True_Pos, False_Neg)

	specificity_NPR <- calculate_specificity(True_Neg, False_Pos)

	positive_predictive_value_PPV <- calculate_positive_predictive_value_PPV(True_Pos, False_Pos)

	negative_predictive_value_NPV <- calculate_negative_predictive_value_NPV(True_Neg, False_Neg)

	accuracy <- calculate_accuracy(True_Pos, True_Neg, False_Pos, False_Neg)
	next_measuring_df <- data.frame(Sensitivity_TPR=sensitivity_TPR, Specificity_NPR=specificity_NPR,
		PPV=positive_predictive_value_PPV, NPV=negative_predictive_value_NPV, Accuracy=accuracy)
	rownames(next_measuring_df) <- package_name	
	measure_table <- rbind(measure_table, next_measuring_df)

}
	return(measure_table)
}
measure_table <- recopilate_values_table(package_truefalse_counts, deg_in_packages)

measure_table <- cbind(measure_table, FC_column)
measure_table <- cbind(measure_table, FDR_means_df)
print(measure_table)

write.table(measure_table, file="measure_table_noiseq.txt", quote=F, col.names=NA, sep="\t")


