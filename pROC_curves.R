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
#suppressPackageStartupMessages(library(pROC))
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

big_final_table <- read.table(opt$input_file, header=T, row.names=1, sep="\t")
print(head(big_final_table))


# subsetting big_final_table in common DEGs and Union DEGs
#--------------------------------------------------
big_final_table <-  na.omit(big_final_table)
print(head(big_final_table))
print(nrow(big_final_table))
big_final_table <- big_final_table[big_final_table$is_union_genenames != "FILTERED",]
print(head(big_final_table))
print(nrow(big_final_table))
#df[df$value>3.0,]
#stopifnot(1>500)
# Subsetting common results
subset_common_DEG_df <- big_final_table[big_final_table$is_union_genenames == "COMMON_DEG",]
subset_common_REJECTED_df <- big_final_table[big_final_table$is_union_genenames == "REJECTED",]
subset_common_df <- rbind(subset_common_DEG_df, subset_common_REJECTED_df)
print(head(subset_common_df))
print(nrow(subset_common_df))

write.table(big_final_table, file="forplotting.txt", quote=F, col.names=NA, sep="\t")


sign_pval_table <- big_final_table[big_final_table$pval_labeling == "SIGN",]
print(sign_pval_table)
write.table(sign_pval_table, file="pval.txt", quote=F, col.names=NA, sep="\t")
#stopifnot(1>500)

two_DEG_df <- big_final_table[big_final_table$coincid_counter == 2,]
three_pck_DEG_df <- big_final_table[big_final_table$coincid_counter >= 3,]
allpck_DEG_df <- big_final_table[big_final_table$coincid_counter == 4,]


drawing_ROC_curves <- function(big_final_table, sign_pval_table, graphname){ #prediction vector entra dentro de la tabla que sea
	pdf(graphname, height=12, width=10)
		roc(big_final_table$prediction, big_final_table$FDR_DESeq2, smooth=TRUE, plot=TRUE, col="blue", main="Statistical comparison", percent=FALSE, add=FALSE, lwd=1)
		roc(big_final_table$prediction, big_final_table$FDR_edgeR, smooth=TRUE, plot=TRUE, col="green", percent=FALSE, add=TRUE, lwd=1)
		roc(big_final_table$prediction, big_final_table$FDR_limma, smooth=TRUE, plot=TRUE, col="red", percent=FALSE, add=TRUE, lwd=1)
		roc(big_final_table$prediction, big_final_table$FDR_NOISeq, smooth=TRUE, plot=TRUE, col="orange", percent=FALSE, add=TRUE, lwd=1)
		roc(big_final_table$prediction, big_final_table$combined_pvalues, smooth=TRUE, plot=TRUE, col="black", percent=FALSE, add=TRUE, lwd=2)		
 		roc(sign_pval_table$prediction, sign_pval_table$combined_pvalues, smooth=TRUE, plot=TRUE, col="violet", percent=FALSE, add=TRUE, lwd=2)

		AUC_DESeq2 <- round(auc(big_final_table$prediction, big_final_table$FDR_DESeq2),3)
 		AUC_edgeR <- round(auc(big_final_table$prediction, big_final_table$FDR_edgeR),3)
 		AUC_limma <- round(auc(big_final_table$prediction, big_final_table$FDR_limma),3)
 		AUC_NOISeq <- round(auc(big_final_table$prediction, big_final_table$FDR_NOISeq),3)
 		AUC_pval_comb <- round(auc(big_final_table$prediction, big_final_table$combined_pvalues),3)
 		AUC_pval_sign <- round(auc(sign_pval_table$prediction, sign_pval_table$combined_pvalues),3)
		#text(50, 50, labels=paste("p-value =", format.pval(big_final_table$prediction)), adj=c(0, .5))
 		#legend("right", legend=c(AUC_DESeq2, AUC_edgeR, AUC_limma, AUC_NOISeq), col=c("blue", "green", "red", "orange"))
		#text(10, 10, labels=paste("\nAUC =", auc(big_final_table$prediction, big_final_table$FDR_DESeq2)))
		legend("bottomright", legend=c(paste("DESeq2 AUC = ", AUC_DESeq2), 
			paste("edgeR AUC = ",AUC_edgeR), paste("limma AUC = ",AUC_limma), 
			paste("NOISeq AUC = ",AUC_NOISeq), paste("Comb pVal AUC = ",AUC_pval_comb),
			paste("Sign pVal AUC = ",AUC_pval_sign)), 
			col=c("blue", "green", "red", "orange", "black", "violet"), lwd=2, cex=1)
		#legend("bottomright", legend=c("DESeq2", "edgeR", "limma", "NOISeq"), col=c("blue", "green", "red", "orange"), lwd=2)
	dev.off()
}


#drawing_ROC_curves(big_final_table, sign_pval_table, "ROC_Curve_all_2.pdf")
#drawing_ROC_curves(subset_common_df, sign_pval_table, "ROC_Curve_common_2.pdf")


drawing_3and4pack_curves <- function(big_final_table, sign_pval_table, graphname){ #prediction vector entra dentro de la tabla que sea
	pdf(graphname, height=12, width=10)
		roc(three_pck_DEG_df$prediction, three_pck_DEG_df$combined_pvalues, smooth=TRUE, plot=TRUE, main="Statistical comparison", col="blue", percent=FALSE, add=FALSE, lwd=1)
		roc(two_DEG_df$prediction, two_DEG_df$combined_pvalues, smooth=TRUE, plot=TRUE, col="red", percent=FALSE, add=TRUE, lwd=2)
		roc(allpck_DEG_df$prediction, allpck_DEG_df$combined_pvalues, smooth=TRUE, plot=TRUE, col="black", percent=FALSE, add=TRUE, lwd=2)		
 		roc(sign_pval_table$prediction, sign_pval_table$combined_pvalues, smooth=TRUE, plot=TRUE, col="violet", percent=FALSE, add=TRUE, lwd=2)

 		AUC_three_pck <- round(auc(three_pck_DEG_df$prediction, three_pck_DEG_df$combined_pvalues),3)
 		AUC_allpck_comb <- round(auc(allpck_DEG_df$prediction, allpck_DEG_df$combined_pvalues),3)
 		AUC_pval_sign <- round(auc(sign_pval_table$prediction, sign_pval_table$combined_pvalues),3)
 		AUC_two_pck <- round(auc(two_DEG_df$prediction, two_DEG_df$combined_pvalues),3)

		legend("bottomright", legend=c(paste("Three package AUC = ", AUC_three_pck), 
			paste("all pck AUC = ",AUC_allpck_comb), paste("two pck AUC = ",AUC_two_pck),
			paste("Sign pVal AUC = ",AUC_pval_sign)), 
			col=c("blue", "black", "red", "violet"), lwd=2, cex=1)
	dev.off()
}


#drawing_3and4pack_curves(big_final_table, sign_pval_table, "ROC_Curve_imp.pdf")






write.table(big_final_table, file="prediction_all.txt", quote=F, col.names=NA, sep="\t")
#write.table(subset_common_df, file="prediction_common.txt"), quote=F, col.names=NA, sep="\t")




############# Create TPR table ################
 print(head(big_final_table))
 #stopifnot(1>500)


stats <- c(TP_counter=0, TN_counter=0, FP_counter=0, FN_counter=0)
package_truefalse_counts <- list(DESeq2_DEG=stats, edgeR_DEG=stats, limma_DEG=stats, NOISeq_DEG=stats)
print(str(package_truefalse_counts))
#stopifnot(1>500)
 deg_in_packages <- c("DESeq2_DEG", "edgeR_DEG", "limma_DEG", "NOISeq_DEG")

 calculate_truefalse_counts_per_package <- function(big_final_table, deg_in_packages){
 for (i in c(1:nrow(big_final_table))){
 	prediction_label <- as.character(big_final_table[["prediction"]][i])
 	print(prediction_label)
 	is_a_DEG <- as.logical(big_final_table[i, deg_in_packages])
 	print(is_a_DEG)
 	 for (i in c(1:length(deg_in_packages))){
 	 	pck_name <- deg_in_packages[i] #Deseq2
 	 	print(pck_name)
 	 	if ((is_a_DEG[i] == TRUE) & (prediction_label == "DEG")){
 	 		counter <- package_truefalse_counts[[pck_name]][1]
 	 		package_truefalse_counts[[pck_name]][1] <- counter + 1
 	 		print(package_truefalse_counts)
 	 	}
 	 	if ((is_a_DEG[i] == FALSE) & (prediction_label == "NOTDEG")){
		counter <- package_truefalse_counts[[pck_name]][2]
		package_truefalse_counts[[pck_name]][2] <- counter + 1
		print(package_truefalse_counts)
 	 	}
 	 	if ((is_a_DEG[i] == TRUE) & (prediction_label == "NOTDEG")){
		counter <- package_truefalse_counts[[pck_name]][3]
		package_truefalse_counts[[pck_name]][3] <- counter + 1
		print(package_truefalse_counts)
 	 	}
 	 	if ((is_a_DEG[i] == FALSE) & (prediction_label == "DEG")){
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


###### Mean logFC per package



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

write.table(measure_table, file="measure_table.txt", quote=F, col.names=NA, sep="\t")

