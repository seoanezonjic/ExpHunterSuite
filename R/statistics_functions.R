acc <- function(df){
	(df$TP + df$TN) / (df$TP + df$TN + df$FP + df$FN)
}

# PPV
ppv <- function(df){
	df$TP / (df$TP + df$FP)
}

# Recall
recall <- function(df){
	df$TP / (df$TP + df$FN)
}

# SPC
spc <- function(df){
	df$TN / (df$TN + df$FP)
}

# F-measure
fmeasure <- function(df){
	2 * (df$Precision * df$Recall) / (df$Precision + df$Recall)
}