#' 
#' @param df 
#' @keywords installation
#' @export
acc <- function(df){
	(df$TP + df$TN) / (df$TP + df$TN + df$FP + df$FN)
}

#' 
#' @param df 
#' @keywords installation
#' @export
ppv <- function(df){
	df$TP / (df$TP + df$FP)
}

#' 
#' @param df 
#' @keywords installation
#' @export
recall <- function(df){
	df$TP / (df$TP + df$FN)
}

#' 
#' @param df 
#' @keywords installation
#' @export
spc <- function(df){
	df$TN / (df$TN + df$FP)
}

#' 
#' @param df 
#' @keywords installation
#' @export
fmeasure <- function(df){
	2 * (df$Precision * df$Recall) / (df$Precision + df$Recall)
}