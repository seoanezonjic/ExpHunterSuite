#' Calculates accuracy from a measures dataframe
#' @param df dataframe with measures 
#' @return calculated metric
#' @keywords metrics
#' @export
#' @examples
#' df <- data.frame(TP = 4, FP = 3, TN = 2, FN = 1)
#' acc(df)
acc <- function(df){
	(df$TP + df$TN) / (df$TP + df$TN + df$FP + df$FN)
}

#' Calculates Positive Predictive Value from a measures dataframe
#' @param df dataframe with measures 
#' @return calculated metric
#' @keywords metrics
#' @export
#' @examples
#' df <- data.frame(TP = 4, FP = 3, TN = 2, FN = 1)
#' df$Precision <- ppv(df)
ppv <- function(df){
	df$TP / (df$TP + df$FP)
}

#' Calculates Recall from a measures dataframe
#' @param df dataframe with measures 
#' @return calculated metric
#' @keywords metrics
#' @export
#' @examples
#' df <- data.frame(TP = 4, FP = 3, TN = 2, FN = 1)
#' df$Recall <- recall(df)
recall <- function(df){
	df$TP / (df$TP + df$FN)
}

#' Calculates specificity from a measures dataframe
#' @param df dataframe with measures 
#' @return calculated metric
#' @keywords metrics
#' @export
#' @examples
#' df <- data.frame(TP = 4, FP = 3, TN = 2, FN = 1)
#' spc(df)
spc <- function(df){
	df$TN / (df$TN + df$FP)
}

#' Calculates F measure from a measures dataframe
#' @param df dataframe with measures 
#' @return calculated metric
#' @keywords metrics
#' @export
#' @examples
#' df <- data.frame(TP = 4, FP = 3, TN = 2, FN = 1)
#' df$Precision <- ppv(df)
#' df$Recall <- recall(df)
#' fmeasure(df)
fmeasure <- function(df){
	2 * (df$Precision * df$Recall) / (df$Precision + df$Recall)
}