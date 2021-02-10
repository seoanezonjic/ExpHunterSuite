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


#' @importFrom MKinfer boot.t.test
prediction_dist_pval <- function(db_distribution) {
  dist_pval <- data.frame()
  
  for (strategy in unique(db_distribution$strategy)) {
    true_dis <- db_distribution[db_distribution$strategy == strategy & db_distribution$step == "predicted", "score"]
    rand_dis <- db_distribution[db_distribution$strategy == strategy & db_distribution$step == "predicted_random", "score"] 
    if(length(true_dis) > 0 && sum(true_dis) > 0){
      res <- MKinfer::boot.t.test(true_dis, y = rand_dis, alternative = "greater" )
      res <- res$boot.p.value
    } else {
      res <- 1
    }
    tem_pval <- data.frame(strategy = strategy, dist_pvalue = res)
    dist_pval <- rbind(dist_pval, tem_pval)
   }
  return(dist_pval)
}