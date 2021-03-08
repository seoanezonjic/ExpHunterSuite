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
    true_dis <- db_distribution[db_distribution$strategy == strategy & 
        db_distribution$step == "predicted", "score"]
    rand_dis <- db_distribution[db_distribution$strategy == strategy & 
        db_distribution$step == "predicted_random", "score"] 
    if(length(true_dis) > 0 && sum(true_dis) > 0){
      res <- MKinfer::boot.t.test(true_dis, 
        y = rand_dis, alternative = "greater" )
      res <- res$boot.p.value
    } else {
      res <- 1
    }
    tem_pval <- data.frame(strategy = strategy, dist_pvalue = res)
    dist_pval <- rbind(dist_pval, tem_pval)
   }
  return(dist_pval)
}



#' @importFrom stats fisher.test
get_strategies_stats <- function(data , input_cols, reference_cols) {
  pval_table <- matrix(NA, ncol = length(reference_cols), 
    nrow = length(input_cols), dimnames = list(input_cols,reference_cols))
  LR_test_matrix <- matrix(NA, ncol = length(reference_cols), 
    nrow = length(input_cols), dimnames = list(input_cols,reference_cols))
  LR_sub_matrix  <- matrix(NA, ncol = length(reference_cols), 
    nrow = length(input_cols), dimnames = list(input_cols,reference_cols))
  for (icol_name in input_cols){
    for (rcol_name in reference_cols){
      icol <- data[,icol_name]
      rcol <- data[,rcol_name]
      c_matrix <- calc_contingency_matrix(icol, rcol)
      print(paste0(icol_name, " ", rcol_name))
      print(c_matrix)
      ftest <- stats::fisher.test(c_matrix, alternative = "greater")
      pval_table[icol_name,rcol_name] <- ftest$p.value
      LR_test_matrix[icol_name,rcol_name] <- calc_LRplus_test(c_matrix)
      LR_sub_matrix[icol_name,rcol_name] <- calc_LRplus_subject(c_matrix)
      # print(c_matrix)
      # message(ftest$p.value)
      LR_test_matrix[is.na(LR_test_matrix)] <- 0
      LR_sub_matrix[is.na(LR_test_matrix)] <- 0

    }
  }
  return(list(pval_table = pval_table, LR_test_matrix = LR_test_matrix, 
    LR_sub_matrix = LR_sub_matrix))
}

calc_contingency_matrix <- function(experiment, gold_standard){
  tpositives <- sum(experiment & gold_standard)
  fpositives <- sum(experiment & !gold_standard)
  tnegatives <- sum(!experiment & !gold_standard)
  fnegatives <- sum(!experiment & gold_standard)

  c_matrix <- matrix(c(tpositives, fnegatives, fpositives, tnegatives), 
    nrow = 2 , dimnames = list(experiment =c(TRUE,FALSE), 
      gold_standard = c(TRUE,FALSE)))
  return(c_matrix)
}

calc_LRplus_subject <- function(c_matrix) {
  tpositives<- c_matrix[1,1]
  fpositives <- c_matrix[1,2]
  tnegatives <- c_matrix[2,2]
  fnegatives <- c_matrix[2,1]
  LR_plus <- ( tpositives / ( tpositives + fpositives ) ) / 
             ( fnegatives / ( fnegatives + tnegatives ))
  return(LR_plus)
}

calc_LRplus_test <- function(c_matrix) {
  tpositives<- c_matrix[1,1]
  fpositives <- c_matrix[1,2]
  tnegatives <- c_matrix[2,2]
  fnegatives <- c_matrix[2,1]
  LR_plus <- ( tpositives / ( tpositives + fnegatives ) ) / 
             ( fpositives / ( fpositives + tnegatives ))
  return(LR_plus)
}

#' @importFrom stats pchisp
vectorial_fisher_method <- function(pval_table){
    log_pval <- log(pval_table) # Log all final p-values
    log_pval[is.na(log_pval)] <- 0 # any NAs made 0s -> not used combined score 
    
    log_pval <- apply(as.matrix(log_pval), 2, 
      function (pval_column) { # When p values of 0 are given
                               #- these become -Inf when logged:   
                               # This code turn -Inf values 
                               #- to smallest possible values
        min_col_val <- sort(unique(pval_column))[2]
        pval_column[pval_column == -Inf] <- min_col_val
        return(pval_column)
    })

    xi_squared <-  -2 * rowSums(log_pval)
    degrees_freedom <- 2 * ncol(log_pval)
    combined_FDR <- stats::pchisq(xi_squared, degrees_freedom, 
                                lower.tail = FALSE)
    return(combined_FDR)
}


