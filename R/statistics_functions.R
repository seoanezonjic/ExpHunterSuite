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
#' df$ppv <- ppv(df)
#' df$recall <- recall(df)
#' fmeasure(df)
fmeasure <- function(df){
    2 * (df$ppv * df$recall) / (df$ppv + df$recall)
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

#' @importFrom stats sd 
permutations_stats <- function(
permutations = 10,
experiment, 
background, 
sample_size){
  random_dist <- c()
  for(i in seq_len(permutations)){
    sampled_bckg <- rand_sample_bool(background, 
        sample_size = sample_size)
    random_dist <- c(random_dist, sum(sampled_bckg &  
                                      experiment))
  }
  return(data.frame(sdev = stats::sd(random_dist),
                    mean = mean(random_dist)))
}  

#Calc contingency table from 2 boolean vectors
calc_ct <- function(experiment, gold_standard){
  TP <- sum(experiment & gold_standard)
  FP <- sum(experiment & !gold_standard)
  TN <- sum(!experiment & !gold_standard)
  FN <- sum(!experiment & gold_standard)

  cont_table <- data.frame(TP = TP, 
                        FN = FN, 
                        FP = FP, 
                        TN = TN)

  return(cont_table)
}

LRplus_sub <- function(df) {
  (df$TP / (df$TP + df$FP)) / ( df$FN / ( df$FN + df$TN ))
}

LRplus_test <- function(df) {
  ( df$TP / ( df$TP + df$FN ) ) /  ( df$FP / ( df$FP + df$TN ))
}

#' @importFrom stats pchisq
vectorial_fisher_method <- function(pval_table){
    log_pval <- log(pval_table) # Log all final p-values
    log_pval[is.na(log_pval)] <- 0 # any NAs made 0s -> not used combined score 
    if (is.vector(log_pval)){
      log_pval <- t(matrix(log_pval))
      log_pval[log_pval == -Inf] <- 10e-200
    } else {
       log_pval <- apply(as.matrix(log_pval), 2, 
         function (pval_column) { # When p values of 0 are given
                                  #- these become -Inf when logged:   
                                  # This code turn -Inf values 
                                  #- to smallest possible values
           min_col_val <- sort(unique(pval_column))[2]
           pval_column[pval_column == -Inf] <- min_col_val
           return(pval_column)
       })
    }
    xi_squared <-  -2 * rowSums(log_pval)
    degrees_freedom <- 2 * ncol(log_pval)
    combined_FDR <- stats::pchisq(xi_squared, degrees_freedom, 
                                lower.tail = FALSE)
    return(combined_FDR)
}


#' Function to calculate full set of measures from a confusion data frame
#' @param df confusion dataframe 
#' @return a data frame with several measures
#' @importFrom data.table setnames
#' @export
#' @examples
#' df <- data.frame(TP = 4, FP = 3, TN = 2, FN = 1)
#' stats <- get_stats(df)
get_stats <- function(df){
  df$acc <- acc(df)
  df$spc <- spc(df)
  df$ppv <- ppv(df)
  df$recall <- recall(df)
  df$fmeasure <- fmeasure(df)
  
  data.table::setnames(
     df,
     c("acc","ppv",      "recall","spc",         "fmeasure"),
     c("ACC","Precision","Recall","Specificity" ,"FMeasure"))

  return(df)
}

#' Obtain full set of measures from a confusion matrix
#' @param cm confusion matrix 
#' @return a data frame with several measures
#' @export
#' @examples
#' cm <- matrix(c(4,3,2,1), ncol = 2)
#' rownames(cm) <- c(TRUE,FALSE)
#' colnames(cm) <- c(TRUE,FALSE)
#' stats <- get_stats_from_cm(cm)
get_stats_from_cm <- function(cm){
  # >>>>>>>[ Real , Pred ]
  TP <- extract_from_cm(cm,"TRUE","TRUE")
  FN <- extract_from_cm(cm,"TRUE","FALSE")
  TN <- extract_from_cm(cm,"FALSE","FALSE")
  FP <- extract_from_cm(cm,"FALSE","TRUE")
  df <- data.frame(TP = TP, FP = FP, TN = TN, FN = FN)
  df <- get_stats(df)
  return(df)
}

#' Extract value from confusion matrix
#' @param cm confusion matrix
#' @param real slot name
#' @param pred slot name
#' @param default to return if value is not available
#' @return value required or default if it is not available
extract_from_cm <- function(cm,real,pred,default = 0){
  value = tryCatch(
    {
      return(cm[real,pred])
    },
    error = function(cond){
      return(0)
    })
  return(value)
}


#' @importFrom stats phyper 
v.fisher.test <- function(df){
  stats::phyper(
    q = df$TP -1, 
    m = df$TP+ df$FN, 
    n = rowSums(df[,c("TP","FP","TN","FN")])-df$TP-df$FN,
    k = df$TP+df$FP, 
    lower.tail = FALSE)
}

#' @importFrom data.table setnames
v.get_stats <- function(
df, 
selected_stats = c("acc", "ppv", "recall", "spc", "fmeasure", 
                 "LRplus_sub", "LRplus_test", "v.fisher.test")
){
  stats_names <- data.frame(orig = c("acc", "ppv", "recall", "spc", "fmeasure", 
                                 "LRplus_sub", "LRplus_test", "v.fisher.test"),
                            renamed = c("Accuracy", "Precision", "Recall", 
                                 "Specificity", "Fmax", "LRplus_sub", 
                                 "LRplus_test", "Pvalue"))
  stats_names <- stats_names[stats_names$orig %in% selected_stats, ]
  for (stat in stats_names$orig){
    df[,stat] <- get(stat)(df)
  }
  data.table::setnames(df, stats_names$orig,stats_names$renamed)
  return(df)
}