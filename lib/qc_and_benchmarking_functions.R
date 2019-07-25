#############################################
### QC AND BENCHMARKING FUNCTIONS 
#############################################

calculate_sensitivity <- function(True_Pos, False_Neg){
  sensitivity <- True_Pos/(True_Pos + False_Neg)
  return(sensitivity)
}

calculate_specificity <- function(True_Neg, False_Pos){
  specificity <- True_Neg/(False_Pos + True_Neg)
  return(specificity)
} 

calculate_positive_predictive_value_PPV <- function(True_Pos, False_Pos){
  positive_predictive_value_PPV <- True_Pos/(True_Pos + False_Pos)
  return(positive_predictive_value_PPV)
}

calculate_negative_predictive_value_NPV <- function(True_Neg, False_Neg){
  negative_predictive_value_NPV <- True_Neg/(False_Neg + True_Neg)
  return(negative_predictive_value_NPV)
}

calculate_accuracy <- function(True_Pos, True_Neg, False_Pos, False_Neg){
  negative_predictive_value_NPV <- True_Neg/(False_Neg + True_Neg)
  return(negative_predictive_value_NPV)
}
