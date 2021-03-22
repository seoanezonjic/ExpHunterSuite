#! /usr/bin/env Rscript

#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>
#' @importFrom naivebayes naive_bayes

#############################################
### CONFIGURE 
#############################################

options(warn=1)

option_list <- list(
  optparse::make_option(c("-t", "--train"), type="character", default = NULL,
    help="[USE ONE] Train file with column 'Prediction'"),
  optparse::make_option(c("-f", "--folder"), type="character", default = NULL,
    help="[USE ONE] Folder with synthetic data and DEG analysis result"),
  optparse::make_option(c("-p", "--path"), type="character", default = NULL,
    help="[USE ONE] Path with several synthetic folders"), 
  optparse::make_option(c("-m", "--model"), type="character", default = NULL,
    help="[USE ONE] Model RData file. Avoid training process"), 
  optparse::make_option(c("-F", "--formula"), type="character", default = "PL",
    help=paste0("String which indicates if formula must include pvalue columns",
      " (P), logFC columns (L) or both (PL). Default: both")),
  optparse::make_option(c("-c", "--change"), type="character", default = "",
    help=paste0("String which indicates if specific columns must be mutated ",
      "applying: -log to pvalue columns (P), ABS to logFC columns (L) or both",
      " (PL).")),
  optparse::make_option(c("-T", "--test"), type="character", default = NULL,
    help=paste0("[OPTIONAL] Test dataset. If include 'Prediction' column stats",
                " will be calculated")),
  optparse::make_option(c("-s", "--save_session"), type="logical",
    default = FALSE, action = "store_true",
    help="Flag to activa SAVE SESSION mode"),
  optparse::make_option(c("-v", "--verbose"), type="logical",
    default = FALSE, action = "store_true",
    help="Activate verbose mode"),
  optparse::make_option(c("-o", "--outfile"), type="character",
    help=paste0("Output file basenames"))
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

if( Sys.getenv('DEGHUNTER_MODE') == 'DEVELOPMENT' ){
  # Obtain this script directory
  full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile),  
                 error=function(e) # works when using R CMD
                normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', 
                  commandArgs())], '='))[2]))
  main_path_script <- dirname(full.fpath)
  root_path <- file.path(main_path_script, '..', '..')
  # Load custom libraries
  custom_libraries <- c('simulate_treatment_control.R')
  for (lib in custom_libraries){
    source(file.path(root_path, 'R', lib))
  }
}else{
  require('ExpHunterSuite')
  root_path <- find.package('ExpHunterSuite')
}


#############################################
### FUNCTIONS 
#############################################
transformPvals <- function(df){
  pvalCols <- which(grepl("pvalue",colnames(df)))
  invisible(lapply(pvalCols,function(i){
    df[,i] <<- -log(df[,i])
  }))
  return(df)
}

transformLogFCs <- function(df){
  logfcCols <- setdiff(which(grepl("logFC",colnames(train))),
                       which(colnames(train) == "mean_logFCs"))
  invisible(lapply(logfcCols,function(i){
    df[,i] <<- abs(df[,i])
  }))
  return(df)
}

transformSet <- function(df,tr){
  if(grepl("P",tr)){
    df <- transformPvals(df)
  }
  if(grepl("L",tr)){
    df <- transformLogFCs(df)
  }
  return(df)
}


#############################################
### LOAD 
#############################################


# Load train set
if(!is.null(opt$train)){
  if(opt$verbose) message("Loading TRAINING dataset from file")
  train <- read.table(file = opt$train, sep = "\t", header = TRUE)
}else if(!is.null(opt$folder)){
  if(opt$verbose) message("Loading TRAINING dataset from folder")
  train <- transformSet(load_synth_dataset(opt$folder),opt$change)

}else if(!is.null(opt$path)){
  if(opt$verbose) message("Loading TRAINING dataset from paths set")
  target_folder = tryCatch(
    {
      normalizePath(opt$path)
    },
    warning=function(w){
      return(file.path(getwd(),opt$path))
    }
  )
  train <- as.data.frame(do.call(rbind,lapply(tail(list.dirs(target_folder, 
                                                        full.names = TRUE),-1), 
                                              function(folder){
    message(paste0("\t",folder))
    return(transformSet(load_synth_dataset(folder),opt$change))
  })))
}else if(!is.null(opt$model)){
  if(opt$verbose) message("Loading already trained model")
  load(opt$model)
}else{
  stop("No training set given. Finishing ...")
}

# Verbose point
if(opt$verbose) message("Loading TESTING dataset")

# Load test set
testLoaded <- FALSE
if(!is.null(opt$test)){
  test <- read.table(file = opt$test, sep = "\t", header = TRUE)
  test <- transformSet(test,opt$change)
  testLoaded <- TRUE
}

#############################################
### TRAIN AND TEST 
#############################################

if(is.null(opt$model)){
  # Verbose point
  if(opt$verbose) message("Creating formula")

  # Prepare formula
  target_cols <- c()
  if(grepl("L",opt$formula)){
    target_cols <- c(target_cols, setdiff(which(grepl("logFC",colnames(train))),
                                      which(colnames(train) == "mean_logFCs")))
  }
  if(grepl("P",opt$formula)){
    target_cols <- c(target_cols, which(grepl("pvalue",colnames(train))))
  }

  formula_t <-paste0("Prediction ~ ",paste(unique(colnames(train)[target_cols]),
                    collapse = " + "))
  formula_t <- as.formula(formula_t)

  # Verbose point
  if(opt$verbose) message("Training model")
  # Train NB model
  model <- naivebayes::naive_bayes(formula_t, data = train)
}

if(testLoaded){
  # Verbose point
  if(opt$verbose) message("Predicting over TESTING set")

  # Predict over dataset
  prediction <- naivebayes:::predict.naive_bayes(model,test)
  if("Prediction" %in% colnames(test)){
    tab <- table(test$Prediction, prediction, dnn = c(TRUE,FALSE))
    prediction_stats <- get_stats_from_cm(tab) 
  }  
}

#############################################
### EXPORT 
#############################################

# Verbose point
if(opt$verbose) message("Exporting results")

# Export model
if(is.null(opt$model)){
  save(model,file = paste(opt$outfile,"model.RData",sep = "_"))
}
# Export session
if(opt$save_session){
  save.image(paste(opt$outfile,"session.RData",sep = "_"))
}
if(testLoaded){
  # Export prediction
  write.table(cbind(rownames(test),as.character(prediction)),
    paste(opt$outfile,"pred",sep = "_"), 
    col.names = FALSE, row.names = FALSE, quote = FALSE)
  if("Prediction" %in% colnames(test)){
    # Export stats over dataset
    write.table(prediction_stats,paste(opt$outfile,"stats",sep = "_"),
      row.names = FALSE, quote = FALSE)
  }  
}
