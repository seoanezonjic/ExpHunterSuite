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
  optparse::make_option(c("-F", "--formula"), type="character", default = "PL",
    help=paste0("String which indicates if formula must include pvalue columns",
      " (P), logFC columns (L) or both (PL). Default: both")),
  optparse::make_option(c("-T", "--test"), type="character",
    help=paste0("Test dataset. If include 'Prediction' column stats will ",
                "be calculated")),
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
### LOAD 
#############################################


# Load train set
if(!is.null(opt$train)){
  if(opt$verbose) message("Loading TRAINING dataset from file")
  train <- read.table(file = opt$train, sep = "\t", header = TRUE)
}else if(!is.null(opt$folder)){
  if(opt$verbose) message("Loading TRAINING dataset from folder")
  train <- load_synth_dataset(opt$folder)
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
    return(load_synth_dataset(folder))
  })))
}else{
  stop("No training set given. Finishing ...")
}

# Verbose point
if(opt$verbose) message("Loading TESTING dataset")

# Load test set
test <- read.table(file = opt$test, sep = "\t", header = TRUE)

#############################################
### TRAIN AND TEST 
#############################################

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

formula <- paste0("Prediction ~ ",paste(unique(colnames(train)[target_cols]),
                  collapse = " + "))
formula <- as.formula(formula)

# Verbose point
if(opt$verbose) message("Training model")

# Train NB model
model <- naivebayes::naive_bayes(formula, data = train)

# Verbose point
if(opt$verbose) message("Predicting over TESTING set")

# Predict over dataset
prediction <- predict(model, test)
if("Prediction" %in% colnames(test)){
  tab <- table(test$Prediction, prediction, dnn = c(TRUE,FALSE))
  prediction_stats <- get_stats_from_cm(tab) 
}

#############################################
### EXPORT 
#############################################

# Verbose point
if(opt$verbose) message("Exporting results")

# Export prediction
write.table(prediction,paste(opt$outfile,"pred",sep = "_"), 
  col.names = FALSE, row.names = FALSE, quote = FALSE)
# Export model
save(model,file = paste(opt$outfile,"model.RData",sep = "_"))
# Export session
if(opt$save_session){
  save.image(paste(opt$outfile,"session.RData",sep = "_"))
}
if("Prediction" %in% colnames(test)){
  # Export stats over dataset
  write.table(prediction_stats,paste(opt$outfile,"stats",sep = "_"),
    row.names = FALSE, quote = FALSE)
}
