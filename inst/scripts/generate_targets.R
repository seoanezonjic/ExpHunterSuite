#!/usr/bin/env Rscript

######################################################
######## INIT
######################################################

CTRL <- 1
TREAT <- 2

######################################################
######## FUNCTIONS
######################################################

load_list <- function(input) {
  if (file.exists(input)) {
    list <- readLines(input)
  } else {
    list <- strsplit(input, ",")[[1]]
  }
  return(list)
}


load_table <- function(input_file, blacklist = NULL, whitelist = NULL, filter = NULL) {
  sample_table <- read.table(input_file, header= TRUE, quote="", sep="\t")
  if(!is.null(blacklist)){
    if(!any(blacklist %in% sample_table[[1]])) {
      warning("None of the samples in blacklist appear in experiment design")
    }
	  sample_table <- sample_table[!sample_table[[1]] %in% blacklist,]
  }
  if(!is.null(whitelist)){
	  sample_table <- sample_table[sample_table[[1]] %in% whitelist,]
    if(!any(blacklist %in% sample_table[[1]])) {
      stop("None of the samples in whitelist appear in experiment design")
    }
  }
  if(!is.null(filter)){ # Keep records with a specific value in a given column
	  col_name <- filter[1]
	  select_value <- filter[2]
	  sample_table <- sample_table[which(sample_table[[col_name]] == select_value),]
  }
  return(sample_table) # TODO: Adjust the output table for the remaining script. Originally was a nested list
}


parse_target <- function(target_string) {
  parsed_target <- list()
  
  # Split by ">"
  target_parts <- strsplit(target_string, ">")[[1]]
  target_name <- target_parts[1]
  all_features <- strsplit(target_parts[2], ";")[[1]]
  
  for (factor in all_features) {
    # Split by ":"
    feature_split <- strsplit(factor, ":")[[1]]
    feature_name <- feature_split[1]
    features <- strsplit(feature_split[2], ",")[[1]]
    
    # Split each feature by "/"
    features <- unlist(lapply(features, function(feature) {
      strsplit(feature, "/")[[1]]
    }))
    
    parsed_target[[feature_name]] <- features
  }
  
  return(list(target_name = target_name, target = parsed_target))
}


build_target <- function(table, target) {
  samples <- experiment_design$sample
  target_ctrls <- target_treats <- TRUE
  for(factor in names(target)) {
    factor_ctrls <- experiment_design[factor] == target[[factor]][CTRL]
    target_ctrls <- target_ctrls & factor_ctrls
    factor_treats <- experiment_design[factor] == target[[factor]][TREAT]
    target_treats <- target_treats & factor_treats
  }
  ctrl_features <- samples[target_ctrls]
  treat_features <- samples[target_treats]
  new_target <- list(Ctrl = ctrl_features, Treat = treat_features)
  return(new_target)
}


filter_features <- function(features, factor) {
  filtered_ft <- list()
  
  for (ft_name in names(features)) {
    feature <- features[[ft_name]]
    # factor is 0-based index in R (CTRL = 0, TREAT = 1)
    filtered_ft[[ft_name]] <- feature[[factor + 1]]
  }
  
  return(filtered_ft)
}


find_features <- function(features, table) {
  samples_list <- c()
  
  for (sample in names(table)) {
    all_features <- table[[sample]]
    include_sample <- TRUE
    
    for (feature_name in names(features)) {
      ft_values <- features[[feature_name]]
      
      if (!(all_features[[feature_name]] %in% ft_values)) {
        include_sample <- FALSE
        break
      }
    }
    
    if (include_sample) {
      samples_list <- c(samples_list, sample)
    }
  }
  
  return(samples_list)
}


save_target <- function(target_name, treats, output_path, experiment_design, additional_columns) {
  output_file <- file.path(output_path, paste0(target_name, "_target.txt"))
  
  # Prepare header
  header <- "sample\ttreat"
  if (length(additional_columns) > 0) {
    header <- paste(c("sample", "treat", additional_columns), collapse = "\t")
  }
  
  # Open file for writing
  con <- file(output_file, "w")
  writeLines(header, con)
  
  # Write data
  for (treat in names(treats)) {
    samples <- treats[[treat]]
    
    for (sample in samples) {
      if (length(additional_columns) > 0) {
        features_by_sample <- c()
        
        for (additional_feature in additional_columns) {
          ft_match <- experiment_design[experiment_design$sample == sample, additional_feature]
          if (!is.null(ft_match)) {
            features_by_sample <- c(features_by_sample, ft_match)
          }
        }
        
        if (!is.null(sample)) {
          line <- paste(c(sample, treat, features_by_sample), collapse = "\t")
          writeLines(line, con)
        }
      } else {
        if (!is.null(sample)) {
          line <- paste(sample, treat, sep = "\t")
          writeLines(line, con)
        }
      }
    }
  }
  
  close(con)
}


save_aux_options <- function(target_name, output_path, aux_options) {
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  aux_file <- file.path(output_path, paste0(target_name, "_target.aux"))
  con <- file(aux_file, "w")
  writeLines(aux_options, con)
  close(con)
}


parse_filter <- function(string) {
  parts <- strsplit(string, "=")[[1]]
  feature_name <- parts[1]
  features <- strsplit(parts[2], ",")[[1]]
  filter <- c(feature_name, features)
  return(filter)
}


######################################################
######## OPTIONS
######################################################

option_list <- list(
  optparse::make_option(c("-e", "--exp_file"), type = "character", default = NULL,
              help = "Tabulated file which describes current experiment"),
  optparse::make_option(c("-f", "--filter"), type = "character", default = NULL,
              help = "Set filter as string 'FEATURE_NAME=feature'"),
  optparse::make_option(c("-t", "--target"), type = "character", default = NULL,
              help = "String which describes target. EXAMPLE: 'TARGET_A>COLUMN_A:FEAT_CTL1,FEAT_TRT1;TARGET_B>COLUMN_B:FEAT_CTL1/FEAT_CTL2,FEAT_TRT1/FEAT_TRT2'"),
  optparse::make_option("--additional_features", type = "character", default = "",
              help = "String with extra factors separated by commas to be added to target."),
  optparse::make_option("--aux_options", type = "character", default = "",
              help = "String with extra options for DEGenesHunter."),
  optparse::make_option(c("-b", "--blacklist"), type = "character", default = NULL,
              help = "List with samples name to exclude from target. File or comma separated string"),
  optparse::make_option(c("-w", "--whitelist"), type = "character", default = NULL,
              help = "List with samples name to accept from target. File or comma separated string"),
  optparse::make_option(c("-o", "--output_path"), type = "character", default = ".",
              help = "Set the output path")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

opt$aux_options <- gsub("\"", "", opt$aux_options)
opt$aux_options <- gsub("'", "", opt$aux_options)

######################################################
######## MAIN
######################################################

blacklist <- NULL
whitelist <- NULL

if (!is.null(opt$blacklist)) {
  blacklist <- load_list(opt$blacklist)
}

if (!is.null(opt$whitelist)) {
  whitelist <- load_list(opt$whitelist)
}

filter <- NULL
if (!is.null(opt$filter)) {
  filter <- parse_filter(opt$filter)
}

experiment_design <- load_table(opt$exp_file, blacklist, whitelist, filter)
parsed <- parse_target(opt$target)
target_name <- parsed$target_name
target <- parsed$target

target <- build_target(experiment_design, target)
additional_features <- NULL
if (opt$additional_features != "") {
  additional_features <- strsplit(opt$additional_features, ",")[[1]]
}

save_target(target_name, target, opt$output_path, experiment_design, additional_features)

if (opt$aux_options != "") {
  save_aux_options(target_name, opt$output_path, opt$aux_options)
}
