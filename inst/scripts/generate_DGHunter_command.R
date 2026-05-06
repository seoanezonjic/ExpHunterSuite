#!/usr/bin/env Rscript

#####################################################
#################### LIBRARIES ######################
#####################################################

library(optparse)

#####################################################
## METHODS


parse_env_variables <- function(variables, mode) {
  get_env_var <- function(variable){ Sys.getenv(variable) }
  hunter_flag_values <- lapply(names(variables), get_env_var)
  var_flag_value <- paste(variables, hunter_flag_values, sep = " ")
  names(var_flag_value) <- names(variables)
  return(as.list(var_flag_value))
}

fix_option <- function(option) {
  split_option <- unlist(strsplit(option, " "))
  res <- list(paste(split_option[-1], collapse = " "))
  names(res) <- split_option[1]
  return(res)
}

parse_string_command <- function(cmd, mode) {
  split_cmd <- unlist(strsplit(cmd, " "))
  if (mode == 'degenes_Hunter') {
    option_list <- list(
    make_option(c("-p", "--p_val_cutoff"), type = "numeric", 
      callback=function(opt, flag, opt_val, parser){return(paste("de_pvalue", flag, opt_val, sep=' '))}),
    make_option(c("-m", "--modules"), type = "character",
      callback=function(opt, flag, opt_val, parser){return(paste("de_packages", flag, opt_val, sep=' '))})
    )
  } else if (mode == 'functional_Hunter') {
    print("jaja")
  }
  opt <- parse_args(OptionParser(option_list = option_list), args = split_cmd)
  # fixed_options <- fix_option(opt[[1]])
  fixed_opt <- list()
  for(i in seq(1, length(opt) - 1)) {
    new_element <- fix_option(opt[[i]])
    fixed_opt[i] <- new_element
    names(fixed_opt)[i] <- names(new_element)
  }
  return(fixed_opt)
}


generate_command <- function(variables) {
  command <- paste(unlist(variables), collapse = " ")
  return(command)
}

########################
## OPTIONS
########################
option_list <- list(
  optparse::make_option(c("-m", "--mode"), type="character", default=NULL,
    help="Set DEGenesHunter mode. Available options are 'degenes_Hunter' and 'functional_Hunter'.")
)

option_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

########################
## MAIN
########################
# we use matrix data to get a data structure that preserves variable name-flag relation clearly.
de_variables <- list(de_pvalue = "-p", de_packages = "-m")
#   "de_min_pack", "-c",
#   "de_logfc", "-f",
#   "WGCNA_mergecutHeight", "--WGCNA_mergecutHeight",
#   "WGCNA_min_genes_cluster", "--WGCNA_min_genes_cluster",
#   "WGCNA_detectcutHeight", "--WGCNA_detectcutHeight",
#   "WGCNA_deepsplit", "--WGCNA_deepsplit",
#   "min_reads", "-r",
#   "filter_type", "--filter_type",
#   "min_libraries", "-l",
#   "string_features", "-S",
#   "numeric_features", "-N",
#   "target_path", "-t",
#   "query_genes", "-q",
#   "seed", "--seed",
#   "count_var_quantile", "--count_var_quantile",
#   "deseq2_var_quantile", "--deseq2_var_quantile"
#   ), ncol = 2, byrow = TRUE 
# )

# fun_variables <- matrix(c(
#   "fun_remote_mode", "-r",
#   "custom_nomenclature", "-C",
#   "fun_an_type", "-f",
#   "GO_modules", "-G",
#   "fun_an_performance", "-A",
#   "fun_pvalue", "-P",
#   "fun_organism", "-m",
#   "annotation_list", "-a",
#   "universe", "-u",
#   "clean_parentals", "--clean_parentals"
#   ), ncol = 2, byrow = TRUE
# )

if (opt$mode == 'degenes_Hunter') {
  var_pairs <- de_variables
} else if (opt$mode == 'functional_Hunter') {
  var_pairs <- fun_variables
}
variables <- parse_env_variables(var_pairs, opt$mode)
additional_options <- parse_string_command(Sys.getenv("ADD_OPTIONS"), opt$mode)
variables <- modifyList(variables, additional_options)

# additional_options <- Sys.getenv("ADD_OPTIONS")
# if (additional_options != "") {
#   variables <- overwrite_options(additional_options, opt$mode, variables)
# }

# if (opt$mode == 'degenes_Hunter') { # Parse auxiliary file for degenes_Hunter mode
#   target_path_idx <- which(variables[,2] == "target_path")
#   aux_path <- sub("\\.txt$", ".aux", variables[target_path_idx, 3])
#   if(length(aux_path) > 1) {
#     if(file.exists(aux_path)) {
#       aux_content <- readLines(aux_path)
#       variables <- overwrite_options(paste(aux_content, collapse = " "), opt$mode, variables)
#     }
#   }
# }

command <- generate_command(variables)
cat(command)
