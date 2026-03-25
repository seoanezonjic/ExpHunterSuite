#!/usr/bin/env Rscript

#####################################################
## METHODS


parse_env_variables <- function(variables, mode) {
  get_env_var <- function(var_pair){ Sys.getenv(var_pair[1]) }
  hunter_flag_values <- apply(variables, MARGIN = 1, get_env_var)
  save(list = ls(all = TRUE), file = "envir.RData")
  var_flag_value <- cbind(variables, hunter_flag_values)
  var_list <- vector(mode = "list", length = nrow(var_flag_value))
  names(var_list) <- var_flag_value[, 1]
  for(row in seq(nrow(var_flag_value))) {
    var <- var_flag_value[row, ]
    var_list[[row]] <- list(flag = var[2], value = var[3])
  }
  rownames(var_flag_value) <- var_flag_value[, 1]
  return(var_flag_value)
}

overwrite_options <- function(additional_options, mode, var_flag_value){
  new_opts <- parse_string_command(additional_options, mode)
  var_flag_value <- c(var_flag_value, new_opts)
  rownames(var_flag_value) <- var_flag_value[, 1]
  return(var_flag_value)
}

parse_string_command <- function(cmd, mode) {
  options <- list()
  
  if (mode == 'degenes_Hunter') {
    # Parse the command string by splitting on spaces
    cmd_parts <- unlist(strsplit(cmd, "\\s+"))
    
    i <- 1
    while (i <= length(cmd_parts)) {
      arg <- cmd_parts[i]
      
      if (arg == "-C" || arg == "--Control_columns") {
        options[["Control_columns"]] <- c("--Control_columns", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "-T" || arg == "--Treatment_columns") {
        options[["Treatment_columns"]] <- c("--Treatment_columns", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "-r" || arg == "--reads") {
        options[["min_reads"]] <- c("-r", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "-l" || arg == "--minlibraries") {
        options[["min_libraries"]] <- c("-l", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "-F" || arg == "--filter_type") {
        options[["filter_type"]] <- c("--filter_type", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "-p" || arg == "--p_val_cutoff") {
        options[["de_pvalue"]] <- c("-p", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "-f" || arg == "--lfc") {
        options[["de_logfc"]] <- c("-f", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "-m" || arg == "--modules") {
        options[["de_packages"]] <- c("-m", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "-c" || arg == "--minpack_common") {
        options[["de_min_pack"]] <- c("-c", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "-t" || arg == "--target_file") {
        options[["target_path"]] <- c("-t", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "-e" || arg == "--external_DEA_file") {
        options[["external_DEA_file"]] <- c("--external_DEA_file", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "-v" || arg == "--model_variables") {
        options[["de_add_factors"]] <- c("-v", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "-S" || arg == "--string_factors") {
        options[["string_features"]] <- c("-S", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "-N" || arg == "--numeric_factors") {
        options[["numeric_features"]] <- c("-N", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "-b" || arg == "--WGCNA_memory") {
        options[["WGCNA_memory"]] <- c("--WGCNA_memory", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "--WGCNA_norm_method") {
        options[["WGCNA_norm_method"]] <- c("--WGCNA_norm_method", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "--WGCNA_deepsplit") {
        options[["WGCNA_deepsplit"]] <- c("--WGCNA_deepsplit", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "--WGCNA_min_genes_cluster") {
        options[["WGCNA_min_genes_cluster"]] <- c("--WGCNA_min_genes_cluster", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "--WGCNA_detectcutHeight") {
        options[["WGCNA_detectcutHeight"]] <- c("--WGCNA_detectcutHeight", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "--WGCNA_mergecutHeight") {
        options[["WGCNA_mergecutHeight"]] <- c("--WGCNA_mergecutHeight", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "-w" || arg == "--WGCNA_all") {
        options[["WGCNA_ALL"]] <- c("-w", "BOOLEAN")
        i <- i + 1
      } else if (arg == "--WGCNA_blockwiseNetworkType") {
        options[["WGCNA_blockwiseNetworkType"]] <- c("--WGCNA_blockwiseNetworkType", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "--WGCNA_blockwiseTOMType") {
        options[["WGCNA_blockwiseTOMType"]] <- c("--WGCNA_blockwiseTOMType", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "--WGCNA_minCoreKME") {
        options[["WGCNA_minCoreKME"]] <- c("--WGCNA_minCoreKME", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "--WGCNA_minCoreKMESize") {
        options[["WGCNA_minCoreKMESize"]] <- c("--WGCNA_minCoreKMESize", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "--WGCNA_minKMEtoStay") {
        options[["WGCNA_minKMEtoStay"]] <- c("--WGCNA_minKMEtoStay", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "--WGCNA_corType") {
        options[["WGCNA_corType"]] <- c("--WGCNA_corType", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "--multifactorial") {
        options[["multifactorial"]] <- c("--multifactorial", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "-q" || arg == "--query_genes") {
        options[["query_genes"]] <- c("--query_genes", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "--seed") {
        options[["seed"]] <- c("--seed", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "--count_var_quantile") {
        options[["count_var_quantile"]] <- c("--count_var_quantile", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "--deseq2_var_quantile") {
        options[["deseq2_var_quantile"]] <- c("--deseq2_var_quantile", cmd_parts[i + 1])
        i <- i + 2
      } else {
        i <- i + 1
      }
    }
    
  } else if (mode == 'functional_Hunter') {
    cmd_parts <- unlist(strsplit(cmd, "\\s+"))
    
    i <- 1
    while (i <= length(cmd_parts)) {
      arg <- cmd_parts[i]
      
      if (arg == "-m" || arg == "--model_organism") {
        options[["fun_organism"]] <- c("-m", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "-a" || arg == "--annot_file") {
        options[["annotation_list"]] <- c("-a", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "-t" || arg == "--input_gene_id") {
        options[["input_gene_id"]] <- c("-t", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "-f" || arg == "--func_annot_db") {
        options[["func_annot_db"]] <- c("-f", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "-G" || arg == "--GO_subont") {
        options[["GO_subont"]] <- c("-G", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "-C" || arg == "--custom") {
        options[["custom_nomenclature"]] <- c("-C", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "-A" || arg == "--analysis") {
        options[["fun_an_performance"]] <- c("-A", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "-r" || arg == "--remote") {
        options[["fun_remote_mode"]] <- c("-r", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "--clean_parentals") {
        options[["clean_parentals"]] <- c("--clean_parentals", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "-P" || arg == "--pthreshold") {
        options[["pthreshold"]] <- c("-P", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "-Q" || arg == "--qthreshold") {
        options[["qthreshold"]] <- c("-Q", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "--max_genes_plot") {
        options[["max_genes_plot"]] <- c("--max_genes_plot", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "-c" || arg == "--cores") {
        options[["cores"]] <- c("-c", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "-s" || arg == "--task_size") {
        options[["task_size"]] <- c("-s", cmd_parts[i + 1])
        i <- i + 2
      } else if (arg == "-u" || arg == "--universe") {
        options[["universe"]] <- c("-u", cmd_parts[i + 1])
        i <- i + 2
      } else {
        i <- i + 1
      }
    }
  }
  
  return(options)
}


generate_command <- function(variables) {
  command <- ""
  
  for (variable in rownames(variables)) {
    attributes <- variables[variable, ]
    flag <- attributes[2]
    value <- attributes[3]
    
    # Check if we should skip this variable
    if ((is.na(value) || value == "" || value == "FALSE") && length(attributes) > 1) {
      next
    }
    
    if (value == "NULL") {
      value <- ""
    }
    
    if (value == "TRUE") {
      option <- paste0(flag, " ")
    } else {
      option <- paste0(flag, " \"", value, "\" ")
    }
    
    command <- paste0(command, option)
  }
  
  return(command)
}

########################
## OPTIONS
########################
option_list <- list(
  optparse::make_option(c("-m", "--mode"), type="character", default=NULL,
    help="Set DEGenesHunter mode. Available options are 'degenes_Hunter' and 'functional_Hunter'.")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

########################
## MAIN
########################
# we use matrix data to get a data structure that preserves variable name-flag relation clearly.
de_variables <- matrix(c(
  "de_pvalue", "-p",
  "de_packages", "-m",
  "de_min_pack", "-c",
  "de_logfc", "-f",
  "WGCNA_mergecutHeight", "--WGCNA_mergecutHeight",
  "WGCNA_min_genes_cluster", "--WGCNA_min_genes_cluster",
  "WGCNA_detectcutHeight", "--WGCNA_detectcutHeight",
  "WGCNA_deepsplit", "--WGCNA_deepsplit",
  "min_reads", "-r",
  "filter_type", "--filter_type",
  "min_libraries", "-l",
  "string_features", "-S",
  "numeric_features", "-N",
  "target_path", "-t",
  "query_genes", "-q",
  "seed", "--seed",
  "count_var_quantile", "--count_var_quantile",
  "deseq2_var_quantile", "--deseq2_var_quantile"
  ), ncol = 2, byrow = TRUE 
)

fun_variables <- matrix(c(
  "fun_remote_mode", "-r",
  "custom_nomenclature", "-C",
  "fun_an_type", "-f",
  "GO_modules", "-G",
  "fun_an_performance", "-A",
  "fun_pvalue", "-P",
  "fun_organism", "-m",
  "annotation_list", "-a",
  "universe", "-u",
  "clean_parentals", "--clean_parentals"
  ), ncol = 2, byrow = TRUE
)

if (opt$mode == 'degenes_Hunter') {
  var_pairs <- de_variables
} else if (opt$mode == 'functional_Hunter') {
  var_pairs <- fun_variables
}
variables <- parse_env_variables(var_pairs, opt$mode)

additional_options <- Sys.getenv("ADD_OPTIONS")
if (additional_options != "") {
  variables <- overwrite_options(additional_options, opt$mode, variables)
}

if (opt$mode == 'degenes_Hunter') { # Parse auxiliary file for degenes_Hunter mode
  target_path_idx <- which(variables[,2] == "target_path")
  aux_path <- sub("\\.txt$", ".aux", variables[target_path_idx, 3])
  if(length(aux_path) > 1) {
    if(file.exists(aux_path)) {
      aux_content <- readLines(aux_path)
      variables <- overwrite_options(paste(aux_content, collapse = " "), opt$mode, variables)
    }
  }
}

variables <- coerce_variables(variables)

coerce_variables <- function(variables) {
  for(var in names(variables)) {
    variables[[var]]$value <- tryCatch(as.numeric())
  }
}

save.image('Testing.RData')

command <- generate_command(variables)
cat(command)