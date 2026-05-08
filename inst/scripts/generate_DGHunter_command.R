#!/usr/bin/env Rscript

#####################################################
#################### LIBRARIES ######################
#####################################################

library(optparse)

#####################################################
##################### METHODS #######################
#####################################################

generate_callback <- function(name) {
  callback <- function(opt, flag, opt_val, parser) {
    return(paste(name, flag, opt_val, sep = " "))
  }
  return(callback)
}

parse_env_variables <- function(variables) {
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
  if(cmd == "") {
    fixed_opt <- NULL
  } else {
    split_cmd <- unlist(strsplit(cmd, " "))
    if (mode == 'degenes_Hunter') {
      option_list <- list(
      make_option(c("-C", "--Control_columns"), type = "character",
        callback = generate_callback("Control_columns")),
      make_option(c("-T", "--Treatment_columns"), type = "character",
        callback = generate_callback("Treatment_columns")),
      make_option(c("-r", "--reads"), type = "integer",
        callback = generate_callback("min_reads")),
      make_option(c("-l", "--minlibraries"), type = "integer",
        callback = generate_callback("min_libraries")),
      make_option(c("-F", "--filter_type"), type = "character",
        callback = generate_callback("filter_type")),
      make_option(c("-p", "--p_val_cutoff"), type = "numeric", 
        callback = generate_callback("de_pvalue")),
      make_option(c("-f", "--lfc"), type = "numeric",
        callback = generate_callback("de_logfc")),
      make_option(c("-m", "--modules"), type = "character",
        callback = generate_callback("de_packages")),
      make_option(c("-c", "--minpack_common"), type = "integer",
        callback = generate_callback("de_min_pack")),
      make_option(c("-t", "--target_file"), type = "character",
        callback = generate_callback("target_path")),
      make_option(c("-e", "--external_DEA_file"), type = "character",
        callback = generate_callback("external_DEA_file")),
      make_option(c("-v", "--model_variables"), type = "character",
        callback = generate_callback("de_add_factors")),
      make_option(c("-S", "--string_factors"), type = "character",
        callback = generate_callback("string_features")),
      make_option(c("-N", "--numeric_factors"), type = "character",
        callback = generate_callback("numeric_features")),
      make_option(c("-b", "--WGCNA_memory"), type = "numeric",
        callback = generate_callback("WGCNA_memory")),
      make_option("--WGCNA_norm_method", type = "character",
        callback = generate_callback("WGCNA_norm_method")),
      make_option("--WGCNA_deepsplit", type = "integer",
        callback = generate_callback("WGCNA_deepsplit")),
      make_option("--WGCNA_min_genes_cluster", type = "integer",
        callback = generate_callback("WGCNA_min_genes_cluster")),
      make_option("--WGCNA_detectcutHeight", type = "integer",
        callback = generate_callback("WGCNA_detectcutHeight")),
      make_option("--WGCNA_mergecutHeight", type = "integer",
        callback = generate_callback("WGCNA_mergecutHeight")),
      make_option(c("-w", "--WGCNA_all"), type = "logical",
        callback = generate_callback("WGCNA_ALL")),
      make_option("--WGCNA_blockwiseNetworkType", type = "character",
        callback = generate_callback("WGCNA_blockwiseNetworkType")),
      make_option("--WGCNA_blockwiseTOMType", type = "character",
        callback = generate_callback("WGCNA_blockwiseTOMType")),
      make_option("--WGCNA_minCoreKME", type = "integer",
        callback = generate_callback("WGCNA_minCoreKME")),
      make_option("--WGCNA_minCoreKMESize", type = "integer",
        callback = generate_callback("WGCNA_minCoreKMESize")),
      make_option("--WGCNA_minKMEtoStay", type = "integer",
        callback = generate_callback("WGCNA_minKMEtoStay")),
      make_option("--WGCNA_corType", type = "character",
        callback = generate_callback("WGCNA_corType")),
      make_option("--multifactorial", type = "character",
        callback = generate_callback("multifactorial")),
      make_option(c("-q", "--query_genes"), type = "character",
        callback = generate_callback("query_genes")),
      make_option("--seed", type = "integer",
        callback = generate_callback("seed")),
      make_option("--count_var_quantile", type = "numeric",
        callback = generate_callback("count_var_quantile")),
      make_option("--deseq2_var_quantile", type = "numeric",
        callback = generate_callback("deseq2_var_quantile"))
      )
    } else if (mode == 'functional_Hunter') {
      option_list <- list(
      make_option(c("-m", "--model_organism"), type = "character",
        callback = generate_callback("fun_organism")),
      make_option(c("-a", "--annot_file"), type = "character",
        callback = generate_callback("annotation_list")),
      make_option(c("-t", "--input_gene_id"), type = "character",
        callback = generate_callback("input_gene_id")),
      make_option(c("-f", "--func_annot_db"), type = "character",
        callback = generate_callback("func_annot_db")),
      make_option(c("-G", "--GO_subont"), type = "character",
        callback = generate_callback("GO_subont")),
      make_option(c("-C", "--custom"), type = "character",
        callback = generate_callback("custom_nomenclature")),
      make_option(c("-A", "--analysis_type"), type = "character",
        callback = generate_callback("fun_an_performance")),
      make_option(c("-r", "--remote"), type = "character",
        callback = generate_callback("fun_remote_mode")),
      make_option("--clean_parentals", type = "logical",
        callback = generate_callback("clean_parentals")),
      make_option(c("-P", "--pthreshold"), type = "logical",
        callback = generate_callback("pthreshold")),
      make_option(c("-Q", "--qthreshold"), type = "double",
        callback = generate_callback("qthreshold")),
      make_option("--max_genes_plot", type = "double",
        callback = generate_callback("max_genes_plot")), 
      make_option(c("-c", "--cores"), type = "integer",
        callback = generate_callback("cores")),
      make_option(c("-s", "--task_size"), type = "integer",
        callback = generate_callback("task_size")),
      make_option(c("-u", "--universe"), type = "character",
        callback = generate_callback("universe"))
      )
    }
    opt <- parse_args(OptionParser(option_list = option_list), args = split_cmd)
    # fixed_options <- fix_option(opt[[1]])
    fixed_opt <- list()
    for(i in seq(1, length(opt) - 1)) {
      new_element <- fix_option(opt[[i]])
      fixed_opt[i] <- new_element
      names(fixed_opt)[i] <- names(new_element)
    }
  }
  return(fixed_opt)
}


generate_command <- function(variables) {
  command <- character(0)
  for(variable in variables) {
    split_var <- unlist(strsplit(variable, " "))
    if(length(split_var) < 2 | split_var[2] == "FALSE") next
    if(is.null(split_var[2])) split_var[2] <- ""
    if(split_var[2] == "TRUE"){
      option <- split_var[1]
    } else {
      option <- paste(split_var, collapse = " ")
    }
    command <- paste(command, option, sep = " ")
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

option_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

########################
## MAIN
########################
# we use matrix data to get a data structure that preserves variable name-flag relation clearly.
de_variables <- list(de_pvalue = "-p", de_packages = "-m", de_min_pack = "-c", de_logfc = "-f",
                     WGCNA_mergecutHeight = "--WGCNA_mergecutHeight", WGCNA_min_genes_cluster = "--WGCNA_min_genes_cluster",
                     WGCNA_detectcutHeight = "--WGCNA_detectcutHeight", WGCNA_deepsplit = "--WGCNA_deepsplit", min_reads = "-r",
                     filter_type = "--filter_type", min_libraries = "-l", string_features = "-S", numeric_features = "-N",
                     target_path = "-t", query_genes = "-q", seed = "--seed", count_var_quantile = "--count_var_quantile",
                     deseq2_var_quantile = "--deseq2_var_quantile")

fun_variables <- list(fun_remote_mode = "-r", custom_nomenclature = "-C", fun_an_type = "-f", GO_modules = "-G",
                      fun_an_performance = "-A", fun_pvalue = "-P", fun_organism = "-m", annotation_list = "-a",
                      universe = "-u", clean_parentals = "--clean_parentals")

if (opt$mode == 'degenes_Hunter') {
  var_pairs <- de_variables
} else if (opt$mode == 'functional_Hunter') {
  var_pairs <- fun_variables
}
variables <- parse_env_variables(var_pairs)
additional_options <- parse_string_command(Sys.getenv("ADD_OPTIONS"), opt$mode)
if(!is.null(additional_options)) {
  variables <- modifyList(variables, additional_options)
}

if (opt$mode == 'degenes_Hunter') { # Parse auxiliary file for degenes_Hunter mode
  target_path <- strsplit(variables$target_path, " ")[[1]][2]
  aux_path <- sub("\\.txt$", ".aux", target_path)
  if(length(aux_path) > 0) {
    if(file.exists(aux_path)) {
      aux_content <- parse_string_command(readLines(aux_path), mode = opt$mode)
      variables <- modifyList(variables, aux_content)
    }
  }
}

command <- generate_command(variables)
cat(command)
