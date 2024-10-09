################## GENERAL FUNCTIONS ########################

#' Method used to import hidden functions legally for CRAN check
#' @param pkg package
#' @param fun function to be luaded
#' @return function requried
get_unexported_function <- function(pkg, fun){
  return(get(fun, envir = asNamespace(pkg),inherits = FALSE))
} 


check_and_create_dir <- function(folder){
 folder <- file.path(folder)
 if (!dir.exists(folder)){
    dir.create(folder, recursive = TRUE)
 } 
}


handling_errors <- function(a){
  normalized_counts <- NULL
  expres_diff <- NULL
  all_genes_df <- NULL
  return(list(expres_diff, normalized_counts, all_genes_df))
}

standard_error <- function(x) {
  sd(x) / sqrt(length(x))
}

########################################################  
# Functions to generate output files
unite_DEG_pack_results <- function(exp_results, 
                                   p_val_cutoff, 
                                  lfc, 
                                  minpack_common) {
  DEG_pack_columns <- exp_results[['DEG_pack_columns']] 
  all_DE <- exp_results[['all_counts_for_plotting']] 
  all_FDR_names <- exp_results[['all_FDR_names']]
  all_LFC_names <- exp_results[['all_LFC_names']] 
  all_pvalue_names <- exp_results[['all_pvalue_names']]
  final_pvalue_names <- exp_results[['final_pvalue_names']]
  final_logFC_names <- exp_results[['final_logFC_names']]
  final_FDR_names <- exp_results[['final_FDR_names']]
  # Reorder all output by gene id so all equal
  all_DE <- lapply(all_DE, function(x) x[order(row.names(x)),])
  gene_ids <- row.names(all_DE[[1]])
  
  # Initialise output dataframe and add p/FDR/logFC 
  #values from the different packages
  all_DE_df <- data.frame(row.names=gene_ids)
  for(i in seq(length(all_DE))) {
    all_DE_df[final_logFC_names[i]] <- all_DE[[i]][, all_LFC_names[i]]
    all_DE_df[final_FDR_names[i]] <- all_DE[[i]][, all_FDR_names[i]]
    all_DE_df[final_pvalue_names[i]] <- all_DE[[i]][, all_pvalue_names[i]]
  }

  # Add TRUE/FALSE for each DEG package
  for(i in seq(length(all_DE))) {
    all_DE_df[DEG_pack_columns[i]] <- all_DE_df[, 
                              final_FDR_names[i]] < p_val_cutoff & 
                              abs(all_DE_df[, final_logFC_names[i]]) >= lfc
    all_DE_df[DEG_pack_columns[i]][
    is.na(all_DE_df[DEG_pack_columns[i]])] <- FALSE
    # Check: if no DE genes for package, give warning
    if(sum(all_DE_df[, DEG_pack_columns[i]]) == 0) 
      warning(paste("No significant", DEG_pack_columns[i], "found"))
  }
  # Get DEG_counts
  all_DE_df["DEG_counts"] <- rowSums(all_DE_df[DEG_pack_columns], na.rm = TRUE)
 
  # Calc and add combined FDR-values
  all_DE_df$combined_FDR <- vectorial_fisher_method(all_DE_df[final_FDR_names])

  # Reorder by combined FDR value
  all_DE_df <- all_DE_df[order(all_DE_df[,"combined_FDR"]), ]
  
  # Label as significant or not using combined FDR values
  all_DE_df[, "FDR_labeling"] <- ifelse(
                                   all_DE_df[, "combined_FDR"] < p_val_cutoff, 
                                   "SIGN", 
                                   "NOTSIGN")

  # Calculate average fold changes
  all_DE_df[, "mean_logFCs"] <- rowMeans(all_DE_df[final_logFC_names])

  # Add PREVALENT_DEG tag if as many as minpack_common; 
  # POSSIBLE_DEG if less but > 0; NOT_DEG if == 0
  genes_tag <- ifelse(all_DE_df[, "DEG_counts"] >= minpack_common, 
                      yes = "PREVALENT_DEG",
                      no = ifelse(all_DE_df[, "DEG_counts"] > 0, 
                                  yes = "POSSIBLE_DEG", 
                                  no = "NOT_DEG")
  )
  all_DE_df[, "genes_tag"] <- genes_tag

  return(all_DE_df)
}

add_filtered_genes <- function(all_DE_df, raw) {
  filtered_genes <- row.names(raw)[! row.names(raw) %in% row.names(all_DE_df)]
  filtered_df <- as.data.frame(matrix(NA, 
                                      ncol = ncol(all_DE_df), 
                                      nrow = length(filtered_genes)))
  colnames(filtered_df) <- colnames(all_DE_df)
  row.names(filtered_df) <- filtered_genes
  filtered_df[, "genes_tag"] <- rep("FILTERED_OUT", dim(filtered_df)[1])
  all_DE_df <- rbind(all_DE_df, filtered_df)

  return(all_DE_df)
}

debug_point <- function(file, message = "Debug point",envir = NULL){
    if(!dir.exists(dirname(file))){
      dir.create(dirname(file), recursive = TRUE)
    }
    if(is.null(envir)){
      current_envir <- environment() # current environment
      envir <- parent.env(current_envir)
    }

    envir$debug_message <- message
    on.exit(file.remove(file))
    save(list = names(envir), file = file, 
         envir = envir, precheck = FALSE)      
    on.exit()
}

#' @importFrom ggplot2 ggplot aes_string geom_bar theme element_text
#' @importFrom grDevices pdf dev.off
#' @importFrom utils write.table
save_times <- function(time_control, 
                       output="times_control.txt", 
                       plot_name = "time_control.pdf"){
 spent_times <- list()
    invisible(lapply(seq(2,length(time_control)), function(time_control_i){
      spent_times[[names(time_control)[time_control_i]]] <<- as.numeric(
                                      unlist(time_control[time_control_i]) - 
                                      unlist(time_control[time_control_i - 1]))
    }))
    spent_times_df <- do.call(rbind.data.frame, spent_times)
    colnames(spent_times_df) <- c("time")
    spent_times_df$control <- names(spent_times)
    pp <- ggplot2::ggplot(spent_times_df, ggplot2::aes_string(x = "control", 
                                                              y = "time")) +
          ggplot2::geom_bar(stat = "identity") +
          ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 25, 
                                                           hjust=1))
    grDevices::pdf(file.path(dirname(output), plot_name))    
      plot(pp)
    grDevices::dev.off()
    utils::write.table(spent_times_df, 
                       file=output, 
                       quote=FALSE, 
                       row.names=FALSE, 
                       sep="\t")
}

#' Parallelize function execution in list with multiple elements
#' @param X R list
#' @param FUNC Function to execute in each element of the list
#' @param workers Number of process to use in parallel execution
#' @param task_size Number of list elements to be executed on each 
#' BiocParallel task
#' @param ... Additional arguments to BiocParallel or to FUNC
#' @keywords performance
#' @return list
#' @importFrom utils str
#' @importFrom BiocParallel MulticoreParam bptry bplapply bpok
parallel_list <- function(X, FUNC, workers=2, task_size=1, ...){
    log <- FALSE
    main_log_path <- file.path(getwd(), 'bcplogs')
    log_path <- NA_character_
    workers <- workers - 1 # Reserve one core for main execution process
    if( workers == 0) workers <- 1
    if(workers > 1){
      timestamp <- as.integer(Sys.time())
      log_path <- file.path(main_log_path, as.character(timestamp))
      if(file.exists(log_path)){
        timestamp = timestamp + 1
        log_path <- file.path(main_log_path, as.character(timestamp))
      }
      log <- TRUE
      dir.create(log_path, recursive = TRUE)
    }
    param <- BiocParallel::MulticoreParam( 
      workers, tasks = ceiling(length(X)/task_size), stop.on.error = TRUE,
      log = log, threshold = "INFO", logdir = log_path
    )
    message(paste('items:', length(X), 'task_size:', task_size))
    print(param)    
    res <- BiocParallel::bptry(
      BiocParallel::bplapply(X, FUNC, BPPARAM = param, ...)
    )
    exec_status <- BiocParallel::bpok(res)
    fails <- which( exec_status == FALSE)
    message(paste('exec_status => exec items:', 
                  length(exec_status), 
                  'fails:', 
                  length(fails)))
    if(length(fails) > 0 ){
      message(tail(attr(res[[fails[1]]], "traceback")))
      stop(paste('Parallel execution has failed at item', fails[1],
                 'and a total of', length(fails) , 'items have failed.'))
    } else {
      unlink(main_log_path, recursive = TRUE)
    }

    return(res)
}



rand_sample_bool <- function(bool, sample_size){
#takes a TRUE FALSE vector and keep only sample_size TRUEs
    sampled_bool <- rep(FALSE, length(bool)) 
    random_i <- sample(which(bool == TRUE), 
                        size = sample_size, 
                        replace = FALSE)
    sampled_bool[random_i] <- TRUE
    return(sampled_bool)
}

check_and_quit <- function(object){
  str(object)
  q()
}


#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom annotatr annotate_regions build_annotations
annotate_genomic_ranges <-function(intervals_df, genome){
 g_regs <- GenomicRanges::makeGRangesFromDataFrame(intervals_df, keep.extra.columns = TRUE)
 annots <- annotatr::build_annotations(genome = genome, annotations = paste0(genome,"_basicgenes"))
 annotated_g_regs <- annotatr::annotate_regions(regions = g_regs,annotations = annots, ignore.strand = FALSE)
 annotated_g_regs <- data.frame(annotated_g_regs)
 annotated_g_regs <- annotated_g_regs[!is.na(annotated_g_regs$annot.symbol),]
 return(annotated_g_regs)
}    

split_str <- function(string, split = NULL) {
  if (is.null(split))
    stop("split character must be specifyed")
  if (is.null(string)) return(NULL)  
  if (nchar(string) <= 1) return(NULL)  

  splitted_str <- unlist(strsplit(string, split = split))
  return(splitted_str)
}

merge_factors <- function(table, target, factors) {
  if (length(factors) == 0 ) {
    return(table)
  }  

  table <- merge(table, target[,factors, drop = FALSE], by = "row.names", all = TRUE)
    row.names(table) <-  table[,"Row.names"]
    table[,"Row.names"] <- NULL 
    return(table)
}

remove_pattern_from_df <- function(dataframe, rgx){
  df_parsed <- as.data.frame(lapply(dataframe,
                    remove_text_from_column, 
                    rgx = rgx))
  return(df_parsed)
}

remove_text_from_column <- function(column, rgx) {
  parsed_col <- gsub(rgx, "", column)
  return(parsed_col)
}

name_column <- function(column, row_names){
  names(column) <- row_names
  column
}

name_all_columns <- function(data_frame) {
   lapply(data_frame, name_column, row_names = rownames(data_frame))
}