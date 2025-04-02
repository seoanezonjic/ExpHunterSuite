#' Main DEgenes Hunter Function
#'
#' This function allows you to perform the DEG analysis with different 
#' algorithms.
#' @param raw Dataframe with raw counts per sample and feature
#' @param target set with control and treatment information
#' @param external_DEA_data external DEA set 
#' @param count_var_quantile quantile variance threshold
#' @param output_files otput files path
#' @param reads set of reads
#' @param minlibraries minimum of libraries to use a set of data
#' @param filter_type filter type to ba applied. Allowed: "separate" or "global"
#' @param p_val_cutoff p-value threshold
#' @param lfc minimum logFold Change
#' @param modules modules to be executed. Allowed: DESeq2 (D), NOISeq (N), 
#' Limma (L), EdgeR (E), WGCNA (W)
#' @param minpack_common minimum of pack that must be significant to tag 
#' a gene as significant
#' @param model_variables custom model
#' @param pseudocounts boolean, activate if the input contains pseudocounts
#' @param numerics_as_factors transform numeric values to factors. Default: TRUE
#' @param string_factors string factors for WGCNA
#' @param numeric_factors numeric factors for WGCNA
#' @param WGCNA_memory see WGCNA package
#' @param WGCNA_norm_method see WGCNA package
#' @param WGCNA_deepsplit see WGCNA package
#' @param WGCNA_min_genes_cluster see WGCNA package
#' @param WGCNA_detectcutHeight see WGCNA package
#' @param WGCNA_mergecutHeight see WGCNA package
#' @param WGCNA_all see WGCNA package
#' @param WGCNA_blockwiseNetworkType see WGCNA package
#' @param WGCNA_blockwiseTOMType see WGCNA package
#' @param WGCNA_minCoreKME see WGCNA package
#' @param WGCNA_minCoreKMESize see WGCNA package
#' @param WGCNA_minKMEtoStay see WGCNA package
#' @param WGCNA_corType see WGCNA package
#' @param multifactorial specify interaction/effect when multifactorial design
#' @param library_sizes NULL or a dataframe with sample names and library sizes
#' @return expression analysis result object with studies performed
#' @keywords method
#' @importFrom rlang .data
#' @export
#' @examples
#' data(toc)
#' data(target)
#' degh_out <- main_degenes_Hunter(raw=toc, target=target, modules="D")

main_degenes_Hunter <- function(
    raw = NULL,
    pseudocounts = FALSE,
    target = NULL,
    count_var_quantile = 0,
    external_DEA_data = NULL,
    output_files = getwd(),
    reads = 2,
    minlibraries = 2,
    filter_type = "separate",
    p_val_cutoff = 0.05,
    lfc = 1,
    modules ="DELNW",
    minpack_common = 4,
    model_variables = "",
    numerics_as_factors = TRUE,
    string_factors = NULL,
    numeric_factors = NULL,
    WGCNA_memory = 5000,
    WGCNA_norm_method = "DESeq2",
    WGCNA_deepsplit = 2,
    WGCNA_min_genes_cluster = 20,
    WGCNA_detectcutHeight = 0.995,
    WGCNA_mergecutHeight = 0.25,
    WGCNA_all = FALSE,
    WGCNA_blockwiseNetworkType = "signed",
    WGCNA_blockwiseTOMType = "signed",
    WGCNA_minCoreKME = 0.5,
    WGCNA_minCoreKMESize = NULL,
    WGCNA_minKMEtoStay = 0.2,
    WGCNA_corType = "pearson",
    multifactorial = "",
    library_sizes = NULL
  ){
    modified_input_args <- check_input_main_degenes_Hunter(raw, 
      minlibraries, pseudocounts,reads, external_DEA_data, modules, model_variables,
      active_modules, WGCNA_all, minpack_common, target, 
      string_factors, numeric_factors, multifactorial)
    modules <- modified_input_args[['modules']]
    active_modules <- modified_input_args[['active_modules']]
    minpack_common <- modified_input_args[['minpack_common']]
    raw <- modified_input_args[['raw']]
    reads <- modified_input_args[['reads']]
    target <- modified_input_args[['target']]
    final_main_params <- list(
      'minpack_common' = minpack_common,
      'p_val_cutoff' = p_val_cutoff,
      'lfc' = lfc,
      'modules' = modules,
      'active_modules' = active_modules
    )

    ############################################################
    ##                         I/O DATA                       ##
    ############################################################


    # Infer replicates and group index from target
    index_control_cols <- as.character(target$sample[target$treat == "Ctrl"])
    index_treatmn_cols <- as.character(target$sample[target$treat == "Treat"])
    replicatesC <- length(index_control_cols)
    replicatesT <- length(index_treatmn_cols)
    design_vector <- c(rep("C", replicatesC), rep("T", replicatesT))
    sample_groups <- data.frame(class = design_vector, 
                                name = c(index_control_cols, 
                                         index_treatmn_cols))
    # Check if there are enough replicates for specified method
    if((replicatesC < 2) | (replicatesT < 2)) 
       stop(paste0('At least two replicates per class (i.e. treatment and',
                   ' control) are required\n'))

    
    numeric_factors <- split_str(numeric_factors, ",")
    if (sum(numeric_factors == "") >= 1) numeric_factors <- NULL #esto esta para controlar que no haya elementos vacios

    string_factors <- split_str(string_factors, ",")

    string_factors <-  c("treat", string_factors)

    
    if(!is.null(target) & grepl("W", modules)) {
      target_numeric_factors <- build_design_for_WGCNA(target, 
           numeric_factors=numeric_factors)
      target_string_factors <- build_design_for_WGCNA(target, 
          string_factors=string_factors)
    }

    if(numerics_as_factors == TRUE) { 
      # Now coerce the targets to factors for the multifactorial analysis
       target <-  data.frame(sapply(target, as.factor), 
         stringsAsFactors=TRUE)

    }
    model_formula_text <- prepare_model_text(model_variables, 
                                             multifactorial)
   
    if(grepl(":(interaction|effect)", multifactorial)){
      target <- prepare_target_for_multifactorial(target, multifactorial)
    }
    # Prepare count table for analysis
    #Indexing selected columns from input count dataframe   
    raw <- raw[c(index_control_cols,index_treatmn_cols)]
    raw[is.na(raw)] <- 0 # Substitute NA values

    filtered_data <- filter_count(reads, minlibraries, raw, filter_type, 
                               index_control_cols, index_treatmn_cols, target)
    raw_filter <- filtered_data[["raw"]]
    cpm_table <- filtered_data[["cpm_table"]]

    var_filter <- filter_by_variance(raw_filter, 
                                     q_filter = count_var_quantile, 
                                     target = target)
    default_norm <- list(default = 
                      as.data.frame(var_filter[["deseq2_normalized_counts"]]))
    raw_filter <- var_filter[["fil_count_mtrx"]]

    

    #computing PCA for all_genes
    full_pca <- compute_pca(pca_data = default_norm$default,
                            target = target,
                            string_factors = string_factors, 
                            numeric_factors = numeric_factors)
    PCA_res <- list("all_genes" = full_pca)
    ############################################################
    ##             PERFORM EXPRESION ANALYSIS                 ##
    ############################################################
    check_and_create_dir(output_files)


    exp_results <- perform_expression_analysis(modules, replicatesC, 
                     replicatesT, raw_filter, p_val_cutoff, target, 
                     model_formula_text, external_DEA_data, 
                     multifactorial)
    exp_results[["all_data_normalized"]] <- c(default_norm,
                                              exp_results[["all_data_normalized"]])

    #################################################################
    ##                       CORRELATION ANALYSIS                   ##
    ##################################################################

    combinations_WGCNA <- NULL
    if(grepl("W", modules)) { # CASE W: WGCNA
      cat('Correlation analysis is performed with WGCNA\n')
      path <- file.path(output_files, "Results_WGCNA")
      dir.create(path)
      all_data_normalized <- exp_results[['all_data_normalized']]
      if(WGCNA_norm_method %in% names(all_data_normalized)) {
        WGCNA_input <- all_data_normalized[[WGCNA_norm_method]]
      } else if(WGCNA_norm_method== "none") {
        WGCNA_input <- raw_filter
      } else {
        warning(paste0("To run WGCNA, you must also run the method chosen for",
          " normalization in the --modules flag, or specify none"))
      }

      combinations_WGCNA <- perform_WGCNA_combinations(WGCNA_all, WGCNA_input, 
        index_treatmn_cols, index_control_cols, path, target_numeric_factors, 
        target_string_factors, WGCNA_memory, WGCNA_deepsplit, 
        WGCNA_detectcutHeight, WGCNA_mergecutHeight, WGCNA_min_genes_cluster, 
        WGCNA_blockwiseNetworkType, WGCNA_blockwiseTOMType, WGCNA_minCoreKME, 
        WGCNA_minCoreKMESize, WGCNA_minKMEtoStay, corType = WGCNA_corType)

      if(length(combinations_WGCNA[['WGCNA_all']]) == 1) {
        warning("Something went wrong with WGCNA on the full dataset")
        modules <- gsub("W", "", modules)
      }
    }

    if(grepl("X", modules)) { # CASE X: diffcoexp
      cat('Correlation analysis is performed with diffcoexp\n')
      path <- file.path(output_files, "Results_diffcoexp")
      check_and_create_dir(path)

      results_diffcoexp <- analysis_diffcoexp(data = raw_filter,
                                           path = path,
                                           target = target)
    }


    #################################################################
    ##                    BUILD MAIN RESULT TABLE                  ##
    #################################################################
    mean_expression_cpm <- rowMeans(raw)
    min_expression_cpm <- matrixStats::rowMins(as.matrix(raw))
    max_expression_cpm <- matrixStats::rowMaxs(as.matrix(raw))
    cpm_stats <- data.frame(mean_expression_cpm, min_expression_cpm, max_expression_cpm)

    DE_all_genes <- NULL
    DEG_pca <- NULL
    if (any(grepl("[DENL]",modules))){
      DE_all_genes <- unite_DEG_pack_results(exp_results, p_val_cutoff, 
                                             lfc, minpack_common)

      #computing PCA for PREVALENT DEG
      prevalent_degs <- rownames(DE_all_genes[DE_all_genes$genes_tag == "PREVALENT_DEG",])
      if (length(prevalent_degs) > 2) {
        pca_deg_data <- default_norm$default
        pca_deg_data <- pca_deg_data[rownames(pca_deg_data) %in% prevalent_degs,]
        DEG_pca <- compute_pca(pca_data = pca_deg_data,
                              target = target,
                              string_factors = string_factors, 
                              numeric_factors = numeric_factors)
      }
    }
    PCA_res[["DEGs"]] <- DEG_pca

    if(grepl("W", modules)) { # Check WGCNA was run and returned proper results
      DE_all_genes <- merge(by.x=0, by.y="ENSEMBL_ID", x= DE_all_genes, 
          y =combinations_WGCNA[['WGCNA_all']][['gene_cluster_info']])
      rownames(DE_all_genes) <- DE_all_genes$Row.names
      DE_all_genes$Row.names <- NULL
    }

    # Add the filtered genes back
    DE_all_genes <- add_filtered_genes(DE_all_genes, raw)
    DE_all_genes <- merge(DE_all_genes, cpm_stats, by=0, sort=FALSE)
    DE_all_genes <- transform(DE_all_genes, row.names=Row.names, Row.names=NULL)
    
    correlation_metrics <- NULL
    if(grepl("P", modules)) { # CASE P: PCIT, TODO: RESTORE FUNCTION, PEDRO 
      # TODO : This is not working, variables "DESeq2_counts" 
      # are not being generated inside this function
      all_data_normalized <- exp_results[['all_data_normalized']]$DESeq2
      raw <- raw[c(index_control_cols,index_treatmn_cols)]
      correlation_metrics <- analysis_diff_correlation(
        DE_all_genes, 
        all_data_normalized, 
        all_data_normalized[index_control_cols], 
        all_data_normalized[index_treatmn_cols], 
        PCIT_filter=FALSE) # The analysis is very very slow, so we disable it
      DE_all_genes <- transform(merge(DE_all_genes, correlation_metrics, by.x=0, by.y=0), row.names=Row.names, Row.names=NULL)
    }

    if(grepl("X", modules)) { #results_diffcoexp
      aux <- row.names(DE_all_genes)
      DE_all_genes$DCG <- aux %in% 
         results_diffcoexp$DCGs$Gene[results_diffcoexp$DCGs$q < 1]
      DE_all_genes$DCL <- aux %in% results_diffcoexp$DCLs$Gene.1 | 
         aux %in% results_diffcoexp$DCLs$Gene.2
    }

   


    coverage_df <- get_counts(cnts_mtx=raw, library_sizes=library_sizes)
    mean_counts_df <- get_mean_counts(cnts_mtx=raw, cpm_table=cpm_table, reads=reads, minlibraries=minlibraries)
    exp_genes_df <- get_gene_stats(cpm_table=cpm_table, reads=reads)    

 

    final_results <- list()
    final_results[['cpm_table']] <- cpm_table
    final_results[['raw_filter']] <- raw_filter
    final_results[['sample_groups']] <- sample_groups
    final_results[['DE_all_genes']] <- DE_all_genes
    final_results[['index_control_cols']] <- index_control_cols
    final_results[['index_treatmn_cols']] <- index_treatmn_cols
    final_results[['design_vector']] <- design_vector
    final_results[['replicatesC']] <- replicatesC
    final_results[['replicatesT']] <- replicatesT
    final_results[['final_main_params']] <- final_main_params
    final_results[["var_filter"]] <- var_filter
    final_results[["coverage_df"]] <- coverage_df
    final_results[["mean_counts_df"]] <- mean_counts_df
    final_results[["exp_genes_df"]] <- exp_genes_df
    final_results[["target"]] <- target 
    final_results[["numeric_factors"]] <- numeric_factors
    final_results[["string_factors"]] <- string_factors
    final_results[["PCA_res"]] <- PCA_res
    final_results[["library_sizes"]] <- library_sizes
    
    if(!is.null(combinations_WGCNA)){
      final_results <- c(final_results, combinations_WGCNA)
    }
    final_results <- c(final_results, exp_results)
    return(final_results)


}


check_input_main_degenes_Hunter <- function(raw, 
                                            minlibraries, 
                                            pseudocounts = FALSE,
                                            reads, 
                                            external_DEA_data, 
                                            modules, 
                                            model_variables, 
                                            active_modules, 
                                            WGCNA_all, 
                                            minpack_common, 
                                            target, 
                                            string_factors, 
                                            numeric_factors,
                                            multifactorial){

    if (pseudocounts) {
      raw <- round(raw)
    }
    if (minlibraries < 1){
      stop(cat(paste0("Minimum library number to check minimum read counts",
        " cannot be less than 1.\nIf you want to avoid filtering, set",
        " --reads to 0.")))
    }

    if (reads == 0){
      cat(paste0("Minimum reads option is set to zero. ",
        "Raw count table will not be filtered."))
    }

    if (!is.null(external_DEA_data)) {
      message("External DEA dataframe given.")
      if(is.null(target)){
        warning(paste0("External DEA dataframe given but no target",
          " - creating a dummy one"))
        target <- data.frame(samples=c("case_dummy_1", "case_dummy_2", 
                                       "control_dummy_1", "control_dummy2"), 
                             treat=c("Treat","Treat", "Ctrl", "Ctrl"))
      }
      if(is.null(raw)) {
        warning(paste0("External DEA dataframe given but no count data",
          " - creating a dummy one"))
        raw <- as.data.frame(matrix(sample(1:100, nrow(external_DEA_data)*nrow(target),replace=TRUE), 
                                    nrow=nrow(external_DEA_data), 
                                    ncol=nrow(target), 
                                  dimnames = list(row.names(external_DEA_data), 
                                                  target$sample)))
        reads <- 0
      }
    }

    if (model_variables != "" & grepl("N", modules) & nchar(modules) == 1) {
      stop(cat(paste0("You cannot run an experimental design that uses the",
        " model_variables option with NOISeq only.")))
    }
    if (model_variables != "" & grepl("N", modules) & nchar(modules) > 1) {
      warning(paste0("NOISeq will not be run as you have an experimental",
        " design that uses the model_variables option."))
      modules <- gsub("N", "", modules)
    }
    if(model_variables != "" & multifactorial != ""){
      stop("Cannot provide both a --multifactorial design and additional vars",
        "via --model_variables")
    }

    if(grepl('P', modules) && !grepl('D', modules)){
      modules <- paste0(modules, 'D')
    }

    active_modules <- nchar(modules)
    user_modules <- unlist(strsplit(modules, '')) # TODO. Maybe use this variable for checks perfomed before
    active_modules <- active_modules - sum(grepl("[WPX]", user_modules)) # Remove from module count modules that are not involved in differential expresion
    if(sum(target$treat == "Ctrl") < 3 | sum(target$treat == "Treat") < 3)
          active_modules <- active_modules - sum(grepl("[LN]", user_modules))
    if(minpack_common > active_modules){
      minpack_common <- active_modules
      warning(paste0("The number of active modules is lower than the thresold",
        " for tag PREVALENT DEG. The thresold is set to the",
        " number of active modules."))
    }


    # If no -t check there is no -v value
    if( is.null(target) & model_variables != "") {
      stop(cat(paste0("You should not include a -v value if you do not",
        " have a target table file.")))
    }

    # If factors are specified but WGCNA not selected, throw a warning.
   
    return(list(modules=modules, 
                active_modules=active_modules, 
                minpack_common=minpack_common, 
                raw=raw, reads=reads, target=target))
}

#' @importFrom edgeR cpm
filter_count <- function(reads, 
                         minlibraries, 
                         raw, 
                         filter_type, 
                         index_control_cols, 
                         index_treatmn_cols,
                         target){
    # Prepare filtered set
    cpm_table <- edgeR::cpm(raw)
    if(reads != 0){
      if (filter_type == "separate") {
        # genes with cpm greater than --reads value for 
        # at least --minlibrariess samples for either case or control samples
        to_keep_control <- rowSums(cpm_table[, index_control_cols] > 
                                reads) >= minlibraries
        to_keep_treatment <- rowSums(cpm_table[, index_treatmn_cols] >
                                  reads) >= minlibraries
        # Keep if at least minlibraries in either of them
        keep_cpm <- to_keep_control | to_keep_treatment 
        raw <- raw[keep_cpm,] # Filter out count data frame
      } else if (filter_type == "global") {
        # genes with cpm greater than --reads value for
        # at least --minlibrariess samples
        keep_cpm <- rowSums(cpm_table >= reads) >= minlibraries 
        raw <- raw[keep_cpm,] # Filter out count data frame
      } else if (grepl("^combined", filter_type) == TRUE) {
        split_filter <- strsplit(filter_type, ":")[[1]]
        filt_factors <- strsplit(split_filter[2], ",")[[1]]
        if(! all(filt_factors %in% colnames(target)))
          stop("Filter factors not found in target")
        combs_list <- split(target, target[,filt_factors])
        samples_per_combo <- lapply(combs_list, function(x) as.vector(x$sample))
        cpm_tab_per_combo <- sapply(samples_per_combo, function(x) cpm_table[, x])
        combos_passing_filter <- sapply(cpm_tab_per_combo, 
          function(x) rowSums(x >= reads) >= minlibraries)
        keep_cpm <- apply(combos_passing_filter, 1, any)
        raw <- raw[keep_cpm,] # Filter out count data frame
      } else {
        warning("Unrecognized minimum read filter type. No filter will be used")
        raw <- raw
      }
    }
    return(list(raw=raw, cpm_table=cpm_table))
}

#' @importFrom stats quantile formula
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq counts
filter_by_variance <- function(count_matrix, q_filter, target){

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_matrix,
                                  colData = target,
                                  design = stats::formula("~ treat"))
  dds <- DESeq2::DESeq(dds)
  normalized_counts <- DESeq2::counts(dds, normalized=TRUE)
  variances <- matrixStats::rowVars(normalized_counts)
  threshold <- stats::quantile(variances, q_filter)
  fil_count_mtrx <- count_matrix[variances >= threshold,]
  return(list(fil_count_mtrx = fil_count_mtrx, 
    variance_dis = variances, 
     thr = q_filter,
     deseq2_normalized_counts = normalized_counts))
}

prepare_model_text <- function(model_variables, 
                               multifactorial){
    # Prepare model text
    if(model_variables != "") {
      model_variables_unlist <- unlist(strsplit(model_variables, ","))
      model_formula_text <- paste("~", paste(model_variables_unlist, 
                                  "+", collapse=" "), "treat")
    } else if (multifactorial != "") {
      mf_text <- split_mf_text(multifactorial)
      mf_factor_A <- mf_text[[1]]
      mf_factor_B <- mf_text[[2]]
      model_formula_text <- paste0("~ ", mf_factor_B, " + ", 
                         mf_factor_A, " + ", mf_factor_B,":", mf_factor_A)
    } else {
      model_formula_text <- "~ treat"
    }
    cat("Model for gene expression analysis is:", model_formula_text, "\n")
    return(model_formula_text)
}

split_mf_text <- function(multifactorial) {
      factors_contrast <- strsplit(multifactorial,":")[[1]]
      factors <- strsplit(factors_contrast[1], ",")[[1]]
      mf_factorA <- factors[1]
      mf_factorB <- factors[2]
      contrast_varA_varB <- strsplit(factors_contrast[2],",")[[1]]
      mf_contrast <- contrast_varA_varB[1]
      mf_varA <- contrast_varA_varB[2]
      mf_varB <- contrast_varA_varB[3]
      return(list(mf_factorA=mf_factorA, mf_factorB=mf_factorB, 
        mf_contrast=mf_contrast, mf_varA=mf_varA, mf_varB=mf_varB))
}

prepare_target_for_multifactorial <- function(target, multifactorial) {
  mf_text <- split_mf_text(multifactorial)
  target[, mf_text[["mf_factorA"]]] <- stats::relevel(
                                             target[,mf_text[["mf_factorA"]]],
                                             mf_text[["mf_varA"]])
  target[, mf_text[["mf_factorB"]]] <- stats::relevel(
                                             target[,mf_text[["mf_factorB"]]],
                                             mf_text[["mf_varB"]])
  return(target)
}

get_counts <- function(cnts_mtx, library_sizes)
{
    if (!is.null(library_sizes)){

        total_counts <- library_sizes[, c("sample","initial_total_sequences")]
        # Total reads might have been counted without taking into account ExpHunterSuite blacklist,
        # which would lead to errors. This next line removes blacklisted samples from total reads table.
        # Also, we have no guarantee they are sorted the same way, so we do it ourselves.
        total_counts <- total_counts[match(colnames(cnts_mtx), total_counts$sample), ]
        # sample_ID column no longer needed
        total_counts <- total_counts$initial_total_sequences
        gene_counts <- colSums(cnts_mtx)
        counted_frac <- gene_counts/total_counts

    } else {
        total_counts <- colSums(cnts_mtx)
        counted_frac <- NULL
    }

    coverage_df <- data.frame(sample_ID = colnames(cnts_mtx),
                                total_counts = total_counts)
    coverage_df$counted_frac <- counted_frac

    coverage_df <- coverage_df[order(coverage_df$total_counts),]
    coverage_df$sample_rank <- seq(1, nrow(coverage_df))

    return(coverage_df)
}

get_mean_counts <- function(cnts_mtx, cpm_table, reads, minlibraries)
{
    passed_filter <- rowSums(cpm_table >= reads) >= minlibraries
    min_1_read <- rowSums(cnts_mtx >= 1) >= minlibraries
    min_10_reads <- rowSums(cnts_mtx >= 10) >= minlibraries

    # Create a vector with the strictest filter passed
    all <- data.frame(counts=rowMeans(cnts_mtx),
                      filter=rep("all", nrow(cnts_mtx)))
    min_1_read <- data.frame(counts=rowMeans(cnts_mtx[min_1_read,]),
                              filter=rep("min_1_read", nrow(cnts_mtx[min_1_read,]))) 
    min_10_reads <- data.frame(counts=rowMeans(cnts_mtx[min_10_reads,]),
                              filter=rep("min_10_reads", nrow(cnts_mtx[min_10_reads,]))) 
    passed_filter <- data.frame(counts=rowMeans(cnts_mtx[passed_filter,]),
                                filter=rep("passed_filter", nrow(cnts_mtx[passed_filter,])))

    res <- rbind(all,min_1_read,min_10_reads,passed_filter)
    return(res)
}

get_gene_stats <- function(cpm_table, reads)
{
    cutoff_passed_matrix <- cpm_table >= reads
    cutoff_passed_matrix <- cutoff_passed_matrix[rowSums(cutoff_passed_matrix) > 0, ]
    exp_genes_df <- data.frame(sample_ID = colnames(cutoff_passed_matrix),
                                expressed_genes = colSums(cutoff_passed_matrix))
    exp_genes_df <- exp_genes_df[order(exp_genes_df$expressed_genes),]
    exp_genes_df$count_rank <- seq(1,nrow(exp_genes_df))
    cutoff_passed_matrix <- cutoff_passed_matrix[, exp_genes_df$sample_ID]
    union_vec = inters_vec <- cutoff_passed_matrix[,1]
    for (sample in 1:ncol(cutoff_passed_matrix))
    {
        union_vec <- union_vec | cutoff_passed_matrix[,sample]
        inters_vec <- inters_vec & cutoff_passed_matrix[,sample]
        exp_genes_df$union_expressed_genes[sample] <- sum(union_vec)
        exp_genes_df$inters_expressed_genes[sample] <- sum(inters_vec)
    }
    return(exp_genes_df)
}