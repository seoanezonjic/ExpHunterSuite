#' Main DEgenes Hunter Function
#'
#' This function allows you to perform the DEG analysis with different algorithms.
#' @param raw Dataframe with raw counts per sample and feature
#' @param target set with control and treatment information
#' @param external_DEA_data external DEA set 
#' @param output_files otput files path
#' @param reads set of reads
#' @param minlibraries minimum of libraries to use a set of data
#' @param filter_type filter type to ba applied. Allowed: "separate" or "global"
#' @param p_val_cutoff p-value threshold
#' @param lfc minimum logFold Change
#' @param modules modules to be executed. Allowed: DESeq2 (D), NOISeq (N), Limma (L), EdgeR (E), WGCNA (W)
#' @param minpack_common minimum of pack that must be significant to tag a gene as significant
#' @param model_variables custom model
#' @param numerics_as_factors transform numeric values to factors. Default: TRUE
#' @param custom_model boolean of usage of custom model. Default: FALSE
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
#' @return expression analysis result object with studies performed
#' @keywords method
#' @export
#' @examples
#' data(toc)
#' data(target)
#' de_out <- main_degenes_Hunter(raw=toc, target=target, modules="D")

main_degenes_Hunter <- function(
    raw = NULL,
    target = NULL,
    external_DEA_data = NULL, # TODO unused variable, remove from interface and affected scripts
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
    custom_model = FALSE,
    string_factors = "",
    numeric_factors = "",
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
    WGCNA_minKMEtoStay = 0.2
  ){


    modified_input_args <- check_input_main_degenes_Hunter(minlibraries, reads, external_DEA_data, modules, model_variables, active_modules, WGCNA_all, minpack_common, target, custom_model, string_factors, numeric_factors)
    modules <- modified_input_args[['modules']]
    active_modules <- modified_input_args[['active_modules']]
    minpack_common <- modified_input_args[['minpack_common']]
    final_main_params <- list(
      'minpack_common' = minpack_common,
      'p_val_cutoff' = p_val_cutoff,
      'lfc' = lfc,
      'modules' = modules
    )

    ############################################################
    ##                         I/O DATA                       ##
    ############################################################

    # Infer replicates and group index from target
    # index_control_cols <- as.character(subset(target, treat == "Ctrl", select = sample, drop=TRUE))
    index_control_cols <- as.character(target$sample[target$treat == "Ctrl"])
    # index_treatmn_cols <- as.character(subset(target, treat == "Treat", select = sample, drop=TRUE))
    index_treatmn_cols <- as.character(target$sample[target$treat == "Treat"])
    replicatesC <- length(index_control_cols)
    replicatesT <- length(index_treatmn_cols)
    design_vector <- c(rep("C", replicatesC), rep("T", replicatesT))
    sample_groups <- data.frame(class = design_vector, name = c(index_control_cols, index_treatmn_cols))   

    # Check if there are enough replicates for specified method
    if((replicatesC < 2) | (replicatesT < 2)) stop('At least two replicates per class (i.e. treatment and control) are required\n')

    if(!is.null(target) & grepl("W", modules)) {
      target_numeric_factors <- build_design_for_WGCNA(target, numeric_factors=numeric_factors)
      target_string_factors <- build_design_for_WGCNA(target, string_factors=string_factors)
    }

    if(numerics_as_factors == TRUE) { # Now coerce the targets to factors for the multifactorial analysis
      target <-  data.frame(sapply(target, as.factor), stringsAsFactors=TRUE)
    }

    model_formula_text <- prepare_model_text(model_variables, custom_model)
   
    # Prepare count table for analysis
    raw <- raw[c(index_control_cols,index_treatmn_cols)] #Indexing selected columns from input count dataframe   
    raw[is.na(raw)] <- 0 # Substitute NA values

    raw_filter <- filter_count(reads, minlibraries, raw, filter_type, index_control_cols, index_treatmn_cols)
    
    ############################################################
    ##             PERFORM EXPRESION ANALYSIS                 ##
    ############################################################
    dir.create(output_files)
   
    exp_results <- perform_expression_analysis(modules, replicatesC, replicatesT, raw_filter, p_val_cutoff, target, model_formula_text)

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
        warning("To run WGCNA, you must also run the method chosen for normalization in the --modules flag, or specify none")
      }

      combinations_WGCNA <- perform_WGCNA_combinations(WGCNA_all, WGCNA_input, index_treatmn_cols, index_control_cols, path, 
          target_numeric_factors, target_string_factors, WGCNA_memory, WGCNA_deepsplit, WGCNA_detectcutHeight, WGCNA_mergecutHeight, 
          WGCNA_min_genes_cluster, WGCNA_blockwiseNetworkType, WGCNA_blockwiseTOMType, WGCNA_minCoreKME, WGCNA_minCoreKMESize, WGCNA_minKMEtoStay)

      if(length(combinations_WGCNA[['WGCNA_all']]) == 1) {
        warning("Something went wrong with WGCNA on the full dataset")
        modules <- gsub("W", "", modules)
      }
    }

    if(grepl("X", modules)) { # CASE X: diffcoexp
      cat('Correlation analysis is performed with diffcoexp\n')
      path <- file.path(output_files, "Results_diffcoexp")
      dir.create(path)

      results_diffcoexp <- analysis_diffcoexp(data = raw_filter,
                                           path = path,
                                           target = target)
    }


    #################################################################################
    ##                       BUILD MAIN RESULT TABLE                  ##
    #################################################################################
    DE_all_genes <- unite_DEG_pack_results(exp_results, p_val_cutoff, lfc, minpack_common)
    
    # if(grepl("P", modules)) { # CASE P: PCIT, TODO: RESTORE FUNCTION, PEDRO 
    #   # TODO : This is not working, variables "DESeq2_counts" are not being generated inside this function
    #   metrics_pcit <- analysis_diff_correlation(DE_all_genes, DESeq2_counts, DESeq2_counts_control, DESeq2_counts_treatment, PCIT_filter=FALSE)
    #   DE_all_genes <- transform(merge(DE_all_genes, metrics_pcit, by.x=0, by.y=0), row.names=Row.names, Row.names=NULL)
    # }
    
    if(grepl("W", modules)) { # Check WGCNA was run and it returned proper results
      DE_all_genes <- transform(merge(DE_all_genes, combinations_WGCNA[['WGCNA_all']][['gene_cluster_info']], by.x=0, by.y="ENSEMBL_ID"), row.names=Row.names, Row.names=NULL)
    }

    if(grepl("X", modules)) { #results_diffcoexp    
      DE_all_genes$DCG <- row.names(DE_all_genes) %in% results_diffcoexp$DCGs$Gene[results_diffcoexp$DCGs$q < 1]
      DE_all_genes$DCL <- row.names(DE_all_genes) %in% results_diffcoexp$DCLs$Gene.1 | row.names(DE_all_genes) %in% results_diffcoexp$DCLs$Gene.2
    }
    
    DE_all_genes <- add_filtered_genes(DE_all_genes, raw) # Add the filtered genes back

    final_results <- list()
    final_results[['raw_filter']] <- raw_filter
    final_results[['sample_groups']] <- sample_groups
    final_results[['DE_all_genes']] <- DE_all_genes
    final_results[['index_control_cols']] <- index_control_cols
    final_results[['index_treatmn_cols']] <- index_treatmn_cols
    final_results[['design_vector']] <- design_vector
    final_results[['replicatesC']] <- replicatesC
    final_results[['replicatesT']] <- replicatesT
    final_results[['final_main_params']] <- final_main_params
    if(!is.null(combinations_WGCNA)){final_results <- c(final_results, combinations_WGCNA)}
    return(c(final_results, exp_results))
}



check_input_main_degenes_Hunter <- function(minlibraries, reads, external_DEA_data, modules, model_variables, active_modules, WGCNA_all, minpack_common, target, custom_model, string_factors, numeric_factors){
    if (minlibraries < 1){
      stop(cat("Minimum library number to check minimum read counts cannot be less than 1.\nIf you want to avoid filtering, set --reads to 0."))
    }

    if (reads == 0){
      cat("Minimum reads option is set to zero. Raw count table will not be filtered.")
    }

    if (!is.null(external_DEA_data)) {
      modules <- paste0(modules, "F")
      warning("External DEA dataframe given. Note that if the count table corresponding to the DEA is non-integer, many DE detection packages will fail.")
    }

    if (model_variables != "" & grepl("N", modules) & nchar(modules) == 1) {
      stop(cat("You cannot run an experimental design that uses the model_variables option with NOISeq only."))
    }
    if (model_variables != "" & grepl("N", modules) & nchar(modules) > 1) {
      warning("NOISeq will not be run as you have an experimental design that uses the model_variables option.")
      modules <- gsub("N", "", modules)
    }

    active_modules <- nchar(modules)
    if(grepl("W", modules)) {
      active_modules <- active_modules - 1
    }
    if(grepl("P", modules)) {
      if (WGCNA_all == FALSE) {
        message("WGCNA only for controls an treatments is not activated (--WGCNA_all option) and is needed for PCIT analysis. Disabling PCIT execution")
        modules <- gsub("P", "", modules)
      } else {
        active_modules <- active_modules - 1
      }
    }
    if(grepl("X", modules)) {
      active_modules <- active_modules - 1
    }

    if(minpack_common > active_modules){
      minpack_common <- active_modules
      warning("The number of active modules is lower than the thresold for tag PREVALENT DEG. The thresold is set to the number of active modules.")
    }


    # If no -t check there is no -v value
    if( is.null(target) & model_variables != "") {
      stop(cat("You should not include a -v value if you do not have a target table file."))
    }
    if(custom_model == TRUE & model_variables == "") {
      stop(cat("If you wish to use a custom model you must provide a value for the model_variables option."))
    }
    # If factors are specified but WGCNA not selected, throw a warning.
    if((string_factors != "" | numeric_factors != "") & (!grepl("W", modules) | is.null(target))) {
      warning("If you wish to use factors for the correlation analysis you must also run WGCNA and include a target file. The -S and -N options will be ignored")
    }
    return(list(modules=modules, active_modules=active_modules, minpack_common=minpack_common))
}

#' @importFrom edgeR cpm
filter_count <- function(reads, minlibraries, raw, filter_type, index_control_cols, index_treatmn_cols){
    # Prepare filtered set
    if(reads != 0){
      if (filter_type == "separate") {
        # genes with cpm greater than --reads value for at least --minlibrariess samples for either case or control samples
        to_keep_control <- rowSums(edgeR::cpm(raw[index_control_cols]) > reads) >= minlibraries
        to_keep_treatment <- rowSums(edgeR::cpm(raw[index_treatmn_cols]) > reads) >= minlibraries
        keep_cpm <- to_keep_control | to_keep_treatment # Keep if at least minlibraries in either of them
        raw <- raw[keep_cpm,] # Filter out count data frame
      } else if (filter_type == "global") {
        keep_cpm <- rowSums(edgeR::cpm(raw) > reads) >= minlibraries # genes with cpm greater than --reads value for at least --minlibrariess samples
        raw <- raw[keep_cpm,] # Filter out count data frame
      } else {
        warning("Unrecognized minimum read filter type. No filter will be used")
        raw <- raw
      }
    }
    return(raw)
}

prepare_model_text <- function(model_variables, custom_model=FALSE){
    # Prepare model text
    if(model_variables != "") {
      if(custom_model == TRUE) {
        model_formula_text <- model_variables
      } else {
        model_variables_unlist <- unlist(strsplit(model_variables, ","))
        model_formula_text <- paste("~", paste(model_variables_unlist, "+", collapse=" "), "treat")
      }
    } else {
      model_formula_text <- "~ treat"
    }
    cat("Model for gene expression analysis is:", model_formula_text, "\n")
    return(model_formula_text)
}