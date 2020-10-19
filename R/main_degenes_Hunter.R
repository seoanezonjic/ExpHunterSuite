#! /usr/bin/env Rscript

main_degenes_Hunter <- function(
    input_file,
    Control_columns,
    Treatment_columns,
    reads,
    minlibraries,
    filter_type,
    output_files,
    p_val_cutoff,
    lfc,
    modules,
    minpack_common,
    target_file,
    external_DEA_file,
    model_variables,
    numerics_as_factors,
    custom_model,
    string_factors,
    numeric_factors,
    WGCNA_memory,
    WGCNA_norm_method,
    WGCNA_deepsplit,
    WGCNA_min_genes_cluster,
    WGCNA_detectcutHeight,
    WGCNA_mergecutHeight,
    WGCNA_all,
    WGCNA_blockwiseNetworkType,
    WGCNA_blockwiseTOMType,
    debug,
    Debug,
    opt,
    template_folder=file.path(find.package('DEgenesHunter'), 'templates')
  ){
    ############################################################
    ##                       CHECK INPUTS                     ##
    ############################################################

    # Check inputs not allowed values
    if (is.null(input_file)){
      stop(cat("No file with RNA-seq counts is provided.\nPlease use -i to submit file"))
    }
    if (minlibraries < 1){
      stop(cat("Minimum library number to check minimum read counts cannot be less than 1.\nIf you want to avoid filtering, set --reads to 0."))
    }
    if (reads == 0){
      cat("Minimum reads option is set to zero. Raw count table will not be filtered.")
    }

    if (!is.null(external_DEA_file)) {
      modules <- paste0(modules, "F")
      warning("External DEA file given. Note that if the count file corresponding to the DEA is non-integer, many DE detection packages will fail.")
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
        message("WGCNA only for controls an tratments is not activated (--WGCNA_all option) and is needed for PCIT analysis. Disabling PCIT execution")
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

    # Check either C and T columns or target file.
    if( (is.null(Treatment_columns) | is.null(Control_columns)) & is.null(target_file)) {
      stop(cat("You must include either the names of the control and treatment columns or a target file with a treat column."))
    }
    # In the case of -C/-T AND -t target - give a warning
    if( (!is.null(Treatment_columns) | !is.null(Control_columns)) & !is.null(target_file)) {
      warning("You have included at least one -C/-T option as well as a -t option for a target file. The target file will take precedence for assigning samples labels as treatment or control.")
    }
    # If no -t check there is no -v value
    if( is.null(target_file) & model_variables != "") {
      stop(cat("You should not include a -v value if you do not have a target table file."))
    }
    if(custom_model == TRUE & model_variables == "") {
      stop(cat("If you wish to use a custom model you must provide a value for the model_variables option."))
    }
    # If factors are specified but WGCNA not selected, throw a warning.
    if((string_factors != "" | numeric_factors != "") & (!grepl("W", modules) | is.null(target_file))) {
      warning("If you wish to use factors for the correlation analysis you must also run WGCNA and include a target file. The -S and -N options will be ignored")
    }

    if(!is.null(Debug)){
      debug <- TRUE
      debug_file <- file.path(Debug)
    }

    if(debug){
      # Define only once
      if(is.null(Debug)){
        debug_file <- file.path(paths$root, debug_files, paste(c("FHunter_Debug_Session_",format(Sys.Date(),format = "%Y%m%d"),".RData"),collapse = ""))
      }
      debug_dir <- dirname(debug_file)
      dir.create(debug_dir, recursive = T)
      debug_dir <- normalizePath(debug_dir)
      # Store session
      time_control <- list(start = Sys.time())
      debug_point(debug_file,"Start point")
    }

    ############################################################
    ##                         I/O DATA                       ##
    ############################################################

    # Create output folder 
    dir.create(output_files)

    # Load target file if it exists, otherwise use the -C and -T flags. Note target takes precedence over target.
    if(! is.null(target_file)) {
      target <- read.table(target_file, header=TRUE, sep="\t", check.names = FALSE)
      # Check there is a column named treat
      if(! "treat" %in% colnames(target)) {
        stop(cat("No column named treat in the target file.\nPlease resubmit"))
      }
      index_control_cols <- as.character(subset(target, treat == "Ctrl", select = sample, drop=TRUE))
      index_treatmn_cols <- as.character(subset(target, treat == "Treat", select = sample, drop=TRUE))
      replicatesC <- length(index_control_cols)
      replicatesT <- length(index_treatmn_cols)
    } else {
      index_control_cols <- unlist(strsplit(Control_columns, ",")) 
      index_treatmn_cols <- unlist(strsplit(Treatment_columns, ","))
      replicatesC <- length(index_control_cols)
      replicatesT <- length(index_treatmn_cols)
      # Create the target data frame needed in the calls to the DE detection methods
      target <- data.frame(sample=c(index_control_cols, index_treatmn_cols), 
        treat=c(rep("Ctrl", length(index_control_cols)), rep("Treat", length(index_treatmn_cols))))
    }

    # FOR WGCNA: Check that the appropriate factor columns can be found in the target file and makes a data frame with the specified factor
    if(exists("target") & grepl("W", modules)) {
      if(string_factors == "") {
        string_factors <- "treat"
      } else {
        string_factors <- paste0("treat,", string_factors)
      }
      string_factors <- unlist(strsplit(string_factors, ","))
      string_factors <- unique(string_factors) # In case we duplicate treat

      if(! all(string_factors %in% colnames(target))) {
        warning("Some factors specified with the --string_factors option cannot be found in the target file.")
      }
      string_factors_index <- colnames(target) %in% string_factors 

      target_string_factors <- target[string_factors_index]
      target_string_factors <-  data.frame(sapply(target_string_factors, as.factor), stringsAsFactors=TRUE)

      if(numeric_factors != "") {
        numeric_factors <- unlist(strsplit(numeric_factors, ","))
        if(! all(numeric_factors %in% colnames(target))) {
          warning("Some factors specified with the --numeric_factors option cannot be found in the target file.")
        }

        numeric_factors_index <- colnames(target) %in% numeric_factors
        if(TRUE %in% numeric_factors_index) {
          target_numeric_factors <- target[numeric_factors_index]
          # Ensure the factors are numeric
          #invisible(lapply(seq(ncol(target_numeric_factors)), function(i){target_numeric_factors[,i] <<- as.numeric(target_numeric_factors[,i])}))
        } else {
          stop(cat("Factors specified with the --numeric_factors option cannot be found in the target file.\nPlease resubmit."))
        }
      } else {
        target_numeric_factors <- ""
      }
    }

    # Now coerce the targets to factors for the multifactorial analysis
    if(numerics_as_factors == TRUE) {
      target <-  data.frame(sapply(target, as.factor), stringsAsFactors=TRUE)
    }

    replicatesC <- length(index_control_cols)
    replicatesT <- length(index_treatmn_cols)

    # Check if there are enough replicates for specified method
    if((replicatesC < 2) | (replicatesT < 2)) stop('At least two replicates per class (i.e. treatment and control) are required\n')

    # Load raw count data
    raw <- read.table(input_file, header=TRUE, row.names=1, sep="\t")
    raw <- raw[c(index_control_cols,index_treatmn_cols)] #Indexing selected columns from input count file

    # Substitute NA values
    raw[is.na(raw)] <- 0

    # Open the external data file for pre-analysed deg data
    if (grepl("F",modules)) {
      external_DEA_data <- read.table(external_DEA_file, header=TRUE, sep="\t", row.names=1)
    }


    ############################################################
    ##                          FILTER                        ##
    ############################################################
    # Prepare filtered set
    if(reads != 0){
      if (filter_type == "separate") {
        # genes with cpm greater than --reads value for at least --minlibrariess samples for either case or control samples
        to_keep_control <- rowSums(edgeR::cpm(raw[index_control_cols]) > reads) >= minlibraries
        to_keep_treatment <- rowSums(edgeR::cpm(raw[index_treatmn_cols]) > reads) >= minlibraries
        keep_cpm <- to_keep_control | to_keep_treatment # Keep if at least minlibraries in either of them
        raw_filter <- raw[keep_cpm,] # Filter out count data frame
      } else if (filter_type == "global") {
        keep_cpm <- rowSums(edgeR::cpm(raw) > reads) >= minlibraries # genes with cpm greater than --reads value for at least --minlibrariess samples
        raw_filter <- raw[keep_cpm,] # Filter out count data frame
      } else {
        warning("Unrecognized minimum read filter type. No filter will be used")
        raw_filter <- raw
      }
    } else { # NO FILTER
      raw_filter <- raw
    }
    # Defining contrast - Experimental design
    design_vector <- c(rep("C", replicatesC), rep("T", replicatesT))


    ############################################################
    ##                      MODEL TEXT                        ##
    ############################################################

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


    ####################### DEBUG POINT #############################
    if(debug){
      time_control$initialize <- Sys.time()
      debug_point(debug_file,"Initialize performed")
    }
    #################################################################


    ############################################################
    ##                       D.EXP ANALYSIS                   ##
    ############################################################

    # Prepare results containers
    all_data_normalized     <- list()
    all_counts_for_plotting <- list()
    all_package_results     <- list()
    package_objects         <- list()
    all_FDR_names           <- c()
    all_LFC_names           <- c()
    all_pvalue_names        <- c()
    final_logFC_names       <- c()
    final_FDR_names         <- c()
    final_pvalue_names      <- c()
    DEG_pack_columns        <- c()
    lfc <- lfc

    #####
    ################## CASE D: DESeq2
    #####
    if(grepl("D",modules)){
      if(replicatesC >= 2 & replicatesT >= 2) {
        # Verbose
        cat('Gene expression analysis is performed with DESeq2.\n')
        dir.create(file.path(output_files, "Results_DESeq2"))
        # Calculate results
        results <- analysis_DESeq2(data   = raw_filter,
                                   p_val_cutoff = p_val_cutoff,
                                   target = target,
                                   model_formula_text = model_formula_text)
        # Store results
        all_data_normalized[['DESeq2']] <- results[[1]]
        all_counts_for_plotting[['DESeq2']] <- results[[2]]
        package_objects[['DESeq2']] <- results[[3]]

        # Result Plot Visualization
        if (!is.null(all_counts_for_plotting[['DESeq2']])){
          all_FDR_names      <- c(all_FDR_names, 'padj')
          all_LFC_names      <- c(all_LFC_names, 'log2FoldChange')
          all_pvalue_names   <- c(all_pvalue_names, 'pvalue')
          final_pvalue_names <- c(final_pvalue_names, 'pvalue_DESeq2')
          final_logFC_names  <- c(final_logFC_names, 'logFC_DESeq2')
          final_FDR_names    <- c(final_FDR_names, 'FDR_DESeq2')
          DEG_pack_columns   <- c(DEG_pack_columns, 'DESeq2_DEG')
        } 
      } else {
      warning("DESeq2 will not be performed due to too few replicates")
      }
    }

    #####
    ################## CASE E: edgeR (AT LEAST 2 REPLICLATES)
    #####
    if(grepl("E", modules)){ 
      if(replicatesC >= 2 & replicatesT >= 2) {
        # Verbose point
        cat('Gene expression analysis is performed with edgeR.\n')
        path <- file.path(output_files, "Results_edgeR")
        dir.create(path)
        # Calculate results
        results <- analysis_edgeR(data   = raw_filter,
                                  target = target,
                                  model_formula_text = model_formula_text)
        # Store results
        all_data_normalized[['edgeR']] <- results[[1]]
        all_counts_for_plotting[['edgeR']] <- results[[2]]
        package_objects[['edgeR']] <- results[[3]]

        # Result Plot Visualization
        if (!is.null(all_counts_for_plotting[['edgeR']])){
          all_FDR_names      <- c(all_FDR_names, 'FDR')
          all_LFC_names      <- c(all_LFC_names, 'logFC')
          all_pvalue_names   <- c(all_pvalue_names, 'PValue')
          final_pvalue_names <- c(final_pvalue_names, 'pvalue_edgeR')
          final_logFC_names  <- c(final_logFC_names, 'logFC_edgeR')
          final_FDR_names    <- c(final_FDR_names, 'FDR_edgeR')
          DEG_pack_columns   <- c(DEG_pack_columns, 'edgeR_DEG')
        }
      } else {
      warning("edgeR will not be performed due to too few replicates")
      }
    }

    #####
    ################## CASE L: limma (AT LEAST 3 REPLICATES)
    #####
    if(grepl("L", modules)){ 
      if(replicatesC >= 3 & replicatesT >= 3) {
        # Verbose
        cat('Gene expression analysis is performed with limma.\n')
        dir.create(file.path(output_files, "Results_limma"))

        # Calculate results
        results <- analysis_limma(data   = raw_filter,
                                  target = target,
                                  model_formula_text = model_formula_text)
        # Store results
        all_data_normalized[['limma']]     <- results[[1]]
        all_counts_for_plotting[['limma']] <- results[[2]]

        # Result Plot Visualization
        if (!is.null(all_counts_for_plotting[['limma']])){
          all_FDR_names      <- c(all_FDR_names, 'adj.P.Val')
          all_LFC_names      <- c(all_LFC_names, 'logFC')
          all_pvalue_names   <- c(all_pvalue_names, 'P.Value')
          final_pvalue_names <- c(final_pvalue_names, 'pvalue_limma')
          final_logFC_names  <- c(final_logFC_names, 'logFC_limma')
          final_FDR_names    <- c(final_FDR_names, 'FDR_limma')
          DEG_pack_columns   <- c(DEG_pack_columns, 'limma_DEG')
        }
      } else {
        warning("limma will not be performed due to too few replicates")
      }
    }

    #####
    ################## CASE N: NOISeq (AT LEAST 3 REPLICLATES)
    #####
    if(grepl("N", modules)){ 
        if(replicatesC >= 3 & replicatesT >= 3) {

        # Verbose
        cat("Gene expression analysis is performed with NOISeqBIO function within NOISeq.\n")
        path <- file.path(output_files, "Results_NOISeq")
        dir.create(path)
        # Calculate results
        results <- analysis_NOISeq(data    = raw_filter, 
                                   target = target)
        # Store results
        all_data_normalized[['NOISeq']]     <- results[[1]]
        all_counts_for_plotting[['NOISeq']] <- results[[2]]
        package_objects[['NOISeq']] <- results[[3]]

        #Result Plot Visualization
        if (!is.null(all_counts_for_plotting[['NOISeq']])){
          all_FDR_names      <- c(all_FDR_names, 'adj.p')
          all_LFC_names      <- c(all_LFC_names, 'log2FC')
          all_pvalue_names   <- c(all_pvalue_names, 'prob')
          final_pvalue_names <- c(final_pvalue_names, 'pvalue_NOISeq')
          final_logFC_names  <- c(final_logFC_names, 'logFC_NOISeq')
          final_FDR_names    <- c(final_FDR_names, 'FDR_NOISeq')
          DEG_pack_columns   <- c(DEG_pack_columns, 'NOISeq_DEG')
        }
      } else {
        warning("NOISeq will not be performed due to too few replicates")
      }
    }

    #####
    ################## CASE F: External DEG file
    #####
    if(grepl("F",modules)){
      cat('External DEA data to be processed \n')
      dir.create(file.path(output_files, "Results_external_DEA"))
      # Calculate results
      results <- analysis_external_DEA(raw_filter, external_DEA_data)
      # Store results
      all_data_normalized[['external_DEA']] <- results[[1]]
      all_counts_for_plotting[['external_DEA']] <- results[[2]]
      package_objects[['external_DEA']] <- results[[3]]

      # Result Plot Visualization
      if (!is.null(all_counts_for_plotting[['external_DEA']])){
        all_pvalue_names   <- c(all_pvalue_names, names(external_DEA_data)[1])
        all_LFC_names      <- c(all_LFC_names, names(external_DEA_data)[2])
        all_FDR_names      <- c(all_FDR_names, names(external_DEA_data)[3])
        final_pvalue_names <- c(final_pvalue_names, 'pvalue_external_DEA')
        final_logFC_names  <- c(final_logFC_names, 'logFC_external_DEA')
        final_FDR_names    <- c(final_FDR_names, 'FDR_external_DEA')
        DEG_pack_columns   <- c(DEG_pack_columns, 'external_DEA_DEG')
      }

    }


    ####################### DEBUG POINT #############################
    if(debug){
      time_control$dea_packages <- Sys.time()
      debug_point(debug_file,"DEA analysis performed")
    }
    ##################################################################

    #################################################################
    ##                       CORRELATION ANALYSIS                   ##
    ##################################################################
    #####
    ################## CASE W: WGCNA
    #####

    if(grepl("W", modules)) {
      cat('Correlation analysis is performed with WGCNA\n')
      path <- file.path(output_files, "Results_WGCNA")
      dir.create(path)

      if(WGCNA_norm_method %in% names(all_data_normalized)) {
        WGCNA_input <- all_data_normalized[[WGCNA_norm_method]]
      } else if(WGCNA_norm_method== "none") {
        WGCNA_input <- raw_filter
      } else {
        warning("To run WGCNA, you must also run the method chosen for normalization in the --modules flag, or specify none")
      }

      WGCNA_input_treatment <- WGCNA_input[, index_treatmn_cols]
      WGCNA_input_control <- WGCNA_input[, index_control_cols]

      if(WGCNA_all == TRUE) {
        WGCNA_treatment_path <- file.path(path, "Treatment_only_data")
        dir.create(WGCNA_treatment_path)
        cat('Performing WGCNA correlation analysis for treated samples\n')
        results_WGCNA_treatment <- analysis_WGCNA(data=WGCNA_input_treatment,
                       path=WGCNA_treatment_path,
                       target_numeric_factors=target_numeric_factors,
                       target_string_factors=target_string_factors,
                       WGCNA_memory=WGCNA_memory,
                       WGCNA_deepsplit=WGCNA_deepsplit,
                       WGCNA_detectcutHeight=WGCNA_detectcutHeight,
                       WGCNA_mergecutHeight=WGCNA_mergecutHeight,
                       WGCNA_min_genes_cluster=WGCNA_min_genes_cluster,
                       cor_only=TRUE, 
                       blockwiseNetworkType = WGCNA_blockwiseNetworkType, 
                       blockwiseTOMType = WGCNA_blockwiseTOMType
        )

        WGCNA_control_path <- file.path(path, "Control_only_data")
        dir.create(WGCNA_control_path)

        cat('Performing WGCNA correlation analysis for control samples\n')
        results_WGCNA_control <- analysis_WGCNA(data=WGCNA_input_control,
                       path=WGCNA_control_path,
                       target_numeric_factors=target_numeric_factors,
                       target_string_factors=target_string_factors,
                       WGCNA_memory=WGCNA_memory,
                       WGCNA_deepsplit=WGCNA_deepsplit,
                       WGCNA_detectcutHeight=WGCNA_detectcutHeight,
                       WGCNA_mergecutHeight=WGCNA_mergecutHeight,
                       WGCNA_min_genes_cluster=WGCNA_min_genes_cluster,                    
                       cor_only=TRUE, 
                       blockwiseNetworkType = WGCNA_blockwiseNetworkType, 
                       blockwiseTOMType = WGCNA_blockwiseTOMType
        )
      }
        # Need to improve the control, probably by removing PCIT
        # if(results_WGCNA_treatment == "NO_POWER_VALUE" | results_WGCNA_control == "NO_POWER_VALUE") {
        #   warning("WGCNA was unable to generate a suitable power value for at least one of the partial datasets")
        # }
      cat('Performing WGCNA correlation analysis for all samples\n')
      results_WGCNA <- analysis_WGCNA(data=WGCNA_input,
                                     path=path,
                                     target_numeric_factors=target_numeric_factors,
                                     target_string_factors=target_string_factors,
                                     WGCNA_memory=WGCNA_memory,
                                     WGCNA_deepsplit=WGCNA_deepsplit,
                                     WGCNA_detectcutHeight=WGCNA_detectcutHeight,
                                     WGCNA_mergecutHeight=WGCNA_mergecutHeight,
                                     WGCNA_min_genes_cluster=WGCNA_min_genes_cluster,
                                     cor_only=FALSE, 
                                     blockwiseNetworkType = WGCNA_blockwiseNetworkType, 
                                     blockwiseTOMType = WGCNA_blockwiseTOMType
      )
      if(length(results_WGCNA) == 1) {
        warning("Something went wrong with WGCNA on the full dataset")
        modules <- gsub("W", "", modules)
      }
    ####################### DEBUG POINT #############################
      if(debug){
        time_control$wgcna <- Sys.time()
        debug_point(debug_file,"WGCNA analysis performed")
      }
    #################################################################
    }

    #####
    ################## CASE X: diffcoexp
    #####
    if(grepl("X", modules)) {
      cat('Correlation analysis is performed with diffcoexp\n')
      path <- file.path(output_files, "Results_diffcoexp")
      dir.create(path)

      results_diffcoexp <- analysis_diffcoexp(data = raw_filter,
                                           path = path,
                                           target = target)
    ####################### DEBUG POINT #############################
      if(debug){
        time_control$diffcoexp <- Sys.time()
        debug_point(debug_file,"Diffcoexp analysis performed")
    #################################################################
      }
    }

    #################################################################################
    ##                       EXPORT FINAL RESULTS AND OTHER FILES                  ##
    #################################################################################

    # Write filtered count data
    write.table(raw_filter, file=file.path(output_files, "filtered_count_data.txt"), quote=FALSE, col.names=NA, sep="\t")

    # Write Counts data and all DE for each  design to file as 2 column data frame
    write.table(data.frame(class = design_vector, name = c(index_control_cols, index_treatmn_cols)), 
      file=file.path(output_files, "control_treatment.txt"), row.names=FALSE, quote=FALSE, sep="\t")
    write_df_list_as_tables(all_data_normalized, prefix = 'Normalized_counts_', root = output_files)
    write_df_list_as_tables(all_counts_for_plotting, prefix = 'allgenes_', root = output_files)

    # Write main results file, normalized counts per package and DE results per package
    DE_all_genes <- unite_DEG_pack_results(DEG_pack_columns, all_counts_for_plotting, all_FDR_names, all_LFC_names, all_pvalue_names, final_pvalue_names, 
      final_logFC_names, final_FDR_names, p_val_cutoff, lfc, minpack_common)

    #####
    ################## CASE P: PCIT
    #####

    if(grepl("P", modules)) {
      metrics_pcit <- analysis_diff_correlation(DE_all_genes, DESeq2_counts, DESeq2_counts_control, DESeq2_counts_treatment, PCIT_filter=FALSE)
      DE_all_genes <- transform(merge(DE_all_genes, metrics_pcit, by.x=0, by.y=0), row.names=Row.names, Row.names=NULL)
    }

    #################################################################################

    # Check WGCNA was run and it returned proper results
    if(grepl("W", modules)) {
      DE_all_genes <- transform(merge(DE_all_genes, results_WGCNA[['gene_cluster_info']], by.x=0, by.y="ENSEMBL_ID"), row.names=Row.names, Row.names=NULL)
    }

    if(grepl("X", modules)) {
      #results_diffcoexp
      DE_all_genes$DCG <- row.names(DE_all_genes) %in% results_diffcoexp$DCGs$Gene[results_diffcoexp$DCGs$q < 1]
      DE_all_genes$DCL <- row.names(DE_all_genes) %in% results_diffcoexp$DCLs$Gene.1 | row.names(DE_all_genes) %in% results_diffcoexp$DCLs$Gene.2
    }
    # Add the filtered genes back
    DE_all_genes <- add_filtered_genes(DE_all_genes, raw)

    # New structure - row names are now actually row names
    dir.create(file.path(output_files, "Common_results"))
    write.table(DE_all_genes, file=file.path(output_files, "Common_results", "hunter_results_table.txt"), quote=FALSE, row.names=TRUE, sep="\t")
    ####################### DEBUG POINT #############################
    if(debug){
      time_control$write_output_tables <- Sys.time()
      debug_point(debug_file,"Full analysis performed")
    }
    #################################################################


    ############################################################
    ##                    GENERATE REPORT                     ##
    ############################################################
    outf <- file.path(normalizePath(output_files),"DEG_report.html")
    rmarkdown::render(file.path(template_folder, 'main_report.Rmd'), 
                      output_file = outf, intermediates_dir = output_files)


    if(debug){
    ####################### DEBUG POINT #############################
        time_control$render_main_report <- Sys.time()
        debug_point(debug_file,"Report printed")
    #################################################################
        save_times(time_control, output = file.path(debug_dir, "degenes_hunter_time_control.txt"), plot_name = "DH_times_control.pdf")
    }

}