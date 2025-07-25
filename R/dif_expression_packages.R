################################################
### FUNCTIONS: DIFFERENCIAL EXPRESSION PACKAGES
################################################
# Manager function
#-----------------------------------------------
perform_expression_analysis <- function(modules, 
                                        replicatesC, 
                                        replicatesT, 
                                        raw_filter, 
                                        p_val_cutoff, 
                                        target, 
                                        model_formula_text, 
                                        external_DEA_data=NULL,
                                        multifactorial){
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

    #####
    ################## CASE D: DESeq2
    #####
    if(grepl("D",modules)){
      if(replicatesC >= 2 & replicatesT >= 2) {
        # Verbose
        message('Gene expression analysis is performed with DESeq2.\n')
        # Calculate results
        results <- analysis_DESeq2(data   = raw_filter,
                                   p_val_cutoff = p_val_cutoff,
                                   target = target,
                                   model_formula_text = model_formula_text,
                                   multifactorial = multifactorial)
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
        message('Gene expression analysis is performed with edgeR.\n')
        # Calculate results
        results <- analysis_edgeR(data   = raw_filter,
                                  target = target,
                                  model_formula_text = model_formula_text,
                                  multifactorial = multifactorial)
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
        message('Gene expression analysis is performed with limma.\n')

        # Calculate results
        results <- analysis_limma(data   = raw_filter,
                                  target = target,
                                  model_formula_text = model_formula_text,
                                  multifactorial = multifactorial)
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
        message("Gene expression analysis is performed with NOISeqBIO",
                   " function within NOISeq.")
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
      message('External DEA data to be processed')
      # Calculate results
      results <- analysis_external_DEA(raw_filter, external_DEA_data)
      # Store results
      all_data_normalized[['externalDEA']] <- results[[1]]
      all_counts_for_plotting[['externalDEA']] <- results[[2]]
      package_objects[['externalDEA']] <- results[[3]]

      # Result Plot Visualization

      if (!is.null(all_counts_for_plotting[['externalDEA']])){
        all_pvalue_names   <- c(all_pvalue_names, "P.Value")
        all_LFC_names      <- c(all_LFC_names, "logFC")
        all_FDR_names      <- c(all_FDR_names, "adj.P.Val")
        final_pvalue_names <- c(final_pvalue_names, 'pvalue_externalDEA')
        final_logFC_names  <- c(final_logFC_names, 'logFC_externalDEA')
        final_FDR_names    <- c(final_FDR_names, 'FDR_externalDEA')
        DEG_pack_columns   <- c(DEG_pack_columns, 'externalDEA_DEG')
      }

    }

    exp_results <- list()
    exp_results[['all_data_normalized']]     <- all_data_normalized
    exp_results[['all_counts_for_plotting']] <- all_counts_for_plotting
    exp_results[['all_package_results']]     <- all_package_results
    exp_results[['package_objects']]         <- package_objects
    exp_results[['all_FDR_names']]           <- all_FDR_names
    exp_results[['all_LFC_names']]           <- all_LFC_names
    exp_results[['all_pvalue_names']]        <- all_pvalue_names
    exp_results[['final_logFC_names']]       <- final_logFC_names
    exp_results[['final_FDR_names']]         <- final_FDR_names
    exp_results[['final_pvalue_names']]      <- final_pvalue_names
    exp_results[['DEG_pack_columns']]        <- DEG_pack_columns

    return(exp_results)
}


# DESeq2
#-----------------------------------------------
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results counts
#' @importFrom stats formula
analysis_DESeq2 <- function(data, p_val_cutoff, target, model_formula_text,
  multifactorial){
    if(grepl(":nested", multifactorial)){
      mf_text <- split_mf_text(multifactorial)
      factor_table <- table(target[, mf_text$mf_factorA])
      ordered_factors <- names(factor_table)[order(factor_table, decreasing=TRUE)]
      # Following necessary to perform design matrix trick in DESeq2 Vignette:
      # Group-specific condition effects, individuals nested within groups
      ordering_factor <- ordered(target[, mf_text$mf_factorA], 
        levels=ordered_factors)
      # Reorder table by group then by patient (paired samples, i.e. factor B)
      target <- target[order(ordering_factor,  target[, mf_text$mf_factorB]), ]
      bigger_grouping <- target[, mf_text$mf_factorA] == ordered_factors[1]
      smaller_grouping <- target[, mf_text$mf_factorA] == ordered_factors[2]
      target[smaller_grouping, mf_text$mf_factorB] <- 
                  target[bigger_grouping, mf_text$mf_factorB][seq(1, sum(smaller_grouping))]
      target <- target[order(as.integer(row.names(target))), ]

      model_formula_text <- paste0("~ ", mf_text$mf_factorA, " + ", 
                         mf_text$mf_factorA, ":", mf_text$mf_factorB, " + ", mf_text$mf_factorA,":treat")
      m1 <- model.matrix(stats::formula(model_formula_text), target)

      all.zero <- apply(m1, 2, function(x) all(x==0))
      idx <- which(all.zero);  
      m1 <- m1[,-idx]
  
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = data,
            colData = target,
            design = ~1)
      dds <- DESeq2::DESeq(dds, full=m1)

    } else {
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = data,
                                  colData = target,
                                  design = stats::formula(model_formula_text))
      dds <- DESeq2::DESeq(dds)
    }

    if(multifactorial != "") {
      mf_options <- get_mf_DE_options(package_name="DESeq2", package_object=dds,
                                      multifactorial=multifactorial, 
                                      target=target)
      all_DESeq2_genes <- do.call(DESeq2::results, 
                                c(object=dds, alpha = p_val_cutoff, mf_options))
    } else {
      all_DESeq2_genes <- DESeq2::results(dds, alpha = p_val_cutoff)
    }
    normalized_counts <- as.data.frame(DESeq2::counts(dds, normalized=TRUE))
    return(list(normalized_counts, 
           as.data.frame(all_DESeq2_genes), 
           list(de_deseq2=all_DESeq2_genes, DESeq2_dataset=dds)))
}

#' @importFrom edgeR DGEList calcNormFactors estimateDisp exactTest glmQLFit 
#' glmQLFTest topTags cpm
#' @importFrom stats formula model.matrix
analysis_edgeR <- function(data, target, model_formula_text, multifactorial){
    # This is the default model_formula_text if nothing else is given.
    if(model_formula_text == "~ treat") {
    # Building edgeR object
        d_edgeR <- edgeR::DGEList(counts=data, group=as.character(target$treat))
        d_edgeR <- edgeR::calcNormFactors(d_edgeR)
        d_edgeR <- edgeR::estimateDisp(d_edgeR)
        d_edgeR_DE <- edgeR::exactTest(d_edgeR, dispersion = "auto", 
                                   pair=c("Ctrl", "Treat"))
    } else {
        d_edgeR <- edgeR::DGEList(counts=data)
        model_design <- stats::model.matrix(
                            stats::formula(paste(model_formula_text)), 
                            target)
        d_edgeR <- edgeR::estimateDisp(d_edgeR, model_design)
        d_edgeR_DE <- edgeR::glmQLFit(d_edgeR, model_design)

      if(multifactorial != "") {
        mf_options <- get_mf_DE_options(package_name="edgeR", 
                                        package_object=d_edgeR_DE,
                                        multifactorial=multifactorial,
                                        target=target)
        d_edgeR_DE <- do.call(edgeR::glmQLFTest, 
                              list(glmfit=d_edgeR_DE, unlist(mf_options)))
      } else {
        d_edgeR_DE <- edgeR::glmQLFTest(d_edgeR_DE)
      }
    }
    all_genes_df <- edgeR::topTags(d_edgeR_DE, 
                                   nrow(d_edgeR), sort.by="none")$table
    dennorm <- edgeR::cpm(d_edgeR, log=TRUE)
    # Getting normalized counts
    normalized_counts <- edgeR::cpm(d_edgeR)
    # Return calculated info
    return(list(normalized_counts, all_genes_df, list(d_edgeR=d_edgeR, 
                                                    dennorm=dennorm)))
}

#' @importFrom limma voom eBayes topTable lmFit
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom stats formula model.matrix
analysis_limma <- function(data, target, model_formula_text, multifactorial){
    model_design <- stats::model.matrix(stats::formula(
                                         paste(model_formula_text)), 
                                      target)

    # Calculating differential expression
    DGE_List <- edgeR::DGEList(counts=data) # Building object (DGEList)
    # Calculation of the normalization factor
  DGE_List <- edgeR::calcNormFactors(DGE_List)  
  # Converting the read counts to log2-cpm
    log2_cpm <- limma::voom(DGE_List, model_design) 
    fit <- limma::lmFit(log2_cpm, model_design)
    fit2 <- limma::eBayes(fit)
    if(multifactorial != "") {
      mf_options <- get_mf_DE_options(package_name="limma", 
                                      package_object=model_design,
                                      multifactorial=multifactorial,
                                      target=target)
      todos_limma <- do.call(limma::topTable, list(fit=fit2, adjust.method="BH",
                                number=nrow(fit2), unlist(mf_options)))
    } else {
          todos_limma <- limma::topTable(fit2, adjust.method="BH",
                                         number=nrow(fit2), coef="treatTreat")
    }
    normalized_counts <- log2_cpm$E
  # Return calculated info
  return(list(normalized_counts, todos_limma))
}

#' @importFrom NOISeq readData noiseqbio tmm
#' @importFrom Biobase assayData
analysis_NOISeq <- function(data, target){
    group <- as.character(target$treat)
    # Experimental design
    myfactors <- data.frame(Tissue = rev(group), TissueRun = rev(group))
    # Calculation of differential expression
    mydata <- NOISeq::readData(data, myfactors) # Building NOISeq Object (eset)

    all_NOISeq <- NOISeq::noiseqbio(mydata, 
                                  k = 0.5, 
                                  norm = "tmm", 
                                  factor="Tissue", 
                                  lc = 1, 
                                  r = 50, 
                                  adj = 1.5, 
                                  plot = FALSE,    
                                  a0per = 0.9, 
                                  random.seed = 12345, 
                                  filter = 1, 
                                  cv.cutoff = 500, 
                                  cpm = 1)

    normalized_counts <- NOISeq::tmm(Biobase::assayData(mydata)$exprs, 
                               long = 1000, lc = 1) # Getting normalized counts
    expres_diff_all <- as.data.frame(all_NOISeq@results)

    prob_all <- expres_diff_all$prob
    # Inverse
    prob_all <- 1 - prob_all
    expres_diff_all <- data.frame(expres_diff_all, adj.p = prob_all)

    # Return calculated info
    return(list(normalized_counts, expres_diff_all, all_NOISeq))
}

get_mf_DE_options <- function(
  package_name, 
  package_object, 
  multifactorial, 
  target) {
  mf_text <- split_mf_text(multifactorial)
  
  if(package_name == "DESeq2") {
    if(mf_text[["mf_contrast"]] == "interaction") {
        mf_opt <- list(name=DESeq2::resultsNames(package_object)[4])
    } else if (mf_text[["mf_contrast"]] == "effect") {

        mf_opt <- DESeq2::resultsNames(package_object)[3]
        mf_opt <- list(contrast=c(
          mf_text[["mf_factorA"]],
          levels(target[, mf_text[["mf_factorA"]]])[2],
           mf_text[["mf_varA"]]
        ))
    } else if (mf_text[["mf_contrast"]] == "nested_int") {
        target_factor <- unique(target[,mf_text[["mf_factorA"]]])
        mf_opt <- list(contrast=list(paste0(mf_text[["mf_factorA"]], mf_text[["mf_varB"]], ".treatTreat"), 
        paste0(mf_text[["mf_factorA"]], 
          target_factor[! target_factor %in% mf_text[["mf_varB"]]], ".treatTreat" ))
        )
    } else if (mf_text[["mf_contrast"]] == "nested_effect") {
      target_factor <- unique(target[,mf_text[["mf_factorA"]]])
      mf_opt <- list(name=paste0(mf_text[["mf_factorA"]],  mf_text[["mf_varB"]], ".treatTreat" ))
    } 
  } else if (package_name %in% c("edgeR","limma")) { 
    if(mf_text[["mf_contrast"]] == "interaction") {
      mf_opt <- list(name=colnames(package_object)[4])
    } else if (mf_text[["mf_contrast"]] == "effect") {
      mf_opt <- list(name=paste0(
          mf_text[["mf_factorA"]],
          levels(target[, mf_text[["mf_factorA"]]])[2]
      ))
    } else if (mf_text[["mf_contrast"]] == "nested") {
      stop("Nested experiment design cannot be used with package other than DESeq2")
    }
  }
  if(! exists("mf_opt")) stop("Check multifactorial arguments flag valid")
  return(mf_opt)
}

# Dummy function to convert external, already processed/analysed data into 
# an object suitable for the big table function
analysis_external_DEA <- function(raw_data, dea_data) {
    return(list(raw_data, dea_data, dea_data))
}
