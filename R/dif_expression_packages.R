################################################
### FUNCTIONS: DIFFERENCIAL EXPRESSION PACKAGES
################################################
# Manager function
#-----------------------------------------------
perform_expression_analysis <- function(modules, replicatesC, replicatesT, raw_filter, p_val_cutoff, target, model_formula_text, external_DEA_data=NULL){
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
        cat('Gene expression analysis is performed with DESeq2.\n')
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
      # Calculate results
      results <- analysis_external_DEA(raw_filter, external_DEA_data)
      # Store results
      all_data_normalized[['externalDEA']] <- results[[1]]
      all_counts_for_plotting[['externalDEA']] <- results[[2]]
      package_objects[['externalDEA']] <- results[[3]]

      # Result Plot Visualization
      if (!is.null(all_counts_for_plotting[['externalDEA']])){
        all_pvalue_names   <- c(all_pvalue_names, names(external_DEA_data)[1])
        all_LFC_names      <- c(all_LFC_names, names(external_DEA_data)[2])
        all_FDR_names      <- c(all_FDR_names, names(external_DEA_data)[3])
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
analysis_DESeq2 <- function(data, p_val_cutoff, target, model_formula_text){
	dds <- DESeq2::DESeqDataSetFromMatrix(countData = data,
								  colData = target,
                                  design = stats::formula(model_formula_text))
	dds <- DESeq2::DESeq(dds)
	all_DESeq2_genes <- DESeq2::results(dds, alpha = p_val_cutoff)
	normalized_counts <- as.data.frame(DESeq2::counts(dds, normalized=TRUE)) # Getting normalized values
	return(list(normalized_counts, as.data.frame(all_DESeq2_genes), list(de_deseq2=all_DESeq2_genes, DESeq2_dataset=dds)))
}

#' @importFrom edgeR DGEList calcNormFactors estimateDisp exactTest glmQLFit glmQLFTest topTags cpm
#' @importFrom stats formula model.matrix
analysis_edgeR <- function(data, target, model_formula_text){
	# This is the default model_formula_text if nothing else is given.
	if(model_formula_text == "~ treat") {
		d_edgeR <- edgeR::DGEList(counts=data, group=as.character(target$treat)) # Building edgeR object  
		d_edgeR <- edgeR::calcNormFactors(d_edgeR)
		d_edgeR <- edgeR::estimateDisp(d_edgeR)
		d_edgeR_DE <- edgeR::exactTest(d_edgeR, dispersion = "auto", pair=c("Ctrl", "Treat"))
	} else {
		d_edgeR <- edgeR::DGEList(counts=data)
		model_design <- stats::model.matrix(stats::formula(paste("~", model_formula_text)), target)
		d_edgeR <- edgeR::estimateDisp(d_edgeR, model_design)
		d_edgeR_DE <- edgeR::glmQLFit(d_edgeR, model_design)
		d_edgeR_DE <- edgeR::glmQLFTest(d_edgeR_DE)
	}
	all_genes_df <- edgeR::topTags(d_edgeR_DE, nrow(d_edgeR), sort.by="none")$table
	dennorm <- edgeR::cpm(d_edgeR, log=TRUE)
	# Getting normalized counts
	normalized_counts <- edgeR::cpm(d_edgeR)
	# Return calculated info
	return(list(normalized_counts, all_genes_df, list(d_edgeR=d_edgeR, dennorm=dennorm)))
}

#' @importFrom limma voom eBayes topTable lmFit
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom stats formula model.matrix
analysis_limma <- function(data, target, model_formula_text){
	model_design <- stats::model.matrix(stats::formula(paste("~", model_formula_text)), target)
	# Calculating differential expression
	DGE_List <- edgeR::DGEList(counts=data) # Building object (DGEList)
	DGE_List <- edgeR::calcNormFactors(DGE_List) # Calculation of the normalization factor 
	log2_cpm <- limma::voom(DGE_List, model_design) # Converting the read counts to log2-cpm
	fit <- limma::lmFit(log2_cpm, model_design)
	fit2 <- limma::eBayes(fit)
	todos_limma <- limma::topTable(fit2, adjust.method="BH", number=nrow(fit2), coef="treatTreat")
	normalized_counts <- log2_cpm$E
    # Return calculated info
    return(list(normalized_counts, todos_limma))
}

#' @importFrom NOISeq readData noiseqbio tmm
#' @importFrom Biobase assayData
analysis_NOISeq <- function(data, target){
	group <- as.character(target$treat)
	# Experimental design
	myfactors = data.frame(Tissue = rev(group), TissueRun = rev(group))
	# Calculation of differential expression
	mydata <- NOISeq::readData(data, myfactors) # Building the NOISeq Object (eset)

	all_NOISeq <- NOISeq::noiseqbio(mydata, k = 0.5, norm = "tmm", factor="Tissue", lc = 1, r = 50, adj = 1.5, plot = FALSE,
	a0per = 0.9, random.seed = 12345, filter = 1, cv.cutoff = 500, cpm = 1)

	normalized_counts <- NOISeq::tmm(Biobase::assayData(mydata)$exprs, long = 1000, lc = 1) # Getting normalized counts
	expres_diff_all <- as.data.frame(all_NOISeq@results)

	prob_all <- expres_diff_all$prob
	# Inverse
	prob_all <- 1 - prob_all
	expres_diff_all <- data.frame(expres_diff_all, adj.p = prob_all)

	# Return calculated info
	return(list(normalized_counts, expres_diff_all, all_NOISeq))
}

# Dummy function to convert external, already processed/analysed data into an object suitable for the big table function
analysis_external_DEA <- function(raw_data, dea_data) {
	return(list(raw_data, dea_data, dea_data))
}
