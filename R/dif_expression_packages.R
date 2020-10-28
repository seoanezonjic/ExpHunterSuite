################################################
### FUNCTIONS: DIFFERENCIAL EXPRESSION PACKAGES
################################################

# DESeq2
#-----------------------------------------------
analysis_DESeq2 <- function(data, p_val_cutoff, target, model_formula_text){
	dds <- DESeq2::DESeqDataSetFromMatrix(countData = data,
								  colData = target,
                                  design = formula(model_formula_text))
	dds <- DESeq2::DESeq(dds)
	all_DESeq2_genes <- DESeq2::results(dds, alpha = p_val_cutoff)
	normalized_counts <- as.data.frame(DESeq2::counts(dds, normalized=TRUE)) # Getting normalized values
	return(list(normalized_counts, as.data.frame(all_DESeq2_genes), list(de_deseq2=all_DESeq2_genes, DESeq2_dataset=dds)))
}

analysis_edgeR <- function(data, target, model_formula_text){
	# This is the default model_formula_text if nothing else is given.
	if(model_formula_text == "~ treat") {
		d_edgeR <- edgeR::DGEList(counts=data, group=as.character(target$treat)) # Building edgeR object  
		d_edgeR <- edgeR::calcNormFactors(d_edgeR)
		d_edgeR <- edgeR::estimateDisp(d_edgeR)
		d_edgeR_DE <- edgeR::exactTest(d_edgeR, dispersion = "auto", pair=c("Ctrl", "Treat"))
	} else {
		d_edgeR <- edgeR::DGEList(counts=data)
		model_design <- model.matrix(formula(paste("~", model_formula_text)), target)
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

analysis_limma <- function(data, target, model_formula_text){
	model_design <- model.matrix(formula(paste("~", model_formula_text)), target)
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