################################################
### FUNCTIONS: DIFFERENCIAL EXPRESSION PACKAGES
################################################

# DESeq2
#-----------------------------------------------
analysis_DESeq2 <- function(data, p_val_cutoff, groups, linked_samples){
	require(DESeq2) # DESeqDataSetFromMatrix
	# Experimental design
	coldat <- DataFrame(grp=factor(groups), each=1) # Experimental design object
	# Add extra column to design if paired samples
	if(linked_samples == TRUE) {
		number_individuals <- length(design_vector)/2
		coldat <- data.frame(coldat, samp_id = factor(rep(1:number_individuals,2)))
		dds <- DESeqDataSetFromMatrix(data, colData=coldat, design =~ samp_id + grp) # proper DESeq2 object
	} else {
		dds <- DESeqDataSetFromMatrix(data, colData=coldat, design =~ grp) # proper DESeq2 object
	}
	dds <- DESeq(dds)
	# Calculating differential expression 	
	all_DESeq2_genes <- results(dds, alpha = p_val_cutoff)
	normalized_counts <- as.data.frame(counts(dds, normalized=TRUE)) # Getting normalized values
	
	return(list(normalized_counts, as.data.frame(all_DESeq2_genes), list(de_deseq2=all_DESeq2_genes, DESeq2_dataset=dds)))
}

analysis_edgeR <- function(data, path, groups, linked_samples){
	require(edgeR)
	d_edgeR <- DGEList(counts=data, group=groups) # Building edgeR object  

	# Obtain necessary data
	cols <- as.numeric(d_edgeR$samples$group)+2
	denraw <- cpm(d_edgeR, log=TRUE) # Counts Per Million

	if(linked_samples == TRUE) {
		number_individuals <- length(design_vector)/2
		samp_id <- factor(rep(1:number_individuals,2))
		Treat <- factor(groups, levels=c("C", "T"))
		design <- model.matrix(~samp_id+Treat)
		d_edgeR <- estimateDisp(d_edgeR, design)
		d_edgeR_DE <- glmQLFit(d_edgeR, design)
		d_edgeR_DE <- glmQLFTest(d_edgeR_DE)
	} else {
		d_edgeR <- calcNormFactors(d_edgeR)
		d_edgeR <- estimateCommonDisp(d_edgeR) # Estimate dispersions (common dispersion, then tagwise dispersion)
		d_edgeR <- estimateTagwiseDisp(d_edgeR)
		d_edgeR_DE <- exactTest(d_edgeR, dispersion = "auto", pair=unique(groups))
	}
	all_genes_df <- topTags(d_edgeR_DE, nrow(d_edgeR), sort.by="p.value")$table

	# Calculate Counts Per Million Mapped Reads
	dennorm <- cpm(d_edgeR, log=TRUE)
	# Getting normalized counts
	normalized_counts <- cpm(d_edgeR)

	# Return calculated info
	return(list(normalized_counts, all_genes_df, list(d_edgeR=d_edgeR, dennorm=dennorm)))
}

analysis_limma <- function(data, groups, linked_samples){
	require(edgeR)
	require(limma)
	# Experimental design
	groups <- factor(groups, levels = c("C", "T"))
	# Calculating differential expression
	DGE_List <- DGEList(counts=data, group=groups) # Building object (DGEList)
	DGE_List <- calcNormFactors(DGE_List) # Calculation of the normalization factor	


	if(linked_samples == TRUE) {
		number_individuals <- length(design_vector)/2
		samp_id <- factor(rep(1:number_individuals,2))
 		design <- model.matrix(~samp_id + groups)
 		log2_cpm <- voom(DGE_List, design) # Converting the read counts to log2-cpm
 		fit <- lmFit(log2_cpm, design)
 		fit2 <- eBayes(fit)
 		todos_limma <- topTable(fit2, adjust.method="BH", number=nrow(fit2), coef="groupsT")

	} else {
		design <- model.matrix(~0+ groups)
		log2_cpm <- voom(DGE_List, design) # Converting the read counts to log2-cpm
		normalized_counts <- log2_cpm$E
		# Fit and make contrast
		fit <- lmFit(log2_cpm, design)
		cont_matrix <- makeContrasts(c("groupsT-groupsC"), levels=design)
		fit2 <- contrasts.fit(fit, cont_matrix)
		fit2 <- eBayes(fit2)
		todos_limma <- topTable(fit2, adjust.method="BH", number=nrow(fit2))
	}
	normalized_counts <- log2_cpm$E
	# Return calculated info
	return(list(normalized_counts, todos_limma))
}

analysis_NOISeq <- function(data, num_controls, num_treatmnts, groups, lfc, path){
	require(NOISeq)
	# Experimental design
	groups_val <- c(rep("A", num_controls), rep("B", num_treatmnts)) 
	myfactors = data.frame(Tissue = rev(groups_val), TissueRun = rev(groups))
	# Calculation of differential expression
	mydata <- readData(data, myfactors) # Building the NOISeq Object (eset)

	all_NOISeq <- noiseqbio(mydata, k = 0.5, norm = "tmm", factor="Tissue", lc = 1, r = 50, adj = 1.5, plot = FALSE,
	a0per = 0.9, random.seed = 12345, filter = 1, cv.cutoff = 500, cpm = 1)

	normalized_counts <- tmm(assayData(mydata)$exprs, long = 1000, lc = 1) # Getting normalized counts
	expres_diff_all <- as.data.frame(all_NOISeq@results)

	prob_all <- expres_diff_all$prob
	# Inverse
	prob_all <- 1 - prob_all
	expres_diff_all <- data.frame(expres_diff_all, adj.p = prob_all)

	# Return calculated info
	return(list(normalized_counts, expres_diff_all, all_NOISeq))
}
