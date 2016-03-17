################################################
### FUNCTIONS: DIFFERENCIAL EXPRESSION PACKAGES
################################################

# DESeq2
#-----------------------------------------------
exp_dif_Deseq2 <- function(raw_filter, coldat){
  dds <- DESeqDataSetFromMatrix(raw_filter, colData=coldat, design =~ grp) # proper DESeq2 object
  dds <- DESeq(dds) #differential expression analysis of the DESeq2 object
  return(dds) # extract results
}

# NOISeq
#-----------------------------------------------
exp_dif_NOISeq <- function(mydata){
	mynoiseq = noiseqbio(mydata, k = 0.5, norm = "tmm", factor="Tissue", lc = 1, r = 50, adj = 1.5, plot = FALSE,
	a0per = 0.9, random.seed = 12345, filter = 1, cv.cutoff = 500, cpm = 1) #funciona en sistema de colas
 
	pdf(file.path(paths[['Results_NOISeq']], "ExpressionPlot.pdf"))
		DE.plot(mynoiseq, q = opt$q_value, graphic = "MD")
    dev.off()
  return(mynoiseq)
}

analysis_DESeq2 <- function(raw_filter, replicatesC, replicatesT, opt, lfc, paths){
	# Experimental design
	coldat = DataFrame(grp=factor(design_vector), each=1) # Experimental design object

	# Calculating differential expression 
	all_genes_df <- exp_dif_Deseq2(raw_filter, coldat)

	expres_diff <- results(all_genes_df)
	normalized_counts <- as.data.frame(counts(all_genes_df, normalized=TRUE)) # Getting normalized values
	if ((replicatesC > 1)&(replicatesT > 1)){ # Filtering DEGs by adjusted p-value only when there are replicates available (if not, only a descriptive analysis is performed)
	  expres_diff <- filter_gene_expression(expres_diff, opt$p_val_cutoff, lfc, "padj", "log2FoldChange")
	}

	expres_diff <- expres_diff[order(rownames(expres_diff)),]	
	return(list(expres_diff, normalized_counts, all_genes_df))
}


analysis_edgeR <- function(raw_filter, replicatesC, replicatesT, opt, lfc, paths){
	# Experimental design
	levels <- c("C", "T")
	groups <- c(rep("C", replicatesC), rep("T", replicatesT)) # design vector
	print(groups)
	# Calculating differential expression
	d_edgeR <- DGEList(counts=raw_filter, group=groups) # Building edgeR object  
	
	d_edgeR <- calcNormFactors(d_edgeR) # Calculation of normalization factor
	normalized_counts <- 1e6 * cpm(d_edgeR) # Getting normalized counts
	d_edgeR <- estimateCommonDisp(d_edgeR) # Estimate dispersions (common dispersion, then tagwise dispersion)
	d_edgeR <- estimateTagwiseDisp(d_edgeR)
	d_edgeR <- exactTest(d_edgeR, dispersion = "auto", pair=levels)

	# Extracts the top DE tags in a data frame for a given pair of groups, ranked by p-value or absolute log-fold change
	all_genes_df <- topTags(d_edgeR, n= nrow(d_edgeR$table), sort.by="p.value")$table # Sorting by p-value
	expres_diff <- filter_gene_expression(all_genes_df, opt$p_val_cutoff, lfc, "FDR", "logFC")

	expres_diff <- expres_diff[order(rownames(expres_diff)),]

	return(list(expres_diff, normalized_counts, all_genes_df))
}


analysis_limma <- function(raw_filter, replicatesC, replicatesT, opt, lfc, paths){
	# Experimental design
	groups <- c(rep("A", replicatesC), rep("B", replicatesT)) # design vector limma
	# Calculating differential expression
	DGE_List <- DGEList(counts=raw_filter, group=groups) # Building object (DGEList)
	DGE_List <- calcNormFactors(DGE_List) # Calculation of the normalization factor	
	normalized_counts <- as.matrix(DGE_List) # Getting normalized counts 
	design <- model.matrix(~0+ groups)
	log2_cpm <- voom(DGE_List, design, plot=TRUE) # Converting the read counts to log2-cpm
	fit <- lmFit(log2_cpm, design)
	cont_matrix <- makeContrasts(c("groupsB-groupsA"), levels=design)
	fit2 <- contrasts.fit(fit, cont_matrix)
	fit2 <- eBayes(fit2)
	results <- decideTests(fit2, p.value = opt$p_val_cutoff, lfc = lfc)
    todos_limma <- topTable(fit2, adjust.method="BH", number=nrow(fit2))
	expres_diff <- topTable(fit2, adjust.method = "BH", number=nrow(fit2), p.value= opt$p_val_cutoff, lfc = lfc)
	expres_diff <- expres_diff[order(rownames(expres_diff)),]
	return(list(expres_diff, normalized_counts, todos_limma))
}


analysis_NOISeq <- function(raw_filter, replicatesC, replicatesT, opt, paths){
	# Experimental design
	groups <- c(rep("A", replicatesC), rep("B", replicatesT)) 
	myfactors = data.frame(Tissue = rev(groups), TissueRun = rev(design_vector))
	print(myfactors)
	# Calculation of differential expression
	mydata <- readData(raw_filter, myfactors) # Building the NOISeq Object (eset)
	print(mydata)
	all_NOISeq <- exp_dif_NOISeq(mydata)
	
	normalized_counts <- tmm(assayData(mydata)$exprs, long = 1000, lc = 1) # Getting normalized counts
	expres_diff <- degenes(all_NOISeq, q = opt$q_value, M = NULL)

	#Calculate FDR from p
	prob <- sort(expres_diff$prob, decreasing=TRUE,na.last=TRUE)
	fdr <- cumsum(1 - prob)/1:length(prob) # calculating adjusted p-value
	expres_diff <- data.frame(expres_diff, adj.p = fdr, threshold=prob)
	expres_diff <- expres_diff[order(rownames(expres_diff)),]
	print(head(expres_diff))
	return(list(expres_diff, normalized_counts, all_NOISeq))
}