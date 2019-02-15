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
	a0per = 0.9, random.seed = 12345, filter = 1, cv.cutoff = 500, cpm = 1)
 
	pdf(file.path(paths[['Results_NOISeq']], "ExpressionPlot_NOISeq.pdf"))
		DE.plot(mynoiseq, q = opt$q_value, graphic = "MD", cex.lab=1.4, cex.axis=1.4)
    dev.off()
  return(mynoiseq)
}

analysis_DESeq2 <- function(raw_filter, replicatesC, replicatesT, opt, lfc, paths){

	# Experimental design
	coldat = DataFrame(grp=factor(design_vector), each=1) # Experimental design object
	print(coldat)
	# Calculating differential expression 
	all_genes_DESeq2_object <- exp_dif_Deseq2(raw_filter, coldat)

	select <- order(rowMeans(counts(all_genes_DESeq2_object, normalized=TRUE)), decreasing=TRUE)[1:20]

	all_DESeq2_genes <- results(all_genes_DESeq2_object)

	normalized_counts <- as.data.frame(counts(all_genes_DESeq2_object, normalized=TRUE)) # Getting normalized values
	
	if ((replicatesC > 1)&(replicatesT > 1)){ # Filtering DEGs by adjusted p-value only when there are replicates available (if not, only a descriptive analysis is performed)
	  expres_diff <- filter_gene_expression(all_DESeq2_genes, opt$p_val_cutoff, lfc, "padj", "log2FoldChange")
	}
	return(list(expres_diff, normalized_counts, all_DESeq2_genes))
}


analysis_edgeR <- function(raw_filter, replicatesC, replicatesT, opt, lfc, paths){
	# Experimental design
	levels <- c("C", "T")
	groups <- c(rep("C", replicatesC), rep("T", replicatesT)) # design vector
	# Calculating differential expression
	d_edgeR <- DGEList(counts=raw_filter, group=groups) # Building edgeR object  

	cols <- as.numeric(d_edgeR$samples$group)+2

	denraw <- cpm(d_edgeR, log=TRUE)

	pdf(file.path(paths[['Results_edgeR']], "group_dendrogram_single.pdf"), w=11, h=8.5)
	  rawt <- t(denraw)
	  hc <- hclust(dist(rawt), "single")
	  plot(hc)
	dev.off()


	pdf(file.path(paths[['Results_edgeR']], "MDSplot.pdf"))
		plotMDS(d_edgeR, col=cols, main="MDS Plot: Treatment colours")
	dev.off()

	d_edgeR <- calcNormFactors(d_edgeR) # Calculation of normalization factor
	d_edgeR <- estimateCommonDisp(d_edgeR) # Estimate dispersions (common dispersion, then tagwise dispersion)
	d_edgeR <- estimateTagwiseDisp(d_edgeR)

	dennorm <- cpm(d_edgeR, log=TRUE)

	pdf(file.path(paths[['Results_edgeR']], "group_dendrogram_norm_average.pdf"), w=11, h=8.5)
  		dent_norm <- t(dennorm)
  		hc_norm <- hclust(dist(dent_norm), "ave")
  		plot(hc_norm)
	dev.off()

	pdf(file.path(paths[['Results_edgeR']], "group_dendrogram_norm_single.pdf"), w=11, h=8.5)
  		dent_norm <- t(dennorm)
  		hc_norm <- hclust(dist(dent_norm), "single")
  		plot(hc_norm)
	dev.off()

	pdf(file.path(paths[['Results_edgeR']], "MDSplot_norm.pdf"))
		plotMDS(d_edgeR, col=cols, main="MDS Plot: Treatment colours - normalized data")
	dev.off()

	normalized_counts <- cpm(d_edgeR) # Getting normalized counts
	
	d_edgeR <- exactTest(d_edgeR, dispersion = "auto", pair=levels)

	# Extracts the top DE tags in a data frame for a given pair of groups, ranked by p-value or absolute log-fold change
	all_genes_df <- topTags(d_edgeR, n= nrow(d_edgeR$table), sort.by="p.value")$table
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

	design <- model.matrix(~0+ groups)
	log2_cpm <- voom(DGE_List, design, plot=TRUE) # Converting the read counts to log2-cpm
	normalized_counts <- log2_cpm$E
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

	all_NOISeq <- exp_dif_NOISeq(mydata)

	normalized_counts <- tmm(assayData(mydata)$exprs, long = 1000, lc = 1) # Getting normalized counts
	expres_diff <- degenes(all_NOISeq, q = opt$q_value, M = NULL)
	expres_diff_all <- as.data.frame(all_NOISeq@results)

	prob_all <- expres_diff_all$prob
	for (i in c(1:nrow(expres_diff_all))){
		prob_all[i] <- 1-(prob_all[i])
	}

	expres_diff_all <- data.frame(expres_diff_all, adj.p = prob_all)

	prob <- expres_diff$prob
	for (i in c(1:nrow(expres_diff))){
		prob[i] <- 1-(prob[i])
	}
	expres_diff <- data.frame(expres_diff, adj.p = prob)
	expres_diff <- expres_diff[order(rownames(expres_diff)),]
	return(list(expres_diff, normalized_counts, expres_diff_all))
}
