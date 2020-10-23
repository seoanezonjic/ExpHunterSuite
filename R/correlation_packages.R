analysis_diffcoexp <- function(data, path, target) {

	dds <- DESeqDataSetFromMatrix(countData = data, colData = target["treat"], design = ~ treat)
	rld <- rlog(dds, blind=FALSE)
	mat_counts <- assay(rld)
	treat <- mat_counts[,target$treat=="Treat"]
	control <- mat_counts[,target$treat=="Ctrl"]

	allowWGCNAThreads()
	res <- diffcoexp(exprs.1 = control, exprs.2 = treat, r.method = "spearman", q.dcgth=1)

	write.table(res$DCLs, file=file.path(path, "DCLs.txt"), sep="\t", quote=FALSE, row.names=FALSE)
	write.table(res$DCGs, file=file.path(path, "DCGs.txt"), sep="\t", quote=FALSE, row.names=FALSE)
	return(res)
}

analysis_WGCNA <- function(data, path, target_numeric_factors, target_string_factors, WGCNA_memory, WGCNA_deepsplit, WGCNA_detectcutHeight, WGCNA_mergecutHeight, WGCNA_min_genes_cluster, cor_only, blockwiseNetworkType, blockwiseTOMType) {

	data <- t(data)#[, 1:500]
	nSamples <- nrow(data)

	####################################################################
	## THRESHOLDING - EFFECTS OF BETA ON TOPOLOGY AND AUTO SELECTION
	####################################################################
	assignInNamespace(x="..minNSamples", value=3, ns="WGCNA") #Overwrite harcoded limit from 4 to 3
	powers <- c(c(1:10), seq(from = 12, to=30, by=2))

 	sft <- WGCNA::pickSoftThreshold(data, powerVector = powers, verbose = 5)

	# Calculate Power automatically
	sft_mfs_r2 <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
	min_pow_ind <- which(sft_mfs_r2 > 0.9)[1] # First time the value passes 0.9

	in_advance <- 2
	if(is.na(min_pow_ind)) {
		# Plan B to calculate power
		for(i in 1:(length(sft_mfs_r2)-in_advance)) {
			vals <- sft_mfs_r2[i:(i+in_advance)]
			diffs <- abs(diff(vals))
			# Condition to meet
			if(sum(diffs < 0.01) == in_advance & sft_mfs_r2[i] > 0.8) {
			  min_pow_ind <- i
	   		  break
			}
		}
	}

	if(is.na(min_pow_ind)) {
		warning("Could not obtain a valid power (beta) value for WGCNA, so the default of 30 will be used - proceed with caution")
		min_pow_ind <- length(sft_mfs_r2) # assumes 30 will be the largest testable power
	}
	pow <- sft$fitIndices[min_pow_ind, "Power"]
	cor <- WGCNA::cor # TO CORRECT FUNCTION OVERRIDE DUE TO OTHER PACKAGES

	pdf(file.path(path, "thresholding.pdf"))
		dev.control(displaylist="enable")
		par(mfrow = c(1,2));
		cex1 = 0.9;
		# Scale-free topology fit index as a function of the soft-thresholding power
		plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
		     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
		     main = paste("Scale independence"));
		text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
		     labels=powers,cex=cex1,col="red");
		# this line corresponds to using an R^2 cut-off of h
		abline(h=0.90, col="red")
		abline(h=0.80, col="red", lty="dashed")
		abline(v=pow, col="black", lty="dotted")

		# Mean connectivity as a function of the soft-thresholding power
		plot(sft$fitIndices[,1], sft$fitIndices[,5],
		     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
		     main = paste("Mean connectivity"))
		text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
		abline(v=pow, col="black", lty="dotted")

		power_threshold_effects <- recordPlot()
	dev.off()

	####################################################################
	## CLUSTER SAMPLES TO GENERATE MODULES
	####################################################################

	# Have to improve this as generating fork.
	tom_file_base <- "generatedTOM"
	if(length(list.files(path=path, pattern=paste0(tom_file_base, ".*\\.RData"))) > 0) {
		loadTOM_TF <- TRUE
		saveTOM_TF <- FALSE
	} else {
		loadTOM_TF <- FALSE
		saveTOM_TF <- TRUE
	}

	net <- WGCNA::blockwiseModules(data, power = pow,
	maxBlockSize = WGCNA_memory, # Increase to memory limit in order to obtain more realistic results
	minModuleSize = WGCNA_min_genes_cluster,
	deepSplit = WGCNA_deepsplit, detectCutHeight = WGCNA_detectcutHeight,
	reassignThreshold = 0, mergeCutHeight = WGCNA_mergecutHeight,
	numericLabels = TRUE, pamRespectsDendro = FALSE,
	loadTOM = loadTOM_TF,	
	saveTOM = saveTOM_TF,
	saveTOMFileBase = file.path(path, tom_file_base), 
	verbose = 5, 
	networkType = blockwiseNetworkType,
	TOMType = blockwiseTOMType)

	moduleColors = WGCNA::labels2colors(net$colors)

	# Plot the dendrogram and the module colors underneath
	pdf(file.path(path, 'wcgnaModules.pdf'))
		WGCNA::plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
	dev.off()

	MEs <- WGCNA::moduleEigengenes(data, net$colors)$eigengenes

	gene_module_cor <- as.data.frame(cor(data, MEs, use = "p"));
	gene_module_cor_p <- as.data.frame(WGCNA::corPvalueStudent(as.matrix(gene_module_cor), nSamples));
	colnames(gene_module_cor_p) = colnames(gene_module_cor) <- gsub("ME", "Cluster_", colnames(gene_module_cor_p) )

	write.table(gene_module_cor, file=file.path(path, "gene_MM.txt"), sep="\t", quote=FALSE)
	write.table(gene_module_cor_p, file=file.path(path, "gene_MM_p_val.txt"), sep="\t", quote=FALSE)
	write.table(MEs, file=file.path(path, "eigen_values_per_samples.txt"), sep="\t", quote=FALSE)

	if(cor_only == TRUE) {
		return("cor_only")
	}

	####################################################################
	## PROCESS TRAIT DATA
	####################################################################

	trait <- data.frame(row.names=rownames(data))

	if(is.data.frame(target_numeric_factors)) {
		trait <- data.frame(trait, target_numeric_factors)
	}
	if(is.data.frame(target_string_factors)) {
		# Code to convert the string factors to numeric (1 vs. 0 for each category). NOTE YOU NEED A MINIMUM COUNT OF 2 FOR EACH CATEGORY IN THE FACTORS
		binarized_string_factors <- lapply(names(target_string_factors), function(factor_name) {
  			string_factor <- target_string_factors[[factor_name]]
  			binarized_table <- WGCNA::binarizeCategoricalVariable(string_factor,
  								minCount = 1,
								includePairwise = FALSE,
								includeLevelVsAll = TRUE)
			colnames(binarized_table) <- paste(factor_name, levels(string_factor),sep="_")
  			return(binarized_table)
		})
		trait <- data.frame(trait, binarized_string_factors)
	}

	###################################################################
	## CLUSTER DATA
	###################################################################

	####################################################################
	## CLUSTER MODULES WITH TRAITS AND PRODUCE HEATMAP
	####################################################################

	MEs_colors <- MEs <- WGCNA::moduleEigengenes(data, net$colors)$eigengenes
	ME_numeric <- as.numeric(gsub("ME", "", colnames(MEs)))
	colnames(MEs_colors) <- paste0("ME", WGCNA::labels2colors(ME_numeric))

	moduleTraitCor = cor(MEs_colors, trait, use = "p")
	moduleTraitPvalue = WGCNA::corPvalueStudent(moduleTraitCor, nSamples)
 
	textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
	                    signif(moduleTraitPvalue, 1), ")", sep = "")
	dim(textMatrix) = dim(moduleTraitCor)
	colnames(moduleTraitCor) <- names(trait)
	rownames(moduleTraitCor) <- names(MEs)
	moduleTraitCor_df <- as.data.frame(as.table(moduleTraitCor))
	colnames(textMatrix) <- names(trait)
	rownames(textMatrix) <- names(MEs)
	textMatrix_df <- as.data.frame(as.table(textMatrix))
	colnames(moduleTraitCor_df) <- c("Module","Trait","Correlation")
	moduleTraitCor_df$Text_correlation <- textMatrix_df[["Freq"]]
	####################################################################
	# Produce trait AND module correlations and p-values
	####################################################################
	ME_numeric_traits <- WGCNA::orderMEs(cbind(MEs, trait))

	trait_and_module_cor <- cor(ME_numeric_traits, use = "p")
	#write.table(trait_and_module_cor, file=file.path(path, "trait_and_module_cor.txt"), sep="\t", quote=FALSE)
	trait_and_module_cor_p <- WGCNA::corPvalueStudent(trait_and_module_cor, nSamples)
	#write.table(trait_and_module_cor_p, file=file.path(path, "trait_and_module_cor_p.txt"), sep="\t", quote=FALSE)
	trait_module_cor_val <- as.data.frame(as.table(trait_and_module_cor))
	trait_module_p_val <- as.data.frame(as.table(trait_and_module_cor_p))

	trait_module_all_vals <- data.frame(trait_module_cor_val, trait_module_p_val[,3])
	names(trait_module_all_vals) <- c("A", "B", "correlation", "p-value")
	write.table(trait_module_all_vals, file=file.path(path, "trait_and_module_all_vals.txt"), sep="\t", quote=FALSE, row.names=FALSE)
	# write.table(trait_and_module_cor, file=file.path(path, "trait_module_all_vals2.txt"), sep="\t", quote=FALSE)

	####################################################################
	## Produce Tables
	####################################################################
	# Report tables: 
	# Genes per module 
	gene_module_cor <- as.data.frame(cor(data, MEs, use = "p"));
	gene_module_cor_p <- as.data.frame(WGCNA::corPvalueStudent(as.matrix(gene_module_cor), nSamples));
	colnames(gene_module_cor_p) = colnames(gene_module_cor) <- gsub("ME", "Cluster_", colnames(gene_module_cor_p) )
	# Genes per trait
	gene_trait_cor <- as.data.frame(cor(data, trait, use = "p"));
	gene_trait_cor_p <- as.data.frame(WGCNA::corPvalueStudent(as.matrix(gene_trait_cor), nSamples));
	# Module per trait (also produced above for the plot - should give smae results.)
	module_trait_cor = cor(MEs, trait, use = "p")
	module_trait_cor_p <- WGCNA::corPvalueStudent(module_trait_cor, nSamples)
	row.names(module_trait_cor) = row.names(module_trait_cor_p) <- gsub("ME", "Cluster_", row.names(module_trait_cor_p))

	write.table(trait, file=file.path(path, "sample_trait.txt"), sep="\t", quote=FALSE)
	write.table(gene_trait_cor, file=file.path(path, "gene_trait.txt"), sep="\t", quote=FALSE)
	write.table(gene_trait_cor_p, file=file.path(path, "gene_trait_p_val.txt"), sep="\t", quote=FALSE)
	write.table(module_trait_cor, file=file.path(path, "module_trait.txt"), sep="\t", quote=FALSE)
	write.table(module_trait_cor_p, file=file.path(path, "module_trait_p_val.txt"), sep="\t", quote=FALSE)

	cluster_ID<- net$colors
	Cluster_MM <- sapply(names(cluster_ID), function(x) gene_module_cor[x, paste0("Cluster_", cluster_ID[x])]) 
	Cluster_MM_pval <- sapply(names(cluster_ID), function(x) gene_module_cor_p[x, paste0("Cluster_", cluster_ID[x])]) # Get the MM value

	# Plot cluster ID vs. ID of cluster with lowest MM p-value for each gene
	MM_Cluster_ID <- apply(gene_module_cor_p, 1, function(x) which(x == min(x))) - 1 # Cluster ID with minimum MM p value
	cluster_vs_MM <- ggplot2::ggplot(as.data.frame(cbind(cluster_ID, MM_Cluster_ID)), ggplot2::aes(x = cluster_ID ,y = MM_Cluster_ID)) + ggplot2::geom_count() + 
		ggplot2::scale_x_continuous("Cluster ID", labels = function(x) paste0("Cluster_", x)) + 
		ggplot2::scale_y_continuous("Cluster with highest MM value", labels = function(x) paste0("Cluster_", x))

	# NEW - get the p value for the cluster id, not the lowest p-value
	# cluster_cor_pval <- sapply(names(cluster_ID), function(x) gene_module_cor_p[x,cluster_ID[x]+1])

	return(list(
			gene_cluster_info = data.frame(ENSEMBL_ID = names(cluster_ID), Cluster_ID = as.numeric(cluster_ID),
				 Cluster_MM = Cluster_MM, Cluster_MM_pval = Cluster_MM_pval),
			package_objects=list(gene_module_cor = gene_module_cor, gene_module_cor_p = gene_module_cor_p,
				gene_trait_cor = gene_trait_cor, gene_trait_cor_p = gene_trait_cor_p, 
				module_trait_cor = module_trait_cor, module_trait_cor_p = module_trait_cor_p),
				plot_objects=list(sorted_colours = WGCNA::labels2colors(ME_numeric),
				trait_vs_module = moduleTraitCor_df,
				trait_and_module = ME_numeric_traits,
				power_threshold_effects = power_threshold_effects,
				cluster_vs_MM = cluster_vs_MM)
			)
		)
}