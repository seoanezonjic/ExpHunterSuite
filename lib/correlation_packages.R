analysis_WGCNA <- function(data, path, target_numeric_factors, target_string_factors, WGCNA_memory, cor_only) {
	data <- t(data)#[, 1:500]

	####################################################################
	## THRESHOLDING - EFFECTS OF BETA ON TOPOLOGY AND AUTO SELECTION
	####################################################################

	powers <- c(c(1:10), seq(from = 12, to=30, by=2))
 	sft <- pickSoftThreshold(data, powerVector = powers, verbose = 5)

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
		abline(h=0.90,col="red")
		# Mean connectivity as a function of the soft-thresholding power
		plot(sft$fitIndices[,1], sft$fitIndices[,5],
		     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
		     main = paste("Mean connectivity"))
		text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
		power_threshold_effects <- recordPlot()
	dev.off()

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

	net <- blockwiseModules(data, power = pow,
	maxBlockSize = WGCNA_memory, # Increase to memory limit in order to obtain more realistic results
	TOMType = "unsigned", minModuleSize = 30,
	reassignThreshold = 0, mergeCutHeight = 0.25,
	numericLabels = TRUE, pamRespectsDendro = FALSE,
	loadTOM = loadTOM_TF,	
	saveTOM = saveTOM_TF,
	saveTOMFileBase = file.path(path, tom_file_base), 
	verbose = 5)

	moduleColors = labels2colors(net$colors)

	# Plot the dendrogram and the module colors underneath
	pdf(file.path(path, 'wcgnaModules.pdf'))
		plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
	dev.off()

	MEs <- moduleEigengenes(data, net$colors)$eigengenes

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
  			binarized_table <- binarizeCategoricalVariable(string_factor,
  								minCount = 2,
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
	sampleTree2 = hclust(dist(data), method = "average")

	# Convert traits to a color representation: white means low, red means high, grey means missing entry
	traitColors = numbers2colors(trait, signed = FALSE);

	# Plots the sample dendrogram and the colors underneath.
	pdf(file.path(path, 'count_trait_data_comparation.pdf'))
		plotDendroAndColors(sampleTree2, traitColors,
			groupLabels = names(trait), 
			main = "Sample dendrogram and trait heatmap")
	dev.off()


	####################################################################
	## CLUSTER MODULES WITH TRAITS AND PRODUCE HEATMAP
	####################################################################

	#save(data, net, trait, moduleColors, file="~/test.RData")


	nSamples <- nrow(data)
	MEs_colors <- MEs <- moduleEigengenes(data, net$colors)$eigengenes
	ME_numeric <- as.numeric(gsub("ME", "", colnames(MEs)))
	colnames(MEs_colors) <- paste0("ME", labels2colors(ME_numeric))

	moduleTraitCor = cor(MEs_colors, trait, use = "p")
	moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
 
	textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
	                    signif(moduleTraitPvalue, 1), ")", sep = "")
	dim(textMatrix) = dim(moduleTraitCor)
	# Display the correlation values within a heatmap plot
	pdf(file.path(path, 'heatmap.pdf'))
		dev.control(displaylist="enable")
		par(mar = c(6, 10, 3, 3));
		labeledHeatmap(Matrix = moduleTraitCor,
	               xLabels = names(trait),
	               yLabels = names(MEs_colors),
	               ySymbols = names(MEs),
	               colorLabels = FALSE,
	               colors = blueWhiteRed(50),
	               textMatrix = textMatrix,
	               setStdMargins = FALSE,
	               cex.text = 0.5,
	               zlim = c(-1,1),
	               main = paste("Module-trait relationships"))
		trait_vs_module_heatmap <- recordPlot()
	dev.off()

	####################################################################
	# Produce dendogram and heatmap
	####################################################################

	#ME_s <- moduleEigengenes(data, labels2colors(net$colors))$eigengenes
	ME_traits = orderMEs(cbind(MEs_colors, trait))

	pdf(file.path(path, 'eigengenes_dendogram.pdf'))
		dev.control(displaylist="enable")
		plotEigengeneNetworks(ME_traits, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
		trait_and_module_dendogram <- recordPlot()
	dev.off()
	
	pdf(file.path(path, 'eigengenes_heatmap.pdf'))
		dev.control(displaylist="enable")
		plotEigengeneNetworks(ME_traits, "Eigengene adjacency heatmap", marHeatmap = c(12,12,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
		trait_and_module_heatmap <- recordPlot()
	dev.off()

	####################################################################
	## Produce Tables
	####################################################################
	# Report tables: 
	# Genes per module 
	gene_module_cor <- as.data.frame(cor(data, MEs, use = "p"));
	gene_module_cor_p <- as.data.frame(corPvalueStudent(as.matrix(gene_module_cor), nSamples));
	colnames(gene_module_cor_p) = colnames(gene_module_cor) <- gsub("ME", "Cluster_", colnames(gene_module_cor_p) )
	# Genes per trait
	gene_trait_cor <- as.data.frame(cor(data, trait, use = "p"));
	gene_trait_cor_p <- as.data.frame(corPvalueStudent(as.matrix(gene_trait_cor), nSamples));
	# Module per trait (also produced above for the plot - should give smae results.)
	module_trait_cor = cor(MEs, trait, use = "p")
	module_trait_cor_p <- corPvalueStudent(module_trait_cor, nSamples)
	row.names(module_trait_cor) = row.names(module_trait_cor_p) <- gsub("ME", "Cluster_", row.names(module_trait_cor_p))
	# Extra text to add to big table
	MM_Cluster_ID <- apply(gene_module_cor_p, 1, function(x) which(x == min(x))) - 1 # Get order the change to -1 as we are 0 indexed for the Cluster IDs
	MM_min_cor_pval <- apply(gene_module_cor_p, 1, function(x) min(x))
	cluster_ID<- net$colors

	return(list(
			gene_cluster_info=data.frame(ENSEMBL_ID = names(cluster_ID), Cluster_ID = as.numeric(cluster_ID),
				MM_min_cor_pval=MM_min_cor_pval, MM_Cluster_ID=MM_Cluster_ID),
			package_objects=list(gene_module_cor=gene_module_cor, gene_module_cor_p=gene_module_cor_p,
				gene_trait_cor=gene_trait_cor, gene_trait_cor_p=gene_trait_cor_p, 
				module_trait_cor=module_trait_cor, module_trait_cor_p=module_trait_cor_p),
			plot_objects=list(trait_vs_module_heatmap = trait_vs_module_heatmap, 
				trait_and_module_dendogram = trait_and_module_dendogram, 
				trait_and_module_heatmap = trait_and_module_heatmap,
				power_threshold_effects = power_threshold_effects)
			)
		)
}