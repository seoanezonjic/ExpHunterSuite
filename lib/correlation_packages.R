analysis_WGCNA <- function(data, path, target_numeric_factors, target_string_factors, WGCNA_memory, cor_only) {
	data <- t(data)#[, 1:500]

	####################################################################
	## THRESHOLDING - EFFECTS OF BETA ON TOPOLOGY AND AUTO SELECTION
	####################################################################

	powers <- c(c(1:10), seq(from = 12, to=30, by=2))
 	sft <- pickSoftThreshold(data, powerVector = powers, verbose = 5)

	pdf(file.path(path, "thresholding.pdf"))
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
		warning("Could not obtain a valid power (beta) value for WGCNA")
		return("NO_POWER_VALUE")
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

	nSamples = nrow(data);
	MEs = moduleEigengenes(data, net$colors)$eigengenes
	MEs_colors = moduleEigengenes(data, moduleColors)$eigengenes

	moduleTraitCor = cor(MEs_colors, trait, use = "p")
	moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
 
	textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
	                    signif(moduleTraitPvalue, 1), ")", sep = "")
	dim(textMatrix) = dim(moduleTraitCor)
	# Display the correlation values within a heatmap plot
	pdf(file.path(path, 'heatmap.pdf'))
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
	# Module per trait
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
				module_trait_cor=module_trait_cor, module_trait_cor_p=module_trait_cor_p)
			)
		)
}

analysis_diff_correlation <- function(data, control, treat, PCIT_filter) {
	# DE and PIF calculation, see "Microarray data processing, normalization and differential expression" in https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000382
	# DE => Diff exp, PIF =>Phenotypic Impact Factor
	sum_control <- apply(control,1,sum)
	sum_treat <- apply(treat,1,sum)
	metrics <-data.frame(de=(sum_control - sum_treat)/10)
	average <- apply(data, 1, mean) # Ai in https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000382
	metrics$pif <- average * metrics$de


	control <- cor(t(control)) # Compute correlations
	treat <- cor(t(treat)) # Compute correlations

	if(PCIT_filter == TRUE) {
		control <- do_pcit(control)
		treat <- do_pcit(treat)
	}

        #diferential_wiring: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000382
        dw <- control - treat
        #data$rif <- get_rif(data, dw) # Regulatory Impact Factor (RIF):  https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000382
        metrics$rif <-get_rif(metrics, dw) # Regulatory Impact Factor (RIF):  https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000382

	#Differential hubbing (differential conectivity): https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000382
	metrics$ctrl_cn <- apply(control,1,function(x) length(which(abs(x) > 0.8)))
	metrics$treat_cn <- apply(treat,1,function(x) length(which(abs(x) > 0.8)))
	metrics$diff_cn <- metrics$ctrl_cn - metrics$treat_cn
        # As Horvath (WGCNA): DiffK (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1998880/pdf/335_2007_Article_9043.pdf)
        k_control <- metrics$ctrl_cn/max(metrics$ctrl_cn)
        k_treat <- metrics$treat_cn/max(metrics$treat_cn)
        metrics$diffK <- k_control - k_treat
	return(metrics)
}

# Regulatory Impact Factor (RIF):  https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000382
compute_rif <- function(g, data, dw){ #eq 4
	#adapted: It will be computed for all genes and not only TF. Also, all genes DW will be take into account instead of only significant DE genes
	# When implemented significant RNAseq DEgenes change to DE only.
	# Also, implement eq 5 and perform the mean between two eqs. A posterior paper, saids that the two measures are very different.
	rif <- sum(data$pif * dw[,g] ^2)
	return(rif)
}

get_rif <- function(data, dw){
	genes <- rownames(data)
	rif <- unlist(lapply(genes, compute_rif, data=data, dw=dw))
}

do_pcit <- function(cor_mat){
        result <- pcit(cor_mat) # Check correlations by PCIT
        signif <- idx(result) # Get significative correlations (vector with linear index of matrix, not pairs of coordinates)
        nonsignif <- idxInvert(nrow(cor_mat), signif) # Get NON significative correlations (vector with linear index of matrix, not pairs of coordinates)
        n_elements <- nrow(cor_mat)
        # plotCorCoeff(cor_mat, list("PCIT Meaningful" = signif), col=c("red")) #TODO make conditional, I think that is important for reports
        cor_mat[nonsignif] <- 0
        # Warning: number of edges and signif length are NOT equal because the diagonal and one matrix triangle are discarded
        #message(paste('Significative corr', (length(signif) - n_elements) / 2, sep="\t"))
        #message(paste('NON significative', length(nonsignif) / 2, sep="\t"))
        return(cor_mat)
}
