analysis_WGCNA <- function(data, path) {

	data <- t(data)#[, 1:500]

	powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
	sft = pickSoftThreshold(data, powerVector = powers, verbose = 5)
	pdf(file.path(path, "thresholding.pdf"))
		# Plot the results:
		# sizeGrWindow(9, 5)
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
	min_pow_ind <- min(which(sft_mfs_r2 > 0.9))
	if(min_pow_ind == Inf) {
		warning("Could not obtain a valid power (beta) value for WGCNA")
		return("NO_POWER_VALUE")
	}
	pow <- sft$fitIndices[min_pow_ind, "Power"]
	# pow <- 12
	cor <- WGCNA::cor # TO CORRECT FUNCTION OVERRIDE DUE TO OTHER PACKAGES

	net <- blockwiseModules(data, power = pow,
			maxBlockSize = 5000, # Increase to memory limit in order to obtain more realistic results
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = file.path(path,"generatedTOM"), 
                       verbose = 5)

	save(net, file=file.path(path,"blockwise_real_output.RData"))
	## PLot
	###############################################################3
	mergedColors = labels2colors(net$colors)
	# Plot the dendrogram and the module colors underneath
	pdf(file.path(path, 'wcgnaModules.pdf'))
		plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
	dev.off()

	moduleLabels = net$colors

	moduleColors = labels2colors(net$colors)
	MEs = net$MEs;
	geneTree = net$dendrograms[[1]];

	module_membership <- net$colors
	return(data.frame(ENSEMBL_ID = names(module_membership), Cluster_ID = as.numeric(module_membership)))
}