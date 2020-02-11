analysis_diff_correlation <- function(data, control, treat, PCIT_filter){
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
