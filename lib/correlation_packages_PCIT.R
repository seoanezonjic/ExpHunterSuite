analysis_diff_correlation <- function(all_genes_stats, data, control, treat, PCIT_filter){
	# DE and PIF calculation, see "Microarray data processing, normalization and differential expression" in https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000382
	# DE => Diff exp, PIF =>Phenotypic Impact Factor
	#sum_control <- apply(control,1,sum)
	#sum_treat <- apply(treat,1,sum)
	#metrics <-data.frame(de=(sum_control - sum_treat)/10)
	nonZeroVal <- 10^-5
	if(sum(data==0) > 0) data[data == 0] <- replicate(sum(data == 0), jitter(nonZeroVal))
	if(sum(control==0) > 0) control[control == 0] <- replicate(sum(control == 0), jitter(nonZeroVal))
	if(sum(treat==0) > 0) treat[treat == 0] <- replicate(sum(treat == 0), jitter(nonZeroVal))
	metrics <- all_genes_stats['logFC_DESeq2']
	names(metrics) <- 'de'
	data <- log(data, 2)
	control <- log(control, 2)
	treat <- log(treat, 2)

	average <- apply(data, 1, mean) # Ai in https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000382
	treat_average <- apply(treat, 1, mean) # Ai in https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000382
	control_average <- apply(control, 1, mean) # Ai in https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000382
	metrics$average <- average
	metrics$treat_average <- treat_average
	metrics$control_average <- control_average
	metrics$pif <- average * metrics$de

	control <- stats::cor(t(control)) # Compute correlations
	treat <- stats::cor(t(treat)) # Compute correlations

	if(PCIT_filter == TRUE) {
		control <- do_pcit(control)
		treat <- do_pcit(treat)
	}

	#diferential_wiring: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000382
	dw <- control - treat
	#data$rif <- get_rif(data, dw) # Regulatory Impact Factor (RIF):  https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000382
	prevalent_gene_stats <- all_genes_stats[all_genes_stats$genes_tag == 'PREVALENT_DEG', ]
	if(nrow(prevalent_gene_stats) == 0){
		message('Using POSSIBLE_DEG instead PREVALENT_DEG')
		prevalent_gene_stats <- all_genes_stats[all_genes_stats$genes_tag == 'POSSIBLE_DEG', ]
	}
	metrics$rif1 <- scale(get_rif(1, metrics, prevalent_gene_stats, dw, control, treat)) # Regulatory Impact Factor (RIF):  https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000382
	metrics$rif2 <- scale(get_rif(2, metrics, prevalent_gene_stats, dw, control, treat)) # Regulatory Impact Factor (RIF):  https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000382
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

get_rif <- function(rif_version, metrics, deg_stats, dw, control_r, treat_r){
	genes <- row.names(metrics)
	rif <- NULL
	if(rif_version == 1){
		rif <- unlist(lapply(genes, compute_rif1, data=metrics, deg_stats=deg_stats, dw=dw))
	}else if(rif_version == 2){
		rif <- unlist(lapply(genes, compute_rif2, data=metrics, deg_stats=deg_stats, control_r=control_r, treat_r=treat_r))
	}
	return(rif)
}

# Regulatory Impact Factor (RIF):  https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000382
compute_rif1 <- function(g, data, deg_stats, dw){ #eq 4
	#adapted: It will be computed for all genes and not only TF. Also, all genes DW will be take into account instead of only significant DE genes
	# When implemented significant RNAseq DEgenes change to DE only.
	# Also, implement eq 5 and perform the mean between two eqs. 
	# A posterior paper, saids that the two measures are very different.
	deg_pif <- data[row.names(deg_stats),"pif" ]
	rif <- sum(deg_pif * dw[row.names(deg_stats), g]^2)/length(deg_pif)
	return(rif)
}

compute_rif2 <- function(g, data, deg_stats, control_r, treat_r){ #eq 5
	#adapted: It will be computed for all genes and not only TF. Also, all genes DW will be take into account instead of only significant DE genes
	# When implemented significant RNAseq DEgenes change to DE only.
	# Also, implement eq 5 and perform the mean between two eqs. 
	# A posterior paper, saids that the two measures are very different.
	rif_treat <- (data[row.names(deg_stats),"treat_average" ] * treat_r[row.names(deg_stats), g])^2
	rif_control <- (data[row.names(deg_stats),"control_average" ] * control_r[row.names(deg_stats), g])^2
	rif <- sum(rif_control - rif_treat)/length(rif_treat)
	return(rif)
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
