merge_dim_tables <- function(dim_data_simp){
	quali_data <- data.frame()
	quanti_data <- data.frame()
	cat_data <- data.frame()
	dim_data_simp$call <- NULL
	for(dimension in names(dim_data_simp)){
		
		dim_data <- dim_data_simp[[dimension]]
		if (!is.null(dim_data$quali)) {

			quali_data_dim <- as.data.frame(dim_data$quali)
			quali_data_dim$dimension <- dimension
			quali_data_dim$factor <- rownames(quali_data_dim)
			quali_data <- rbind(quali_data, quali_data_dim)	
        }
     if (!is.null(dim_data$quanti)){
			quanti_data_dim <-  as.data.frame(dim_data$quanti)
			quanti_data_dim$dimension <- dimension
			quanti_data_dim$factor <- rownames(quanti_data_dim)
			quanti_data <- rbind(quanti_data, quanti_data_dim)	
		}
		if (!is.null(dim_data$category)){
			cat_data_dim <-  as.data.frame(dim_data$category)
			cat_data_dim$dimension <- dimension
			cat_data_dim$factor <- rownames(cat_data_dim)
			cat_data <- rbind(cat_data, cat_data_dim)	
		}

	}
	return(list(quantitative = quanti_data,
		qualitative = quali_data,
		qual_category = cat_data))
}

#' @importFrom FactoInvestigate dimRestrict eigenRef
get_PCA_dimensions <- function(pca_obj, min_dimensions = 2, time = "10s") {
    ref = FactoInvestigate::eigenRef(pca_obj, time = time, parallel=FALSE) # to avoid use parallel computation that greedy takes all cpu cores
    rand = c(ref$inertia[1], diff(ref$inertia)) * 100
    keep_dimensions <- FactoInvestigate::dimRestrict(pca_obj, rand = rand)
    if(keep_dimensions < min_dimensions){
      keep_dimensions <- min_dimensions
      message('Significant axis are less than 2. The first two axis will be selected to continue the analysis')
    }
    return(keep_dimensions)
}

#' @importFrom FactoMineR PCA dimdesc HCPC
compute_pca <- function(pca_data, 
						target = NULL, 
						string_factors = NULL, 
						numeric_factors = NULL, 
						transpose = TRUE,
						add_samples = NULL,
						min_dimensions = 2) {

	if (transpose) 
		pca_data <- as.data.frame(t(pca_data))
	 
	if (!is.null(target)) {
		rownames(target) <- as.character(target$sample)

		if (!is.null(numeric_factors)){
			pca_data <- merge_factors(pca_data, target, numeric_factors)

		}

		if (!is.null(string_factors)) {
		  pca_data <- merge_factors(pca_data, target, string_factors)
		}
	} 

	raw_pca_data <- pca_data[!rownames(pca_data) %in% add_samples,
							 ! colnames(pca_data) %in% c(numeric_factors, string_factors)]

	std_pca <- FactoMineR::PCA(raw_pca_data, scale.unit=TRUE, 
	  							graph = FALSE)                                                     
	dim_to_keep <- get_PCA_dimensions(std_pca, min_dimensions = min_dimensions)

	
	if (!is.null(add_samples)) {
		add_samples_idx <- match(add_samples, rownames(pca_data))
	} else {
		add_samples_idx <- NULL
	}

	pca_res <- FactoMineR::PCA(pca_data, ncp = dim_to_keep, 
										scale.unit=TRUE, 
										graph = FALSE,
										quanti.sup = numeric_factors, 
										quali.sup=string_factors,
										ind.sup = add_samples_idx)	
 	dim_data <- FactoMineR::dimdesc(pca_res, axes=seq(1, dim_to_keep))
    dim_data_merged <- merge_dim_tables(dim_data)

  res.hcpc <- FactoMineR::HCPC(pca_res, graph = FALSE, nb.clust = -1)

	return(list(pca_data = pca_res,
							dim_to_keep = dim_to_keep,
							dim_data = dim_data,
							dim_data_merged = dim_data_merged,
							res.hcpc = res.hcpc))
}


#' @importFrom FactoMineR CA dimdesc HCPC
compute_mca <- function(mca_data, 
						target = NULL, 
						string_factors = NULL, 
						numeric_factors = NULL, 
						transpose = TRUE,
						add_samples = NULL,
						min_dimensions = 2) {

	if (transpose) 
		mca_data <- as.data.frame(t(mca_data))
	 
	if (!is.null(target)) {
		rownames(target) <- as.character(target$sample)

		if (!is.null(numeric_factors)){
			mca_data <- merge_factors(mca_data, target, numeric_factors)

		}

		if (!is.null(string_factors)) {
		  mca_data <- merge_factors(mca_data, target, string_factors)
		}
	} 

	raw_mca_data <- mca_data[!rownames(mca_data) %in% add_samples,
							 ! colnames(mca_data) %in% c(numeric_factors, string_factors)]

	std_mca <- FactoMineR::MCA(raw_mca_data, 	graph = FALSE)                                                     
	dim_to_keep <- get_PCA_dimensions(std_mca, min_dimensions = min_dimensions)

	
	if (!is.null(add_samples)) {
		add_samples_idx <- match(add_samples, rownames(mca_data))
	} else {
		add_samples_idx <- NULL
	}

	mca_res <- FactoMineR::MCA(mca_data, ncp = dim_to_keep, 
										graph = FALSE,
										quanti.sup = numeric_factors, 
										quali.sup=string_factors,
										ind.sup = add_samples_idx)	
 	dim_data <- FactoMineR::dimdesc(mca_res, axes=seq(1, dim_to_keep))
    dim_data_merged <- merge_dim_tables(dim_data)

  res.hcpc <- FactoMineR::HCPC(mca_res, graph = FALSE, nb.clust = -1)

	return(list(pca_data = mca_res,
							dim_to_keep = dim_to_keep,
							dim_data = dim_data,
							dim_data_merged = dim_data_merged,
							res.hcpc = res.hcpc))
}



#' @importFrom FactoMineR CA dimdesc HCPC
compute_mfa <- function(input_data_list, 
						table_types = c("n"),
						target = NULL, 
						string_factors = NULL, 
						numeric_factors = NULL, 
						transpose = TRUE,
						add_samples = NULL,
						min_dimensions = 2) {

	table_lengths <- sapply(input_data_list, ncol)
	input_data <- Reduce(function(d1,d2) merge(d1, d2, all=T, by = "row.names"), input_data_list)

	rownames(input_data) <- input_data$Row.names
	input_data$Row.names <- NULL

	if (transpose) 
		input_data <- as.data.frame(t(input_data))
	 
	if (!is.null(target)) {
		rownames(target) <- as.character(target$sample)

		if (!is.null(numeric_factors)){
			input_data <- merge_factors(input_data, target, numeric_factors)

		}

		if (!is.null(string_factors)) {
		  input_data <- merge_factors(input_data, target, string_factors)
		}
	} 

	raw_input_data <- input_data[!rownames(input_data) %in% add_samples,
							 ! colnames(input_data) %in% c(numeric_factors, string_factors)]

	std_mfa <- FactoMineR::MCA(raw_input_data, 	graph = FALSE)                                                     
	dim_to_keep <- get_PCA_dimensions(std_mfa, min_dimensions = min_dimensions)

	
	if (!is.null(add_samples)) {
		add_samples_idx <- match(add_samples, rownames(mca_data))
	} else {
		add_samples_idx <- NULL
	}

	mca_res <- FactoMineR::MCA(mca_data, ncp = dim_to_keep, 
										graph = FALSE,
										quanti.sup = numeric_factors, 
										quali.sup=string_factors,
										ind.sup = add_samples_idx)	
 	dim_data <- FactoMineR::dimdesc(mca_res, axes=seq(1, dim_to_keep))
    dim_data_merged <- merge_dim_tables(dim_data)

  res.hcpc <- FactoMineR::HCPC(mca_res, graph = FALSE, nb.clust = -1)

	return(list(pca_data = mca_res,
							dim_to_keep = dim_to_keep,
							dim_data = dim_data,
							dim_data_merged = dim_data_merged,
							res.hcpc = res.hcpc))
}

get_cluster_string_assoc <- function(res.hcpc, string_factors){

	all_factor_clusters <- data.frame()
	for (add_factor in string_factors) {
		clusters_ct <- ct_from_qual(res.hcpc$data.clust[,c("clust", add_factor)])
		clusters_ct$sub <- paste(add_factor,clusters_ct$sub, sep = ":")
		colnames(clusters_ct)[match(c("ref","sub"),colnames(clusters_ct))] <- c("Cluster","Category")
		clusters_ct$fisher.test.pval <- v.fisher.test(clusters_ct)
		clusters_ct$FDR <- p.adjust(clusters_ct$fisher.test.pval, method = "BH")
		all_factor_clusters <- rbind(all_factor_clusters, clusters_ct)
	}
	return(all_factor_clusters)
}

parse_eigenvectors <- function(eigenvectors, eig_abs_cutoff = NULL) {

	parsed_eigvec <- list()
	for (Dim in names(eigenvectors)) {
		pDim <- gsub("Dim","PC",Dim)
		eigen_vec <- eigenvectors[[Dim]]
		
		if (!is.null(eig_abs_cutoff))
			eigen_vec <- eigen_vec[abs(eigen_vec) > eig_abs_cutoff]
			# eigen_vec[abs(eigen_vec) < eig_abs_cutoff] <- 0

		parsed_eigvec[[paste("positive", pDim, sep = ".")]] <- eigen_vec
		parsed_eigvec[[paste("negative", pDim, sep = ".")]] <- eigen_vec * -1
	}
	return(parsed_eigvec)
}

get_and_parse_pca_eigenvectors <- 	function(pca_res, cor_pval_cutoff = 0.05, eig_abs_cutoff = NULL){
	dimnames(pca_res$pca_data$svd$V) <- dimnames(pca_res$pca_data$var$cor)
	eigenvectors <- name_all_columns(as.data.frame(pca_res$pca_data$svd$V))
	correlations <- pca_res$pca_data$var$cor
	n <- nrow(pca_res$pca_data$ind$cos2)
	cor_sig <- cor_pval(correlations, n)
	eigenvectors <- filter_eigenvectors(eigenvectors, cor_sig, pval_cutoff = cor_pval_cutoff)
	eigenvectors <- parse_eigenvectors(eigenvectors, eig_abs_cutoff = eig_abs_cutoff)
	return(eigenvectors)

}

get_and_parse_clusters <- function(pca_res){
	clusters <- pca_res$res.hcpc$call$X
	dims <- colnames(clusters)[grepl("Dim.", colnames(clusters))]
	weights <- clusters |> dplyr::group_by(clust) |> dplyr::summarise_at(dims, mean)
	dimnames(pca_res$pca_data$svd$V) <- dimnames(pca_res$pca_data$var$cor)
	clust_names <- paste("Cluster",weights$clust )
	weights <- as.data.frame(dplyr::select(weights, -clust))
	rownames(weights) <- clust_names
	weighted_eigen <- apply(weights,1, get_clust_contributions, 
	                                   eigenvectors = pca_res$pca_data$svd$V)
	colnames(weighted_eigen) <- rownames(weights)
	weighted_eigen <- as.data.frame(weighted_eigen)
	weighted_eigen <- lapply(weighted_eigen, function(col) setNames(col, rownames(weighted_eigen)))
	return(weighted_eigen)
}

get_clust_contributions <- function(weigths, eigenvectors)  {
	weighted_eig <- weigths * eigenvectors
	rowSums(weighted_eig)
}


filter_eigenvectors <- function(eigenvectors, cor_sig, pval_cutoff) {
 	eigenvec_fil <- lapply(names(eigenvectors), function(dimension) {
		eigenvec <- eigenvectors[[dimension]]
		significance <- cor_sig[,dimension]
	  eigenvec[significance <= pval_cutoff] 
	})
	names(eigenvec_fil) <- names(eigenvectors)
	return(eigenvec_fil)
}