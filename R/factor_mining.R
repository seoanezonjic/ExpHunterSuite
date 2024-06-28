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
						add_samples = NULL) {

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
	dim_to_keep <- get_PCA_dimensions(std_pca)

	
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


