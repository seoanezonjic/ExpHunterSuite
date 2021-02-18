################### FUNCTIONS

#' @importFrom utils read.table
#' @importFrom stringr str_remove
load_DEGH_information <- function(execution_path){ 
	execution <- list()

	execution[['DH_results']] <- utils::read.table(file.path(execution_path, "Common_results/hunter_results_table.txt"), header=TRUE, row.names=1, sep="\t")
	execution[['DH_results']]$gene_name <- rownames(execution[['DH_results']])
	rownames(execution[['DH_results']]) <- NULL
	execution$DH_results <- execution$DH_results[execution$DH_results$gene_name != "id",]

	execution[['normalized_counts']] <- as.matrix(utils::read.table(file.path(execution_path, "Results_DESeq2/Normalized_counts_DESeq2.txt"), header=TRUE, row.names=1, sep="\t"))


	execution[['Eigengene']] <- as.matrix(utils::read.table(file.path(execution_path, "Results_WGCNA/eigen_values_per_samples.txt"), header=TRUE, row.names=1, sep="\t"))
	colnames(execution[["Eigengene"]]) <- stringr::str_remove(colnames(execution[["Eigengene"]]), "ME")
 	
	execution[['normalized_counts']] <- t(execution[['normalized_counts']]) 

	execution[['hub_1']] <- get_hub_genes_by_MM(execution[['normalized_counts']], execution[['DH_results']])
	return(execution)
}

# translate_genes <- function(){
	# if (!is.null(translation_table)) {
	# # 	#translate DH_results
	# 	temp_translation <- as.vector(translation_table[match(execution[['DH_results']]$gene_name, translation_table$input_ID), "mature_ID"])
	# 	if (!only_known) {
	# 		for (element in seq(1:length(temp_translation))){ 
	# 			if(is.na(temp_translation[element])) {
	# 				temp_translation[element] <- execution[["DH_results"]][element, "gene_name"]
	# 			}
	# 		}
	# 	}
	# 	execution[["DH_results"]]$gene_name <- temp_translation
	# 	execution[["DH_results"]] <- execution[["DH_results"]][!is.na(execution[["DH_results"]]$gene_name), ]

	# 	#translate normalized counts
	# 	temp_translation <- as.vector(translation_table[match(rownames(execution[['normalized_counts']]), translation_table$input_ID), "mature_ID"])
	# 	if (!only_known) {
	# 		for (element in seq(1:length(temp_translation))){ 
	# 			if(is.na(temp_translation[element])) {
	# 				temp_translation[element] <- rownames(execution[['normalized_counts']])[element]
	# 			}
	# 		}
	# 	}

	# 	rownames(execution[['normalized_counts']]) <- temp_translation
	# 	execution[['normalized_counts']] <- execution[['normalized_counts']][!is.na(rownames(execution[['normalized_counts']])), ]
		
	# } 	

# }

filter_DEGH_data <- function(DGH_data, MM_cutoff){
	# save(list = ls(all.names = TRUE), file = "/mnt/scratch/users/bio_267_uma/josecordoba/NGS_projects/LaforaRNAseq/target_miRNA/test.RData", envir = environment())

	DGH_results <- DGH_data$DH_results 
	all_degs <- DGH_results$genes_tag %in% c("PREVALENT_DEG", "POSSIBLE_DEG")
	modules_with_DEGS <- unique(DGH_results[all_degs & DGH_results$Cluster_MM >= MM_cutoff, "Cluster_ID"])

   	candidate_not_deg <- DGH_results$Cluster_ID %in% modules_with_DEGS &
   							 DGH_results$genes_tag == "NOT_DEG" &
   							  DGH_results$Cluster_MM >= MM_cutoff

	candidate_not_deg <- DGH_results[candidate_not_deg, "gene_name"]

	relevant_genes <- c(candidate_not_deg, DGH_results[all_degs, "gene_name"]) 

	DGH_data$DH_results$relevant_genes <- DGH_results$gene_name %in% relevant_genes
	
	DGH_data$normalized_counts <- DGH_data$normalized_counts[, colnames(DGH_data$normalized_counts) %in% relevant_genes ]


	DGH_data$Eigengene_0 <- DGH_data$Eigengene[,"0", drop = FALSE]

	DGH_data$Eigengene <- DGH_data$Eigengene[,colnames(DGH_data$Eigengene) %in% modules_with_DEGS]
	DGH_data$hub_1 <- DGH_data$hub_1[,colnames(DGH_data$hub_1) %in% modules_with_DEGS]

	return(DGH_data)
}


summarize_multimir <- function(
	multimir_table, 
	selected_predicted_databases = c("diana_microt", "elmmo", "microcosm", "miranda","mirdb", "pictar", "pita", "targetscan"),
	filter_unmaintained = 0
	) {
	# save(multimir_table,selected_predicted_databases, file = file.path("/mnt/scratch/users/bio_267_uma/josecordoba/NGS_projects/pmm2_belen/target_miRNA_2020", "test1.RData"))

	multimir_table <- as.data.frame(multimir_table)
	pred_dbs <- c("diana_microt", "elmmo", "microcosm", "miranda","mirdb", "pictar", "pita", "targetscan")

	dbs <- list(
		"validated_c" = c("mirecords", "mirtarbase", "tarbase"),
		"predicted_c" = pred_dbs[pred_dbs %in% unlist(selected_predicted_databases)]
		)
	# for (db_type in names(dbs)) {
	multimir_table[["validated_c"]] <- rowSums(multimir_table[,dbs[["validated_c"]]], na.rm = TRUE)
	# }
	for (pred_db in dbs$predicted_c) {
		multimir_table[is.na(multimir_table[,pred_db]), pred_db] <- 0
	}
	multimir_table[["predicted_c"]] <- rowSums(multimir_table[,dbs[["predicted_c"]]] > 0, na.rm = TRUE)

	multimir_table$score <- (rowSums(multimir_table[,dbs$predicted_c], na.rm = TRUE) / multimir_table$predicted_c) 
	multimir_table[is.na(multimir_table$score), "score"] <- 0


	gc()
	if ( filter_unmaintained > 0) {
		unmaintained_dbs <- c("diana_microt", "elmmo", "microcosm", "miranda", "pictar", "pita")
		aux_multimir_table <- multimir_table[,unmaintained_dbs[unmaintained_dbs %in% unlist(selected_predicted_databases) ]]
		# save(list = ls(all.names = TRUE), file = "/mnt/scratch/users/bio_267_uma/josecordoba/NGS_projects/pmm2_belen/target_miRNA_2020/test.RData", envir = environment())
		# q()
		aux_multimir_table$unmaintained_count <- rowSums(aux_multimir_table)
		multimir_table <- multimir_table[aux_multimir_table$unmaintained_count >= filter_unmaintained |
											 multimir_table$predicted_c > 0,]
	}

	return(multimir_table)
}

parse_approaches <- function(approaches){

	strat_equivalence <- list(
		"d" = "normalized_counts",
		"E" = "Eigengene",
		"h" = "hub_1")
	parsed_approaches <- c()
	for (approach in approaches) {
		splitted_approach <- unlist(strsplit(approach, ""))
		RNA_profile <- strat_equivalence[[splitted_approach[1]]]
		miRNA_profile <- strat_equivalence[[splitted_approach[2]]]
		parsed_approaches <- c(parsed_approaches, paste0(RNA_profile, "_RNA_vs_miRNA_", miRNA_profile))
	}
	return(parsed_approaches)
}


# filter_correlations

# tag_strategy <- function(all_pairs_info, tag){
# 	all_pairs_info$strategy <- rep("zero", nrow(all_pairs_info))
# 	return(all_pairs_info)
# }

# translate_column <- function(column, translation_table) { #translation_table, a dataframe with translated_ids in first column and ids_to_translate in second column
# 	column <- translation_table[match(column, translation_table[,2]), 1]
# 	return(column)
# }


perform_correlations <- function(strategy = "normalized_counts_RNA_vs_miRNA_normalized_counts", RNAseq, miRNAseq, corrected_positions, cor_cutoff = 0, pval_cutoff = 0.05){ #correct_positions is a mirna_RNA pairs vector

	
	if (strategy == "DEGs_vs_DEMs") {
		DEGS <- miRNAseq$DH_results[miRNAseq$DH_results$genes_tag %in% c("PREVALENT_DEG", "POSSIBLE_DEG"), "gene_name"]
		DEMS <- miRNAseq$DH_results[miRNAseq$DH_results$genes_tag %in% c("PREVALENT_DEG", "POSSIBLE_DEG"), "gene_name"]
		all_pairs_info <- expand.grid(DEGS, DEMS, stringsAsFactors= FALSE, KEEP.OUT.ATTRS = FALSE)
		colnames(all_pairs_info) <- c("RNAseq", "miRNAseq")
		# all_pairs_info$RNAseq <- as.character(all_pairs_info$RNAseq)
		# all_pairs_info$miRNAseq <- as.character(all_pairs_info$miRNAseq)
		all_pairs_info$correlation <- -1
		all_pairs_info$pval <- 0
		
	} else {


		strat_description <- unlist(strsplit(strategy, "_RNA_vs_miRNA_"))

		RNA_profiles <- RNAseq[[strat_description[1]]]
		miRNA_profiles <- miRNAseq[[strat_description[2]]]
		# print(strategy)
		print(str(RNA_profiles))
		print(str(miRNA_profiles))

		# save(miRNA_profiles, file = "/mnt/scratch/users/bio_267_uma/josecordoba/NGS_projects/LaforaRNAseq/target_miRNA/test.RData")
		if (ncol(RNA_profiles) == 0 || ncol(miRNA_profiles) == 0) {
			return(data.frame(NULL))
		}
		all_pairs_info <- correlate_profiles(RNA_profiles, miRNA_profiles)
		all_pairs_info <- data.table::as.data.table(all_pairs_info)

		if (strat_description[1] != "normalized_counts"){
			# message("expand RNA")
			relevant_genes <- RNAseq$DH_results$relevant_genes
			if (strat_description[1] == "Eigengene_0") {
				relevant_genes <- rep(TRUE, length(relevant_genes))
			}
			all_pairs_info <- expand_module(all_pairs_info = all_pairs_info, tag = "RNAseq", DH_results = RNAseq$DH_results[relevant_genes, ])
		}
		
		if (strat_description[2] != "normalized_counts"){
			# message("expand miRNA")
			relevant_genes <- miRNAseq$DH_results$relevant_genes
			if (strat_description[2] == "Eigengene_0") {
				relevant_genes <- rep(TRUE, length(relevant_genes))
			}

			all_pairs_info <- expand_module(all_pairs_info = all_pairs_info, tag = "miRNAseq", DH_results = miRNAseq$DH_results[relevant_genes, ])
		}
	}

	all_pairs_info$correlated_pairs <- all_pairs_info$correlation <= cor_cutoff & all_pairs_info$pval < pval_cutoff
	print(sum(all_pairs_info$correlated_pairs))
	if (!is.null(corrected_positions)){
		all_pairs <- paste0(all_pairs_info$RNAseq, "_", all_pairs_info$miRNAseq)

		all_pairs_info$pair_n <- match(all_pairs, corrected_positions) 
		# all_pairs_info[is.na(all_pairs_info$all_permutations), "all_permutations"] <- FALSE
		# all_pairs_info[is.na(all_pairs_info$correlated_pairs), "correlated_pairs"] <- FALSE
		# rownames(all_pairs_info) <- 1:nrow(all_pairs_info)	
		all_pairs_info[,c("RNAseq", "miRNAseq")] <- NULL
	}
	return(all_pairs_info)
}


#' @importFrom WGCNA corPvalueStudent
#' @importFrom stats cor
correlate_profiles <- function(RNA_profiles, miRNA_profiles) {
	# str(RNA_profiles)
	# str(miRNA_profiles)

	correlations <- WGCNA::cor(RNA_profiles, miRNA_profiles)
	
	nSamples <- ncol(RNA_profiles) 
	cor_pval <- as.data.frame(WGCNA::corPvalueStudent(as.matrix(correlations), nSamples))

	names_matrix <- as.matrix(correlations)  
	all_pairs_info <- as.data.frame(as.table(correlations)) 
	colnames(all_pairs_info) <- c("RNAseq", "miRNAseq", "correlation")
	
	all_pairs_info$correlation <- round(all_pairs_info$correlation, 2)

	all_pairs_info$pval <- as.data.frame(as.table(as.matrix(cor_pval)))$Freq

	all_pairs_info$miRNAseq <- as.character(all_pairs_info$miRNAseq)
	all_pairs_info$RNAseq <- as.character(all_pairs_info$RNAseq)
	return(all_pairs_info)
}



#' @importFrom dplyr desc between row_number filter arrange group_by
get_hub_genes_by_MM <- function(normalized_counts, hunter_results, top = 1){
	"%>%" <- magrittr::"%>%"
	# hunter_results <- hunter_results[hunter_results$relevant_genes,]

	hub_genes <- hunter_results %>% 
				dplyr::filter(Cluster_MM_pval <= 0.05) %>% 
				dplyr::arrange(dplyr::desc(Cluster_MM)) %>% 
				dplyr::group_by(Cluster_ID) %>% dplyr::filter(dplyr::between(dplyr::row_number(), 1, top))

	hub_genes_profile <- as.matrix(normalized_counts[,colnames(normalized_counts) %in% hub_genes$gene_name])
	
	colnames(hub_genes_profile) <- as.data.frame(hub_genes)[match(colnames(hub_genes_profile), hub_genes$gene_name), "Cluster_ID"]
	return(hub_genes_profile)
}


#' @importFrom data.table as.data.table merge.data.table  
#' @importFrom dplyr select
expand_module <- function(all_pairs_info, tag, DH_results){
		"%>%" <- magrittr::"%>%"
		# save(list = ls(all.names = TRUE), file = "/mnt/scratch/users/bio_267_uma/josecordoba/NGS_projects/pmm2_belen/target_miRNA_2020/test_3.RData", envir = environment())
		# q()

		mod_tag <- paste0(tag, "_mod")
		names(all_pairs_info)[names(all_pairs_info)== tag] <- mod_tag

		partial_expanded <- DH_results %>% dplyr::select(Cluster_ID, gene_name)
		partial_expanded <- data.table::as.data.table(partial_expanded)
		colnames(partial_expanded) <- c("module", tag)
		# print(head(RNAseq))
		partial_expanded$module <- as.character(partial_expanded$module)
		all_pairs_info <- data.table::as.data.table(all_pairs_info)
		all_pairs_info <- data.table::merge.data.table(x = all_pairs_info, y = partial_expanded, by.x = mod_tag, by.y = "module", allow.cartesian  = TRUE)
		all_pairs_info[,mod_tag] <- NULL
		return(all_pairs_info)
}

#' @importFrom data.table as.data.table merge.data.table
#' @importFrom dplyr select
add_multimir_info <- function(all_pairs_info, multimir_summary = NULL) {
	all_pairs_info <- data.table::as.data.table(all_pairs_info)
		
	if (!is.null(multimir_summary)){
		multimir_summary <- data.table::as.data.table(multimir_summary)
		all_pairs_info <- data.table::merge.data.table(all_pairs_info, multimir_summary, by.x = c("RNAseq","miRNAseq"), by.y =c("target_ensembl","mature_mirna_acc"), all.x = TRUE, suffixes = c("", ""), no.dups = TRUE)
	}
	all_pairs_info <- as.data.frame(all_pairs_info)
	for (pred_db in c("predicted_c", "validated_c","score","diana_microt", "elmmo", "microcosm", "miranda","mirdb", "pictar", "pita", "targetscan")) {
		all_pairs_info[is.na(all_pairs_info[,pred_db]), pred_db] <- 0

	}
	for (pred_db in c("mirecords", "mirtarbase", "tarbase")) {
		all_pairs_info[is.na(all_pairs_info[,pred_db]), pred_db] <- FALSE
	}
	all_pairs_info$predicted <- all_pairs_info$predicted_c > 0
	all_pairs_info$validated <- all_pairs_info$validated_c > 0
	all_pairs_info$pred_and_val <- all_pairs_info$predicted & all_pairs_info$validated
	all_pairs_info$multimir <- all_pairs_info$predicted | all_pairs_info$validated
	gc()
	return(all_pairs_info)
}

#' @importFrom data.table as.data.table rbindlist
get_db_scores <- function(all_db_info){
		all_db_info <- as.data.frame(all_db_info)
		all_pred_dbs <- c("diana_microt", "elmmo", "microcosm", "miranda","mirdb", "pictar", "pita", "targetscan")
		strats_r_scores <- lapply(all_pred_dbs, function(pred_db){
			pred_db_score <- data.table::data.table(database = pred_db, r_scores = all_db_info[, pred_db], significant = all_db_info$correlated_pairs)
			# print("1")
			pred_db_score <- pred_db_score[pred_db_score$r_scores > 0, ]
			# print("2")
			return(pred_db_score) 
		})
		strats_r_scores <- data.table::rbindlist(strats_r_scores)
		return(as.data.frame(strats_r_scores))
}


get_prediction_stats <- function(all_pairs_info, all_possible_pairs){
		strategy_correlated_pairs <- rep(FALSE, nrow(all_possible_pairs))
		strategy_correlated_pairs[all_pairs_info$pair_n[all_pairs_info$correlated_pairs]] <- TRUE
		# save(list = ls(all.names = TRUE), file = file.path("/mnt/scratch/users/bio_267_uma/josecordoba/NGS_projects/pmm2_belen/target_miRNA_2020", "test.RData"), envir = environment())
		# q()

		# str(all_pairs_info$pair_n)
		# print(paste0(length(strategy_correlated_pairs), " teeeeest ", nrow(all_possible_pairs)))

		correlated_in_multimir <- strategy_correlated_pairs & all_possible_pairs$possible_positives
		
		predicted_pairs <- all_possible_pairs$predicted_c > 0
		validated_pairs <- all_possible_pairs$validated_c > 0
		both_pairs <- predicted_pairs & validated_pairs
		# save(list = ls(all.names = TRUE), file = "/mnt/scratch/users/bio_267_uma/josecordoba/NGS_projects/pmm2_belen/target_miRNA_2020/test.RData", envir = environment())
		# q()
		## Get precision recall F1 stats
		stats <- data.frame(stringsAsFactors = FALSE,
							predicted = pred_stats(correlated_in_multimir, predicted_pairs),
							validated = pred_stats(correlated_in_multimir, validated_pairs),
							both = pred_stats(correlated_in_multimir, both_pairs)
							)

		rownames(stats) <- c("precision", "recall", "F1")
		stats <- as.data.frame(t(as.matrix(stats)))
		stats$gold_standard <- rownames(stats)
		return(stats)
}


get_module_genes <- function(hunter_results, module_name, column_name){
	module_genes <- hunter_results[["gene_name"]][hunter_results[["Cluster_ID"]] == module_name]

	module_genes <- data.frame(module = rep(module_name, length(module_genes)), genes = module_genes)
	colnames(module_genes) <- c("module", column_name)
	return(module_genes)
}

#' @importFrom data.table data.table rbindlist as.data.table
#' @importFrom stats sd 
add_randoms <- function(background = NULL, filters_summary = NULL, permutations = 10, all_strategies, db_distribution = data.table(stringsAsFactors = FALSE)){
	background <- data.table::as.data.table(background)
	# str(filters_summary)
	db_distribution <- list("0" = db_distribution)
	for (strategy_name in unique(filters_summary$strategy)) {
		random_dist <- list()
		message("debug_1")	
		# sig_pairs <- rep(FALSE, nrow(background))
		# sig_pairs[as.numeric(all_strategies[all_strategies$strategy == strategy_name, "pair_n"])] <- TRUE
		# strat_background <- background[!sig_pairs,]
		sig_pairs_count <- filters_summary[filters_summary$strategy == strategy_name & filters_summary$type == "known_miRNAs", "pairs"]
		# print(sig_pairs_count)
		# save(list = ls(all.names = TRUE), file = "/mnt/scratch/users/bio_267_uma/josecordoba/NGS_projects/pmm2_belen/target_miRNA_2020/test.RData", envir = environment())
		# q()
		for( i in 1:permutations){ 
			message("\tdebug1_1")
			random_indices <- sample(nrow(background), size = sig_pairs_count, replace = FALSE)	
			random_set <- as.data.frame(background[random_indices,])
			random_summary <- data.table::data.table(multiMiR = sum(random_set$predicted_c > 0 | random_set$validated_c > 0),
											predicted = sum(random_set$predicted_c > 0),
											validated = sum(random_set$validated_c > 0),
											both= sum(random_set$predicted_c > 0 & random_set$validated_c > 0)
										)
			random_dist[[i]] <- random_summary
			
			random_score_dist <- data.table::data.table(strategy = strategy_name,
											step = "predicted_random",
											score = random_set[random_set$predicted_c > 0, "score"] )
			
			db_distribution <- c(db_distribution, list(random_score_dist))
			
		}
		random_dist <- as.data.frame(data.table::rbindlist(random_dist))
		
		message("debug_2")

		random_summary <- lapply(c("multiMiR", "predicted", "validated", "both"), function(type){
			type_dist <- random_dist[, type]
			type_pairs <- filters_summary[filters_summary$strategy == strategy_name & filters_summary$type == type, "pairs"]
			data.frame(type = paste0(type, "_random"),
						pairs = mean(type_dist),
						sdev = stats::sd(type_dist),
						strategy = strategy_name,
						p_val = calc_pval(type_pairs, type_dist),
						quantile = calc_quantile(type_pairs, type_dist)
						)
		})
		
		message("debug_4")
		

		filters_summary <- as.data.frame(data.table::rbindlist(c(list(filters_summary), random_summary), use.names=TRUE))	
	}
	db_distribution <- data.table::rbindlist(db_distribution, use.names=TRUE)
	gc()
	random_obj <- list("filters_summary" = filters_summary,
							"db_distribution" = db_distribution)
	return(random_obj)
}


#' @importFrom stats sd pnorm
calc_pval <- function(value, distribution){
	d_mean <- mean(distribution)
	d_sd <- stats::sd(distribution)
	pval <- 1 - stats::pnorm(value, mean = d_mean, sd = d_sd)
	return(pval)
}

#' @importFrom stats ecdf
calc_quantile <- function(value, distribution){
	Fn <-  stats::ecdf(distribution)
	quantile <- 1 - Fn(value)
	return(quantile)
}

pred_stats <- function(predicted, condition){ #takes 2 true/false vectors and returns another vector c(precision, recall, f measure)
	# start_time <- Sys.time()	
	tpositives <- sum(predicted & condition)
	fnegatives <- sum(!predicted & condition)
	fpositives <- sum(predicted & !condition)
	precision <- tpositives / (tpositives + fpositives)
	recall <- tpositives / (tpositives + fnegatives)
	fmeasure <- calc_f_measure(precision, recall)
	# print(Sys.time() - start_time)
	return(c(precision, recall, fmeasure))
}


calc_f_measure <- function(precision, recall) { #precision and recall must be numeric vectors or numbers
	f_measure <- ( 2 * precision * recall) / (precision + recall)
	return(f_measure) 
}

#' @importFrom data.table as.data.table merge.data.table
get_entrez_symbol_translations <- function(ensembl_ids, organism_info){

	input_to_entrezgene <- id_translation_orgdb(input_ids = ensembl_ids,
												 organism_db = organism_info$Bioconductor_DB[1],
												 org_var_name = organism_info$Bioconductor_VarName[1])

	colnames(input_to_entrezgene) <- c("ensembl_gene_id", "entrezgene") # Fix names


	input_to_symbol <- id_translation_orgdb(input_ids = input_to_entrezgene$entrezgene,
										 		organism_db = organism_info$Bioconductor_DB[1],
												org_var_name = organism_info$Bioconductor_VarName_SYMBOL[1])
	colnames(input_to_symbol) <- c("entrezgene", "Symbol")
	input_to_entrezgene <- data.table::as.data.table(input_to_entrezgene)
	input_to_symbol <- data.table::as.data.table(input_to_symbol)
	gene_id_translation <- data.table::merge.data.table(input_to_entrezgene, input_to_symbol, by = "entrezgene", all.x = TRUE)
	return(gene_id_translation)
}
