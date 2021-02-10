#' @importFrom utils read.table 
#' @importFrom miRBaseConverter getAllMiRNAs miRNA_AccessionToName
#' @importFrom data.table as.data.table rbindlist
#' @importFrom dplyr summarize group_by
#' @importFrom reshape2 melt
#' @importFrom rmarkdown render
miRNA_RNAseq_analysis <- function(
	RNAseq_folder,
	miRNAseq_folder,
	output_files,
	approaches,
	organism,
	multimir_db,
	p_val_cutoff,
	corr_cutoff,
	permutations, 
	module_membership_cutoff,
	report, #this option can be set to default value
	translation_file,
	databases,
	translate_ensembl = FALSE,
	filter_unmaintained,
	mc_cores = 1,
	template_folder = file.path(find.package('ExpHunterSuite'), "templates"),
	organism_table_path = file.path(find.package('ExpHunterSuite'), "inst", "external_data", "organism_table.txt") 
	){
	#create output folder
	"%>%" <- magrittr::"%>%"


		
	dir.create(output_files, recursive = TRUE)

	#parse strategies
	approaches <- c("normalized_counts_RNA_vs_miRNA_normalized_counts", parse_approaches(unlist(strsplit(approaches, ","))))

	#translate miRNA ID
	# translation_miRNA <- translation_file
	# if (!is.null(translation_miRNA)) {
	# 	translation_miRNA <- utils::read.table(translation_miRNA, sep = "\t")
	# 	colnames(translation_miRNA) <- c("mature_ID", "input_ID")
	# }
	# Prepare for RNA ID translation
	organism_info <- utils::read.table(organism_table_path, header = TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE, fill = NA) 
	organism_info <- subset(organism_info, organism_info$KeggCode %in% organism)  


	#Prepare multiMiR
	# if (!is.null(multimir_db)){
		load(multimir_db)
		selected_databases <- strsplit(opt$databases, ",")
		# save(selected_databases, file = file.path(output_files, "test.RData"))
		# q()
			multimir_summary <- summarize_multimir(multimir_summary, selected_predicted_databases = selected_databases, filter_unmaintained= filter_unmaintained)
		message("multiMiR database has been parsed and summarized")
	# }
	# save(multimir_summary, file = file.path("/mnt/scratch/users/bio_267_uma/josecordoba/NGS_projects/pmm2_belen/target_miRNA_2020", "test2.RData"))
	# q()

	#Load and prepare DGH data
	RNAseq <- load_DEGH_information(RNAseq_folder)
	miRNAseq <- load_DEGH_information(miRNAseq_folder)

	#prepare all possible pairs
	strategies <- list()


	strategies$all_possible_pairs <- expand.grid(colnames(RNAseq$normalized_counts), colnames(miRNAseq$normalized_counts))
	names(strategies$all_possible_pairs) <- c("RNAseq", "miRNAseq")
	strategies$all_possible_pairs <- add_multimir_info(strategies$all_possible_pairs, multimir_summary = multimir_summary)
	

	#Filter DGH data
	RNAseq <- filter_DEGH_data(RNAseq, module_membership_cutoff)
	miRNAseq <- filter_DEGH_data(miRNAseq, module_membership_cutoff)
	# save(RNAseq, miRNAseq, file = "/mnt/scratch/users/bio_267_uma/josecordoba/NGS_projects/LaforaRNAseq/target_miRNA/test.RData")
	# q()

	# strategies$all_possible_pairs[ ,c("correlation", "pval", "all_permutations", "correlated_pairs")] <- NULL

	# #all_possible_pairs: Is a scaffold data.table with all posible permutations between RNA and miRNA expressed in samples. It also include multimiR info. The pairs of each strategy must be reordered according this table.
	# strategies$all_possible_pairs <- perform_correlations(strategy = "normalized_counts_RNA_vs_miRNA_normalized_counts", 
	# 																RNAseq = RNAseq, 
	# 																miRNAseq = miRNAseq, 
	# 																corrected_positions = NULL, 
	# 																multimir_summary = multimir_summary, 
	# 																cor_cutoff = corr_cutoff, 
	# 																pval_cutoff = p_val_cutoff
	# 																)
	# strategies$normalized_counts_RNA_vs_miRNA_normalized_counts <- strategies$all_possible_pairs[ ,c("correlation", "pval", "all_permutations", "correlated_pairs")] 
	# message(paste0("\nnormalized_counts_RNA_vs_miRNA_normalized_counts"," strategy has been performed"))

	##### DEBUG_CONTROL A
		time_control <- Sys.time()

	gc()
	corrected_positions <- paste0(strategies$all_possible_pairs$RNAseq, "_", strategies$all_possible_pairs$miRNAseq)
	names(corrected_positions) <- seq(length(corrected_positions))
	selected_strategies <- lapply(approaches, function(approach){ #this is not a for loop because it can be parallelized with mclapply
		actual_strategy <- perform_correlations(strategy = approach, 
												RNAseq = RNAseq, 
												miRNAseq = miRNAseq, 
												corrected_positions = corrected_positions, 
												cor_cutoff = corr_cutoff, 
												pval_cutoff = p_val_cutoff)
		message(paste0("\n", approach, " strategy has been performed"))
			gc()
		return(actual_strategy)
	})
	names(selected_strategies) <- approaches
	corrected_positions <- NULL


	strategies <- c(strategies, selected_strategies)
	# approaches <- c(approaches, "normalized_counts_RNA_vs_miRNA_normalized_counts")
	# save(strategies, file = file.path(output_files, "test4.RData"))
	# q()
	all_known_miRNAs <- miRBaseConverter::getAllMiRNAs(version = "v22", type = "all", species = organism)$Accession
	strategies$all_possible_pairs$known_miRNA <- strategies$all_possible_pairs$miRNAseq %in% all_known_miRNAs

	#possible_positives: Pairs with RNA and miRNA included in multiMiR. 
	strategies$all_possible_pairs$possible_positives <- strategies$all_possible_pairs$known_miRNA & 
																strategies$all_possible_pairs$miRNAseq %in% multimir_summary$mature_mirna_acc &
																  strategies$all_possible_pairs$RNAseq %in% multimir_summary$target_ensembl 
	multimir_summary <- NULL



	# q()
	all_strategies <- list()
	prec_recall <- list()
	score_filter <- list()


	strategies$all_possible_pairs[is.na(strategies$all_possible_pairs$predicted_c), "predicted_c"] <- 0 
	strategies$all_possible_pairs[is.na(strategies$all_possible_pairs$validated_c), "validated_c"] <- 0

	message("PARSING STRATEGIES")

	for (strategy in approaches) {
		# message(strategy)
		# message(paste0(sum(strategies[[strategy]]$correlated_pairs), " significant pairs and ", sum(strategies$all_possible_pairs$known_miRNA), " pairs in background"))

		if (sum(strategies[[strategy]]$correlated_pairs) == 0 || sum(strategies$all_possible_pairs[strategies[[strategy]]$pair_n, "known_miRNA"]) == 0 ){
			strategies[[strategy]] <- NULL
			message("No significant pairs can be used for multimir comparison in strategy:")
			# message(strategy)
			next 
		}
		all_pairs_info <- data.table::as.data.table(strategies[[strategy]])

		### Preparing precision, recall & F1 for correlated pairs
		# message("Getting prediction stats")
		# str(all_pairs_info)
		prec_recall[[strategy]] <- get_prediction_stats(all_pairs_info = all_pairs_info, 
													all_possible_pairs = strategies$all_possible_pairs)
		
		### get Ranking Scores 
		# message("Getting ranking scores")
		all_db_info <- strategies$all_possible_pairs[all_pairs_info$pair_n, c("diana_microt", "elmmo", "microcosm", "miranda","mirdb", "pictar", "pita", "targetscan")]
		all_db_info$correlated_pairs <- all_pairs_info$correlated_pairs
		# score_filter[[strategy]]$significant <- all_pairs_info$correlated[all_pairs_info$pair_n]
		
		score_filter[[strategy]] <- get_db_scores(all_db_info = all_db_info)
		# all_pairs_info$pair_n <- rownames(all_pairs_info)
		all_pairs_info$score <- strategies$all_possible_pairs[all_pairs_info$pair_n, "score"]
		all_strategies[[strategy]] <- data.table::as.data.table(all_pairs_info[all_pairs_info$correlated_pairs,])

	}

	message("MERGING RESULTS")

	# save(prec_recall, file = file.path(output_files, "test.RData"))
	# q()
	all_strategies <- as.data.frame(data.table::rbindlist(all_strategies, use.names=TRUE, idcol= "strategy"))
	prec_recall <- as.data.frame(data.table::rbindlist(prec_recall, use.names = TRUE, idcol = "strategy"))
	all_strategies <- cbind(all_strategies, strategies$all_possible_pairs[as.numeric(all_strategies$pair_n), c("known_miRNA", "predicted_c", "validated_c")])

	# ##### DEBUG_CONTROL B
	# time_control <- Sys.time() - time_control
	# debug_point(file.path(output_files, "Debug_mirna_new.RData"),message = "Debug point",envir = environment())
	# q()	

	filters_summary <- all_strategies %>%
						dplyr::group_by(strategy) %>%
						dplyr::summarize(known_miRNAs = sum(known_miRNA == TRUE),
									novel_miRNAs = sum(known_miRNA == FALSE),
									predicted = sum(predicted_c > 0),
									validated = sum(validated_c > 0), 
									both  = sum(predicted_c > 0 & validated_c > 0))
	filters_summary <- as.data.frame(filters_summary, stringsAsFactors = FALSE)
	filters_summary <- reshape2::melt(filters_summary)
	names(filters_summary) <- c("strategy", "type", "pairs")
	filters_summary$type <- as.character(filters_summary$type)
	filters_summary$sdev <- NA
	filters_summary$p_val <- NA
	filters_summary$quantile <- NA


	#generate and merge randoms

	db_distribution <- data.table::data.table(strategy = all_strategies[all_strategies$predicted_c > 0, "strategy"],
									step = "predicted",
									score = all_strategies[all_strategies$predicted_c > 0, "score"])
	message("ADDING RANDOMS")
	random_obj <- add_randoms(background = strategies$all_possible_pairs, 
									filters_summary = filters_summary,
									db_distribution = db_distribution,
									permutations = permutations,
									all_strategies = all_strategies)

	filters_summary <- random_obj[["filters_summary"]]
	db_distribution <- as.data.frame(random_obj[["db_distribution"]])


	random_obj <- NULL
	filters_summary$type <- factor(filters_summary$type, levels=c("novel_miRNAs","known_miRNAs","predicted", "predicted_random", "validated","validated_random", "both", "both_random"))

# save(list = ls(all.names = TRUE), file = file.path(output_files, "test2.RData"), envir = environment())

	# save(all_strategies, file = file.path(output_files, "test2.RData"))
	################# GET IDS TARNSALTIONS

	mirna_names <- unique(strategies$all_possible_pairs$miRNAseq)
	mirna_names <- mirna_names[grepl("MIMAT", mirna_names)]
	mirna_names <- miRBaseConverter::miRNA_AccessionToName(mirna_names, targetVersion = "v22")
	gene_id_translation <- NULL
	if (! organism_info$Bioconductor_VarName_SYMBOL[1] == "" && translate_ensembl) {
		gene_id_translation <- get_entrez_symbol_translations(ensembl_ids = unique(strategies$all_possible_pairs$RNAseq), 
															organism_info = organism_info)
	} 


	miRNA_cor_results <- list(
		all_strategies = all_strategies,
		filters_summary = filters_summary,
		miRNAseq = miRNAseq,
		RNAseq = RNAseq,
		strategies = strategies,
		approaches = approaches,
		db_distribution = db_distribution,
		prec_recall = prec_recall,
		RNAseq_folder = RNAseq_folder,
		miRNAseq_folder = miRNAseq_folder,
		corr_cutoff = corr_cutoff,
		mirna_names = mirna_names,
		gene_id_translation = gene_id_translation,
		# ref_strategy = ref_strategy,
		score_filter =score_filter
		 
	)

	return(miRNA_cor_results)
}

