#' @importFrom utils read.table 
#' @importFrom miRBaseConverter getAllMiRNAs miRNA_AccessionToName
#' @importFrom data.table as.data.table rbindlist
#' @importFrom dplyr summarize group_by
#' @importFrom reshape2 melt
#' @importFrom rmarkdown render
miRNA_RNAseq_analysis <- function(
	RNAseq_folder,
	miRNAseq_folder,
	deg_tag,
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
	translate_ensembl = FALSE,
	mc_cores = 1,
	template_folder = file.path(find.package('DEgenesHunter'), "templates"),
	organism_table_path = file.path(find.package('DEgenesHunter'), "inst", "external_data", "organism_table.txt") 
	){
#create output folder
dir.create(output_files, recursive = T)

#parse strategies
approaches <- parse_approaches(unlist(strsplit(approaches, ",")))

#tranbslate miRNA ID
translation_miRNA <- translation_file
if (!is.null(translation_miRNA)) {
	translation_miRNA <- utils::read.table(translation_miRNA, sep = "\t")
	colnames(translation_miRNA) <- c("mature_ID", "input_ID")
}
#load and prepare files
RNAseq <- load_DEGH_information(RNAseq_folder, MM_cutoff = module_membership_cutoff)

miRNAseq <- load_DEGH_information(miRNAseq_folder, translation_table =  translation_miRNA, MM_cutoff = module_membership_cutoff)
organism_info <- utils::read.table(organism_table_path, header = TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE, fill = NA) 
organism_info <- subset(organism_info, organism_info$KeggCode %in% organism)  

#summarize multimir
if (!is.null(multimir_db)){
	multimir_path <- file.path(multimir_db, paste0("parsed_", organism, ".RData"))
	load(multimir_path)
	multimir_summary <- summarize_multimir(multimir_summary)
	message("multiMiR database has been parsed and summarized")
}

#perform dd strategy and prepare main object
strategies <- list()
#all_possible_permutations: Is a scaffold data.table with all posible permutations between RNA and miRNA expressed in samples. It also include multimiR info. The pairs of each strategy must be reordered according this table.
strategies$all_possible_permutations <- perform_correlations(strategy = "normalized_counts_RNA_vs_miRNA_normalized_counts", 
																RNAseq = RNAseq, 
																miRNAseq = miRNAseq, 
																corrected_positions = NULL, 
																multimir_summary = multimir_summary, 
																cor_cutoff = corr_cutoff, 
																pval_cutoff = p_val_cutoff
																)
strategies$normalized_counts_RNA_vs_miRNA_normalized_counts <- strategies$all_possible_permutations[ ,c("correlation", "pval", "all_permutations", "correlated_pairs")] 
strategies$all_possible_permutations[ ,c("correlation", "pval", "all_permutations", "correlated_pairs")] <- NULL
print(paste0("dd"," strategy has finished"))

print(gc())

selected_strategies <- lapply(approaches, function(approach){ #this is not a for loop because it can be parallelized with mclapply
	actual_strategy <- perform_correlations(strategy = approach, 
											RNAseq = RNAseq, 
											miRNAseq = miRNAseq, 
											corrected_positions = paste0(strategies$all_possible_permutations$RNAseq, "_", strategies$all_possible_permutations$miRNAseq), 
											multimir_summary = NULL, 
											cor_cutoff = corr_cutoff, 
											pval_cutoff = p_val_cutoff)
	print(paste0(approach, " strategy has finished"))
	gc()
	return(actual_strategy)
})
names(selected_strategies) <- approaches

strategies <- c(strategies, selected_strategies)
approaches <- c(approaches, "normalized_counts_RNA_vs_miRNA_normalized_counts")
ref_strategy <- "normalized_counts_RNA_vs_miRNA_normalized_counts"

all_known_miRNAs <- miRBaseConverter::getAllMiRNAs(version = "v22", type = "all", species = organism)$Accession
strategies$all_possible_permutations$known_miRNA <- strategies$all_possible_permutations$miRNAseq %in% all_known_miRNAs

#possible_positives: Pairs with RNA and miRNA included in multiMiR. 
strategies$all_possible_permutations$possible_positives <- strategies$all_possible_permutations$known_miRNA & 
															strategies$all_possible_permutations$miRNAseq %in% multimir_summary$mature_mirna_acc &
															  strategies$all_possible_permutations$RNAseq %in% multimir_summary$target_ensembl 


# q()
all_strategies <- list()
prec_recall <- list()
score_filter <- list()


strategies$all_possible_permutations[is.na(strategies$all_possible_permutations$predicted_c), "predicted_c"] <- 0 
strategies$all_possible_permutations[is.na(strategies$all_possible_permutations$validated_c), "validated_c"] <- 0

for (strategy in approaches) {
	print(strategy)

	
	if (!sum(strategies[[strategy]]$correlated_pairs & strategies$all_possible_permutations$known_miRNA) > 0 ){
		strategies[[strategy]] <- NULL
		message("No significant pairs can be used for multimir comparison in strategy:")
		message(print(strategy))
		next 
	}
	all_pairs_info <- data.table::as.data.table(strategies[[strategy]])

	### Preparing precision, recall & F1 for correlated pairs
	print("Getting prediction stats")
	prec_recall[[strategy]] <- get_prediction_stats(strategy_correlated_pairs = all_pairs_info$correlated_pairs, 
												all_possible_permutations = strategies$all_possible_permutations)
	
	### get Ranking Scores 
	print("Getting ranking scores")

	score_filter[[strategy]] <- get_db_scores(all_db_info = strategies$all_possible_permutations[rownames(all_pairs_info), c("diana_microt", "elmmo", "microcosm", "miranda","mirdb", "pictar", "pita", "targetscan")])

	all_pairs_info$pair <- rownames(all_pairs_info)
	all_pairs_info$score <- strategies$all_possible_permutations[as.numeric(all_pairs_info$pair), "score"]
	all_strategies[[strategy]] <- data.table::as.data.table(all_pairs_info[all_pairs_info$correlated_pairs,])


}



all_strategies <- as.data.frame(data.table::rbindlist(all_strategies, use.names=TRUE, idcol= "strategy"))
prec_recall <- as.data.frame(data.table::rbindlist(prec_recall, use.names = TRUE, idcol = "strategy"))

all_strategies <- cbind(all_strategies, strategies$all_possible_permutations[as.numeric(all_strategies$pair), c("known_miRNA", "predicted_c", "validated_c")])


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

random_obj <- add_randoms(background = strategies$all_possible_permutations, 
								filters_summary = filters_summary,
								db_distribution = db_distribution,
								permutations = permutations,
								all_strategies = all_strategies)

filters_summary <- random_obj[["filters_summary"]]
db_distribution <- random_obj[["db_distribution"]]
 

random_obj <- NULL
filters_summary$type <- factor(filters_summary$type, levels=c("novel_miRNAs","known_miRNAs","predicted", "predicted_random", "validated","validated_random", "both", "both_random"))


################# GET IDS TARNSALTIONS

mirna_names <- unique(strategies$all_possible_permutations$miRNAseq)
mirna_names <- mirna_names[grepl("MIMAT", mirna_names)]
mirna_names <- miRBaseConverter::miRNA_AccessionToName(mirna_names, targetVersion = "v22")
gene_id_translation <- NULL
if (! organism_info$Bioconductor_VarName_SYMBOL[1] == "" && translate_ensembl) {
	gene_id_translation <- get_entrez_symbol_translations(ensembl_ids = unique(strategies$all_possible_permutations$RNAseq), 
														organism_info = organism_info)
} 

save(strategies, all_strategies , file = file.path(output_files, "test.RData"))

rmarkdown::render(file.path(template_folder, 'miRNA_RNA.Rmd'), 
                  output_file = file.path(output_files, report), intermediates_dir = output_files)
}

