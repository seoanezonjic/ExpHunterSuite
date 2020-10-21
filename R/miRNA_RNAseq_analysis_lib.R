miRNA_RNAseq_analysis <- function(
	RNAseq_folder,
	miRNAseq_folder,
	deg_tag,
	output_files,
	aproaches,
	organism,
	multimir_db,
	only_known,
	p_val_cutoff,
	corr_cutoff,
	module_membership_cutoff,
	report,
	translation_file,
	template_folder = file.path(find.package('DEgenesHunter'), "templates"),
	root_path = file.path(find.package('DEgenesHunter')) #esta variable es provisional
	){

time_control <- Sys.time()
aproaches <- unlist(strsplit(aproaches, ","))
print("packages")
print(mem_used())
translation_miRNA <- translation_file
if (!is.null(translation_miRNA)) {
	translation_miRNA <- read.table(translation_miRNA, sep = "\t")
	colnames(translation_miRNA) <- c("mature_ID", "input_ID")
}

RNAseq <- load_DEGH_information(RNAseq_folder)
miRNAseq <- load_DEGH_information(miRNAseq_folder, translation_table =  translation_miRNA)
organism_info <- read.table(file.path(root_path, 'R', "organism_table.txt"), header = TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE, fill = NA) # organism_table.txt hay que cambiarlo de sitio
organism_info <- subset(organism_info, organism_info$KeggCode %in% organism)  




print("files")
print(mem_used())

if (!is.null(multimir_db)){
	multimir_path <- file.path(multimir_db, paste0("parsed_", organism, ".RData"))
	load(multimir_path)
	multimir_summary <- summarize_multimir(multimir_summary)
	message("multiMiR database has been parsed and summarized")
	print(gc())
print(mem_used())
}

strategies <- list()
strategies[["dd"]] <- correlate_profiles(RNAseq[["normalized_counts"]], miRNAseq[["normalized_counts"]])
strategies[["dd"]]$plot_obj <- corr2plot(strategies[["dd"]])
strategies[["dd"]][["name"]] <- "counts_RNA_vs_counts_miRNA" 
strategies[["dd"]]$plot_obj <- add_multimir_info(strategies[["dd"]]$plot_obj, expand_miRNA = NULL , expand_RNA = NULL, multimir_summary = multimir_summary)
print(paste0("dd"," strategy has finished"))

print(gc())
print(mem_used())

if("EE" %in% aproaches){
# Aproach 1: Anticorrelation between RNAseq modules Eigengene and miRNAseq modules Eigengene
# ###### RNAseq modules as rows and miRNA modules as columns
	strategies[["EE"]] <- correlate_profiles(RNAseq[["Eigengene"]], miRNAseq[["Eigengene"]])
	strategies[["EE"]]$plot_obj <- corr2plot(strategies[["EE"]])
	strategies[["EE"]][["name"]] <- "Eigen_RNA_v_Eigen_miRNA" 
	strategies[["EE"]]$plot_obj <- add_multimir_info(strategies[["EE"]]$plot_obj, expand_miRNA = miRNAseq[['DH_results']] , expand_RNA = RNAseq[['DH_results']], multimir_summary = multimir_summary)
	print(paste0("EE"," strategy has finished"))

	print(gc())
	print(mem_used()) 
}

if("Eh" %in% aproaches){
# # Aproach 2: Anticorrelation between RNAseq modules Eigengene and miRNAseq hub genes 
# ###### RNAseq modules as rows and miRNA hub genes as columns
	miRNAseq$hub_1 <- get_hub_genes_by_MM(miRNAseq[["normalized_counts"]], miRNAseq[["DH_results"]], top = 1)
	strategies[["Eh"]] <- correlate_profiles(RNAseq[["Eigengene"]], miRNAseq[["hub_1"]])
	strategies[["Eh"]]$plot_obj <- corr2plot(strategies[["Eh"]])
	strategies[["Eh"]][["name"]] <- "Eigen_RNA_v_Hub1_miRNA" 
	strategies[["Eh"]]$plot_obj <- add_multimir_info(strategies[["Eh"]]$plot_obj, expand_miRNA = miRNAseq[['DH_results']], expand_RNA = RNAseq[['DH_results']], multimir_summary = multimir_summary) 
	print(paste0("Eh"," strategy has finished"))

	print(gc())
	print(mem_used()) 
}

if("Ed" %in% aproaches){
# Aproach 3: Anticorrelation between RNAseq modules Eigengene and miRNAseq expression profiles
###### RNAseq modules as rows and miRNA genes as columns
	strategies[["Ed"]] <- correlate_profiles(RNAseq[["Eigengene"]], miRNAseq[["normalized_counts"]])
	strategies[["Ed"]]$plot_obj <- corr2plot(strategies[["Ed"]])
	strategies[["Ed"]][["name"]] <- "Eigen_RNA_v_counts_miRNA" 
	strategies[["Ed"]]$plot_obj <- add_multimir_info(strategies[["Ed"]]$plot_obj, expand_miRNA = NULL , expand_RNA = RNAseq[['DH_results']], multimir_summary = multimir_summary) 
	print(paste0("Ed"," strategy has finished"))

	print(gc())
	print(mem_used()) 

}

if("hd" %in% aproaches){
# Aproach 4: Anticorrelation between RNAseq hub genes and miRNAseq expression profiles
	RNAseq$hub_1 <- get_hub_genes_by_MM(RNAseq[["normalized_counts"]], RNAseq[["DH_results"]], top = 1)
	strategies[["hd"]] <- correlate_profiles(RNAseq[["hub_1"]], miRNAseq[["normalized_counts"]])
	strategies[["hd"]]$plot_obj <- corr2plot(strategies[["hd"]])
	strategies[["hd"]][["name"]] <- "Hub1_RNA_v_counts_miRNA" 
	strategies[["hd"]]$plot_obj <- add_multimir_info(strategies[["hd"]]$plot_obj, expand_miRNA = NULL , expand_RNA = RNAseq[['DH_results']], multimir_summary = multimir_summary) 
	print(paste0("hd"," strategy has finished"))

	print(gc())
	print(mem_used()) 
}

if("hE" %in% aproaches){
# Aproach 5: Correlation between RNAseq hub genes and miRNAseq modules Eigengene
	RNAseq$hub_1 <- get_hub_genes_by_MM(RNAseq[["normalized_counts"]], RNAseq[["DH_results"]], top = 1)
	strategies[["hE"]] <- correlate_profiles(RNAseq[["hub_1"]], miRNAseq[["Eigengene"]])
	strategies[["hE"]]$plot_obj <- corr2plot(strategies[["hE"]])
	strategies[["hE"]][["name"]] <- "Hub1_RNA_v_Eigen_miRNA" 
	strategies[["hE"]]$plot_obj <- add_multimir_info(strategies[["hE"]]$plot_obj, expand_miRNA = miRNAseq[['DH_results']] , expand_RNA = RNAseq[['DH_results']], multimir_summary = multimir_summary) 
	print(paste0("hE"," strategy has finished"))
	
	print(gc())
	print(mem_used()) 
}


all_known_miRNAs <- miRBaseConverter::getAllMiRNAs(version = "v22", type = "all", species = organism)$Accession
all_strategies <- list()
prec_recall <- list()
print("TEST")
# save.image(file = file.path(output_files, "debug.RData"))
for (strategy in names(strategies)) {
	print(strategy)

	strategy_name <- strategies[[strategy]][["name"]]
	strategies[[strategy]]$plot_obj[is.na(strategies[[strategy]]$plot_obj$predicted_c), "predicted_c"] <- 0
	strategies[[strategy]]$plot_obj[is.na(strategies[[strategy]]$plot_obj$validated_c), "validated_c"] <- 0
	plot_obj <- as.data.table(strategies[[strategy]]$plot_obj)
	print("TEST2")

	significant_pairs <- (plot_obj$correlation <= corr_cutoff & plot_obj$pval <= p_val_cutoff )
	print(paste0("before filtering ", sum(significant_pairs), " pairs"))
	significant_pairs <- significant_pairs & plot_obj$miRNAseq %in% all_known_miRNAs

	
	if (!sum(significant_pairs) > 0 ){
		strategies[[strategy]] <- NULL
		next 
	}
	predicted_pairs <- plot_obj$predicted_c > 0
	validated_pairs <- plot_obj$validated_c > 0
	both_pairs <- predicted_pairs & validated_pairs
	print(paste0("before filtering ", sum(significant_pairs), " pairs"))
	significant_pairs <- significant_pairs & strategies$dd$plot_obj$miRNAseq %in% multimir_summary$mature_mirna_acc & strategies$dd$plot_obj$RNAseq %in% multimir_summary$target_ensembl
	print(paste0("after filtering ", sum(significant_pairs), " pairs"))

	if (sum(significant_pairs) > 0 ){
		 
		stats <- data.frame(stringsAsFactors = FALSE,
							predicted = pred_stats(significant_pairs, predicted_pairs),
							validated = pred_stats(significant_pairs, validated_pairs),
							both = pred_stats(significant_pairs, both_pairs)
							)
	} else {
		message("No significant pairs can be used for multimir comparison in strategy:")
		message(print(strategy))
	}

	print("TEST3")

	rownames(stats) <- c("precision", "recall", "F1")

	stats <- as.data.frame(t(as.matrix(stats)))
	stats$gold_standard <- rownames(stats)

	all_strategies[[strategy_name]] <- as.data.table(plot_obj)[significant_pairs,]
	prec_recall[[strategy_name]] <- stats
	# str(all_strategies[[strategy_name]])
	# str(stats)
	print("TEST4")

	gc()
}

print("TEST")


all_strategies <- as.data.frame(rbindlist(all_strategies, use.names=TRUE, idcol= "strategy"))
prec_recall <- as.data.frame(rbindlist(prec_recall, use.names = TRUE, idcol = "strategy"))



all_strategies$known_miRNA <- all_strategies$miRNAseq %in% all_known_miRNAs 
print("TEST")

filters_summary <- all_strategies %>%
					group_by(strategy) %>%
					summarize(known_miRNAs = sum(known_miRNA == TRUE),
								novel_miRNAs = sum(known_miRNA == FALSE),
								predicted = sum(predicted_c > 0),
								validated = sum(validated_c > 0), 
								both = sum(predicted_c > 0 & validated_c > 0))
filters_summary <- as.data.frame(filters_summary, stringsAsFactors = FALSE)
filters_summary <- reshape2::melt(filters_summary)
names(filters_summary) <- c("strategy", "type", "pairs")
filters_summary$type <- as.character(filters_summary$type)
filters_summary$sdev <- rep(NA, nrow(filters_summary))
filters_summary$p_val <- rep(NA, nrow(filters_summary))
filters_summary$quantile <- rep(NA, nrow(filters_summary))


#generate and merge randoms
# str(filters_summary)
print("TEST")
background_pairs <- strategies$dd$plot_obj[strategies$dd$plot_obj$miRNAseq %in% all_known_miRNAs & strategies$dd$plot_obj$miRNAseq %in% multimir_summary$mature_mirna_acc & strategies$dd$plot_obj$RNAseq %in% multimir_summary$target_ensembl,]
print(paste0("background before: ", nrow(strategies$dd$plot_obj), " background after: ", nrow(background_pairs)))

filters_summary <- add_randoms(background = background_pairs, 
								filters_summary = filters_summary)

filters_summary$type <- factor(filters_summary$type, levels=c("novel_miRNAs","known_miRNAs","predicted", "predicted_random", "validated","validated_random", "both", "both_random"))
print(filters_summary)

################ PREPARE ROC DATA

# roc_data <- data.frame(counts_RNA_vs_counts_miRNA = rep(NA, nrow(strategies$dd$plot_obj)), stringsAsFactors = FALSE)
# rownames(roc_data) <- paste0(strategies$dd$plot_obj$RNAseq, strategies$dd$plot_obj$miRNAseq)
# for (strategy in names(strategies)){
# 	strategy_pairs <- paste0(strategies[[strategy]]$plot_obj$RNAseq, strategies[[strategy]]$plot_obj$miRNAseq)
# 	ordered_correlations <- strategies[[strategy]]$plot_obj[match(rownames(roc_data), strategy_pairs),"correlation"]
# 	scaled_score <- 1 - ordered_correlations
# 	scaled_score[is.na(scaled_score)] <- 0
# 	roc_data[[strategies[[strategy]]$name]] <- scaled_score
# 	# scaled_score
# }

# rownames(roc_data) <- NULL
# roc_data$predicted <- as.numeric(strategies$dd$plot_obj$predicted_c > 0)
# roc_data$validated <- as.numeric(strategies$dd$plot_obj$validated_c > 0)

# save(, roc_data, file = "/mnt/home/users/bio_267_uma/josecordoba/test/test_miRNA-RNA/test.RData")
# q()
################# GET IDS TARNSALTIONS
mirna_names <- unique(all_strategies$miRNAseq)
mirna_names <- mirna_names[grepl("MIMAT", mirna_names)]
mirna_names <- miRNA_AccessionToName(mirna_names, targetVersion = "v22")

gene_id_translation <- NULL
if (! organism_info$Bioconductor_VarName_SYMBOL[1] == "") {

	input_to_entrezgene <- id_translation_orgdb(input_ids = unique(all_strategies$RNAseq),
												 organism_db = organism_info$Bioconductor_DB[1],
												 org_var_name = organism_info$Bioconductor_VarName[1])

	colnames(input_to_entrezgene) <- c("ensembl_gene_id", "entrezgene") # Fix names


	input_to_symbol <- id_translation_orgdb(input_ids = input_to_entrezgene$entrezgene,
										 		organism_db = organism_info$Bioconductor_DB[1],
												org_var_name = organism_info$Bioconductor_VarName_SYMBOL[1])
	colnames(input_to_symbol) <- c("entrezgene", "Symbol")
	input_to_entrezgene <- as.data.table(input_to_entrezgene)
	input_to_symbol <- as.data.table(input_to_symbol)
	gene_id_translation <- merge(input_to_entrezgene, input_to_symbol, by = "entrezgene", all.x = TRUE)

}

rmarkdown::render(file.path(template_folder, 'miRNA_RNA.Rmd'), 
                  output_file = file.path(output_files, report), intermediates_dir = output_files)

}
