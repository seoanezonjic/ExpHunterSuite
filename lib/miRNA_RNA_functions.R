################### FUNCTIONS

load_DEGH_information <- function(execution_path, translation_table = NULL, only_known = FALSE){
	execution <- list()
	execution[['normalized_counts']] <- as.matrix(read.table(file.path(execution_path, "Results_DESeq2/Normalized_counts_DESeq2.txt"), header=TRUE, row.names=1, sep="\t"))

	execution[['DH_results']] <- read.table(file.path(execution_path, "Common_results/hunter_results_table.txt"), header=TRUE, row.names=1, sep="\t")
	execution[['DH_results']]$gene_name <- rownames(execution[['DH_results']])
	rownames(execution[['DH_results']]) <- NULL
	execution[["DH_results"]] <-  execution[["DH_results"]][execution[["DH_results"]]$gene_name != "id",] #Esto es un control porque en algun momento del flujo se ha insertado una columna en la matriz de conteos que no debe
	
	execution[['Eigengene']] <- read.table(file.path(execution_path, "Results_WGCNA/eigen_values_per_samples.txt"), header=TRUE, row.names=1, sep="\t")
	execution[['Eigengene']] <- execution[['Eigengene']] #TODO => Proponer el cambio de output del DGH
	names(execution[["Eigengene"]]) <- str_remove(names(execution[["Eigengene"]]), "ME")
	execution[["Eigengene"]] <- as.matrix(execution[["Eigengene"]][,names(execution[["Eigengene"]]) != 0])
 
	if (!is.null(translation_table)) {
	# 	#translate DH_results
		temp_translation <- as.vector(translation_table[match(execution[['DH_results']]$gene_name, translation_table$input_ID), "mature_ID"])
		if (!only_known) {
			for (element in seq(1:length(temp_translation))){ 
				if(is.na(temp_translation[element])) {
					temp_translation[element] <- execution[["DH_results"]][element, "gene_name"]
				}
			}
		}
		execution[["DH_results"]]$gene_name <- temp_translation
		execution[["DH_results"]] <- execution[["DH_results"]][!is.na(execution[["DH_results"]]$gene_name), ]

		#translate normalized counts
		temp_translation <- as.vector(translation_table[match(rownames(execution[['normalized_counts']]), translation_table$input_ID), "mature_ID"])
		if (!only_known) {
			for (element in seq(1:length(temp_translation))){ 
				if(is.na(temp_translation[element])) {
					temp_translation[element] <- rownames(execution[['normalized_counts']])[element]
				}
			}
		}

		rownames(execution[['normalized_counts']]) <- temp_translation
		execution[['normalized_counts']] <- execution[['normalized_counts']][!is.na(rownames(execution[['normalized_counts']])), ]
		
	} 	
	execution[['normalized_counts']] <- t(execution[['normalized_counts']])
	return(execution)
}



# filter_correlations
get_hub_genes_by_MM <- function(normalized_counts, hunter_results, top = 1){
	hub_genes <- hunter_results %>% 
				filter(Cluster_MM_pval <= 0.05) %>% 
				arrange(desc(Cluster_MM)) %>% 
				group_by(Cluster_ID) %>% filter(between(row_number(), 1, top))

	hub_genes_profile <- as.matrix(normalized_counts[,colnames(normalized_counts) %in% hub_genes$gene_name])
	
	colnames(hub_genes_profile) <- as.data.frame(hub_genes)[match(colnames(hub_genes_profile), hub_genes$gene_name), "Cluster_ID"]
	return(hub_genes_profile)
}

tag_strategy <- function(plot_obj, tag){
	plot_obj$strategy <- rep("zero", nrow(plot_obj))
	return(plot_obj)
}

translate_column <- function(column, translation_table) { #translation_table, a dataframe with translated_ids in first column and ids_to_translate in second column
	column <- translation_table[match(column, translation_table[,2]), 1]
	return(column)
}

correlate_profiles <- function(profiles_A, profiles_B) {
	
	correlations <- cor(profiles_A, profiles_B)
	
	nSamples <- ncol(profiles_A) 
	cor_pval <- as.data.frame(corPvalueStudent(as.matrix(correlations), nSamples))

	correlated_profiles <- list(
		corr = correlations,
		pval = cor_pval
	)
	return(correlated_profiles)
}


corr2plot <- function(correlated_profiles) {
	names_matrix <- as.matrix(correlated_profiles$corr)  
	plot_df <- as.data.frame(as.table(correlated_profiles$corr)) 
	colnames(plot_df) <- c("RNAseq", "miRNAseq", "correlation")
	
	plot_df$correlation <- round(plot_df$correlation, 2)

	plot_df$pval <- as.data.frame(as.table(as.matrix(correlated_profiles$pval)))$Freq

	plot_df$miRNAseq <- as.character(plot_df$miRNAseq)
	plot_df$RNAseq <- as.character(plot_df$RNAseq)
	return(plot_df)
}

summarize_multimir <- function(multimir_table) {
	multimir_table <- as.data.frame(multimir_table)

	dbs <- list(
		validated_c = c("mirecords", "mirtarbase", "tarbase"),
		predicted_c = c("diana_microt", "elmmo", "microcosm", "miranda","mirdb", "pictar", "pita", "targetscan")

		)
	for (db_type in names(dbs)) {
		multimir_table[[db_type]] <- rowSums(multimir_table[,dbs[[db_type]]])
	}
	return(multimir_table)
}



add_multimir_info <- function(plot_obj, expand_miRNA = NULL , expand_RNA = NULL, multimir_summary = NULL) {
	plot_obj <- as.data.table(plot_obj)
	if (!is.null(expand_RNA)) {
		print("expand RNA")
		names(plot_obj)[names(plot_obj)=="RNAseq"] <- "RNAseq_mod"

		RNAseq <- expand_RNA  %>% select(Cluster_ID, gene_name)
		RNAseq <- as.data.table(RNAseq)
		colnames(RNAseq) <- c("module", "RNAseq")
		print(head(RNAseq))
		RNAseq$module <- as.character(RNAseq$module)

		plot_obj <- merge(x = plot_obj, y = RNAseq, by.x = "RNAseq_mod", by.y = "module", by =.EACHI, allow.cartesian  = TRUE)
		plot_obj$RNAseq_mod <- NULL
		
	}
 #HASTA AQUI BIEN

	if (!is.null(expand_miRNA)) {
		# print("test2")
		print("expand RNA")
		names(plot_obj)[names(plot_obj) == "miRNAseq"] <- "miRNAseq_mod"
		
		miRNAseq <- expand_miRNA %>% select(Cluster_ID, gene_name) # get_module_genes(expand_miRNA, miRNAseq, "miRNAseq")
		miRNAseq <- as.data.table(miRNAseq)
		colnames(miRNAseq) <- c("module", "miRNAseq")
		miRNAseq$module <- as.character(miRNAseq$module)

		plot_obj <- merge(x = plot_obj, y = miRNAseq, by.x = "miRNAseq_mod", by.y = "module", by =.EACHI, allow.cartesian  = TRUE)
		plot_obj$miRNAseq_mod <- NULL
	}
	
	if (!is.null(multimir_summary)){
		multimir_summary <- as.data.table(multimir_summary)
		plot_obj <- merge(plot_obj, multimir_summary, by.x = c("RNAseq","miRNAseq"), by.y =c("target_ensembl","mature_mirna_acc"), all.x = TRUE, suffixes = c("", ""), no.dups = TRUE)
	}
	gc()
	return(as.data.frame(plot_obj))
}

expand_module <- function(module_genes, plot_obj, column_name){ #column_name is the name of column to expand
	gene_corr_list <- lapply(1:nrow(plot_obj), function(i) {

	 	 	this_row <- plot_obj[i,]
	 	 	module_genes2 <- module_genes[module_genes$module %in% this_row[[column_name]], column_name]
	 	 	this_row[[column_name]] <- NULL
	  		expanded_data <- data.table(this_row, module_genes2)
	  		names(expanded_data)[names(expanded_data) == "module_genes2"] <- column_name
	  		return(as.data.frame(expanded_data))
		})
	gene_corr_list <- rbindlist(gene_corr_list)

	return(gene_corr_list)
}


get_module_genes <- function(hunter_results, module_name, column_name){
	module_genes <- hunter_results[["gene_name"]][hunter_results[["Cluster_ID"]] == module_name]

	module_genes <- data.frame(module = rep(module_name, length(module_genes)), genes = module_genes)
	colnames(module_genes) <- c("module", column_name)
	return(module_genes)
}