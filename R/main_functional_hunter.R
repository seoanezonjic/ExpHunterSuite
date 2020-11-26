#' This function allows you to perform the Functional analysis with different enrichment corpus.
#' @param hunter_results list with DEG analysis results
#' @param model_organism organisms which genes are being studied
#' @param annot_table (OPTIONAL) annotation table to translate gene IDs
#' @param organisms_table (OPTIONAL) configuration table for given organism. Use see get_organism_table()
#' @param input_gene_id (OPTIONAL) type og gene IDs given. Allowed: ENSEMBL (E), entrezgene (e), TAIR/Arabidopsis (T), Gene Names (G).
#' @param func_annot_db (OPTIONAL) modules to be enriched
#' @param GO_subont (OPTIONAL) specific GO subontologies to be enriched
#' @param custom (OPTIONAL) list of dataframes with GMT corpus
#' @param analysis_type (OPTIONAL) enrichment type to be performed. Allowed: ORA (o) and GSEA (g)
#' @param remote (OPTIONAL) remote allowed actions. Are: Biomart query (b) and KEGG enrichment (k)
#' @param save_query
#' @param pthreshold pvalue threshold
#' @param qthreshold qvalue threshold
#' @param cores cores to be used if parallel features are going to be used. Default: 1
#' @param output_files output folder
#' @param fc_colname main logFC colname (into hunter_results dataframe)
#' @keywords 
#' @export
#' @examples
functional_hunter <- function(
	hunter_results,
	model_organism,
	annot_table = NULL,
	organisms_table = get_organism_table(),
	input_gene_id = "E",
	func_annot_db = "gKR",
	GO_subont = "BMC",
	custom = NULL,
	analysis_type = "o", #g
	remote = "",
	save_query = FALSE,
	pthreshold = 0.1,
	qthreshold = 0.2,
	cores = 1,
	output_files = "results",
	fc_colname = "mean_logFCs"
	){

	func_results <- list()
	main_params <- list(model_organism = model_organism,
						input_gene_id = input_gene_id,
						func_annot_db = func_annot_db,
						GO_subont = GO_subont,
						analysis_type = analysis_type,
						remote = remote,
						save_query = save_query,
						pthreshold = pthreshold,
						qthreshold = qthreshold,
						cores = cores,
						fc_colname = fc_colname)
	func_results[["final_main_params"]] <- main_params

	# Check organism selected
	if (any(is.null(model_organism), !model_organism %in% rownames(organisms_table))) {
		stop('Model organism has not been indicated or is not available. Please indicate the model organism.')
	}else{ # ORGANISM AVAILABLE --> LOAD
		current_organism_info <- subset(organisms_table, rownames(organisms_table) %in% model_organism)  
	}

	# experiments <- read.table(file.path(input_hunter_folder, "control_treatment.txt"), sep = "\t", quote = "", header = TRUE, stringsAsFactors = FALSE)
	experiments <- hunter_results$sample_groups
	experiments$class[experiments$class == "C"] <- "Control"
	experiments$class[experiments$class == "T"] <- "Treatment"
	sample_classes <- apply(experiments, 1, function(x) paste0("* [", x[1], "] ", x[2]))

	## Load DEGenesHunter results
	## DEGH_results => DEgenes Hunter results table and modifications
	DEGH_results <- hunter_results$DE_all_genes
	DEGH_results <- DEGH_results[DEGH_results$genes_tag != "FILTERED_OUT", ]

	# Temporary bodge to enable funcitonal hunter to be run on externally processed data and DESeq2 not run
	# TODO : This is not working
	if(! grepl("DESeq2", names(DEGH_results)) && grepl("external_DEA", names(DEGH_results))) {
		names(DEGH_results) <- gsub("external_DEA", "DESeq2", names(DEGH_results))
		# This is particularly horrible - like this to ensure it matches EXACTLY the file loaded if WGCNA
		external_DEA_folder <- file.path(input_hunter_folder, "Results_external_DEA")
		DESeq2_dummy_folder <- file.path(input_hunter_folder, "Results_DESeq2")
		dir.create(DESeq2_dummy_folder)
		file.copy(file.path(external_DEA_folder, "Normalized_counts_external_DEA.txt"),
				  file.path(DESeq2_dummy_folder, "Normalized_counts_DESeq2.txt"))
	}

	## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
	##                                                                                                                   ##
	##                                         CONFIGURE ENRICHMENTS                                                     ## 
	##                                                                                                                   ##
	## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

	# Check remote actions
	remote_actions <- list(biomart = grepl("b", remote),
	              		   kegg    = grepl("k", remote))

	# Check which enrichments are gonna be performed
	# flags => A list of enricgment analysis options
	flags <- list(GO    = grepl("G", func_annot_db),
	              KEGG  = grepl("K", func_annot_db),
	              GO_cp = grepl("g", func_annot_db),
	              REACT = grepl("R", func_annot_db),
	              GSEA  = grepl("g", analysis_type),
	              ORA   = grepl("o", analysis_type),
	              WGCNA = "Cluster_ID" %in% colnames(DEGH_results))


	need_translate <- TRUE
	# Prepare and check given IDs 
	if(input_gene_id == "e"){
		input_gene_id <- 'entrezgene'
		need_translate <- FALSE
		keytypes <- "ENTREZID"
		input_to_entrezgene <- cbind(row.names(DEGH_results), row.names(DEGH_results))
		colnames(input_to_entrezgene) <- c("ensembl_gene_id", "entrezgene")
	}else if(input_gene_id == "E"){
	  input_gene_id <- 'ensembl_gene_id'
	  keytypes <- "ENTREZID"
	}else if(input_gene_id == "T"){
	  input_gene_id <- ''
	  keytypes <- "TAIR"
	}else if(input_gene_id == "G"){
	  input_gene_id <- ''
	  keytypes <- "GENENAME"
	}else{
	  stop(paste0("Given ID type (", input_gene_id, ") is still not allowed. Please check -t option."))
	}


	if(flags$GO_cp && current_organism_info$Bioconductor_VarName[1] == ""){
		flags$GO_cp <- FALSE
		message("Specified organism is not allowed to be used with GO (clusterProfiler) module. Please check your IDs table")
	}

	if(flags$KEGG && current_organism_info$KeggCode[1] == ""){
		flags$KEGG <- FALSE
		warning("Specified organism is not allowed to be used with KEGG module. Please check your IDs table")
	}

	if(flags$REACT){
		if (keytypes %in% c("TAIR","GENENAME")) {
			flags$REACT <- FALSE
			warning("Reactome module can not be used with GENENAME identifiers")
		} else if (current_organism_info$Reactome_ID[1] == "" || (keytypes != "ENTREZID")) {
			flags$REACT <- FALSE
			warning("Specified organism is not allowed to be used with Reactome module. Please check your IDs table")
		} else {
			require(ReactomePA)
		}
	}

	if(!any(flags$GP_cp, flags$KEGG, flags$REACT)){
		if(!is.null(custom)){
			flags$ORA <- FALSE
			flags$GSEA <- FALSE 
		}
	}

	# Check which subontologies are gonna be performed
	# GO_subontologies => vector with subontologies to be perfromed
	if(any(flags$GO, flags$GO_cp)){
		GO_subontologies <- c()
		if(grepl("M", GO_subont)){
			GO_subontologies <- c(GO_subontologies,"MF")
		}
		if(grepl("B", GO_subont)){
			GO_subontologies <- c(GO_subontologies,"BP")
		}
		if(grepl("C", GO_subont)){
			GO_subontologies <- c(GO_subontologies,"CC")
		}
	}


	## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
	##                                                                                                                   ##
	##                                         LOAD COMPLEMENTARY FILES                                                  ## 
	##                                                                                                                   ##
	## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
	#Load annotation gene table for translation ID
	if(!is.null(annot_table)){
	  	DEGH_results$input_IDs <- translate_id(rownames(DEGH_results), annot_table)
	} else {
		DEGH_results$input_IDs <- rownames(DEGH_results)
	}
	added_cols <- c("input_IDs")

	# Verbose point 
	aux <- table(DEGH_results$genes_tag)
	for(i in seq_along(aux)) {
		message(paste(names(aux)[i],aux[i]))
	}


	## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
	##                                                                                                                   ##
	##                                            Translate IDs                                                          ## 
	##                                                                                                                   ##
	## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
	if (need_translate) {
		# valid_genes => A vector containing all unique genes that were not filtered out. These genes are translated according annot_table given with -a
		valid_genes <- unique(DEGH_results$input_IDs[!is.na(DEGH_results$input_IDs)])
		if (remote_actions$biomart) { # REMOTE MODE
			############# BIOMART QUERY  #############

			# Prepare organism info
			organism_mart <- as.character(current_organism_info[,"Mart"])
			organism_dataset <- as.character(current_organism_info[,"Dataset"])
			organism_host <- as.character(current_organism_info[,"biomaRt_Host"])
			organism_attr <- c(
			  input_gene_id, 
			  as.character(current_organism_info[,"Attribute_GOs"]),
			  as.character(current_organism_info[,"Attribute_entrez"])
			)

			# Obtain references from biomart
			input_to_entrezgene <- obtain_info_from_biomaRt(orthologues = as.character(valid_genes),
			                                            id_type = input_gene_id, 
			                                            mart = organism_mart,
			                                            dataset = organism_dataset, 
			                                            host = organism_host, 
			                                            attr = organism_attr) 

		} else if (current_organism_info$Bioconductor_DB[1] != "") { # LOCAL MODE
			# Check genetical items to be used
			if (input_gene_id == 'ensembl_gene_id') {
				input_to_entrezgene <- id_translation_orgdb(input_ids = valid_genes,
													 organism_db = current_organism_info$Bioconductor_DB[1],
													 org_var_name = current_organism_info$Bioconductor_VarName[1])
				colnames(input_to_entrezgene) <- c("ensembl_gene_id", "entrezgene") # Fix names
			}
		} else {
			stop("Specified organism not available in LOCAL mode. Please try REMOTE mode")
		}

		if (save_query) { ## TODO => ESTO HAY QUE DARLE UNA FUNCIONALIDAD O QUITARLO
		  saveRDS(input_to_entrezgene, file=file.path("query_results_temp"))
		} 
		
		################# ADD ENTREZ IDS TO INPUT FILE #############
		DEGH_results$entrezgene <- input_to_entrezgene[match(DEGH_results$input_IDs, input_to_entrezgene$ensembl_gene_id), "entrezgene"]
		added_cols <- c(added_cols, "entrezgene")
		message(paste(input_gene_id, "IDs have been translated to entrezgene"))
	} else {
		################# ADD ENTREZ IDS TO INPUT FILE #############
		DEGH_results$entrezgene <- DEGH_results$input_IDs
	}

	################# ADD GENE SYMBOLS TO INPUT FILE #############
	# Currently only runs if we are local, have an org db available and a symbol correspondence. The ENTREZ bit should become universal even if no SYMBOL available
	if ( ! remote_actions$biomart &&
		! current_organism_info$Bioconductor_VarName_SYMBOL[1] == "") {

		input_to_symbol <- id_translation_orgdb(input_ids = unique(DEGH_results$entrezgene[!is.na(DEGH_results$entrezgene)]),
											 		organism_db = current_organism_info$Bioconductor_DB[1],
													org_var_name = current_organism_info$Bioconductor_VarName_SYMBOL[1])
		colnames(input_to_symbol) <- c("entrezgene", "Symbol")
		
		DEGH_results$Symbol <- input_to_symbol[match(DEGH_results$entrezgene, input_to_symbol$entrezgene), "Symbol"]
		added_cols <- c(added_cols, "Symbol")
		message("entrezgene IDs have been translated to SYMBOLS")
	}



	## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
	##                                                                                                                   ##
	##                                     PREPARE DATA FOR ENRICHMENTS                                                  ## 
	##                                                                                                                   ##
	## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

	# Obtain prevalent items
	# Obtain significant sbusets
	#  > likely_degs_df : ONLY FOR SUBSETTING: subset of DEGenes Hunter annotation file with PREVALENT_DEG flag
	#  > likely_degs : USED FOR topGO. genes identifiers included into likely_degs_df
	#  > likely_degs_entrez : USED FOR ORA ENRICHMENTS:  list of entrez gene present into likely_degs_df AND input_to_entrezgene with filtered count data
	#  > union_DEGs_df : subset of DEGenes Hunter annotation file with POSSIBLE_DEG flag or PREVALENT_DEG flag
	#  > union_DEGs : genes identifiers included into union_DEGs_df
	#  > union_annot_DEGs_df : reference table subset with identifiers included into union_DEGs

	###############################
	## topGO & ORA
	likely_degs_df <- subset(DEGH_results, genes_tag == "PREVALENT_DEG", !is.na(input_IDs))
	likely_degs_entrez <- unique(likely_degs_df$entrezgene)
	# likely_degs <- unique(likely_degs_df$input_IDs)

	# Verbose point
	message(paste("IDs used in ORA:",length(likely_degs_entrez)))
	## TODO => ESTARIA BIEN REFLEJAR ESTA INFORMACION EN EL REPORT
	"%>%" <- magrittr::"%>%"
	union_DEGs_df <- subset(DEGH_results, genes_tag %in% c("POSSIBLE_DEG", "PREVALENT_DEG"))
	union_DEGs <- union_DEGs_df[!is.na(union_DEGs_df$input_IDs), "input_IDs"] %>% unique
	union_annot_DEGs_df <- subset(input_to_entrezgene, input_to_entrezgene[,1] %in% union_DEGs)

	###############################
	## GSEA
	# Calculate geneList
	#geneList: A NAMED VECTOR OF logFC; USED IN GSEA AND PLOTS
	geneList <- DEGH_results[!is.na(DEGH_results$entrezgene),  fc_colname]
	names(geneList) <- DEGH_results[!is.na(DEGH_results$entrezgene), "entrezgene"]
	# Sort FC
	geneList <- sort(geneList, decreasing = TRUE)

	#############################################
	### WGCNA MODULES
	if (flags$WGCNA) {
		cls <- unique(DEGH_results$Cluster_ID)
		# DELETE GREY MODULE
		if (any(c(0,"grey") %in% cls)) {
			cls <- cls[!cls %in% c(0,"grey")]
		} else {
			warning("Module Zero/Grey not found")
		}

		clgenes <- lapply(cls,function(cl) { # Find
			unique(DEGH_results$entrezgene[which(DEGH_results$Cluster_ID == cl)])
		}) 
		names(clgenes) <- cls

		if (flags$ORA) enrichments_ORA_expanded <- list()

		if (flags$GSEA) {
			enrichments_GSEA_expanded <- list()
			cl_gene_fc <- lapply(cls, function(cl) {
				cl_results <- DEGH_results[!is.na(DEGH_results$entrezgene) & DEGH_results$Cluster_ID == cl,]			
				gene_fc <- as.vector(cl_results[,fc_colname])
				names(gene_fc) <- cl_results$entrezgene
				gene_fc <- sort(gene_fc, decreasing = TRUE)

				return(gene_fc)
			})
			names(cl_gene_fc) <- as.character(cls)
		}
	}

	message("Data prepared for enrichments")



	## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
	##                                                                                                                   ##
	##                                           PERFORM ENRICHMENTS                                                     ## 
	##                                                                                                                   ##
	## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

	if (flags$GO) {
		pos_logFC_likely_degs <- subset(likely_degs_df, likely_degs_df[fc_colname] > 0, !is.na(entrezgene) )$entrezgene
		neg_logFC_likely_degs <- subset(likely_degs_df, likely_degs_df[fc_colname] < 0, !is.na(entrezgene))$entrezgene
		pos_logFC_union_DEGs <- subset(union_DEGs_df, union_DEGs_df[fc_colname] > 0, !is.na(entrezgene))$entrezgene
		neg_logFC_union_DEGs <- subset(union_DEGs_df, union_DEGs_df[fc_colname] < 0, !is.na(entrezgene))$entrezgene
		# Generate output
		# Check execution mode
		if (remote_actions$biomart) { # REMOTE MODE
			# Prepare necessary info
			go_attr_name <- as.character(current_organism_info[,"Attribute_GOs"])
			# Launch topGO analysis
			invisible(lapply(GO_subontologies,function(mod) {
				# Common
# TODO : this methods create files (should be exported out of this main)
				perform_topGO(go_attr_name, likely_degs_entrez, union_annot_DEGs_df, mod, paste("GOgraph_preval_",mod,".pdf",sep=""),input_gene_id, output_files)
				perform_topGO(go_attr_name, pos_logFC_likely_degs, union_annot_DEGs_df, mod, paste("GOgraph_preval_overex_",mod,".pdf",sep=""),input_gene_id, output_files)
				perform_topGO(go_attr_name, neg_logFC_likely_degs, union_annot_DEGs_df, mod, paste("GOgraph_preval_underex_",mod,".pdf",sep=""),input_gene_id, output_files)
				# Union
				perform_topGO(go_attr_name, union_DEGs, input_to_entrezgene, mod, paste("GOgraph_allpos_",mod,".pdf",sep=""),input_gene_id, output_files)
				perform_topGO(go_attr_name, pos_logFC_union_DEGs, input_to_entrezgene, mod, paste("GOgraph_allpos_overex_",mod,".pdf",sep=""),input_gene_id, output_files)
				perform_topGO(go_attr_name, neg_logFC_union_DEGs, input_to_entrezgene, mod, paste("GOgraph_allpos_underex_",mod,".pdf",sep=""),input_gene_id, output_files)    				
			}))
		} else { # LOCAL MODE

			reference_ids_union <- unique(DEGH_results$entrezgene)
			reference_ids_common <- unique(DEGH_results$entrezgene)

			# Launch topGO analysis
			invisible(lapply(GO_subontologies,function(mod) {
				# Common
# TODO : this methods create files (should be exported out of this main)
				perform_topGO_local(likely_degs_entrez, reference_ids_common, mod, file.path(output_files, paste("GO_preval",mod,sep="_")),current_organism_info$Bioconductor_DB[1])
				perform_topGO_local(pos_logFC_likely_degs, reference_ids_common, mod, file.path(output_files, paste("GO_preval_overex",mod,sep="_")),current_organism_info$Bioconductor_DB[1])
				perform_topGO_local(neg_logFC_likely_degs, reference_ids_common,mod, file.path(output_files, paste("GO_preval_underex",mod,sep="_")),current_organism_info$Bioconductor_DB[1])
				# # Union
				perform_topGO_local(union_DEGs, reference_ids_union, mod, file.path(output_files, paste("GO_allpos",mod,sep="_")),current_organism_info$Bioconductor_DB[1])
				perform_topGO_local(pos_logFC_union_DEGs, reference_ids_union, mod, file.path(output_files, paste("GO_allpos_overex",mod,sep="_")),current_organism_info$Bioconductor_DB[1])
				perform_topGO_local(neg_logFC_union_DEGs, reference_ids_union, mod, file.path(output_files, paste("GO_allpos_underex",mod,sep="_")),current_organism_info$Bioconductor_DB[1])    
			}))
		} # END LOCAL/REMOTE IF
		message("topGO analysis finished")
	}



	#############################################
	### GO ENRICHMENT (clusterProfiler)
	if (flags$GO_cp) {
		message("Performing GO enrichments")
		###########
		### ORA ENRICHMENTS
		if (flags$ORA) {
			enrich_go <- lapply(GO_subontologies,function(mod) {
				enrich <- enrichment_ORA(genes = likely_degs_entrez, organism = current_organism_info$Bioconductor_DB[1],keyType = keytypes,pvalueCutoff = pthreshold,pAdjustMethod = "BH", ont = paste0("GO_",mod), qvalueCutoff = qthreshold)
				return(enrich)
			})
			# Add names
			names(enrich_go) <- GO_subontologies
			func_results$GO_ORA <- enrich_go
			#PERFORM WGCNA MODULES ANALYSIS
			if (flags$WGCNA) {
				for (subont in GO_subontologies) {
					subont_name <- paste0("GO_", subont)
					enrichments_ORA_expanded[[subont_name]] <- enrichment_clusters_ORA(genes = clgenes,
		                                 organism = current_organism_info$Bioconductor_DB[1],
		                                 keyType = keytypes,
		                                 pvalueCutoff = pthreshold,
		                                 pAdjustMethod = "BH",
		                                 ont = subont_name,
		                                 qvalueCutoff = qthreshold,
		                                 useInternal = FALSE,
		                                 mc.cores = cores)
				}
			}
		}

		###########
		### GSEA ENRICHMENTS
		if (flags$GSEA) {
			# Enrich
			enrich_go_gsea <- lapply(GO_subontologies,function(mod) {
				enrich <- enrich_GSEA(geneList = geneList,organism = current_organism_info$Bioconductor_DB[1],keyType = keytypes,pvalueCutoff = pthreshold,pAdjustMethod = "BH",ont = paste0("GO_",mod))
				return(enrich)
			})
			# Add names
			names(enrich_go_gsea) <- GO_subontologies
			func_results$GO_GSEA <- enrich_go_gsea
			#PERFORM WGCNA MODULES ANALYSIS
			if (flags$WGCNA) {
				for (subont in GO_subontologies) {
					subont_name <- paste0("GO_", subont)
					enrichments_GSEA_expanded[[subont_name]] <- perform_GSEA_clusters(all_clusters = cl_gene_fc,
										organism = current_organism_info$Bioconductor_DB[1],
										keyType = keytypes,
										pvalueCutoff = pthreshold,
										pAdjustMethod = "BH",
										ont = subont_name, 
										useInternal = FALSE)
				}
			}
		}
		message("clusterProfiler GO analysis finished")
	}


	#############################################
	### KEGG ENRICHMENT
	if (flags$KEGG) {
		message("Performing KEGG enrichments")
		###########
		### ORA ENRICHMENTS
		if (flags$ORA) {
			# Enrich
			enrich_ora <- enrichment_ORA(genes = likely_degs_entrez,organism = current_organism_info$KeggCode[1],keyType = "kegg",pvalueCutoff = pthreshold,pAdjustMethod = "BH",ont = "KEGG",useInternal = !remote_actions$kegg, qvalueCutoff = qthreshold)
			func_results$KEGG_ORA <- enrich_ora
			#PERFORM WGCNA MODULES ANALYSIS
			if (flags$WGCNA) {
				enrichments_ORA_expanded[["KEGG"]] <- enrichment_clusters_ORA(genes = clgenes,
	                                 organism = current_organism_info$KeggCode[1],
	                                 keyType = "kegg",
	                                 pvalueCutoff = pthreshold,
	                                 pAdjustMethod = "BH",
	                                 ont = "KEGG",
	                                 qvalueCutoff = qthreshold,
	                                 useInternal = !remote_actions$kegg,
	                                 mc.cores = cores)
			}
		}

		###########
		### GSEA ENRICHMENTS
		if (flags$GSEA) {
			enrich_gsea <- enrich_GSEA(geneList = geneList,organism = current_organism_info$KeggCode[1],pvalueCutoff = pthreshold,ont = "KEGG",useInternal = !remote_actions$kegg)
			func_results$KEGG_GSEA <- enrich_gsea
			#PERFORM WGCNA MODULES ANALYSIS
			if (flags$WGCNA) {
				enrichments_GSEA_expanded[["KEGG"]] <- perform_GSEA_clusters(all_clusters = cl_gene_fc,
									organism = current_organism_info$KeggCode[1],
									keyType = "ENTREZID",
									pvalueCutoff = pthreshold,
									pAdjustMethod = "BH",
									ont = "KEGG", 
									useInternal = !remote_actions$kegg)
			}
		}
		message("clusterProfiler KEGG analysis finished")
	}


	#############################################
	### REACTOME ENRICHMENT

	if (flags$REACT) {
		message("Performing Reactome enrichments")
		
		###########
		### ORA ENRICHMENTS
		if (flags$ORA) {
			# Make enrichment (ORA)
			enrich_react <- enrichment_ORA(genes = likely_degs_entrez,organism = current_organism_info$Reactome_ID[1],keyType = "ENTREZID",pvalueCutoff = pthreshold,pAdjustMethod = "BH",ont = "REACT", qvalueCutoff = qthreshold)		
			func_results$REACT_ORA <- enrich_react
			#PERFORM WGCNA MODULES ANALYSIS
			if (flags$WGCNA) {
				enrichments_ORA_expanded[["REACT"]] <- enrichment_clusters_ORA(genes = clgenes,
	                                 organism = current_organism_info$Reactome_ID[1],
	                                 keyType = "ENTREZID",
	                                 pvalueCutoff = pthreshold,
	                                 pAdjustMethod = "BH",
	                                 ont = "REACT",
	                                 qvalueCutoff = qthreshold,
	                                 useInternal = FALSE,
	                                 mc.cores = cores)
			}
		}

		###########
		### GSEA ENRICHMENTS
		if (flags$GSEA) {
			enrich_react_gsea <- enrich_GSEA(geneList = geneList, organism = current_organism_info$Reactome_ID[1], pvalueCutoff = pthreshold, pAdjustMethod = "BH", ont = "REACT")
			func_results$REACT_GSEA <- enrich_react_gsea

			if (flags$WGCNA) {
				enrichments_GSEA_expanded[["REACT"]] <- perform_GSEA_clusters(all_clusters = cl_gene_fc,
									organism = current_organism_info$Reactome_ID[1],
									keyType = "ENTREZID",
									pvalueCutoff = pthreshold,
									pAdjustMethod = "BH",
									ont = "REACT", 
									useInternal = FALSE)
			}
		}
		message("clusterProfiler REACTOME analysis finished")
	}

	#############################################
	### CUSTOM ENRICHMENT
	if (!is.null(custom)) {
		message("Performing CUSTOM enrichments")
		###########
		### CUSTOM ENRICHMENTS
		custom_enrichments <- enrich_all_customs(custom_sets = custom, 
												 p_val_threshold = pthreshold, 
												 genes = likely_degs_entrez)
		func_results$CUSTOM <- custom_enrichments
		if (flags$WGCNA){
			custom_cls_ORA_expanded <- lapply(custom, function(gmt){
				enrich_clusters_with_gmt(custom_set = gmt, 
										genes_in_modules = clgenes, 
										p_val_threshold = pthreshold, 
										cores = cores)
			})
		}
		message("CUSTOM enrichments finished")
	} else {
		custom_enrichments <- NULL
	}

	## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
	##                                                                                                                   ##
	##                                               EXPORT DATA                                                         ## 
	##                                                                                                                   ##
	## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
	if (flags$WGCNA) {
		if (flags$ORA) {
			enrichments_ORA <- lapply(enrichments_ORA_expanded, clusterProfiler::merge_result)
			enrichments_ORA <- lapply(enrichments_ORA, function(res){
				if(nrow(res) > 0) res <- catched_pairwise_termsim(res)
				return(res)
			})
			func_results$WGCNA_ORA <- enrichments_ORA
			func_results$WGCNA_ORA_expanded <- enrichments_ORA_expanded
		}

		if (flags$GSEA) {
			enrichments_GSEA <- lapply(enrichments_GSEA_expanded, clusterProfiler::merge_result)
			func_results$WGCNA_GSEA <- enrichments_GSEA
			func_results$WGCNA_GSEA_expanded <- enrichments_GSEA_expanded
		}

		if (!is.null(custom)){
			custom_cls_ORA <- lapply(custom_cls_ORA_expanded, clusterProfiler::merge_result)
			custom_cls_ORA <- lapply(custom_cls_ORA, function(res){
				if(nrow(res) > 0) res <- catched_pairwise_termsim(res)
				return(res)
			})
			func_results$WGCNA_CUSTOM <- custom_cls_ORA
			func_results$WGCNA_CUSTOM_expanded <-custom_cls_ORA_expanded
		}
	}

	DEGH_results <- DEGH_results[c(added_cols, setdiff(names(DEGH_results), added_cols))] # Reorder columns so annotated first
	DEGH_results <- DEGH_results[order(DEGH_results[,"combined_FDR"]),] # Reorder rows by combined FDR
	func_results$DEGH_results_annot <- DEGH_results


	func_results$flags <- flags
	return(func_results)
}


#' Table with information abaut all organism available
#' @param file to be loaded. Default: internal organism table
#' @return organism table
#' @keywords 
#' @export
get_organism_table <- function(file = file.path(find.package('DEgenesHunter'), "external_data", "organism_table.txt")){
	return(read.table(file, header = TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE, fill = NA))
}
