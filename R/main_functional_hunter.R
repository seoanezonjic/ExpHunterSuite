functional_hunter <- function(
	input_hunter_folder,
	model_organism,
	annot_file,
	organisms_table = file.path(find.package('DEgenesHunter'), "R", "organism_table.txt"),
	template_folder = file.path(find.package('DEgenesHunter'), "templates"),
	List_organisms = FALSE,
	input_gene_id = "E",
	func_annot_db = "gKR",
	GO_subont = "BMC",
	custom = NULL,
	analysis_type = c("go"),
	remote = "",
	save_query = FALSE,
	pthreshold = 0.1,
	qthreshold = 0.2,
	debug_file = NULL,
	cores = 1,
	output_files = "results",
	fc_colname = "mean_logFCs"
	){

	if(!file.exists(input_hunter_folder)) 
		stop("No input degenes_Hunter folder")

	############ CREATE DEBUG FOLDERS #########
	if(!is.null(debug_file)) {
		debug <- TRUE
		debug_dir <- dirname(debug_file)
		dir.create(debug_dir, recursive = T)
		debug_dir <- normalizePath(debug_dir)
		# Store session
	####################### DEBUG POINT #############################
		time_control <- list(start = Sys.time())
		debug_point(debug_file, "Start point", environment())
	#################################################################
	}



	## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
	##                                                                                                                   ##
	##                                       LOAD AND CHECK MAIN FILES                                                   ## 
	##                                                                                                                   ##
	## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

	## Load available organisms
	## all_organisms_info => Table with information abaut all organism available
	all_organisms_info <- read.table(organisms_table, header = TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE, fill = NA)

	# Check organism selected
	if (List_organisms){
	  print(as.character(rownames(all_organisms_info)))
	  stop('Check this list and choose one model species.')
	} else if (any(is.null(model_organism), !model_organism %in% rownames(all_organisms_info))) {
	  stop('Model organism has not been indicated or is not available. Please indicate the model organism using parameter -m. Use -L to list all model organisms available')
	} else { # ORGANISM AVAILABLE --> LOAD
	  current_organism_info <- subset(all_organisms_info, rownames(all_organisms_info) %in% model_organism)  
	}

	# Load DEGenesHunter config
	# DEGenesHunter_expression_opt => DEGenesHunter options 
	DEGenesHunter_expression_opt <- read.table(file.path(input_hunter_folder, "opt_input_values.txt"), header = FALSE, stringsAsFactors = FALSE, sep = "\t")
	degh_exp_threshold <- as.numeric(DEGenesHunter_expression_opt[which(DEGenesHunter_expression_opt[,1] == "p_val_cutoff"), 2])

	experiments <- read.table(file.path(input_hunter_folder, "control_treatment.txt"), sep = "\t", quote = "", header = TRUE, stringsAsFactors = FALSE)
	experiments$class[experiments$class == "C"] <- "Control"
	experiments$class[experiments$class == "T"] <- "Treatment"
	sample_classes <- apply(experiments, 1, function(x) paste0("* [", x[1], "] ", x[2]))

	## Load DEGenesHunter results
	## DEGH_results => DEgenes Hunter results table and modifications
	DEGH_results <- read.table(file.path(input_hunter_folder, "Common_results", "hunter_results_table.txt"), header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE)
	DEGH_results <- DEGH_results[DEGH_results$genes_tag != "FILTERED_OUT", ]

	# Temporary bodge to enable funcitonal hunter to be run on externally processed data and DESeq2 not run
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
	if(!is.null(annot_file)){
		annot_table <- read.table(annot_file, header=FALSE, row.names=NULL, sep="\t", stringsAsFactors = FALSE, quote = "") # Load
	  	DEGH_results$input_IDs <- translate_id(rownames(DEGH_results), annot_table)
	} else {
		DEGH_results$input_IDs <- rownames(DEGH_results)
	}
	added_cols <- c("input_IDs")

	# Load Normalized correlation data
	if (flags$WGCNA) {
		####
		# LOAD NORMALIZED COUNTS
		norm_counts <- as.matrix(read.table(file.path(input_hunter_folder, "Results_DESeq2", "Normalized_counts_DESeq2.txt"), header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE))
		scaled_counts <- scale_data_matrix(data_matrix = norm_counts, transpose = TRUE)
		scaled_counts_table <- as.data.frame(as.table(scaled_counts))
		colnames(scaled_counts_table) <- c("Gene","Sample","Count")
			
		####
		# LOAD WGCNA clusters representative profiles with samples
		cl_eigvalues <- as.matrix(read.table(file.path(input_hunter_folder, "Results_WGCNA", "eigen_values_per_samples.txt"), header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE))
		cl_eigvalues <- as.data.frame(as.table(cl_eigvalues),stringsAsFactors = FALSE)
		colnames(cl_eigvalues) <- c("Sample","Cluster_ID","Count") 
		cl_eigvalues_gnorm <- cl_eigvalues
		cl_eigvalues_gnorm$Count <- (cl_eigvalues_gnorm$Count + 1) / 2 
		
		####
		# LOAD WGCNA - PVal (Cluster - Trait)
		wgcna_pval_cl_trait <- as.matrix(read.table(file.path(input_hunter_folder, "Results_WGCNA", "module_trait_p_val.txt"), header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE))
		wgcna_corr_cl_trait <- as.matrix(read.table(file.path(input_hunter_folder, "Results_WGCNA", "module_trait.txt"), header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE))
		
		####
		# LOAD WGCNA - Correlation (Sample - Trait)
		wgcna_count_sample_trait <- as.matrix(read.table(file.path(input_hunter_folder, "Results_WGCNA", "sample_trait.txt"), header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE))
		wgcna_count_sample_trait <- scale_data_matrix(wgcna_count_sample_trait)
	}


	###############################
	## LOAD CUSTOM GMT 	
	all_custom_gmt <- NULL
	if (!is.null(custom)) {
		custom_files <- unlist(strsplit(custom, ","))
		all_custom_gmt <- lapply(custom_files, load_and_parse_gmt)
		names(all_custom_gmt) <- custom_files
	}

	# Verbose point 
	aux <- table(DEGH_results$genes_tag)
	for(i in seq_along(aux)) {
		message(paste(names(aux)[i],aux[i]))
	}


	####################### DEBUG POINT #############################
	if (debug) {
		time_control$load_data <- Sys.time()
		debug_point(debug_file, "Files have been loaded and parsed", environment())
	}
	message("Files have been loaded and parsed")
	#################################################################



	## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
	##                                                                                                                   ##
	##                                            Translate IDs                                                          ## 
	##                                                                                                                   ##
	## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
	if (need_translate) {
		# valid_genes => A vector containing all unique genes that were not filtered out. These genes are translated according annot_file given with -a
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
				# colnames(input_to_entrezgene) <- c("ensembl_gene_id", current_organism_info[,"Attribute_entrez"])
				colnames(input_to_entrezgene) <- c("ensembl_gene_id", "entrezgene") # Fix names
			}
		} else {
		# print("it works")
			stop("Specified organism not available in LOCAL mode. Please try REMOTE mode")
		}

		if (save_query) { ## TODO => ESTO HAY QUE DARLE UNA FUNCIONALIDAD O QUITARLO
		  saveRDS(input_to_entrezgene, file=file.path("query_results_temp"))
		} 
		
		################# ADD ENTREZ IDS TO INPUT FILE #############
		DEGH_results$entrezgene <- input_to_entrezgene[match(DEGH_results$input_IDs, input_to_entrezgene$ensembl_gene_id), "entrezgene"]
		added_cols <- c(added_cols, "entrezgene")
		message(paste(opt$input_gene_id, "IDs have been translated to entrezgene"))
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

	####################### DEBUG POINT #############################
	if (debug) {
		time_control$translate_id <- Sys.time()
		debug_point(debug_file, "IDs translated", environment())
	}
	#################################################################


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
			
				gene_fc <- as.vector(cl_results$fc_colname)
				names(gene_fc) <- cl_results$entrezgene
				gene_fc <- sort(gene_fc, decreasing = TRUE)

				return(gene_fc)
			})
		}
	}

	####################### DEBUG POINT #############################
	if (debug) {
		time_control$data_preparation <- Sys.time()
		debug_point(debug_file, "Data prepared for enrichments", environment())
	} 
	message("Data prepared for enrichments")
	#################################################################



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
				perform_topGO(go_attr_name, likely_degs_entrez, union_annot_DEGs_df, mod, paste("GOgraph_preval_",mod,".pdf",sep=""),input_gene_id)
				perform_topGO(go_attr_name, pos_logFC_likely_degs, union_annot_DEGs_df, mod, paste("GOgraph_preval_overex_",mod,".pdf",sep=""),input_gene_id)
				perform_topGO(go_attr_name, neg_logFC_likely_degs, union_annot_DEGs_df, mod, paste("GOgraph_preval_underex_",mod,".pdf",sep=""),input_gene_id)
				# Union
				perform_topGO(go_attr_name, union_DEGs, input_to_entrezgene, mod, paste("GOgraph_allpos_",mod,".pdf",sep=""),input_gene_id)
				perform_topGO(go_attr_name, pos_logFC_union_DEGs, input_to_entrezgene, mod, paste("GOgraph_allpos_overex_",mod,".pdf",sep=""),input_gene_id)
				perform_topGO(go_attr_name, neg_logFC_union_DEGs, input_to_entrezgene, mod, paste("GOgraph_allpos_underex_",mod,".pdf",sep=""),input_gene_id)    				
			}))
		} else { # LOCAL MODE

			reference_ids_union <- unique(DEGH_results$entrezgene)
			reference_ids_common <- unique(DEGH_results$entrezgene)

			# Launch topGO analysis
			invisible(lapply(GO_subontologies,function(mod) {
				# Common
				perform_topGO_local(likely_degs_entrez, reference_ids_common, mod, file.path(output_files, paste("GO_preval",mod,sep="_")),current_organism_info$Bioconductor_DB[1])
				perform_topGO_local(pos_logFC_likely_degs, reference_ids_common, mod, file.path(output_files, paste("GO_preval_overex",mod,sep="_")),current_organism_info$Bioconductor_DB[1])
				perform_topGO_local(neg_logFC_likely_degs, reference_ids_common,mod, file.path(output_files, paste("GO_preval_underex",mod,sep="_")),current_organism_info$Bioconductor_DB[1])
				# # Union
				perform_topGO_local(union_DEGs, reference_ids_union, mod, file.path(output_files, paste("GO_allpos",mod,sep="_")),current_organism_info$Bioconductor_DB[1])
				perform_topGO_local(pos_logFC_union_DEGs, reference_ids_union, mod, file.path(output_files, paste("GO_allpos_overex",mod,sep="_")),current_organism_info$Bioconductor_DB[1])
				perform_topGO_local(neg_logFC_union_DEGs, reference_ids_union, mod, file.path(output_files, paste("GO_allpos_underex",mod,sep="_")),current_organism_info$Bioconductor_DB[1])    
			}))
		} # END LOCAL/REMOTE IF
		####################### DEBUG POINT #############################
		if (debug) {
			time_control$TopGO <- Sys.time()
			debug_point(debug_file, "topGO analysis finished", environment())
		}
		message("topGO analysis finished")
		#################################################################
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
			# Write results
			write.table(as.data.frame(do.call(rbind,lapply(enrich_go,function(res) {as.data.frame(res)}))), file=file.path(output_files, "GO_CL_ora"), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")
			
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
				enrich <- enrichment_GSEA(geneList = geneList,organism = current_organism_info$Bioconductor_DB[1],keyType = keytypes,pvalueCutoff = pthreshold,pAdjustMethod = "BH",ont = paste0("GO_",mod))
				return(enrich)
			})
			# Add names
			names(enrich_go_gsea) <- GO_subontologies
			# Write results
			write.table(as.data.frame(do.call(rbind,lapply(enrich_go_gsea,function(res) {as.data.frame(res)}))), file=file.path(output_files, "GO_CL_gsea"), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")	

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
		####################### DEBUG POINT #############################
		if (debug) {
			time_control$clusterProfiler_GO <- Sys.time()
			debug_point(debug_file, "clusterProfiler GO analysis finished", environment())
		}
		message("clusterProfiler GO analysis finished")
		#################################################################
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
			# Write output
			write.table(enrich_ora, file=file.path(output_files, "KEGG_results"), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")
			
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
			enrich_gsea <- enrichment_GSEA(geneList = geneList,organism = current_organism_info$KeggCode[1],pvalueCutoff = pthreshold,ont = "KEGG",useInternal = !remote_actions$kegg)
			# Write output
			write.table(enrich_gsea, file=file.path(output_files, "KEGG_GSEA_results"), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")
			
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
		####################### DEBUG POINT #############################
		if (debug) {
			time_control$clusterProfiler_KEGG <- Sys.time()
			debug_point(debug_file, "clusterProfiler KEGG analysis finished", environment())
		}
		message("clusterProfiler KEGG analysis finished")
		#################################################################
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
			# Write output
			write.table(enrich_react, file=file.path(output_files, "REACT_results"), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")	
			
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
			enrich_react_gsea <- enrichment_GSEA(geneList = geneList, organism = current_organism_info$Reactome_ID[1], pvalueCutoff = pthreshold, pAdjustMethod = "BH", ont = "REACT")
			# Write output
			write.table(enrich_react_gsea, file=file.path(output_files, "REACT_GSEA_results"), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")
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
		####################### DEBUG POINT #############################
		if (debug) {
			time_control$clusterProfiler_REACTOME <- Sys.time()
			debug_point(debug_file, "clusterProfiler REACTOME analysis finished", environment())
		}
		message("clusterProfiler REACTOME analysis finished")
		#################################################################
		#PERFORM WGCNA MODULES ANALYSIS
	}

	#############################################
	### CUSTOM ENRICHMENT
	if (!is.null(custom)) {
		# Obtain custom files
		custom_files <- unlist(strsplit(custom,","))

		###########
		### CUSTOM ENRICHMENTS
		# Per each file, launch enrichment
		custom_enrichments <- enrich_all_customs(custom_files = custom_files, 
												 p_val_threshold = pthreshold, 
												 likely_degs_entrez = likely_degs_entrez)
		if (flags$WGCNA){
			custom_cls_ORA_expanded <- lapply(all_custom_gmt, function(gmt){
				enrich_clusters_with_gmt(gmt, clgenes)
			})
		}
		####################### DEBUG POINT #############################
		if (debug) {
			time_control$TopGO <- Sys.time()
			debug_point(debug_file, "CUSTOM enrichments finished", environment())
		}
		message("CUSTOM enrichments finished")
		#################################################################
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
			enrichments_ORA <- lapply(enrichments_ORA_expanded, merge_result)
			# Write output
			for(enrichment_i in 1:length(enrichments_ORA)) {
				# Concat
				df <- enrichplot:::fortify.compareClusterResult(enrichments_ORA[[enrichment_i]])
				# Write table
				write.table(df, file=file.path(output_files, paste0(names(enrichments_ORA[enrichment_i]),"_cls_ora")), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")
			}
		}

		if (flags$GSEA) {
			enrichments_GSEA <- lapply(enrichments_GSEA_expanded, merge_result)

			for (enrichment_i in 1:length(enrichments_GSEA)) {
				# Concat
				# df <- clusterProfiler:::fortify.compareClusterResult(enrichments_GSEA[[enrichment_i]])
				df <- enrichments_GSEA[[enrichment_i]]@compareClusterResult
				# Write table
				write.table(df, file=file.path(output_files, paste0(names(enrichments_GSEA[enrichment_i]),"_cls_gsea")), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")
			}
		}

		if (!is.null(opt$custom)){
			custom_cls_ORA <- lapply(custom_cls_ORA_expanded, merge_result)
			for(enrichment_i in 1:length(custom_cls_ORA)) {
				# Concat
				df <- enrichplot:::fortify.compareClusterResult(custom_cls_ORA[[enrichment_i]])
				# Write table
				write.table(df, file=file.path(output_files, paste0(basename(names(custom_cls_ORA[enrichment_i])),"_cls_ORA")), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")
			}
		}
	}

	DEGH_results <- DEGH_results[c(added_cols, setdiff(names(DEGH_results), added_cols))] # Reorder columns so annotated first
	DEGH_results <- DEGH_results[order(DEGH_results[,"combined_FDR"]),] # Reorder rows by combined FDR
	write.table(DEGH_results, file=file.path(output_files, "hunter_results_table_annotated.txt"), quote=FALSE, col.names=NA, sep="\t")

	####################### DEBUG POINT #############################
	if (debug) {
		time_control$TopGO <- Sys.time()
		debug_point(debug_file, "Results saved", environment())
	}
	message("Results have been saved")
	#################################################################

	## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
	##                                                                                                                   ##
	##                                                 RENDERING REPORTS                                                 ##                                                     
	##                                                                                                                   ##
	## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

	message("RENDERING REPORT ...")

	############################################################
	##                GENERATE CLUSTER REPORTS                ##
	############################################################
	results_path <- normalizePath(output_files)

	if (flags$WGCNA) { # Clustered
		message("Rendering specific cluster reports")
		invisible(mclapply(cls, function(cl) {
			# Take output name
			aux <- paste0("cl_func_",cl,".html")
			outf_cls_i <- file.path(results_path, aux)
			# Generate report
			rmarkdown::render(file.path(template_folder, 'cl_func_report.Rmd'), output_file = outf_cls_i, intermediates_dir = results_path)
			if (debug) {
				time_control[[paste0("cl_func_",cl)]] <<- Sys.time()
			}
		}, mc.cores = cores
	))

		message("\tRendering clustered report")
		outf_cls <- file.path(results_path, "clusters_func_report.html")
		rmarkdown::render(file.path(template_folder, 'clusters_main_report.Rmd'),output_file = outf_cls, intermediates_dir = results_path)
		if (debug) {
			time_control$render_cluster_main_report <- Sys.time()
		}
	}

	############################################################
	##              GENERATE DEG FUNCTIONAL REPORT            ##
	############################################################
	message("\tRendering regular report")
	outf <- file.path(results_path, "functional_report.html")
	rmarkdown::render(file.path(template_folder, 'functional_report.Rmd'), output_file = outf, intermediates_dir = results_path)

	if (debug) {
	####################### DEBUG POINT #############################
	    time_control$render_main_report <- Sys.time()
	    debug_point(debug_file,"Report printed", environment())
	#################################################################
	    save_times(time_control, output = file.path(debug_dir, "function_hunter_time_control.txt"), plot_name = "FH_times_control.pdf")
	}
}
