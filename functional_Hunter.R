#! /usr/bin/env Rscript
#############################################
############## FUNCTIONAL HUNTER ###########
#############################################

# this is wrapped in a tryCatch. The first expression works when source executes, the
# second expression works when R CMD does it.
full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile),  # works when using source
               error=function(e) # works when using R CMD
              normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', commandArgs())], '='))[2]))
main_path_script <- dirname(full.fpath)

#############################################
### OPTIONS
#############################################
#Loading optparse  
suppressPackageStartupMessages(require(optparse))
# Parse command line
#------------------------------------------------

option_list <- list(
  make_option(c("-i", "--input_hunter_folder"), type="character",
    help="DEgenes Hunter differential expression analysis output folder"), 
  make_option(c("-m", "--model_organism"), type="character",
    help="Species to use for functional enrichment. You can see all available species running functional_Hunter.R -L"),
  make_option(c("-L", "--List_organisms"), action="store_true", type="logical", default=FALSE, 
    help="Print all organisms available and ends the program"),
  make_option(c("-a", "--annot_file"), type="character",
  	help="If the species used does not exist in organism_table.txt, please add a two-column file mapping orthologue ENSEMBL gene ids from a model organism (first column) to the original IDs (second column)."),
  make_option(c("-t", "--input_gene_id"), type="character", default="E",
    help="Input gene ID type. ENSEMBL (E), TAIR/Arabidopsis (T), Gene Names (G). [Default:%default]"),      	
  make_option(c("-f", "--func_annot_db"), type="character", default="gKR",
    help="Functional annotation database and enrichment method(s) to use (topGO: G = GO | clusterProfiler: K = KEGG, g = GO, R = Reactome). [Default=%default]"),
  make_option(c("-G", "--GO_subont"), type="character", default=c("BMC"),
    help="GO sub-ontologies to use for functional analysis (M = Molecular Function, B = Biological Process, C = Celular Component). Default=%default"), # Not Checked
  make_option(c("-C", "--custom"), ,type = "character", default=NULL,
    help="Files with custom functional annotation database (in GMT format) separated by commas (,)"),
  make_option(c("-A", "--analysis_type"), type="character", default=c("go"),
    help="Analysis performance (g = Gene Set Enrichment Analysis, o = Over Representation Analysis). Default=%default"), # Not Checked
  # make_option(c("-K", "--Kegg_organism"), type="character", default=NULL, 
  #   help="Indicate organism to look for in the Kegg database for doing the path enrichment"), # Not Checked
  make_option(c("-r", "--remote"), ,type = "character", default="",
    help="Flags to activate remote query from enrichments and Genes translation. Use (b) to launch biomaRt translation; (k) to use Kegg remote data base"),
  make_option(c("-q", "--save_query"), type="logical", action = "store_true", default=FALSE,
    help="Flag to save biomaRt query."), 
  make_option(c("-P", "--pthreshold"), type="double", default=0.1,
    help="Enrichment p-value threshold. [Default = %default]"),
  make_option(c("-Q", "--qthreshold"), type="double", default=0.2,
    help="Enrichment q-value threshold. [Default = %default]"),
  make_option(c("--debug"), type="logical", default=FALSE, action = "store_true",
    help="Activate debug mode, which stores RData sessions at different points of the pipeline"),
  make_option(c("--Debug"), type="character", default=NULL,
    help="Activate debug mode and uses given filename. File must have '.RData' extension"),
  make_option(c("-c", "--cores"), ,type = "numeric", default=1,
    help="Cores to be used to parallelize clusters enrichments. Default : %default"),
  make_option(c("-o", "--output_files"), type="character", default="results",
    help="Output path. Default=%default")
)
opt <- parse_args(OptionParser(option_list=option_list))

#############################################
### LOAD LIBRARIES
#############################################
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(biomaRt)) 
suppressPackageStartupMessages(require(topGO))
suppressPackageStartupMessages(require(KEGGREST))
suppressPackageStartupMessages(require(stringr))
suppressPackageStartupMessages(require(plyr))
suppressPackageStartupMessages(require(knitr))
suppressPackageStartupMessages(require(clusterProfiler))
suppressPackageStartupMessages(require(dplyr))
# source custom functions
source(file.path(main_path_script, 'lib', 'general_functions.R'))
source(file.path(main_path_script, 'lib', 'functional_analysis_library.R'))
source(file.path(main_path_script, 'lib', 'plotting_functions.R'))

#############################################
### MAIN 
#############################################
# Load Reference-NonRefernce models gene IDs relations
if(! file.exists(opt$input_hunter_folder)) stop("No input degenes_Hunter folder")

# Special IDs
fc_colname <- "mean_logFCs"

# Check remote actions
remote_actions <- list(biomart = grepl("b", opt$remote),
              		   kegg    = grepl("k", opt$remote))

############ CREATE FOLDERS #########3
paths <- list()
dir.create(opt$output_files)
paths$root <-opt$output_files

if(!is.null(opt$Debug)){
	opt$debug <- TRUE
	debug_file <- file.path(opt$Debug)
}

if(opt$debug){
	# Define only once
	if(is.null(opt$Debug)){
		debug_file <- file.path(paths$root, debug_files, paste(c("FHunter_Debug_Session_",format(Sys.Date(),format = "%Y%m%d"),".RData"),collapse = ""))
	}
	debug_dir <- dirname(debug_file)
	dir.create(debug_dir, recursive = T)
	debug_dir <- normalizePath(debug_dir)
	# Store session
####################### DEBUG POINT #############################
	time_control <- list(start = Sys.time())
	debug_point(debug_file,"Start point")
#################################################################
}

#############################################
### LOAD AND PARSE 
#############################################

# Load available organisms
all_organisms_info <- read.table(file.path(main_path_script, "lib", "organism_table.txt"), header = TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE, fill = NA)

# Load DEGH_results
#DEGH_results => DEgenes Hunter results table and modifications
DEGH_results <- read.table(file.path(opt$input_hunter_folder, "Common_results", "hunter_results_table.txt"), header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE)
aux <- which(DEGH_results$genes_tag == "FILTERED_OUT") 
if(length(aux) > 0){
	DEGH_results <- DEGH_results[-aux,]
}
#TODO => ESTO ES RARO; BUSCAR UNA MANERA MEJOR DE HACERLO

####
# LOAD DEGenesHunter expression execution config
DEGenesHunter_expression_opt <- read.table(file.path(opt$input_hunter_folder, "opt_input_values.txt"), header = FALSE, stringsAsFactors = FALSE, sep = "\t")
degh_exp_threshold <- as.numeric(DEGenesHunter_expression_opt[which(DEGenesHunter_expression_opt[,1] == "p_val_cutoff"),2])

experiments <- read.table(file.path(opt$input_hunter_folder, "control_treatment.txt"), sep = "\t", quote = "", header = TRUE, stringsAsFactors = FALSE)
exp_names <- paste("[Control]",experiments[which(experiments[,1] == "C"),2],sep=" ")
exp_names <- c(exp_names,paste("[Treatment]",experiments[which(experiments[,1] == "T"),2],sep=" "))


if(!is.null(opt$annot_file)){
	annot_table <- read.table(opt$annot_file, header=FALSE, row.names=NULL, sep="\t", stringsAsFactors = FALSE, quote = "") # Load
  	DEGH_results$Annot_IDs <- translate_id(rownames(DEGH_results), annot_table)
}else{
	DEGH_results$Annot_IDs <- rownames(DEGH_results)
}
valid_genes <- unique(DEGH_results$Annot_IDs[!is.na(DEGH_results$Annot_IDs)])

# Prepare ID type
if(opt$input_gene_id == "E"){
  opt$input_gene_id <- 'ensembl_gene_id'
  keytypes <- "ENTREZID"
# }else if(opt$input_gene_id == "R"){
#   opt$input_gene_id <- 'refseq_peptide'
#   keytypes <- "ENTREZID"
}else if(opt$input_gene_id == "T"){
  opt$input_gene_id <- ''
  keytypes <- "TAIR"
}else if(opt$input_gene_id == "G"){
  opt$input_gene_id <- ''
  keytypes <- "GENENAME"
}else{
  stop(paste("Given ID type (",opt$input_gene_id,") is still not allowed. Please check -t option.",sep=""))
}

# Check organism selected
if(opt$List_organisms == TRUE){
  print(as.character(rownames(all_organisms_info)))
  stop('Check this list and choose one model species.')

}else if(is.null(opt$model_organism)){
  stop('No model organism has been indicated indicated. Please indicate the model organism using parameter -m. Use -L to list all model organisms available')
# }else if(!opt$organisms %in% rownames(all_organisms_info)){
#   stop("Organism selected is not available. Use -L to display available organims list")
}else{ # ORGANISM AVAILABLE --> LOAD
  current_organism_info <- subset(all_organisms_info, rownames(all_organisms_info) %in% opt$model_organism)  
}


# Check which enrichments are gonna be performed
flags <- list(GO    = grepl("G", opt$func_annot_db),
              KEGG  = grepl("K", opt$func_annot_db),
              GO_cp = grepl("g", opt$func_annot_db),
              REACT = grepl("R", opt$func_annot_db),
              GSEA  = grepl("g", opt$analysis_type),
              ORA   = grepl("o", opt$analysis_type),
              Clustered = "Cluster_ID" %in% colnames(DEGH_results))

# TODO =>  REMOVE
# Special case
# if(exists("annot_table")){
# 	DEGH_results$Annot_IDs <- translate_id(rownames(DEGH_results), annot_table)
# } else {
# 	DEGH_results$Annot_IDs <- rownames(DEGH_results)
# }

# Verbose point
aux <- table(DEGH_results$genes_tag)
for(i in seq_along(aux)) {
	message(paste(names(aux)[i],aux[i]))
}


####################### DEBUG POINT #############################
if(opt$debug){
	time_control$load_data <- Sys.time()
	debug_point(debug_file,"Data loaded")
}
#################################################################

#############################################
### PREPARE AND TRANSFORM DATA
#############################################

if(remote_actions$biomart){ # REMOTE MODE
	############# BIOMART QUERY  #############

	# Prepare organism info
	organism_mart <- as.character(current_organism_info[,"Mart"])
	organism_dataset <- as.character(current_organism_info[,"Dataset"])
	organism_host <- as.character(current_organism_info[,"biomaRt_Host"])
	organism_attr <- c(
	  opt$input_gene_id, 
	  as.character(current_organism_info[,"Attribute_GOs"]),
	  as.character(current_organism_info[,"Attribute_entrez"])
	)

	# Obtain references from biomart
	reference_table <- obtain_info_from_biomaRt(orthologues = as.character(valid_genes),
	                                            id_type = opt$input_gene_id, 
	                                            mart = organism_mart,
	                                            dataset = organism_dataset, 
	                                            host = organism_host, 
	                                            attr = organism_attr) 

}else if(!is.na(current_organism_info$Bioconductor_DB[1])){ # LOCAL MODE
	# Check genetical items to be used
	if(opt$input_gene_id == 'ensembl_gene_id'){
		reference_table <- ensembl_to_entrez(ensembl_ids = valid_genes,
											 organism_db = current_organism_info$Bioconductor_DB[1],
											 organism_var = current_organism_info$Bioconductor_VarName[1])
		# colnames(reference_table) <- c("ensembl_gene_id", current_organism_info[,"Attribute_entrez"])
		colnames(reference_table) <- c("ensembl_gene_id", "entrezgene") # Fix names

	}else if(opt$input_gene_id == 'refseq_peptide'){
		stop(paste("This genes type (",opt$input_gene_id,") is not allowed in local mode yet",sep="")) ##################################### NOT IMPLEMENTED
	}
}else{
	error("Specified organism are not available to be studiend in LOCAL model. Please try REMOTE mode")
}



if(opt$save_query == TRUE){ ## TODO => CACHE IS SAVED BUT NEVER IS LOADED
  saveRDS(reference_table, file=file.path("query_results_temp"))
} 


# Load Normalized correlation data
if(flags$Clustered){
	####
	# LOAD NORMALIZED COUNTS
	norm_counts_raw <- as.matrix(read.table(file.path(opt$input_hunter_folder, "Results_DESeq2", "Normalized_counts_DESeq2.txt"), header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE))
	# Translate genes
	if(exists("reference_table")){
		aux <- unlist(lapply(rownames(norm_counts_raw),function(id){reference_table$entrezgene[which(reference_table$ensembl_gene_id == id)][1]}))
		aux_indx <- which(is.na(aux))
		if(length(aux_indx)>0) aux[aux_indx] <- rownames(norm_counts_raw)[aux_indx]
		rownames(norm_counts_raw) <- aux
	}
	# Normalize
	norm_counts_raw_gnorm <- norm_counts_raw
	invisible(lapply(seq_along(norm_counts_raw[,1]),function(i){
		m <- min(norm_counts_raw[i,])
		M <- max(norm_counts_raw[i,])
		dff <- M - m
		norm_counts_raw_gnorm[i,] <<- (norm_counts_raw[i,] - m) / dff
	}))
	# Modify to plot better later
	norm_counts <- as.data.frame(as.table(norm_counts_raw))
	colnames(norm_counts) <- c("Gene","Sample","Count")
	norm_counts_gnorm <- as.data.frame(as.table(norm_counts_raw_gnorm))
	colnames(norm_counts_gnorm) <- c("Gene","Sample","Count")
	# norm_counts_gnorm <- cbind(norm_counts_gnorm,list(Type = rep("Regular",nrow(norm_counts_gnorm))))

	####
	# LOAD WGCNA clusters representative profiles with samples
	cl_eigvalues <- as.matrix(read.table(file.path(opt$input_hunter_folder, "Results_WGCNA", "eigen_values_per_samples.txt"), header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE))
	cl_eigvalues <- as.data.frame(as.table(cl_eigvalues),stringsAsFactors = FALSE)
	colnames(cl_eigvalues) <- c("Sample","Cluster_ID","Count") 
	cl_eigvalues_gnorm <- cl_eigvalues
	cl_eigvalues_gnorm$Count <- (cl_eigvalues_gnorm$Count + 1) / 2 
	
	####
	# LOAD WGCNA - PVal (Cluster - Trait)
	wgcna_pval_cl_trait <- as.matrix(read.table(file.path(opt$input_hunter_folder, "Results_WGCNA", "module_trait_p_val.txt"), header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE))
	wgcna_corr_cl_trait <- as.matrix(read.table(file.path(opt$input_hunter_folder, "Results_WGCNA", "module_trait.txt"), header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE))
	####
	# LOAD WGCNA - Correlation (Sample - Trait)
	wgcna_count_sample_trait <- as.matrix(read.table(file.path(opt$input_hunter_folder, "Results_WGCNA", "sample_trait.txt"), header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE))
	wgcna_count_sample_trait_gnorm <- as.data.frame(do.call(cbind,lapply(seq(ncol(wgcna_count_sample_trait)),function(j){
		(wgcna_count_sample_trait[,j] - min(wgcna_count_sample_trait[,j],na.rm = TRUE)) / (max(wgcna_count_sample_trait[,j],na.rm = TRUE) - min(wgcna_count_sample_trait[,j],na.rm = TRUE))
	})))
	colnames(wgcna_count_sample_trait_gnorm) <- colnames(wgcna_count_sample_trait)
}

####################### DEBUG POINT #############################
if(opt$debug){
	time_control$data_extended <- Sys.time()
	debug_point(debug_file,"Data extended")
}
#################################################################

################# PREPROCESSING INPUT DATA FOR FUNCTIONAL ANALYSIS #############
# Obtain prevalent items
# Obtain significant sbusets
#  > likely_degs_df : ONLY FOR SUBSETTING: subset of DEGenes Hunter annotation file with PREVALENT_DEG flag
#  > likely_degs : USED FOR topGO. genes identifiers included into likely_degs_df
#  > likely_degs_entrez : USED FOR ORA ENRICHMENTS:  list of entrez gene present into likely_degs_df AND reference_table with filtered count data
#  > union_DEGs_df : subset of DEGenes Hunter annotation file with POSSIBLE_DEG flag or PREVALENT_DEG flag
#  > union_DEGs : genes identifiers included into union_DEGs_df
#  > union_annot_DEGs_df : reference table subset with identifiers included into union_DEGs

likely_degs_df <- subset(DEGH_results, genes_tag == "PREVALENT_DEG", !is.na(Annot_IDs))
likely_degs <- unique(likely_degs_df$Annot_IDs)
likely_degs_entrez <- reference_table[reference_table$ensembl_gene_id %in% likely_degs, "entrezgene"] %>% unique

# Verbose point
message(paste("IDs used to enrich:",length(likely_degs_entrez)))
## TODO => ESTARIA BIEN REFLEJAR ESTA INFORMACION EN EL REPORT

union_DEGs_df <- subset(DEGH_results, genes_tag %in% c("POSSIBLE_DEG", "PREVALENT_DEG"))
union_DEGs <- union_DEGs_df[!is.na(union_DEGs_df$Annot_IDs), "Annot_IDs"] %>% unique
union_annot_DEGs_df <- subset(reference_table, reference_table[,1] %in% union_DEGs)

################# ADD ENTREZ IDS AND GENE SYMBOLS TO INPUT FILE #############
# Currently only runs if we are local, have an org db available and a symbol correspondence. The ENTREZ bit should become universal even if no SYMBOL available

if( ! remote_actions$biomart &
	! current_organism_info$Bioconductor_DB[1] == "" & 
	! current_organism_info$Bioconductor_VarName_SYMBOL[1] == ""){

	DEGH_results_Symbol <- DEGH_results
	if(exists("annot_table")){
		ENSEMBL_IDs <- DEGH_results_Symbol$Annot_IDs
	}else{
		ENSEMBL_IDs <- rownames(DEGH_results_Symbol)
	}

	DEGH_results_Symbol$Ensembl <- ENSEMBL_IDs # Necessary step for merging
	DEGH_results_Symbol <- merge(DEGH_results_Symbol, reference_table, by.x="Ensembl", by.y="ensembl_gene_id", all.x=TRUE)

	reference_table_symbol <- ensembl_to_entrez(ensembl_ids = DEGH_results_Symbol$entrezgene,
										 		organism_db = current_organism_info$Bioconductor_DB[1],
												organism_var = current_organism_info$Bioconductor_VarName_SYMBOL[1])
	colnames(reference_table_symbol) <- c("entrezgene", "Symbol")
	DEGH_results_Symbol <- merge(DEGH_results_Symbol, reference_table_symbol, by.x="entrezgene", by.y="entrezgene", all.x=TRUE)
	first_cols <- c("Ensembl", "entrezgene", "Symbol")

	DEGH_results_Symbol <- DEGH_results_Symbol[c(first_cols, setdiff(names(DEGH_results_Symbol), first_cols))] # Reorder columns so annotated first
	DEGH_results_Symbol <- DEGH_results_Symbol[order(DEGH_results_Symbol[,"combined_FDR"]),] # Reorder rows by combined FDR

	# Write output here for now, until the rest of the workflow can use this object directly
	write.table(DEGH_results_Symbol, file=file.path(paths$root, "DEG_results_annotated.txt"), quote=F, row.names=FALSE, sep="\t")
}

#############################################
### EXPORT DATA
#############################################
write.table(DEGH_results, file=file.path(paths$root, "Annotated_table.txt"), quote=F, col.names=NA, sep="\t")
write.table(reference_table, file=file.path(paths$root, "ENSEMBL2ENTREZ.txt"), quote=F, col.names=NA, sep="\t")


#############################################
### GO ENRICHMENT (topGO)
#############################################
if(flags$GO){ #TODO =>  ESTO HAY QUE TOCARLO LUEGO, POR QUE SE UTILIZA PARA UNO ENSEMBL Y OTRO ENTREZ? 
	# Prepare special subsets to be studied
	# if(exists("annot_table")){
	pos_logFC_likely_degs <- subset(likely_degs_df, likely_degs_df[fc_colname] > 0)$Annot_IDs
	neg_logFC_likely_degs <- subset(likely_degs_df, likely_degs_df[fc_colname] < 0)$Annot_IDs
	pos_logFC_union_DEGs <- subset(union_DEGs_df, union_DEGs_df[fc_colname] > 0)$Annot_IDs
	neg_logFC_union_DEGs <- subset(union_DEGs_df, union_DEGs_df[fc_colname] < 0)$Annot_IDs
	# }else{
	# 	pos_logFC_likely_degs <- rownames(subset(likely_degs_df, likely_degs_df[fc_colname] > 0))
	# 	neg_logFC_likely_degs <- rownames(subset(likely_degs_df, likely_degs_df[fc_colname] < 0))
	# 	pos_logFC_union_DEGs <- rownames(subset(union_DEGs_df, union_DEGs_df[fc_colname] > 0))
	# 	neg_logFC_union_DEGs <- rownames(subset(union_DEGs_df, union_DEGs_df[fc_colname] < 0))
	# }

	# Prepare modules to be loaded
	GO_subontologies <- c()
	if(grepl("M", opt$GO_subont)){
		GO_subontologies <- c(GO_subontologies,"MF")
	}
	if(grepl("B", opt$GO_subont)){
		GO_subontologies <- c(GO_subontologies,"BP")
	}
	if(grepl("C", opt$GO_subont)){
		GO_subontologies <- c(GO_subontologies,"CC")
	}
	if(length(GO_subontologies) == 0){
		warning("No GO sub-ontology have been selected. Please check -G option")
	}

	# Generate output
	if(length(GO_subontologies) > 0){
		# Check execution mode
		if(remote_actions$biomart){ # REMOTE MODE
			# Prepare necessary info
			go_attr_name <- as.character(current_organism_info[,"Attribute_GOs"])
			# Launch GSEA analysis
			invisible(lapply(GO_subontologies,function(mod){
				# Common
				perform_GSEA_analysis(go_attr_name, likely_degs, union_annot_DEGs_df, mod, paste("GOgraph_preval_",mod,".pdf",sep=""),opt$input_gene_id)
				perform_GSEA_analysis(go_attr_name, pos_logFC_likely_degs, union_annot_DEGs_df, mod, paste("GOgraph_preval_overex_",mod,".pdf",sep=""),opt$input_gene_id)
				perform_GSEA_analysis(go_attr_name, neg_logFC_likely_degs, union_annot_DEGs_df, mod, paste("GOgraph_preval_underex_",mod,".pdf",sep=""),opt$input_gene_id)
				# Union
				perform_GSEA_analysis(go_attr_name, union_DEGs, reference_table, mod, paste("GOgraph_allpos_",mod,".pdf",sep=""),opt$input_gene_id)
				perform_GSEA_analysis(go_attr_name, pos_logFC_union_DEGs, reference_table, mod, paste("GOgraph_allpos_overex_",mod,".pdf",sep=""),opt$input_gene_id)
				perform_GSEA_analysis(go_attr_name, neg_logFC_union_DEGs, reference_table, mod, paste("GOgraph_allpos_underex_",mod,".pdf",sep=""),opt$input_gene_id)    				
			}))
		}else{ # LOCAL MODE
			if(opt$input_gene_id == 'ensembl_gene_id'){
				# Transform ENSEMBL ids to Entrez ids
				entrez_likely_degs <- ensembl_to_entrez(likely_degs,current_organism_info$Bioconductor_DB[1],current_organism_info$Bioconductor_VarName[1]) 
				entrez_pos_logFC_likely_degs <- ensembl_to_entrez(pos_logFC_likely_degs,current_organism_info$Bioconductor_DB[1],current_organism_info$Bioconductor_VarName[1])
				entrez_neg_logFC_likely_degs <- ensembl_to_entrez(neg_logFC_likely_degs,current_organism_info$Bioconductor_DB[1],current_organism_info$Bioconductor_VarName[1])
				entrez_union_DEGs <- ensembl_to_entrez(union_DEGs,current_organism_info$Bioconductor_DB[1],current_organism_info$Bioconductor_VarName[1])
				entrez_pos_logFC_union_DEGs <- ensembl_to_entrez(pos_logFC_union_DEGs,current_organism_info$Bioconductor_DB[1],current_organism_info$Bioconductor_VarName[1])
				entrez_neg_logFC_union_DEGs <- ensembl_to_entrez(neg_logFC_union_DEGs,current_organism_info$Bioconductor_DB[1],current_organism_info$Bioconductor_VarName[1])
			}else if(opt$input_gene_id == 'refseq_peptide'){
				stop(paste("This genes type (",opt$input_gene_id,") is not allowed in local mode yet",sep="")) ##################################### NOT IMPLEMENTED
			}else{
				stop("Unchecked biomart filter type")
			}

			reference_ids_union <- unique(reference_table$entrezgene)
			reference_ids_common <- unique(union_annot_DEGs_df$entrezgene)

			# Launch GSEA analysis
			invisible(lapply(GO_subontologies,function(mod){
				# Common
				perform_GSEA_analysis_local(entrez_likely_degs$ENTREZ, reference_ids_common, mod, file.path(paths$root, paste("GO_preval",mod,sep="_")),current_organism_info$Bioconductor_DB[1])
				perform_GSEA_analysis_local(entrez_pos_logFC_likely_degs$ENTREZ, reference_ids_common, mod, file.path(paths$root, paste("GO_preval_overex",mod,sep="_")),current_organism_info$Bioconductor_DB[1])
				perform_GSEA_analysis_local(entrez_neg_logFC_likely_degs$ENTREZ, reference_ids_common,mod, file.path(paths$root, paste("GO_preval_underex",mod,sep="_")),current_organism_info$Bioconductor_DB[1])
				# Union
				perform_GSEA_analysis_local(entrez_union_DEGs$ENTREZ, reference_ids_union, mod, file.path(paths$root, paste("GO_allpos",mod,sep="_")),current_organism_info$Bioconductor_DB[1])
				perform_GSEA_analysis_local(entrez_pos_logFC_union_DEGs$ENTREZ, reference_ids_union, mod, file.path(paths$root, paste("GO_allpos_overex",mod,sep="_")),current_organism_info$Bioconductor_DB[1])
				perform_GSEA_analysis_local(entrez_neg_logFC_union_DEGs$ENTREZ, reference_ids_union, mod, file.path(paths$root, paste("GO_allpos_underex",mod,sep="_")),current_organism_info$Bioconductor_DB[1])    
			}))
		} # END LOCAL/REMOTE IF
	}
}

####################### DEBUG POINT #############################
if(opt$debug){
	time_control$TopGO <- Sys.time()
	debug_point(debug_file,"TopGO executed")
}
#################################################################

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
##                                                                                                                   ##
##                                               NORMALIZED ENRICHMENTS                                              ##                                                     
##                                                                                                                   ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
# Calculate geneList
#geneList: A NAMED VECTOR OF logFC; USED IN GSEA AND PLOTS
# if(exists("annot_table")){
	aux <- subset(reference_table, reference_table[,1] %in% DEGH_results$Annot_IDs)
	geneList <- as.vector(DEGH_results[which(DEGH_results$Annot_IDs %in% aux[,"ensembl_gene_id"]), fc_colname])
	names(geneList) <- DEGH_results$Annot_IDs[which(DEGH_results$Annot_IDs %in% aux[,"ensembl_gene_id"])]
	names(geneList) <- aux[match(names(geneList),aux[,"ensembl_gene_id"]),"entrezgene"]
# }else{
# 	aux <- subset(reference_table, reference_table[,1] %in% rownames(DEGH_results))
# 	geneList <- as.vector(DEGH_results[which(rownames(DEGH_results) %in% aux[,"ensembl_gene_id"]),fc_colname])
# 	names(geneList) <- rownames(DEGH_results)[which(rownames(DEGH_results) %in% aux[,"ensembl_gene_id"])]
# 	names(geneList) <- aux[match(names(geneList),aux[,"ensembl_gene_id"]),"entrezgene"]
# }
# Sort FC
geneList <- sort(geneList, decreasing = TRUE)


# Prepare executions
if(any(unlist(flags[c("GO_cp","KEGG","REACT")]),!is.null(opt$custom))){
	options(stringsAsFactors=FALSE)
	ora_config <- data.frame(Fun = character(),Onto = character(),Organism = character(),KeyType = character(), UseInternal = logical(), stringsAsFactors = FALSE)
	gsea_config <- data.frame(Onto = character(),Organism = character(),KeyType = character(),UseInternal = logical(), stringsAsFactors = FALSE)

	###################
	## GO
	if(flags$GO_cp){
		# Check
		if(is.na(current_organism_info$Bioconductor_DB[1]) | is.na(current_organism_info$Bioconductor_VarName[1])){
			flags$GO_cp <- FALSE
			warning("Specified organism is not allowed to be used with GO (clusterProfiler) module. Please check your IDs table")
		}else{
			# Load necessary packages	
			GO_subontologies <- c()
			if(grepl("M", opt$GO_subont)){
				GO_subontologies <- "MF"
			}
			if(grepl("B", opt$GO_subont)){
				GO_subontologies <- c(GO_subontologies,"BP")
			}
			if(grepl("C", opt$GO_subont)){
				GO_subontologies <- c(GO_subontologies,"CC")
			}
			if(length(GO_subontologies) == 0){
				warning("Any GO sub-ontology have been selected. Use -G input command")
			}
			# Add execution info
			if(flags$Clustered){
				if(flags$ORA){ 
					invisible(lapply(GO_subontologies,function(mod){ora_config <<- rbind(ora_config,list(Fun = "enrichGO",Onto=paste0("GO_",mod),Organism=current_organism_info$Bioconductor_DB[1],KeyType=keytypes,UseInternal=FALSE))}))
				}
				if(flags$GSEA){
					invisible(lapply(GO_subontologies,function(mod){gsea_config <<- rbind(gsea_config,list(Onto=paste0("GO_",mod),Organism=current_organism_info$Bioconductor_DB[1],KeyType=keytypes,UseInternal=FALSE))}))					
				}
			}			
		}
	}


	###################
	## KEGG
	if(flags$KEGG){ 
		if(is.na(current_organism_info$KeggCode[1])){
			flags$KEGG <- FALSE
			warning("Specified organism is not allowed to be used with KEGG module. Please check your IDs table")
		}else{
			if(!remote_actions$kegg){
				require(KEGG.db)
			}
			if(flags$Clustered){
				if(flags$ORA) ora_config <- rbind(ora_config,list(Fun = "enrichKEGG",Onto="KEGG",Organism=current_organism_info$KeggCode[1],KeyType="kegg",UseInternal=!remote_actions$kegg))
				if(flags$GSEA) gsea_config <- rbind(gsea_config,list(Onto="KEGG",Organism=current_organism_info$KeggCode[1],KeyType="ENTREZID",UseInternal=!remote_actions$kegg))
			}
		}
	}

	###################
	## REACTOME
	if(flags$REACT){
		if(keytypes == "GENENAME"){
			flags$REACT <- FALSE
			warning("Reactome module can not be used with GENENAME identifiers")
		}else if(is.na(current_organism_info$Reactome_ID[1]) | (keytypes != "ENTREZID")){
			flags$REACT <- FALSE
			warning("Specified organism is not allowed to be used with Reactome module. Please check your IDs table")
		}else{
			require(ReactomePA)
			if(flags$Clustered){
				if(flags$ORA) ora_config <- rbind(ora_config,list(Fun = "enrichPathway",Onto="REACT",Organism=current_organism_info$Reactome_ID[1],KeyType="ENTREZID",UseInternal=FALSE))				
				if(flags$GSEA) gsea_config <- rbind(gsea_config,list(Onto="REACT",Organism=current_organism_info$Reactome_ID[1],KeyType="ENTREZID",UseInternal=FALSE))				
			}
		}
	}


	###################
	## CUSTOM
	# Prepare custom enrichments
	if(!is.null(opt$custom)) {
		if(nchar(opt$custom)>0) {
			# Obtain custom files
			custom_files <- unlist(strsplit(opt$custom,","))
			# Per each file, launch enrichment
			custom_sets <- lapply(custom_files,function(f){
				# Load info
				c_terms <- unlist(read.table(file = f, sep = "\n", header = FALSE, stringsAsFactors = FALSE))
				# Split all
				c_terms <- as.data.frame(do.call(rbind,lapply(c_terms,function(GeneTerms){
					aux <- unlist(strsplit(GeneTerms,"\t"))
					return(data.frame(Term = rep(aux[1],length(aux)-2),
									  Gene = tail(aux,-2),
									  stringsAsFactors = FALSE))
				})))
				return(c_terms)
			})
			names(custom_sets) <- custom_files
		}else{
			custom_sets <- NULL
		}
	}else{
		custom_sets <- NULL
	}
}


#############################################
### Clustered enrichments
#############################################
if(flags$Clustered){
	cls <- unique(DEGH_results$Cluster_ID)
	# Check
	if(any(c(0,"grey") %in% cls)){
		cls <- cls[!cls %in% c(0,"grey")]
	}else{
		warning("Cluster Zero/Grey not found")
	}

	clgenes <- lapply(cls,function(cl){ # Find
		indx <- which(DEGH_results$Cluster_ID == cl)
		if(exists("annot_table")){
			unique(DEGH_results$Annot_IDs[indx])
		}else{
			unique(rownames(DEGH_results[indx,]))			
		}
	}) 
	clgenes <- lapply(clgenes,function(genes){reference_table$entrezgene[which(reference_table$ensembl_gene_id %in% genes)]}) # Translate
	names(clgenes) <- cls
	if(flags$ORA){
		message("Performing ORA enrichments")
		enrichments_ORA_expanded <- lapply(seq(nrow(ora_config)),function(i){
			# Perform per each cluster
			enr <- enrichment_clusters_ORA(genes = clgenes,organism = ora_config$Organism[i],keyType = ora_config$KeyType[i],pvalueCutoff = opt$pthreshold,pAdjustMethod = "BH",ont = ora_config$Onto[i],qvalueCutoff = opt$qthreshold, useInternal = ora_config$UseInternal[i], mc.cores = opt$cores)
			# enr <- merge_result(enr)
			return(enr)
		})
		names(enrichments_ORA_expanded) <- ora_config$Onto
		enrichments_ORA <- lapply(enrichments_ORA_expanded,function(enrCL){merge_result(enrCL)})
		# Write output
		invisible(lapply(seq_along(enrichments_ORA),function(i){
			# Concat
			df <- enrichplot:::fortify.compareClusterResult(enrichments_ORA[[i]])
			# Write table
			write.table(df, file=file.path(paths$root, paste0(ora_config$Onto[i],"_cls_ora")), quote=F, col.names=TRUE, row.names = FALSE, sep="\t")
		}))

	}
	if(flags$GSEA){
		message("Performing GSEA enrichments")
		if(exists("annot_table")){
			aux <- subset(reference_table, reference_table[,1] %in% DEGH_results$Annot_IDs)
			geneListCL <- lapply(cls, function(cl){
				glist <- as.vector(DEGH_results[which(DEGH_results$Annot_IDs %in% aux[,"ensembl_gene_id"] & DEGH_results$Cluster_ID == cl),fc_colname])
				names(glist) <- DEGH_results$Annot_IDs[which(DEGH_results$Annot_IDs %in% aux[,"ensembl_gene_id"] & DEGH_results$Cluster_ID == cl)]
				names(glist) <- aux[match(names(glist),aux[,"ensembl_gene_id"]),"entrezgene"]
				glist <- sort(glist, decreasing = TRUE)
				return(glist)
			})
		}else{
			aux <- subset(reference_table, reference_table[,1] %in% rownames(DEGH_results))
			geneListCL <- lapply(cls,function(cl){
				glist <- as.vector(DEGH_results[which(rownames(DEGH_results) %in% aux[,"ensembl_gene_id"] & DEGH_results$Cluster_ID == cl),fc_colname])
				names(glist) <- rownames(DEGH_results)[which(rownames(DEGH_results) %in% aux[,"ensembl_gene_id"] & DEGH_results$Cluster_ID == cl)]
				names(glist) <- aux[match(names(glist),aux[,"ensembl_gene_id"]),"entrezgene"]
				glist <- sort(glist, decreasing = TRUE)
				return(glist)
			})
		}
		names(geneListCL) <- cls
		enrichments_GSEA_expanded <- lapply(seq(nrow(gsea_config)),function(i){
			# Perform per each cluster
			enr <- lapply(geneListCL,function(genes){
				# Check
				if(length(genes) <= 0) return(NULL)
				# Enrich
				curr_enr <- enrichment_GSEA(geneList = genes,organism = gsea_config$Organism[i],keyType = gsea_config$KeyType[i],pvalueCutoff = opt$pthreshold,pAdjustMethod = "BH",ont = gsea_config$Onto[i], useInternal = gsea_config$UseInternal[i])
				# Check
				if(is.null(curr_enr)) return(NULL)
				if(nrow(curr_enr) <= 0) return(NULL)
				# Return
				return(curr_enr)
			})
			# # Merge all results
			# enr <- merge_result(enr)
			# Return
			return(enr)
		})
		names(enrichments_GSEA_expanded) <- gsea_config$Onto
		enrichments_GSEA <- lapply(enrichments_GSEA_expanded,function(enrCL){merge_result(enrCL)})
		names(enrichments_GSEA) <- gsea_config$Onto
		# Write output
		invisible(lapply(seq_along(enrichments_GSEA),function(i){
			# Concat
			# df <- clusterProfiler:::fortify.compareClusterResult(enrichments_GSEA[[i]])
			df <- enrichments_GSEA[[i]]@compareClusterResult
			# Write table
			write.table(df, file=file.path(paths$root, paste0(gsea_config$Onto[i],"_cls_gsea")), quote=F, col.names=TRUE, row.names = FALSE, sep="\t")
		}))
	}


	## CUSTOM
	if(!is.null(custom_sets)){
		# Execute ORA per each set
		custom_cls_ORA_expanded <- lapply(seq_along(custom_sets),function(i){
			cs_set <- custom_sets[[i]]
			# Enrich
			cs_enr <- mclapply(clgenes,function(genesset){
				enricher(genesset, pvalueCutoff = opt$pthreshold, TERM2GENE = cs_set)
			},mc.cores = opt$cores)
			names(cs_enr) <- names(clgenes)
			# Return
			return(cs_enr)
		})
		names(custom_cls_ORA_expanded) <- names(custom_sets)
		custom_cls_ORA <- lapply(custom_cls_ORA_expanded,function(enrCL){merge_result(enrCL)})

		invisible(lapply(seq_along(custom_cls_ORA),function(i){
			# Concat
			df <- enrichplot:::fortify.compareClusterResult(custom_cls_ORA[[i]])
			# Write table
			write.table(df, file=file.path(paths$root, paste0(basename(names(custom_sets)[i]),"_cls_ORA")), quote=F, col.names=TRUE, row.names = FALSE, sep="\t")
		}))
	}

####################### DEBUG POINT #############################
	if(opt$debug){
		time_control$cluster_enrichment <- Sys.time()
		debug_point(debug_file,"Clusters enriched")
	}
#################################################################

}



###########################################################################################################################
##                                                                                                                       ##
##                                                        NOT CLUSTERED                                                  ##
##                                                                                                                       ##
###########################################################################################################################

#############################################
### GO ENRICHMENT (clusterProfiler)
#############################################
if(flags$GO_cp){
	message("Performing GO enrichments")
	###########
	### ORA ENRICHMENTS
	if(flags$ORA){
		enrich_go <- lapply(GO_subontologies,function(mod){
			enrich <- enrichment_ORA(genes = likely_degs_entrez,organism = current_organism_info$Bioconductor_DB[1],keyType = keytypes,pvalueCutoff = opt$pthreshold,pAdjustMethod = "BH",ont = paste0("GO_",mod), qvalueCutoff = opt$qthreshold)
			return(enrich)
		})
		# Add names
		names(enrich_go) <- GO_subontologies
		# Write results
		write.table(as.data.frame(do.call(rbind,lapply(enrich_go,function(res){as.data.frame(res)}))), file=file.path(paths$root, "GO_CL_ora"), quote=F, col.names=TRUE, row.names = FALSE, sep="\t")
	}


	### GSEA ENRICHMENTS
	if(flags$GSEA){
		# Enrich
		enrich_go_gsea <- lapply(GO_subontologies,function(mod){
			enrich <- enrichment_GSEA(geneList = geneList,organism = current_organism_info$Bioconductor_DB[1],keyType = keytypes,pvalueCutoff = opt$pthreshold,pAdjustMethod = "BH",ont = paste0("GO_",mod))
			return(enrich)
		})
		# Add names
		names(enrich_go_gsea) <- GO_subontologies
		# Write results
		write.table(as.data.frame(do.call(rbind,lapply(enrich_go_gsea,function(res){as.data.frame(res)}))), file=file.path(paths$root, "GO_CL_gsea"), quote=F, col.names=TRUE, row.names = FALSE, sep="\t")	
	}
}



#############################################
### KEGG ENRICHMENT
#############################################

if(flags$KEGG){
	message("Performing KEGG enrichments")
	if(flags$ORA){
		# Enrich
		enrich_ora <- enrichment_ORA(genes = likely_degs_entrez,organism = current_organism_info$KeggCode[1],keyType = "kegg",pvalueCutoff = opt$pthreshold,pAdjustMethod = "BH",ont = "KEGG",useInternal = !remote_actions$kegg, qvalueCutoff = opt$qthreshold)
		# Write output
		write.table(enrich_ora, file=file.path(paths$root, "KEGG_results"), quote=F, col.names=TRUE, row.names = FALSE, sep="\t")
	}
	# Launch GSEA
	if(flags$GSEA){
		enrich_gsea <- enrichment_GSEA(geneList = geneList,organism = current_organism_info$KeggCode[1],pvalueCutoff = opt$pthreshold,ont = "KEGG",useInternal = !remote_actions$kegg)
		# Write output
		write.table(enrich_gsea, file=file.path(paths$root, "KEGG_GSEA_results"), quote=F, col.names=TRUE, row.names = FALSE, sep="\t")
	}
}


#############################################
### REACTOME ENRICHMENT
#############################################

if(flags$REACT){
	message("Performing Reactome enrichments")
	if(flags$ORA){
		# Make enrichment (ORA)
		enrich_react <- enrichment_ORA(genes = likely_degs_entrez,organism = current_organism_info$Reactome_ID[1],keyType = "ENTREZID",pvalueCutoff = opt$pthreshold,pAdjustMethod = "BH",ont = "REACT", qvalueCutoff = opt$qthreshold)		
		# Write output
		write.table(enrich_react, file=file.path(paths$root, "REACT_results"), quote=F, col.names=TRUE, row.names = FALSE, sep="\t")	
	}

	# Make enrichment GSEA
	if(flags$GSEA){
		enrich_react_gsea <- enrichment_GSEA(geneList = geneList, organism = current_organism_info$Reactome_ID[1], pvalueCutoff = opt$pthreshold, pAdjustMethod = "BH", ont = "REACT")
		# Write output
		write.table(enrich_react_gsea, file=file.path(paths$root, "REACT_GSEA_results"), quote=F, col.names=TRUE, row.names = FALSE, sep="\t")
	}
}

####################### DEBUG POINT #############################
if(opt$debug){
	time_control$regular_enrichments <- Sys.time()
	debug_point(debug_file,"Regular enrichments executed")
}
#################################################################


#############################################
### CUSTOM ENRICHMENT
#############################################
if(!is.null(opt$custom)) {
	if(nchar(opt$custom)>0) {
		# Obtain custom files
		custom_files <- unlist(strsplit(opt$custom,","))
		# Per each file, launch enrichment
		custom_enrichments <- lapply(custom_files,function(f){
			# Load info
			c_terms <- unlist(read.table(file = f, sep = "\n", header = FALSE, stringsAsFactors = FALSE))
			# Split all
			c_terms <- as.data.frame(do.call(rbind,lapply(c_terms,function(GeneTerms){
				aux <- unlist(strsplit(GeneTerms,"\t"))
				return(data.frame(Term = rep(aux[1],length(aux)-2),
								  Gene = tail(aux,-2),
								  stringsAsFactors = FALSE))
			})))
			enr <- enricher(likely_degs_entrez, pvalueCutoff = opt$pthreshold, TERM2GENE = c_terms)
			# Store results
			write.table(enr, file=file.path(paths$root, paste0(basename(f),"_ora_results")), quote=F, col.names=TRUE, row.names = FALSE, sep="\t")			
			# Return
			return(list(File = f,
						Result = enr))
		})
	}else{
		custom_enrichments <- NULL
	}
}else{
	custom_enrichments <- NULL
}

####################### DEBUG POINT #############################
if(opt$debug){
	time_control$custom_enrichments <- Sys.time()
	debug_point(debug_file,"All enrichments executed")
}
#################################################################



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
##                                                                                                                   ##
##                                                 RENDERING REPORTS                                                 ##                                                     
##                                                                                                                   ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##


message("RENDERING REPORT ...")

############################################################
##                    GENERATE REPORT                     ##
############################################################
results_path <- normalizePath(paths$root)

if(flags$Clustered){ # Clustered
	message("Rendering specific cluster reports")
	invisible(lapply(cls,function(cl){
		# Take output name
		aux <- paste0("cl_func_",cl,".html")
		outf_cls_i <- file.path(results_path, aux)
		# Generate report
		rmarkdown::render(file.path(main_path_script, 'templates', 'cl_func_report.Rmd'), output_file = outf_cls_i, intermediates_dir = results_path)
		if(opt$debug){
			time_control[[paste0("cl_func_",cl)]] <<- Sys.time()
		}	
	}))
	message("\tRendering clustered report")
	outf_cls <- file.path(results_path, "clusters_func_report.html")
	rmarkdown::render(file.path(main_path_script, 'templates', 'clusters_main_report.Rmd'),output_file = outf_cls, intermediates_dir = results_path)
	if(opt$debug){
		time_control$render_cluster_main_report <- Sys.time()
	}
}
message("\tRendering regular report")
outf <- file.path(results_path, "functional_report.html")
rmarkdown::render(file.path(main_path_script, 'templates', 'functional_report.Rmd'), output_file = outf, intermediates_dir = results_path)

if(opt$debug){
####################### DEBUG POINT #############################
    time_control$render_main_report <- Sys.time()
    debug_point(debug_file,"Report printed")
#################################################################
    save_times(time_control, output = file.path(debug_dir, "function_hunter_time_control.txt"), plot_name = "FH_times_control.pdf")
}
