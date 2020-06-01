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
 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
##                                                                                                                   ##
##                                                     INITIALIZE                                                    ##                                                     
##                                                                                                                   ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
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
    help="Input gene IDs. Available IDs are: ENSEMBL (E), entrezgene (e), TAIR/Arabidopsis (T), Gene Names (G). [Default:%default]"),      	
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

# Special IDs
fc_colname <- "mean_logFCs"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
##                                                                                                                   ##
##                                              CREATE FOLDERS                                                       ## 
##                                                                                                                   ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
if (! file.exists(opt$input_hunter_folder)) stop("No input degenes_Hunter folder")


############ CREATE OUTPUT FOLDERS #########
paths <- list()
dir.create(opt$output_files)
paths$root <-opt$output_files

############ CREATE DEBUG FOLDERS #########
if (!is.null(opt$Debug)) {
	opt$debug <- TRUE
	debug_file <- file.path(opt$Debug)
}

if (opt$debug) {
	# Define only once
	if (is.null(opt$Debug)) {
		debug_file <- file.path(paths$root, debug_files, paste(c("FHunter_Debug_Session_", format(Sys.Date(), format = "%Y%m%d"), ".RData"), collapse = ""))
	}
	debug_dir <- dirname(debug_file)
	dir.create(debug_dir, recursive = T)
	debug_dir <- normalizePath(debug_dir)
	# Store session
####################### DEBUG POINT #############################
	time_control <- list(start = Sys.time())
	debug_point(debug_file, "Start point")
#################################################################
}



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
##                                                                                                                   ##
##                                       LOAD AND CHECK MAIN FILES                                                   ## 
##                                                                                                                   ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

## Load available organisms
## all_organisms_info => Table with information abaut all organism available
all_organisms_info <- read.table(file.path(main_path_script, "lib", "organism_table.txt"), header = TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE, fill = NA)

# Check organism selected
if (opt$List_organisms == TRUE) {
  print(as.character(rownames(all_organisms_info)))
  stop('Check this list and choose one model species.')

} else if (any(is.null(opt$model_organism), !opt$organism %in% rownames(all_organisms_info))) {
  stop('Model organism has not been indicated or is not available. Please indicate the model organism using parameter -m. Use -L to list all model organisms available')

} else { # ORGANISM AVAILABLE --> LOAD
  current_organism_info <- subset(all_organisms_info, rownames(all_organisms_info) %in% opt$model_organism)  
}

# Load DEGenesHunter config
# DEGenesHunter_expression_opt => DEGenesHunter options 
DEGenesHunter_expression_opt <- read.table(file.path(opt$input_hunter_folder, "opt_input_values.txt"), header = FALSE, stringsAsFactors = FALSE, sep = "\t")
degh_exp_threshold <- as.numeric(DEGenesHunter_expression_opt[which(DEGenesHunter_expression_opt[,1] == "p_val_cutoff"), 2])

experiments <- read.table(file.path(opt$input_hunter_folder, "control_treatment.txt"), sep = "\t", quote = "", header = TRUE, stringsAsFactors = FALSE)
experiments$class[experiments$class == "C"] <- "Control"
experiments$class[experiments$class == "T"] <- "Treatment"
sample_classes <- apply(experiments, 1, function(x) paste0("* [", x[1], "] ", x[2]))

## Load DEGenesHunter results
## DEGH_results => DEgenes Hunter results table and modifications
DEGH_results <- read.table(file.path(opt$input_hunter_folder, "Common_results", "hunter_results_table.txt"), header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE)
DEGH_results <- DEGH_results[DEGH_results$genes_tag != "FILTERED_OUT", ]



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
##                                                                                                                   ##
##                                         CONFIGURE ENRICHMENTS                                                     ## 
##                                                                                                                   ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

# Check remote actions
remote_actions <- list(biomart = grepl("b", opt$remote),
              		   kegg    = grepl("k", opt$remote))

# Check which enrichments are gonna be performed
# flags => A list of enricgment analysis options
flags <- list(GO    = grepl("G", opt$func_annot_db),
              KEGG  = grepl("K", opt$func_annot_db),
              GO_cp = grepl("g", opt$func_annot_db),
              REACT = grepl("R", opt$func_annot_db),
              GSEA  = grepl("g", opt$analysis_type),
              ORA   = grepl("o", opt$analysis_type),
              WGCNA = "Cluster_ID" %in% colnames(DEGH_results))


need_translate <- TRUE
# Prepare and check given IDs 
if (opt$input_gene_id == "e"){
	opt$input_gene_id <- 'entrezgene'
	need_translate <- FALSE

} else if (opt$input_gene_id == "E") {
  opt$input_gene_id <- 'ensembl_gene_id'
  keytypes <- "ENTREZID"

} else if (opt$input_gene_id == "T") {
  opt$input_gene_id <- ''
  keytypes <- "TAIR"

} else if (opt$input_gene_id == "G") {
  opt$input_gene_id <- ''
  keytypes <- "GENENAME"

} else {
  stop(paste0("Given ID type (", opt$input_gene_id, ") is still not allowed. Please check -t option."))

}


if (flags$GO_cp && current_organism_info$Bioconductor_VarName[1] == "") {
		flags$GO_cp <- FALSE
		message("Specified organism is not allowed to be used with GO (clusterProfiler) module. Please check your IDs table")
}

if (flags$KEGG && current_organism_info$KeggCode[1] == "") {
		flags$KEGG <- FALSE
		warning("Specified organism is not allowed to be used with KEGG module. Please check your IDs table")
}

if (flags$REACT) {
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

# Check which subontologies are gonna be performed
# GO_subontologies => vector with subontologies to be perfromed
if (any(flags$GO, flags$GO_cp)) {
	GO_subontologies <- c()
	if (grepl("M", opt$GO_subont)) {
		GO_subontologies <- c(GO_subontologies,"MF")
	}
	if (grepl("B", opt$GO_subont)) {
		GO_subontologies <- c(GO_subontologies,"BP")
	}
	if (grepl("C", opt$GO_subont)) {
		GO_subontologies <- c(GO_subontologies,"CC")
	}
}


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
##                                                                                                                   ##
##                                         LOAD COMPLEMENTARY FILES                                                  ## 
##                                                                                                                   ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
#Load annotation gene table for translation ID
if (!is.null(opt$annot_file)) {
	annot_table <- read.table(opt$annot_file, header=FALSE, row.names=NULL, sep="\t", stringsAsFactors = FALSE, quote = "") # Load
  	DEGH_results$input_IDs <- translate_id(rownames(DEGH_results), annot_table)
} else {
	DEGH_results$input_IDs <- rownames(DEGH_results)
}
added_cols <- c("input_IDs")

# Load Normalized correlation data
if (flags$WGCNA) {

	####
	# LOAD NORMALIZED COUNTS
	norm_counts <- as.matrix(read.table(file.path(opt$input_hunter_folder, "Results_DESeq2", "Normalized_counts_DESeq2.txt"), header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE))
	scaled_counts <- scale_counts_table(norm_counts)
	scaled_counts_table <- as.data.frame(as.table(scaled_counts))
	colnames(scaled_counts_table) <- c("Gene","Sample","Count")
		
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
	wgcna_count_sample_trait <- scale_counts_table(wgcna_count_sample_trait)
}


###############################
## LOAD CUSTOM GMT 	
all_custom_gmt <- NULL
if (!is.null(opt$custom)) {
	custom_files <- unlist(strsplit(opt$custom, ","))
	all_custom_gmt <- lapply(custom_files, load_and_parse_gmt)
	names(all_custom_gmt) <- custom_files
}

# Verbose point 
aux <- table(DEGH_results$genes_tag)
for(i in seq_along(aux)) {
	message(paste(names(aux)[i],aux[i]))
}


####################### DEBUG POINT #############################
if (opt$debug) {
	time_control$load_data <- Sys.time()
	debug_point(debug_file, "Files have been loaded and parsed")
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
		  opt$input_gene_id, 
		  as.character(current_organism_info[,"Attribute_GOs"]),
		  as.character(current_organism_info[,"Attribute_entrez"])
		)

		# Obtain references from biomart
		input_to_entrezgene <- obtain_info_from_biomaRt(orthologues = as.character(valid_genes),
		                                            id_type = opt$input_gene_id, 
		                                            mart = organism_mart,
		                                            dataset = organism_dataset, 
		                                            host = organism_host, 
		                                            attr = organism_attr) 

	} else if (current_organism_info$Bioconductor_DB[1] != "") { # LOCAL MODE
		# Check genetical items to be used
		if (opt$input_gene_id == 'ensembl_gene_id') {
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

	if (opt$save_query == TRUE) { ## TODO => ESTO HAY QUE DARLE UNA FUNCIONALIDAD O QUITARLO
	  saveRDS(input_to_entrezgene, file=file.path("query_results_temp"))
	} 
	
	################# ADD ENTREZ IDS TO INPUT FILE #############
	DEGH_results$entrezgene <- input_to_entrezgene[match(DEGH_results$input_IDs, input_to_entrezgene$ensembl_gene_id), "entrezgene"]
	added_cols <- c(added_cols, "entrezgene")
	message(paste(opt$input_gene_id, "IDs have been translated to entrezgene"))
} else {
	################# ADD ENTREZ IDS TO INPUT FILE #############
	DEGH_results$entrezgene <- DEGH_results$Input_IDs
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
if (opt$debug) {
	time_control$translate_id <- Sys.time()
	debug_point(debug_file, "IDs translated")
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
if (opt$debug) {
	time_control$data_preparation <- Sys.time()
	debug_point(debug_file, "Data prepared for enrichments")
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
			perform_topGO(go_attr_name, likely_degs_entrez, union_annot_DEGs_df, mod, paste("GOgraph_preval_",mod,".pdf",sep=""),opt$input_gene_id)
			perform_topGO(go_attr_name, pos_logFC_likely_degs, union_annot_DEGs_df, mod, paste("GOgraph_preval_overex_",mod,".pdf",sep=""),opt$input_gene_id)
			perform_topGO(go_attr_name, neg_logFC_likely_degs, union_annot_DEGs_df, mod, paste("GOgraph_preval_underex_",mod,".pdf",sep=""),opt$input_gene_id)
			# Union
			perform_topGO(go_attr_name, union_DEGs, input_to_entrezgene, mod, paste("GOgraph_allpos_",mod,".pdf",sep=""),opt$input_gene_id)
			perform_topGO(go_attr_name, pos_logFC_union_DEGs, input_to_entrezgene, mod, paste("GOgraph_allpos_overex_",mod,".pdf",sep=""),opt$input_gene_id)
			perform_topGO(go_attr_name, neg_logFC_union_DEGs, input_to_entrezgene, mod, paste("GOgraph_allpos_underex_",mod,".pdf",sep=""),opt$input_gene_id)    				
		}))
	} else { # LOCAL MODE

		reference_ids_union <- unique(DEGH_results$entrezgene)
		reference_ids_common <- unique(DEGH_results$entrezgene)

		# Launch topGO analysis
		invisible(lapply(GO_subontologies,function(mod) {
			# Common
			perform_topGO_local(likely_degs_entrez, reference_ids_common, mod, file.path(paths$root, paste("GO_preval",mod,sep="_")),current_organism_info$Bioconductor_DB[1])
			perform_topGO_local(pos_logFC_likely_degs, reference_ids_common, mod, file.path(paths$root, paste("GO_preval_overex",mod,sep="_")),current_organism_info$Bioconductor_DB[1])
			perform_topGO_local(neg_logFC_likely_degs, reference_ids_common,mod, file.path(paths$root, paste("GO_preval_underex",mod,sep="_")),current_organism_info$Bioconductor_DB[1])
			# # Union
			perform_topGO_local(union_DEGs, reference_ids_union, mod, file.path(paths$root, paste("GO_allpos",mod,sep="_")),current_organism_info$Bioconductor_DB[1])
			perform_topGO_local(pos_logFC_union_DEGs, reference_ids_union, mod, file.path(paths$root, paste("GO_allpos_overex",mod,sep="_")),current_organism_info$Bioconductor_DB[1])
			perform_topGO_local(neg_logFC_union_DEGs, reference_ids_union, mod, file.path(paths$root, paste("GO_allpos_underex",mod,sep="_")),current_organism_info$Bioconductor_DB[1])    
		}))
	} # END LOCAL/REMOTE IF
	####################### DEBUG POINT #############################
	if (opt$debug) {
		time_control$TopGO <- Sys.time()
		debug_point(debug_file, "topGO analysis finished")
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
			enrich <- enrichment_ORA(genes = likely_degs_entrez, organism = current_organism_info$Bioconductor_DB[1],keyType = keytypes,pvalueCutoff = opt$pthreshold,pAdjustMethod = "BH", ont = paste0("GO_",mod), qvalueCutoff = opt$qthreshold)
			return(enrich)
		})
		# Add names
		names(enrich_go) <- GO_subontologies
		# Write results
		write.table(as.data.frame(do.call(rbind,lapply(enrich_go,function(res) {as.data.frame(res)}))), file=file.path(paths$root, "GO_CL_ora"), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")
		
		#PERFORM WGCNA MODULES ANALYSIS
		if (flags$WGCNA) {
			for (subont in GO_subontologies) {
				subont_name <- paste0("GO_", subont)
				enrichments_ORA_expanded[[subont_name]] <- enrichment_clusters_ORA(genes = clgenes,
	                                 organism = current_organism_info$Bioconductor_DB[1],
	                                 keyType = keytypes,
	                                 pvalueCutoff = opt$pthreshold,
	                                 pAdjustMethod = "BH",
	                                 ont = subont_name,
	                                 qvalueCutoff = opt$qthreshold,
	                                 useInternal = FALSE,
	                                 mc.cores = opt$cores)
			}
		}
	}

	###########
	### GSEA ENRICHMENTS
	if (flags$GSEA) {
		# Enrich
		enrich_go_gsea <- lapply(GO_subontologies,function(mod) {
			enrich <- enrichment_GSEA(geneList = geneList,organism = current_organism_info$Bioconductor_DB[1],keyType = keytypes,pvalueCutoff = opt$pthreshold,pAdjustMethod = "BH",ont = paste0("GO_",mod))
			return(enrich)
		})
		# Add names
		names(enrich_go_gsea) <- GO_subontologies
		# Write results
		write.table(as.data.frame(do.call(rbind,lapply(enrich_go_gsea,function(res) {as.data.frame(res)}))), file=file.path(paths$root, "GO_CL_gsea"), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")	

		#PERFORM WGCNA MODULES ANALYSIS
		if (flags$WGCNA) {
			for (subont in GO_subontologies) {
				subont_name <- paste0("GO_", subont)
				enrichments_GSEA_expanded[[subont_name]] <- perform_GSEA_clusters(all_clusters = cl_gene_fc,
									organism = current_organism_info$Bioconductor_DB[1],
									keyType = keytypes,
									pvalueCutoff = opt$pthreshold,
									pAdjustMethod = "BH",
									ont = subont_name, 
									useInternal = FALSE)
			}
		}
	}
	####################### DEBUG POINT #############################
	if (opt$debug) {
		time_control$clusterProfiler_GO <- Sys.time()
		debug_point(debug_file, "clusterProfiler GO analysis finished")
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
		enrich_ora <- enrichment_ORA(genes = likely_degs_entrez,organism = current_organism_info$KeggCode[1],keyType = "kegg",pvalueCutoff = opt$pthreshold,pAdjustMethod = "BH",ont = "KEGG",useInternal = !remote_actions$kegg, qvalueCutoff = opt$qthreshold)
		# Write output
		write.table(enrich_ora, file=file.path(paths$root, "KEGG_results"), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")
		
		#PERFORM WGCNA MODULES ANALYSIS
		if (flags$WGCNA) {
			enrichments_ORA_expanded[["KEGG"]] <- enrichment_clusters_ORA(genes = clgenes,
                                 organism = current_organism_info$KeggCode[1],
                                 keyType = "kegg",
                                 pvalueCutoff = opt$pthreshold,
                                 pAdjustMethod = "BH",
                                 ont = "KEGG",
                                 qvalueCutoff = opt$qthreshold,
                                 useInternal = !remote_actions$kegg,
                                 mc.cores = opt$cores)
		}
	}

	###########
	### GSEA ENRICHMENTS
	if (flags$GSEA) {
		enrich_gsea <- enrichment_GSEA(geneList = geneList,organism = current_organism_info$KeggCode[1],pvalueCutoff = opt$pthreshold,ont = "KEGG",useInternal = !remote_actions$kegg)
		# Write output
		write.table(enrich_gsea, file=file.path(paths$root, "KEGG_GSEA_results"), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")
		
		#PERFORM WGCNA MODULES ANALYSIS
		if (flags$WGCNA) {
			enrichments_GSEA_expanded[["KEGG"]] <- perform_GSEA_clusters(all_clusters = cl_gene_fc,
								organism = current_organism_info$KeggCode[1],
								keyType = "ENTREZID",
								pvalueCutoff = opt$pthreshold,
								pAdjustMethod = "BH",
								ont = "KEGG", 
								useInternal = !remote_actions$kegg)
		}

	}
	####################### DEBUG POINT #############################
	if (opt$debug) {
		time_control$clusterProfiler_KEGG <- Sys.time()
		debug_point(debug_file, "clusterProfiler KEGG analysis finished")
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
		enrich_react <- enrichment_ORA(genes = likely_degs_entrez,organism = current_organism_info$Reactome_ID[1],keyType = "ENTREZID",pvalueCutoff = opt$pthreshold,pAdjustMethod = "BH",ont = "REACT", qvalueCutoff = opt$qthreshold)		
		# Write output
		write.table(enrich_react, file=file.path(paths$root, "REACT_results"), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")	
		
		#PERFORM WGCNA MODULES ANALYSIS
		if (flags$WGCNA) {
			enrichments_ORA_expanded[["REACT"]] <- enrichment_clusters_ORA(genes = clgenes,
                                 organism = current_organism_info$Reactome_ID[1],
                                 keyType = "ENTREZID",
                                 pvalueCutoff = opt$pthreshold,
                                 pAdjustMethod = "BH",
                                 ont = "REACT",
                                 qvalueCutoff = opt$qthreshold,
                                 useInternal = FALSE,
                                 mc.cores = opt$cores)
		}
	}

	###########
	### GSEA ENRICHMENTS
	if (flags$GSEA) {
		enrich_react_gsea <- enrichment_GSEA(geneList = geneList, organism = current_organism_info$Reactome_ID[1], pvalueCutoff = opt$pthreshold, pAdjustMethod = "BH", ont = "REACT")
		# Write output
		write.table(enrich_react_gsea, file=file.path(paths$root, "REACT_GSEA_results"), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")
		if (flags$WGCNA) {
			enrichments_GSEA_expanded[["REACT"]] <- perform_GSEA_clusters(all_clusters = cl_gene_fc,
								organism = current_organism_info$Reactome_ID[1],
								keyType = "ENTREZID",
								pvalueCutoff = opt$pthreshold,
								pAdjustMethod = "BH",
								ont = "REACT", 
								useInternal = FALSE)
		}
	}
	####################### DEBUG POINT #############################
	if (opt$debug) {
		time_control$clusterProfiler_REACTOME <- Sys.time()
		debug_point(debug_file, "clusterProfiler REACTOME analysis finished")
	}
	message("clusterProfiler REACTOME analysis finished")
	#################################################################
	#PERFORM WGCNA MODULES ANALYSIS
}

#############################################
### CUSTOM ENRICHMENT
if (!is.null(opt$custom)) {
	# Obtain custom files
	custom_files <- unlist(strsplit(opt$custom,","))

	###########
	### CUSTOM ENRICHMENTS
	# Per each file, launch enrichment
	custom_enrichments <- enrich_all_customs(custom_files = custom_files, 
											p_val_threshold = opt$pthreshold, 
												likely_degs_entrez = likely_degs_entrez)
	if (flags$WGCNA){
		custom_cls_ORA_expanded <- lapply(all_custom_gmt, function(gmt){
			enrich_clusters_with_gmt(gmt, clgenes)
		})
	}
	####################### DEBUG POINT #############################
	if (opt$debug) {
		time_control$TopGO <- Sys.time()
		debug_point(debug_file, "CUSTOM enrichments finished")
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
			write.table(df, file=file.path(paths$root, paste0(names(enrichments_ORA[enrichment_i]),"_cls_ora")), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")
		}
	}

	if (flags$GSEA) {
		enrichments_GSEA <- lapply(enrichments_GSEA_expanded, merge_result)

		for (enrichment_i in 1:length(enrichments_GSEA)) {
			# Concat
			# df <- clusterProfiler:::fortify.compareClusterResult(enrichments_GSEA[[enrichment_i]])
			df <- enrichments_GSEA[[enrichment_i]]@compareClusterResult
			# Write table
			write.table(df, file=file.path(paths$root, paste0(names(enrichments_ORA[enrichment_i]),"_cls_gsea")), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")
		}
	}

	if (!is.null(opt$custom)){
		custom_cls_ORA <- lapply(custom_cls_ORA_expanded, merge_result)
		for(enrichment_i in 1:length(custom_cls_ORA)) {
			# Concat
			df <- enrichplot:::fortify.compareClusterResult(custom_cls_ORA[[enrichment_i]])
			# Write table
			write.table(df, file=file.path(paths$root, paste0(names(enrichments_ORA[enrichment_i]),"_cls_ORA")), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")
		}
	}
}

DEGH_results <- DEGH_results[c(added_cols, setdiff(names(DEGH_results), added_cols))] # Reorder columns so annotated first
DEGH_results <- DEGH_results[order(DEGH_results[,"combined_FDR"]),] # Reorder rows by combined FDR
write.table(DEGH_results, file=file.path(paths$root, "hunter_results_table_annotated.txt"), quote=FALSE, col.names=NA, sep="\t")

####################### DEBUG POINT #############################
if (opt$debug) {
	time_control$TopGO <- Sys.time()
	debug_point(debug_file, "Results saved")
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
results_path <- normalizePath(paths$root)

if (flags$WGCNA) { # Clustered
	message("Rendering specific cluster reports")
	invisible(lapply(cls, function(cl) {
		# Take output name
		aux <- paste0("cl_func_",cl,".html")
		outf_cls_i <- file.path(results_path, aux)
		# Generate report
		rmarkdown::render(file.path(main_path_script, 'templates', 'cl_func_report.Rmd'), output_file = outf_cls_i, intermediates_dir = results_path)
		if (opt$debug) {
			time_control[[paste0("cl_func_",cl)]] <<- Sys.time()
		}
	}))

	message("\tRendering clustered report")
	outf_cls <- file.path(results_path, "clusters_func_report.html")
	rmarkdown::render(file.path(main_path_script, 'templates', 'clusters_main_report.Rmd'),output_file = outf_cls, intermediates_dir = results_path)
	if (opt$debug) {
		time_control$render_cluster_main_report <- Sys.time()
	}
}

############################################################
##              GENERATE DEG FUNCTIONAL REPORT            ##
############################################################
message("\tRendering regular report")
outf <- file.path(results_path, "functional_report.html")
rmarkdown::render(file.path(main_path_script, 'templates', 'functional_report.Rmd'), output_file = outf, intermediates_dir = results_path)

if (opt$debug) {
####################### DEBUG POINT #############################
    time_control$render_main_report <- Sys.time()
    debug_point(debug_file,"Report printed")
#################################################################
    save_times(time_control, output = file.path(debug_dir, "function_hunter_time_control.txt"), plot_name = "FH_times_control.pdf")
}
