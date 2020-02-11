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


#Loading libraries  
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(biomaRt)) 
suppressPackageStartupMessages(require(topGO))
suppressPackageStartupMessages(require(KEGGREST))
suppressPackageStartupMessages(require(stringr))
suppressPackageStartupMessages(require(plyr))
suppressPackageStartupMessages(require(knitr))
suppressPackageStartupMessages(require(clusterProfiler))

# Load custom libraries
source(file.path(main_path_script, 'lib', 'general_functions.R'))
source(file.path(main_path_script, 'lib', 'functional_analysis_library.R'))


#############################################
### MAIN 
#############################################

# Parse command line
#------------------------------------------------

option_list <- list(
  make_option(c("-i", "--input_hunter_folder"), type="character",
    help="DEgenes Hunter's differential expression analysis output folder"), 
  make_option(c("-m", "--model_organism"), type="character",
    help="Ortologue Species"),
  make_option(c("-a", "--annot_file"), type="character",
  	help="Two column file with annotations for functional analysis of a non-model organism. First column must be a gene ensembl id or a refseq id from a model organism. This id must be a orthologue of the gene id of the second column, that is the custom id from a non-model organism (whose functional analysis is desired)"),
  make_option(c("-t", "--biomaRt_filter"), type="character", default="E",
    help="ID types. ENSEMBL (E) Refseq_peptide (R), TAIR/Arabidopsis (T), Gene Names (G) Gene SYMBOLS (S). [Default:%default]"),      	
  make_option(c("-f", "--functional_analysis"), type="character", default="GK",
    help="Type of functional analyses to be performed (G = GO [topGO], K = KEGG, g = GO [clusterProfiler], R [Reactome]). [Default=%default]"),
  make_option(c("-G", "--GO_graphs"), type="character", default=c("M"),
    help="Modules to able go enrichments (M = Molecular Function, B = Biological Process, C = Celular Components). By default Default=%default GO cathegory is performed"), # Not Checked
  make_option(c("-A", "--analysis"), type="character", default=c("go"),
    help="Analysis performance (g = Gene Set Enrichment Analysis, o = Over Representation Analysis). By default Default=%default analysis is performed"), # Not Checked
  make_option(c("-K", "--Kegg_organism"), type="character", default=NULL, 
    help="Indicate organism to look for in the Kegg database for doing the path enrichment"), # Not Checked
  make_option(c("-q", "--save_query"), type="logical", action = "store_true", default=FALSE,
    help="Flag to save biomaRt query"),
  make_option(c("-L", "--List_organisms"), action="store_true", type="logical", default=FALSE, 
    help="Print all organisms available at biomaRt table and ends the program"),
  make_option(c("-o", "--output_files"), type="character", default="results",
    help="Output path. Default=%default"),
  make_option(c("-r", "--remote"), ,type = "character", default="",
    help="Flags to activate remote query from enrichments and Genes translation. Use (b) to launch biomaRt translation; (k) to use Kegg remote data base"),
  make_option(c("-C", "--custom"), ,type = "character", default=NULL,
    help="Files with custom nomenclature (in GMT format) separated by commas (,)"),
  make_option(c("-T", "--threshold"), type="double", default=0.1,
    help="Enrichment p-value threshold. [Default = %default]"),
  make_option(c("-Q", "--qthreshold"), type="double", default=0.2,  # Currently only used in CLUSTERED enrichments
    help="Enrichment q-value threshold. [Default = %default]")

)
opt <- parse_args(OptionParser(option_list=option_list))

# Special IDs
fc_colname <- "mean_logFCs"

# Check remote actions
remote_actions <- list(biomart = grepl("b", opt$remote),
              		   kegg    = grepl("k", opt$remote))

#############################################
### LOAD AND PARSE 
#############################################

# Load available organisms
biomaRt_query_info <- read.table(file.path(main_path_script, "lib", "biomaRt_organism_table.txt"), header = TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE, fill = NA)

# Load Reference-NonRefernce models gene IDs relations
if(! file.exists(opt$input_hunter_folder)) stop("No input degenes_Hunter folder")

DEG_annot_table <- read.table(file.path(opt$input_hunter_folder, "Common_results", "hunter_results_table.txt"), header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE)

experiments <- read.table(file.path(opt$input_hunter_folder, "control_treatment.txt"), sep = "\t", quote = "", header = TRUE, stringsAsFactors = FALSE)
exp_names <- paste("[Control]",experiments[which(experiments[,1] == "C"),2],sep=" ")
exp_names <- c(exp_names,paste("[Treatment]",experiments[which(experiments[,1] == "T"),2],sep=" "))

if(!is.null(opt$annot_file)){
  annot_table <- read.table(opt$annot_file, header=FALSE, row.names=NULL, sep="\t", stringsAsFactors = FALSE, quote = "")
  reference_table <- annot_table
}else{
  reference_table <- DEG_annot_table
  reference_table[,1] <- row.names(DEG_annot_table)
}

# Prepare ID type
if(opt$biomaRt_filter == "E"){
  opt$biomaRt_filter <- 'ensembl_gene_id'
  keytypes <- "ENTREZID"
}else if(opt$biomaRt_filter == "R"){
  opt$biomaRt_filter <- 'refseq_peptide'
  keytypes <- "ENTREZID"
}else if(opt$biomaRt_filter == "T"){
  opt$biomaRt_filter <- ''
  keytypes <- "TAIR"
}else if(opt$biomaRt_filter == "G"){
  opt$biomaRt_filter <- ''
  keytypes <- "GENENAME"
}else{
  stop(paste("Given ID type (",opt$biomaRt_filter,") is not allowed.",sep=""))
}

# Check organism selected
if(opt$List_organisms == TRUE){
  print(as.character(rownames(biomaRt_query_info)))
  stop('program ends')
}else if(is.null(opt$model_organism)){
  stop('No model organism indicated. Please indicate the model organism using parameter -m. Use -L to list all model organisms available')
# }else if(!opt$organisms %in% rownames(biomaRt_query_info)){
#   stop("Organism selected is not available. Use -L to display available organims list")
}else{ # ORGANISM AVAILABLE --> LOAD
  biomaRt_organism_info <- subset(biomaRt_query_info, rownames(biomaRt_query_info) %in% opt$model_organism)  
}


# Check which enrichments are gonna be performed
flags <- list(GO    = grepl("G", opt$functional_analysis),
              KEGG  = grepl("K", opt$functional_analysis),
              GO_cp = grepl("g", opt$functional_analysis),
              REACT = grepl("R", opt$functional_analysis),
              GSEA  = grepl("g", opt$analysis),
              ORA   = grepl("o", opt$analysis),
              Clustered = "Cluster_ID" %in% colnames(DEG_annot_table))

# Special case
if(exists("annot_table")){
	DEG_annot_table$Annot_IDs <- unlist(lapply(rownames(DEG_annot_table),function(id){
		# Find model ID for non-model ID
		indx <- which(annot_table[,2] == id)
		if(length(indx) == 0){
			return(id)
		}else{
			return(annot_table[indx[1],1])
		}
	}))
}

# Verbose point
aux <- table(DEG_annot_table$genes_tag)
for(i in seq_along(aux)) 
	message(paste(names(aux)[i],aux[i]))

############ CREATE FOLDERS #########3
paths <- list()
dir.create(opt$output_files)
paths$root <-opt$output_files



#############################################
### PREPARE AND TRANSFORM DATA
#############################################

if(remote_actions$biomart){ # REMOTE MODE
	############# BIOMART QUERY  #############

	# Prepare organism info
	organism_mart <- as.character(biomaRt_organism_info[,"Mart"])
	organism_dataset <- as.character(biomaRt_organism_info[,"Dataset"])
	organism_host <- as.character(biomaRt_organism_info[,"biomaRt_Host"])
	organism_attr <- c(
	  opt$biomaRt_filter, 
	  as.character(biomaRt_organism_info[,"Attribute_GOs"]),
	  as.character(biomaRt_organism_info[,"Attribute_entrez"])
	)

	# Obtain references from biomart
	reference_table <- obtain_info_from_biomaRt(orthologues = as.character(reference_table[,1]),
	                                            id_type = opt$biomaRt_filter, 
	                                            mart = organism_mart,
	                                            dataset = organism_dataset, 
	                                            host = organism_host, 
	                                            attr = organism_attr) 

}else if(!is.na(biomaRt_organism_info$Bioconductor_DB[1])){ # LOCAL MODE
	# Check genetical items to be used
	if(opt$biomaRt_filter == 'ensembl_gene_id'){
		cat("test here\n")
		reference_table <- ensembl_to_entrez(ensembl_ids = reference_table[,1],
											 organism_db = biomaRt_organism_info$Bioconductor_DB[1],
											 organism_var = biomaRt_organism_info$Bioconductor_VarName[1])
		cat("test finished\n")

		colnames(reference_table) <- c("ensembl_gene_id", biomaRt_organism_info[,"Attribute_entrez"])

	}else if(opt$biomaRt_filter == 'refseq_peptide'){
		stop(paste("This genes type (",opt$biomaRt_filter,") is not allowed in local mode yet",sep="")) ##################################### NOT IMPLEMENTED
	}
}else{
	error("Specified organism are not available to be studiend in LOCAL model. Please try REMOTE mode")
}


if(opt$save_query == TRUE){
  saveRDS(reference_table, file=file.path("query_results_temp"))
} 

################# PREPROCESSING INPUT DATA FOR FUNCTIONAL ANALYSIS #############
# Obtain prevalent items
common_DEGs_df <- subset(DEG_annot_table, genes_tag=="PREVALENT_DEG")

# Check
if(nrow(common_DEGs_df) == 0){
  stop('Not detected genes tagged as PREVALENT_DEG. Likely some modules of degenes_Hunter.R failed and not generated output. Please, revise the input data')
} 


# Obtain significant sbusets
#  > common_DEGs_df : subset of DEGenes Hunter annotation file with PREVALENT_DEG flag
#  > common_DEGs : genes identifiers included into common_DEGs_df
#  > common_unique_entrez : list of entrez gene present into common_DEGs_df AND reference_table with filtered count data
#  > union_DEGs_df : subset of DEGenes Hunter annotation file with POSSIBLE_DEG flag or PREVALENT_DEG flag
#  > union_DEGs : genes identifiers included into union_DEGs_df
#  > union_annot_DEGs_df : reference table subset with identifiers included into union_DEGs
if(exists("annot_table")){
	common_DEGs <- unique(common_DEGs_df$Annot_IDs)
}else{
	common_DEGs <- rownames(common_DEGs_df)	
}

common_unique_entrez <- subset(reference_table, reference_table[,1] %in% common_DEGs)
common_unique_entrez <- unique(common_unique_entrez[,biomaRt_organism_info[,"Attribute_entrez"]])

# Verbose point
message(paste("IDs used to enrich:",length(common_unique_entrez)))

union_DEGs_df <- subset(DEG_annot_table, genes_tag=="POSSIBLE_DEG")
union_DEGs_df <- rbind(common_DEGs_df, union_DEGs_df)

if(exists("annot_table")){
	union_DEGs <- unique(union_DEGs_df$Annot_IDs)
}else{
	union_DEGs <- rownames(union_DEGs_df)
}

union_annot_DEGs_df <- subset(reference_table, reference_table[,1] %in% union_DEGs)


################# ADD ENTREZ IDS AND GENE SYMBOLS TO INPUT FILE #############
# Currently only runs if we are local, have an org db available and a symbol correspondence. The ENTREZ bit should become universal even if no SYMBOL available

if(!remote_actions$biomart &
	! biomaRt_organism_info$Bioconductor_DB[1] == "" & 
	! biomaRt_organism_info$Bioconductor_VarName_SYMBOL[1] == "") {

	DEG_annot_table_Symbol <- DEG_annot_table
	if(exists("annot_table")){
		ENSEMBL_IDs <- unique(DEG_annot_table_Symbol$Annot_IDs)
	}else{
		ENSEMBL_IDs <- rownames(DEG_annot_table_Symbol)
	}

	DEG_annot_table_Symbol$Ensembl <- ENSEMBL_IDs # Necessary step for merging
	DEG_annot_table_Symbol <- merge(DEG_annot_table_Symbol, reference_table, by.x="Ensembl", by.y="ensembl_gene_id", all.x=TRUE)

	reference_table_symbol <- ensembl_to_entrez(ensembl_ids = DEG_annot_table_Symbol$entrezgene,
										 		organism_db = biomaRt_organism_info$Bioconductor_DB[1],
												organism_var = biomaRt_organism_info$Bioconductor_VarName_SYMBOL[1])
	colnames(reference_table_symbol) <- c(biomaRt_organism_info[,"Attribute_entrez"], "Symbol")
	DEG_annot_table_Symbol <- merge(DEG_annot_table_Symbol, reference_table_symbol, by.x="entrezgene", by.y="entrezgene", all.x=TRUE)
	first_cols <- c("Ensembl", "entrezgene", "Symbol")

	DEG_annot_table_Symbol <- DEG_annot_table_Symbol[c(first_cols, setdiff(names(DEG_annot_table_Symbol), first_cols))] # Reorder columns so annotated first
	DEG_annot_table_Symbol <- DEG_annot_table_Symbol[order(DEG_annot_table_Symbol[,"combined_FDR"]),] # Reorder rows by combined FDR

	# Write output here for now, until the rest of the workflow can use this object directly
	write.table(DEG_annot_table_Symbol, file=file.path(paths$root, "DEG_results_annotated.txt"), quote=F, row.names=FALSE, sep="\t")
}

#############################################
### EXPORT DATA
#############################################
write.table(DEG_annot_table, file=file.path(paths$root, "Annotated_table.txt"), quote=F, col.names=NA, sep="\t")
write.table(reference_table, file=file.path(paths$root, "ENSEMBL2ENTREZ.txt"), quote=F, col.names=NA, sep="\t")


#############################################
### GO ENRICHMENT (topGO)
#############################################
if(flags$GO){
	# Prepare special subsets to be studied
	if(exists("annot_table")){
		pos_logFC_common_DEGs <- subset(common_DEGs_df, common_DEGs_df[fc_colname] > 0)$Annot_IDs
		neg_logFC_common_DEGs <- subset(common_DEGs_df, common_DEGs_df[fc_colname] < 0)$Annot_IDs
		pos_logFC_union_DEGs <- subset(union_DEGs_df, union_DEGs_df[fc_colname] > 0)$Annot_IDs
		neg_logFC_union_DEGs <- subset(union_DEGs_df, union_DEGs_df[fc_colname] < 0)$Annot_IDs
	}else{
		pos_logFC_common_DEGs <- rownames(subset(common_DEGs_df, common_DEGs_df[fc_colname] > 0))
		neg_logFC_common_DEGs <- rownames(subset(common_DEGs_df, common_DEGs_df[fc_colname] < 0))
		pos_logFC_union_DEGs <- rownames(subset(union_DEGs_df, union_DEGs_df[fc_colname] > 0))
		neg_logFC_union_DEGs <- rownames(subset(union_DEGs_df, union_DEGs_df[fc_colname] < 0))
	}

	# Prepare modules to be loaded
	modules_to_export <- c()
	if(grepl("M", opt$GO_graphs)){
		modules_to_export <- "MF"
	}
	if(grepl("B", opt$GO_graphs)){
		modules_to_export <- c(modules_to_export,"BP")
	}
	if(grepl("C", opt$GO_graphs)){
		modules_to_export <- c(modules_to_export,"CC")
	}
	if(length(modules_to_export) == 0){
		warning("Any GO sub-ontology have been selected. Use -G input command")
	}

	# Generate output
	if(length(modules_to_export) > 0){
		# Check execution mode
		if(remote_actions$biomart){ # REMOTE MODE
			# Prepare necessary info
			go_attr_name <- as.character(biomaRt_organism_info[,"Attribute_GOs"])
			# Launch GSEA analysis
			invisible(lapply(modules_to_export,function(mod){
				# Common
				perform_GSEA_analysis(go_attr_name, common_DEGs, union_annot_DEGs_df, mod, paste("GOgraph_preval_",mod,".pdf",sep=""),opt$biomaRt_filter)
				perform_GSEA_analysis(go_attr_name, pos_logFC_common_DEGs, union_annot_DEGs_df, mod, paste("GOgraph_preval_overex_",mod,".pdf",sep=""),opt$biomaRt_filter)
				perform_GSEA_analysis(go_attr_name, neg_logFC_common_DEGs, union_annot_DEGs_df, mod, paste("GOgraph_preval_underex_",mod,".pdf",sep=""),opt$biomaRt_filter)
				# Union
				perform_GSEA_analysis(go_attr_name, union_DEGs, reference_table, mod, paste("GOgraph_allpos_",mod,".pdf",sep=""),opt$biomaRt_filter)
				perform_GSEA_analysis(go_attr_name, pos_logFC_union_DEGs, reference_table, mod, paste("GOgraph_allpos_overex_",mod,".pdf",sep=""),opt$biomaRt_filter)
				perform_GSEA_analysis(go_attr_name, neg_logFC_union_DEGs, reference_table, mod, paste("GOgraph_allpos_underex_",mod,".pdf",sep=""),opt$biomaRt_filter)    				
			}))
		}else{ # LOCAL MODE
			if(opt$biomaRt_filter == 'ensembl_gene_id'){
				# Transform ENSEMBL ids to Entrez ids
				entrez_common_DEGs <- ensembl_to_entrez(common_DEGs,biomaRt_organism_info$Bioconductor_DB[1],biomaRt_organism_info$Bioconductor_VarName[1]) 
				entrez_pos_logFC_common_DEGs <- ensembl_to_entrez(pos_logFC_common_DEGs,biomaRt_organism_info$Bioconductor_DB[1],biomaRt_organism_info$Bioconductor_VarName[1])
				entrez_neg_logFC_common_DEGs <- ensembl_to_entrez(neg_logFC_common_DEGs,biomaRt_organism_info$Bioconductor_DB[1],biomaRt_organism_info$Bioconductor_VarName[1])
				entrez_union_DEGs <- ensembl_to_entrez(union_DEGs,biomaRt_organism_info$Bioconductor_DB[1],biomaRt_organism_info$Bioconductor_VarName[1])
				entrez_pos_logFC_union_DEGs <- ensembl_to_entrez(pos_logFC_union_DEGs,biomaRt_organism_info$Bioconductor_DB[1],biomaRt_organism_info$Bioconductor_VarName[1])
				entrez_neg_logFC_union_DEGs <- ensembl_to_entrez(neg_logFC_union_DEGs,biomaRt_organism_info$Bioconductor_DB[1],biomaRt_organism_info$Bioconductor_VarName[1])
			}else if(opt$biomaRt_filter == 'refseq_peptide'){
				stop(paste("This genes type (",opt$biomaRt_filter,") is not allowed in local mode yet",sep="")) ##################################### NOT IMPLEMENTED
			}else{
				stop("Unchecked biomart filter type")
			}

			reference_ids_union <- unique(reference_table$entrezgene)
			reference_ids_common <- unique(union_annot_DEGs_df$entrezgene)

			# Launch GSEA analysis
			invisible(lapply(modules_to_export,function(mod){
				# Common
				perform_GSEA_analysis_local(entrez_common_DEGs$ENTREZ, reference_ids_common, mod, file.path(paths$root, paste("GO_preval",mod,sep="_")),biomaRt_organism_info$Bioconductor_DB[1])
				perform_GSEA_analysis_local(entrez_pos_logFC_common_DEGs$ENTREZ, reference_ids_common, mod, file.path(paths$root, paste("GO_preval_overex",mod,sep="_")),biomaRt_organism_info$Bioconductor_DB[1])
				perform_GSEA_analysis_local(entrez_neg_logFC_common_DEGs$ENTREZ, reference_ids_common,mod, file.path(paths$root, paste("GO_preval_underex",mod,sep="_")),biomaRt_organism_info$Bioconductor_DB[1])
				# Union
				perform_GSEA_analysis_local(entrez_union_DEGs$ENTREZ, reference_ids_union, mod, file.path(paths$root, paste("GO_allpos",mod,sep="_")),biomaRt_organism_info$Bioconductor_DB[1])
				perform_GSEA_analysis_local(entrez_pos_logFC_union_DEGs$ENTREZ, reference_ids_union, mod, file.path(paths$root, paste("GO_allpos_overex",mod,sep="_")),biomaRt_organism_info$Bioconductor_DB[1])
				perform_GSEA_analysis_local(entrez_neg_logFC_union_DEGs$ENTREZ, reference_ids_union, mod, file.path(paths$root, paste("GO_allpos_underex",mod,sep="_")),biomaRt_organism_info$Bioconductor_DB[1])    
			}))
		} # END LOCAL/REMOTE IF
	}
}



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
##                                                                                                                   ##
##                                               NORMALIZED ENRICHMENTS                                              ##                                                     
##                                                                                                                   ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
# Calcualte geneList
if(exists("annot_table")){
	aux <- subset(reference_table, reference_table[,1] %in% DEG_annot_table$Annot_IDs)
	geneList <- as.vector(DEG_annot_table[which(DEG_annot_table$Annot_IDs %in% aux[,"ensembl_gene_id"]),fc_colname])
	names(geneList) <- DEG_annot_table$Annot_IDs[which(DEG_annot_table$Annot_IDs %in% aux[,"ensembl_gene_id"])]
	names(geneList) <- aux[match(names(geneList),aux[,"ensembl_gene_id"]),biomaRt_organism_info[,"Attribute_entrez"]]
}else{
	aux <- subset(reference_table, reference_table[,1] %in% rownames(DEG_annot_table))
	geneList <- as.vector(DEG_annot_table[which(rownames(DEG_annot_table) %in% aux[,"ensembl_gene_id"]),fc_colname])
	names(geneList) <- rownames(DEG_annot_table)[which(rownames(DEG_annot_table) %in% aux[,"ensembl_gene_id"])]
	names(geneList) <- aux[match(names(geneList),aux[,"ensembl_gene_id"]),biomaRt_organism_info[,"Attribute_entrez"]]
}
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
		if(is.na(biomaRt_organism_info$Bioconductor_DB[1]) | is.na(biomaRt_organism_info$Bioconductor_VarName[1])){
			flags$GO_cp <- FALSE
			warning("Specified organism is not allowed to be used with GO (clusterProfiler) module. Please check your IDs table")
		}else{
			# Load necessary packages	
			modules_to_export <- c()
			if(grepl("M", opt$GO_graphs)){
				modules_to_export <- "MF"
			}
			if(grepl("B", opt$GO_graphs)){
				modules_to_export <- c(modules_to_export,"BP")
			}
			if(grepl("C", opt$GO_graphs)){
				modules_to_export <- c(modules_to_export,"CC")
			}
			if(length(modules_to_export) == 0){
				warning("Any GO sub-ontology have been selected. Use -G input command")
			}
			# Add execution info
			if(flags$Clustered){
				if(flags$ORA){ 
					invisible(lapply(modules_to_export,function(mod){ora_config <<- rbind(ora_config,list(Fun = "enrichGO",Onto=paste0("GO_",mod),Organism=biomaRt_organism_info$Bioconductor_DB[1],KeyType=keytypes,UseInternal=FALSE))}))
				}
				if(flags$GSEA){
					invisible(lapply(modules_to_export,function(mod){gsea_config <<- rbind(gsea_config,list(Onto=paste0("GO_",mod),Organism=biomaRt_organism_info$Bioconductor_DB[1],KeyType=keytypes,UseInternal=FALSE))}))					
				}
			}			
		}
	}


	###################
	## KEGG
	if(flags$KEGG){ 
		if(is.na(biomaRt_organism_info$KeggCode[1])){
			flags$KEGG <- FALSE
			warning("Specified organism is not allowed to be used with KEGG module. Please check your IDs table")
		}else{
			if(!remote_actions$kegg){
				require(KEGG.db)
			}
			if(flags$Clustered){
				if(flags$ORA) ora_config <- rbind(ora_config,list(Fun = "enrichKEGG",Onto="KEGG",Organism=biomaRt_organism_info$KeggCode[1],KeyType="kegg",UseInternal=!remote_actions$kegg))
				if(flags$GSEA) gsea_config <- rbind(gsea_config,list(Onto="KEGG",Organism=biomaRt_organism_info$KeggCode[1],KeyType="ENTREZID",UseInternal=!remote_actions$kegg))
			}
		}
	}

	###################
	## REACTOME
	if(flags$REACT){
		if(keytypes == "GENENAME"){
			flags$REACT <- FALSE
			warning("Reactome module can not be used with GENENAME identifiers")
		}else if(is.na(biomaRt_organism_info$Reactome_ID[1]) | (keytypes != "ENTREZID")){
			flags$REACT <- FALSE
			warning("Specified organism is not allowed to be used with Reactome module. Please check your IDs table")
		}else{
			require(ReactomePA)
			if(flags$Clustered){
				if(flags$ORA) ora_config <- rbind(ora_config,list(Fun = "enrichPathway",Onto="REACT",Organism=biomaRt_organism_info$Reactome_ID[1],KeyType="ENTREZID",UseInternal=FALSE))				
				if(flags$GSEA) gsea_config <- rbind(gsea_config,list(Onto="REACT",Organism=biomaRt_organism_info$Reactome_ID[1],KeyType="ENTREZID",UseInternal=FALSE))				
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
	cls <- unique(DEG_annot_table$Cluster_ID)
	# Check
	if(any(c(0,"grey") %in% cls)){
		cls <- cls[-(cls %in% c(0,"grey"))]
	}else{
		warning("Cluster Zero/Grey not found")
	}
	clgenes <- lapply(cls,function(cl){unique(rownames(DEG_annot_table[which(DEG_annot_table$Cluster_ID == cl),]))}) # Find
	clgenes <- lapply(clgenes,function(genes){reference_table$entrezgene[which(reference_table$ensembl_gene_id %in% genes)]}) # Translate
	names(clgenes) <- cls
	if(flags$ORA){
		message("Performing ORA enrichments")
		enrichments_ORA <- lapply(seq(nrow(ora_config)),function(i){
			# Perform per each cluster
			enr <- enrichment_clusters_ORA(genes = clgenes,organism = ora_config$Organism[i],keyType = ora_config$KeyType[i],pvalueCutoff = opt$threshold,pAdjustMethod = "BH",ont = ora_config$Onto[i],qvalueCutoff = opt$qthreshold, useInternal = ora_config$UseInternal[i])
		})
		names(enrichments_ORA) <- ora_config$Onto
		# Write output
		invisible(lapply(seq_along(enrichments_ORA),function(i){
			# Concat
			df <- clusterProfiler:::fortify.compareClusterResult(enrichments_ORA[[i]])
			# Write table
			write.table(df, file=file.path(paths$root, paste0(ora_config$Onto[i],"_cls_ora")), quote=F, col.names=TRUE, row.names = FALSE, sep="\t")
		}))

	}
	if(flags$GSEA){
		message("Performing GSEA enrichments")
		if(exists("annot_table")){
			aux <- subset(reference_table, reference_table[,1] %in% DEG_annot_table$Annot_IDs)
			geneListCL <- lapply(cls, function(cl){
				glist <- as.vector(DEG_annot_table[which(DEG_annot_table$Annot_IDs %in% aux[,"ensembl_gene_id"] & DEG_annot_table$Cluster_ID == cl),fc_colname])
				names(glist) <- DEG_annot_table$Annot_IDs[which(DEG_annot_table$Annot_IDs %in% aux[,"ensembl_gene_id"] & DEG_annot_table$Cluster_ID == cl)]
				names(glist) <- aux[match(names(glist),aux[,"ensembl_gene_id"]),biomaRt_organism_info[,"Attribute_entrez"]]
				glist <- sort(glist, decreasing = TRUE)
				return(glist)
			})
		}else{
			aux <- subset(reference_table, reference_table[,1] %in% rownames(DEG_annot_table))
			geneListCL <- lapply(cls,function(cl){
				glist <- as.vector(DEG_annot_table[which(rownames(DEG_annot_table) %in% aux[,"ensembl_gene_id"] & DEG_annot_table$Cluster_ID == cl),fc_colname])
				names(glist) <- rownames(DEG_annot_table)[which(rownames(DEG_annot_table) %in% aux[,"ensembl_gene_id"] & DEG_annot_table$Cluster_ID == cl)]
				names(glist) <- aux[match(names(glist),aux[,"ensembl_gene_id"]),biomaRt_organism_info[,"Attribute_entrez"]]
				glist <- sort(glist, decreasing = TRUE)
				return(glist)
			})
		}
		names(geneListCL) <- cls
		enrichments_GSEA <- lapply(seq(nrow(gsea_config)),function(i){
			# Perform per each cluster
			enr <- lapply(geneListCL,function(genes){
				# Check
				if(length(genes) <= 0) return(NULL)
				# Enrich
				curr_enr <- enrichment_GSEA(geneList = genes,organism = gsea_config$Organism[i],keyType = gsea_config$KeyType[i],pvalueCutoff = opt$threshold,pAdjustMethod = "BH",ont = gsea_config$Onto[i], useInternal = gsea_config$UseInternal[i])
				# Check
				if(is.null(curr_enr)) return(NULL)
				if(nrow(curr_enr) <= 0) return(NULL)
				# Return
				return(curr_enr)
			})
			# Merge all results
			enr <- merge_result(enr)
			# Return
			return(enr)
		})
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
		custom_cls_ORA <- lapply(seq_along(custom_sets),function(i){
			cs_set <- custom_sets[[i]]
			# Enrich
			cs_enr <- lapply(clgenes,function(genesset){
				enricher(genesset, pvalueCutoff = opt$threshold, TERM2GENE = cs_set)
			})
			names(cs_enr) <- names(clgenes)
			# Concat
			cs_enr <- merge_result(cs_enr)
			# Store results
			write.table(cs_enr, file=file.path(paths$root, paste0(names(custom_sets)[i],"_cls_ORA")), quote=F, col.names=TRUE, row.names = FALSE, sep="\t")
			return(cs_enr)
		})
		names(custom_cls_ORA) <- names(custom_sets)
	}
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
		enrich_go <- lapply(modules_to_export,function(mod){
			enrich <- enrichment_ORA(genes = common_unique_entrez,organism = biomaRt_organism_info$Bioconductor_DB[1],keyType = keytypes,pvalueCutoff = opt$threshold,pAdjustMethod = "BH",ont = paste0("GO_",mod), qvalueCutoff = opt$qthreshold)
			return(enrich)
		})
		# Add names
		names(enrich_go) <- modules_to_export
		# Write results
		write.table(as.data.frame(do.call(rbind,lapply(enrich_go,function(res){as.data.frame(res)}))), file=file.path(paths$root, "GO_CL_ora"), quote=F, col.names=TRUE, row.names = FALSE, sep="\t")
	}


	### GSEA ENRICHMENTS
	if(flags$GSEA){
		# Enrich
		enrich_go_gsea <- lapply(modules_to_export,function(mod){
			enrich <- enrichment_GSEA(geneList = geneList,organism = biomaRt_organism_info$Bioconductor_DB[1],keyType = keytypes,pvalueCutoff = opt$threshold,pAdjustMethod = "BH",ont = paste0("GO_",mod))
			return(enrich)
		})
		# Add names
		names(enrich_go_gsea) <- modules_to_export
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
		enrich_ora <- enrichment_ORA(genes = common_unique_entrez,organism = biomaRt_organism_info$KeggCode[1],keyType = "kegg",pvalueCutoff = opt$threshold,pAdjustMethod = "BH",ont = "KEGG",useInternal = !remote_actions$kegg, qvalueCutoff = opt$qthreshold)
		# Write output
		write.table(enrich_ora, file=file.path(paths$root, "KEGG_results"), quote=F, col.names=TRUE, row.names = FALSE, sep="\t")
	}
	# Launch GSEA
	if(flags$GSEA){
		enrich_gsea <- enrichment_GSEA(geneList = geneList,organism = biomaRt_organism_info$KeggCode[1],pvalueCutoff = opt$threshold,ont = "KEGG",useInternal = !remote_actions$kegg)
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
		enrich_react <- enrichment_ORA(genes = common_unique_entrez,organism = biomaRt_organism_info$Reactome_ID[1],keyType = "ENTREZID",pvalueCutoff = opt$threshold,pAdjustMethod = "BH",ont = "REACT", qvalueCutoff = opt$qthreshold)		
		# Write output
		write.table(enrich_react, file=file.path(paths$root, "REACT_results"), quote=F, col.names=TRUE, row.names = FALSE, sep="\t")	
	}

	# Make enrichment GSEA
	if(flags$GSEA){
		enrich_react_gsea <- enrichment_GSEA(geneList = geneList, organism = biomaRt_organism_info$Reactome_ID[1], pvalueCutoff = opt$threshold, pAdjustMethod = "BH", ont = "REACT")
		# Write output
		write.table(enrich_react_gsea, file=file.path(paths$root, "REACT_GSEA_results"), quote=F, col.names=TRUE, row.names = FALSE, sep="\t")
	}
}









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
			enr <- enricher(common_unique_entrez, pvalueCutoff = opt$threshold, TERM2GENE = c_terms)
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


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
##                                                                                                                   ##
##                                                 RENDERING REPORTS                                                 ##                                                     
##                                                                                                                   ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##


############################################################
############################################################
############################################################
# REMOVE THIS WHEN DEVELOP/TEST PHASE ENDS
# save.image("test.RData")

message("RENDERING REPORT ...")

############################################################
##                    GENERATE REPORT                     ##
############################################################
results_path <- paste(normalizePath(dirname(paths$root)),paths$root,sep=.Platform$file.sep)
message("\tRendering regular report")
outf <- paste(dirname(normalizePath(paths$root,"functional_report.html")),"functional_report.html",sep=.Platform$file.sep)

rmarkdown::render(file.path(main_path_script, 'templates', 'functional_report.Rmd'), output_file = outf, intermediates_dir = paths$root)	
if(flags$Clustered){ # Clustered
	message("\tRendering clustered report")
	outf_cls <- paste(dirname(normalizePath(paths$root,"clusters_func_report.html")),"clusters_func_report.html",sep=.Platform$file.sep)
	rmarkdown::render(file.path(main_path_script, 'templates', 'clusters_main_report.Rmd'),output_file = outf_cls, intermediates_dir = paths$root)
}
