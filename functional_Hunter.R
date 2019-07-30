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

# Load custom libraries
source(file.path(main_path_script, 'lib', 'general_functions.R'))
source(file.path(main_path_script, 'lib', 'functional_analysis_library.R'))









#############################################
### MAIN 
#############################################

# Parse command line
#------------------------------------------------

option_list <- list(
  make_option(c("-i", "--input_hunter_file"), type="character", default="hunter_DE_results/Common_results/hunter_results_table.txt",
    help="DEgenes Hunter's differential expression analysis output file"), # Not Checked 
  make_option(c("-c", "--countdata_file"), type="character", default="hunter_DE_results/filtered_count_data.txt",
    help="Filtered count data file"), # Not Checked
  make_option(c("-a", "--annot_file"), type="character",
    help="Two column file with annotations for functional analysis of a non-model organism. First column must be a gene ensembl id or a refseq id from a model organism. This id must be a orthologue of the gene id of the second column, that is the custom id from a non-model organism (whose functional analysis is desired)"),
  make_option(c("-m", "--model_organism"), type="character",
    help="Ortologue Species"),
  make_option(c("-t", "--biomaRt_filter"), type="character", default="E",
    help="ID types. ENSEMBL (E) Refseq_peptide (R), TAIR/Arabidopsis (T), Gene Names (G) Gene SYMBOLS (S). [Default:%default]"),      	
  make_option(c("-f", "--functional_analysis"), type="character", default="GK",
    help="Type of functional analyses to be performed (G = GO [topGO], K = KEGG, g = GO [clusterProfiler], R [Reactome]). [Default=%default]"),
  make_option(c("-G", "--GO_graphs"), type="character", default=c("M"),
    help="Modules to able go enrichments (M = Molecular Function, B = Biological Process, C = Celular Components). By default Default=%default GO cathegory is performed"), # Not Checked
  make_option(c("-K", "--Kegg_organism"), type="character", default=NULL, 
    help="Indicate organism to look for in the Kegg database for doing the path enrichment"), # Not Checked
  make_option(c("-q", "--save_query"), type="logical", default=TRUE,
    help="Save biomaRt query. By default the biomaRt query is saved"),
  make_option(c("-L", "--List_organisms"), action="store_true", type="logical", default=FALSE, 
    help="Print all organisms available at biomaRt table and ends the program"),
  make_option(c("-o", "--output_files"), type="character", default="results",
    help="Output path. Default=%default"),
  make_option(c("-r", "--remote"), ,type = "logical", action="store_true", default=FALSE,
    help="Flag to activate remote query from enrichments and Genes translation"),
  make_option(c("-T", "--threshold"), type="double", default=0.1,
    help="Enrichment p-value threshold. [Default = %default]")
)
opt <- parse_args(OptionParser(option_list=option_list))





# Special IDs
fc_colname <- "mean_logFCs"





#############################################
### LOAD AND PARSE 
#############################################

# Load available organisms
message(main_path_script)
biomaRt_query_info <- read.table(file.path(main_path_script, "lib/biomaRt_organism_table.txt"), header=T, row.names=1, sep="\t", stringsAsFactors = FALSE, fill=NA)

# Load Reference-NonRefernce models gene IDs relations
if(file.exists(opt$countdata_file)){
	count_data <- read.table(opt$countdata_file, header=TRUE, row.names=NULL, sep="\t", stringsAsFactors = FALSE)
	# exp_names <- colnames(count_data)[-1]
	dir <- dirname(opt$countdata_file)
	if(file.exists(file.path(dir,"control_treatment.txt"))){
		experiments <- read.table(file = file.path(dir,"control_treatment.txt"), sep = "\t", quote = "", header = TRUE, stringsAsFactors = FALSE)
		exp_names <- paste("[Control]",experiments[which(experiments[,1] == "C"),2],sep=" ")
		exp_names <- c(exp_names,paste("[Treatment]",experiments[which(experiments[,1] == "T"),2],sep=" "))
	}else{
		exp_names <- "EXPERIMENT NAMES NOT AVAILABLE"		
	}
}else{
	exp_names <- "EXPERIMENT NAMES NOT AVAILABLE"
}


if(!is.null(opt$annot_file)){
  annot_table <- read.table(opt$annot_file, header=FALSE, row.names=NULL, sep="\t", stringsAsFactors = FALSE, quote = "")
  reference_table <- annot_table
}else if(!file.exists(opt$countdata_file)){
    stop('Count file not exists, check the PATH given to the -c flag')
}else{
  reference_table <- count_data
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
              REACT = grepl("R", opt$functional_analysis))


# Load input
DEG_annot_table <- read.table(opt$input_hunter_file, header=T, row.names=1, sep="\t", stringsAsFactors = FALSE)
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



############ CREATE FOLDERS #########3
paths <- list()
dir.create(opt$output_files)
paths$root <-opt$output_files








#############################################
### PREPARE AND TRANSFORM DATA
#############################################

if(opt$remote){ # REMOTE MODE
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
	message(biomaRt_organism_info$Bioconductor_DB[1])
	message(biomaRt_organism_info$Bioconductor_VarName[1])

	# Check genetical items to be used
	if(opt$biomaRt_filter == 'ensembl_gene_id'){
		reference_table <- ensembl_to_entrez(ensembl_ids = reference_table[,1],
											 organism_db = biomaRt_organism_info$Bioconductor_DB[1],
											 organism_var = biomaRt_organism_info$Bioconductor_VarName[1])
		colnames(reference_table) <- c("ensembl_gene_id", biomaRt_organism_info[,"Attribute_entrez"])
	}else if(opt$biomaRt_filter == 'refseq_peptide'){
		stop(paste("This genes type (",opt$biomaRt_filter,") is not allowed in local mode yet",sep="")) ##################################### NOT IMPLEMENTED
	}
}else{
	error("Specified organism are not availvable to be studiend in LOCAL model. Please try REMOTE mode")
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

union_DEGs_df <- subset(DEG_annot_table, genes_tag=="POSSIBLE_DEG")
union_DEGs_df <- rbind(common_DEGs_df, union_DEGs_df)

if(exists("annot_table")){
	union_DEGs <- unique(union_DEGs_df$Annot_IDs)
}else{
	union_DEGs <- rownames(union_DEGs_df)
}

union_annot_DEGs_df <- subset(reference_table, reference_table[,1] %in% union_DEGs)








#############################################
### EXPORT DATA
#############################################
write.table(DEG_annot_table, file=file.path(paths$root, "Annotated_table.txt"), quote=F, col.names=NA, sep="\t")
write.table(reference_table, file=file.path(paths$root, "entrez_Gos.txt"), quote=F, col.names=NA, sep="\t")








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
		if(opt$remote){ # REMOTE MODE
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




#############################################
### GO ENRICHMENT (clusterProfiler)
#############################################

if(flags$GO_cp & !is.na(biomaRt_organism_info$Bioconductor_DB[1]) & !is.na(biomaRt_organism_info$Bioconductor_VarName[1])){
	# Load necessary packages
	require(clusterProfiler)
	
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



	# ORA ENRICHMENTS
	enrich_go <- lapply(modules_to_export,function(mod){
		enrich <-  enrichGO(gene          = common_unique_entrez, #genes,
							OrgDb         = biomaRt_organism_info$Bioconductor_DB[1], #organism,
							keyType       = keytypes, #keyType,
							ont           = mod, # SubOntology
							pvalueCutoff  = opt$threshold, #pvalueCutoff,
							pAdjustMethod = "BH") #qvalueCutoff)
		# enrich <- as.data.frame(enrich)
		return(enrich)
	# })))
	})

	

	### GSEA ENRICHMENTS
	# Obtain target genes
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

	# Enrich
	enrich_go_gsea <- lapply(modules_to_export,function(mod){
		enrich <- gseGO(geneList      = geneList,
						   OrgDb        = biomaRt_organism_info$Bioconductor_DB[1],
						   keyType       = keytypes, #keyType,
						   ont           = mod, # SubOntology
						   pvalueCutoff  = opt$threshold, #pvalueCutoff,
						   pAdjustMethod = "BH")
		return(enrich)
	# })))
	})

	# Add names
	names(enrich_go) <- modules_to_export
	names(enrich_go_gsea) <- modules_to_export

	# Write results
	write.table(as.data.frame(do.call(rbind,lapply(enrich_go,function(res){as.data.frame(res)}))), file=file.path(paths$root, "GO_CL_ora"), quote=F, col.names=TRUE, row.names = FALSE, sep="\t")	
	write.table(as.data.frame(do.call(rbind,lapply(enrich_go_gsea,function(res){as.data.frame(res)}))), file=file.path(paths$root, "GO_CL_gsea"), quote=F, col.names=TRUE, row.names = FALSE, sep="\t")	
}else if(flags$GO_cp){
	warning("Specified organism is not allowed to be used with GO (clusterProfiler) module. Please check your IDs table")
}








#############################################
### KEGG ENRICHMENT
#############################################

if(flags$KEGG & !is.na(biomaRt_organism_info$KeggCode[1])){

	# Load necessary packages
	require(clusterProfiler)
	if(!opt$remote){
		require(KEGG.db)
	}

	# Enrich
	enrich_ora <-  enrichKEGG(gene          = common_unique_entrez, #genes,
							  organism      = biomaRt_organism_info$KeggCode[1], #organism,
							  keyType       = "kegg", #keyType,
							  pvalueCutoff  = opt$threshold, #pvalueCutoff,
							  pAdjustMethod = "BH", #pAdjustMethod,
							  use_internal_data = !opt$remote)

	if(!exists("geneList")){
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
	}

	geneList <- sort(geneList, decreasing = TRUE)


	enrich_gsea <- gseKEGG(geneList     = geneList,
						   organism     = biomaRt_organism_info$KeggCode[1],
						   use_internal_data = !opt$remote,
						   # nPerm        = 1000,
						   # minGSSize    = 120,
						   pvalueCutoff = opt$threshold,
						   verbose      = FALSE)

	# Write output
	write.table(enrich_ora, file=file.path(paths$root, "KEGG_results"), quote=F, col.names=TRUE, row.names = FALSE, sep="\t")
	write.table(enrich_gsea, file=file.path(paths$root, "KEGG_GSEA_results"), quote=F, col.names=TRUE, row.names = FALSE, sep="\t")
}else if(flags$KEGG){
	warning("Specified organism is not allowed to be used with KEGG module. Please check your IDs table")
}








#############################################
### REACTOME ENRICHMENT
#############################################

if(flags$REACT & (!is.na(biomaRt_organism_info$Reactome_ID[1]) & (keytypes == "ENTREZID"))){
	# Load necessary packages
	require(ReactomePA)

	# Make enrichment (ORA)
	enrich_react <- enrichPathway(common_unique_entrez,
								 organism = biomaRt_organism_info$Reactome_ID[1],
								 pAdjustMethod = "BH",
								 # minGSSize = 10,
								 # maxGSSize = 500, 
								 # readable = FALSE,
 								 pvalueCutoff = opt$threshold)


	# Prepare enrichment (GSEA)
	if(!exists("geneList")){
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
	}

	geneList <- sort(geneList, decreasing = TRUE)

	# Make enrichment GSEA
	enrich_react_gsea <- gsePathway(geneList, 
									organism = biomaRt_organism_info$Reactome_ID[1],
									# exponent = 1, 
									# nPerm = 1000,
									# minGSSize = 10, 
									# maxGSSize = 500, 
									pvalueCutoff = opt$threshold,
									pAdjustMethod = "BH")


	# Write output
	write.table(enrich_react, file=file.path(paths$root, "REACT_results"), quote=F, col.names=TRUE, row.names = FALSE, sep="\t")
	write.table(enrich_react_gsea, file=file.path(paths$root, "REACT_GSEA_results"), quote=F, col.names=TRUE, row.names = FALSE, sep="\t")

}else if(flags$REACT & !keytypes == "GENENAME"){
	warning("Specified organism is not allowed to be used with Reactome module. Please check your IDs table")
}else if(flags$REACT){
	warning("Reactome module can not be used with GENENAME identifiers")
}









############################################################
##                    GENERATE REPORT                     ##
############################################################
# message(dirname(normalizePath(paths$root, "report.html")))
results_path <- paste(normalizePath(dirname(paths$root)),paths$root,sep=.Platform$file.sep)
outf <- paste(dirname(normalizePath(paths$root,"functional_report.html")),"functional_report.html",sep=.Platform$file.sep)

rmarkdown::render(file.path(main_path_script, 'templates', 'functional_report.Rmd'), 
                  output_file = outf, intermediates_dir = paths$root)