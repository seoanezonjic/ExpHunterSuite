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
    help="IDtype, 'E' for 'ensembl_gene_id' or 'R' for 'refseq_peptide'. [Default:%default]"),      	
  make_option(c("-f", "--functional_analysis"), type="character", default="GK",
    help="Type of functional analyses to be performed (G = GO, K = KEGG). [Default=%default]"),
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
    help="Flag to activate remote query from enrichments and Genes translation")
)
opt <- parse_args(OptionParser(option_list=option_list))











#############################################
### LOAD AND PARSE 
#############################################

# Load available organisms
biomaRt_query_info <- read.table(file.path(main_path_script, "lib/biomaRt_organism_table.txt"), header=T, row.names=1, sep="\t", stringsAsFactors = FALSE)

# Load Reference-NonRefernce models gene IDs relations
if(!is.null(opt$annot_file)){
  reference_table <- read.table(opt$annot_file, header=FALSE, row.names=NULL, sep="\t", stringsAsFactors = FALSE)
  # reference_table <- annot  
}else if(!file.exists(opt$countdata_file)){
    stop('Count file not exists, check the PATH given to the -c flag')
}else{
  reference_table <- read.table(opt$countdata_file, header=TRUE, row.names=NULL, sep="\t", stringsAsFactors = FALSE)
}


# Prepare ID type
if(opt$biomaRt_filter == "E"){
  opt$biomaRt_filter <- 'ensembl_gene_id'
  #opt$biomaRt_filter <- as.character(biomaRt_organism_info[1,7])
}else if(opt$biomaRt_filter == "R"){
  opt$biomaRt_filter <- 'refseq_peptide'
  #opt$biomaRt_filter <- as.character(biomaRt_organism_info[1,8])
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
flags <- list(GO   = grepl("G", opt$functional_analysis),
              KEGG = grepl("K", opt$functional_analysis))


# Load input
DEG_annot_table <- read.table(opt$input_hunter_file, header=T, row.names=1, sep="\t", stringsAsFactors = FALSE)




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

}else{ # LOCAL MODE
	# Check genetical items to be used
	if(opt$biomaRt_filter == 'ensembl_gene_id'){
		reference_table <- human_ENSEMBL_to_ENTREZ(reference_table[,1])
		colnames(reference_table) <- c("ensembl_gene_id", "entrezgene")
	}else if(opt$biomaRt_filter == 'refseq_peptide'){
		stop(paste("This genes type (",opt$biomaRt_filter,") is not allowed in local mode yet",sep="")) ##################################### NOT IMPLEMENTED
	}
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
common_DEGs <- rownames(common_DEGs_df)
common_unique_entrez <- subset(reference_table, reference_table[,1] %in% common_DEGs)
common_unique_entrez <- unique(common_unique_entrez[,"entrezgene"])

union_DEGs_df <- subset(DEG_annot_table, genes_tag=="POSSIBLE_DEG")
union_DEGs_df <- rbind(common_DEGs_df, union_DEGs_df)
union_DEGs <- rownames(union_DEGs_df)

union_annot_DEGs_df <- subset(reference_table, reference_table[,1] %in% union_DEGs)








#############################################
### EXPORT DATA
#############################################
write.table(DEG_annot_table, file=file.path(paths$root, "Annotated_table.txt"), quote=F, col.names=NA, sep="\t")
write.table(reference_table, file=file.path(paths$root, "entrez_Gos.txt"), quote=F, col.names=NA, sep="\t")








#############################################
### GO ENRICHMENT
#############################################
if(flags$GO){
	# Prepare special subsets to be studied
	pos_logFC_common_DEGs <- rownames(subset(common_DEGs_df, common_DEGs_df["logFC_DESeq2"] > 0))
	neg_logFC_common_DEGs <- rownames(subset(common_DEGs_df, common_DEGs_df["logFC_DESeq2"] < 0))
	pos_logFC_union_DEGs <- rownames(subset(union_DEGs_df, union_DEGs_df["logFC_DESeq2"] > 0))
	neg_logFC_union_DEGs <- rownames(subset(union_DEGs_df, union_DEGs_df["logFC_DESeq2"] < 0))

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
				entrez_common_DEGs <- human_ENSEMBL_to_ENTREZ(common_DEGs) 
				entrez_pos_logFC_common_DEGs <- human_ENSEMBL_to_ENTREZ(pos_logFC_common_DEGs)
				entrez_neg_logFC_common_DEGs <- human_ENSEMBL_to_ENTREZ(neg_logFC_common_DEGs)
				entrez_union_DEGs <- human_ENSEMBL_to_ENTREZ(union_DEGs)
				entrez_pos_logFC_union_DEGs <- human_ENSEMBL_to_ENTREZ(pos_logFC_union_DEGs)
				entrez_neg_logFC_union_DEGs <- human_ENSEMBL_to_ENTREZ(neg_logFC_union_DEGs)
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
				perform_GSEA_analysis_local(entrez_common_DEGs$ENTREZ, reference_ids_common, mod, file.path(paths$root, paste("GO_res_common",mod,sep="_")))
				perform_GSEA_analysis_local(entrez_pos_logFC_common_DEGs$ENTREZ, reference_ids_common, mod, file.path(paths$root, paste("GO_res_common_pos",mod,sep="_")))
				perform_GSEA_analysis_local(entrez_neg_logFC_common_DEGs$ENTREZ, reference_ids_common,mod, file.path(paths$root, paste("GO_res_common_neg",mod,sep="_")))
				# Union
				perform_GSEA_analysis_local(entrez_union_DEGs$ENTREZ, reference_ids_union, mod, file.path(paths$root, paste("GO_res_union",mod,sep="_")))
				perform_GSEA_analysis_local(entrez_pos_logFC_union_DEGs$ENTREZ, reference_ids_union, mod, file.path(paths$root, paste("GO_res_union_pos",mod,sep="_")))
				perform_GSEA_analysis_local(entrez_neg_logFC_union_DEGs$ENTREZ, reference_ids_union, mod, file.path(paths$root, paste("GO_res_union_neg",mod,sep="_")))    
			}))
		} # END LOCAL/REMOTE IF
	}
}










#############################################
### KEGG ENRICHMENT
#############################################

if(flags$KEGG){
	# Load necessary packages
	require(clusterProfiler)
	if(!opt$remote){
		require(KEGG.db)
	}

	# Enrich
	enrich <-  enrichKEGG(gene          = common_unique_entrez, #genes,
						  organism      = "hsa", #organism,
						  keyType       = "kegg", #keyType,
						  pvalueCutoff  = 1, #pvalueCutoff,
						  pAdjustMethod = "BH", #pAdjustMethod,
						  use_internal_data = opt$remote, 
						  qvalueCutoff  = 1) #qvalueCutoff)

	# Write output
	write.table(enrich, file=file.path(paths$root, "KEGG_results"), quote=F, col.names=TRUE, row.names = FALSE, sep="\t")
}







#############################################
### OUTPUT
#############################################

################## CREATE FINAL ANNOTATION TABLE ##########

#### functional
generate_FA_report()
