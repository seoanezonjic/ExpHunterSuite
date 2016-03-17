#! /usr/bin/env Rscript
#############################################
### Rocio Bautista & Isabel Gonzalez , 2014. 
#############################################

# this is wrapped in a tryCatch. The first expression works when source executes, the
# second expression works when R CMD does it.
full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile),  # works when using source
               error=function(e) # works when using R CMD
              normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', commandArgs())], '='))[2]))
main_path_script <- dirname(full.fpath)


#Loading libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(biomaRt)) 
suppressPackageStartupMessages(library(topGO))
suppressPackageStartupMessages(library(KEGGREST))
suppressPackageStartupMessages(library(FGNet))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(GOplot))
suppressPackageStartupMessages(library(KEGGprofile))

# Load custom libraries
source(file.path(main_path_script, 'lib', 'general_functions.R'))
source(file.path(main_path_script, 'lib', 'functional_analysis_functions.R'))

countdata_file <- read.table(file.path(main_path_script, "filtered_count_data.txt"), header=T, row.names=NULL, sep="\t")
biomaRt_query_info <- read.table(file.path(main_path_script, "lib/biomaRt_organism_table.txt"), header=T, row.names=1, sep="\t")

#############################################
### MAIN 
#############################################

### Input/Output (I/O)
#############################################

# Parse command line
#------------------------------------------------

option_list <- list(
	make_option(c("-i", "--input_hunter_file"), type="character", default=NULL,
		help="DEgenes Hunter's differential expression analysis output file"), 
	make_option(c("-a", "--annot_file"), type="character", default=NULL,
		help="File with annotations for functional analysis of a non-model organism"),
  make_option(c("-L", "--List_organisms"), action="store_true", type="logical", default=FALSE, 
      help="List all organisms provided"),
  make_option(c("-m", "--model_organism"), type="character", default=NULL,
    	help="Ortologue Species. Default=%default"),
  make_option(c("-t", "--biomaRt_filter"), type="character", default="E",
      help="IDtype"),       	
	make_option(c("-f", "--functional_analysis"), type="character", default=c("GK"),
		help="Selection of type of functional analyses to be performed (G = GO, K = KEGG).
    By default the following modules Default=%default are performed"),
  make_option(c("-g", "--GO_graphs"), type="character", default=c("M"),
    help="Modules to able/disable (M = Molecular Function, B = Biological Process, C = Celular Components). By default Default=%default GO cathegory is performed"),
  make_option(c("-q", "--save_query"), type="logical", default=FALSE, 
    help="Save biomaRt query. By default the biomaRt query is saved"), 
	make_option(c("-o", "--output_files"), type="character", default="",
    	help="Output path. Default=%default")
)
opt <- parse_args(OptionParser(option_list=option_list))


# Parse input data 
#--------------------------------------------------

DEG_annot_table <- read.table(opt$input_hunter_file, header=T, row.names=1, sep="\t")

if (!is.null(opt$annot_file)){
  annot <- read.table(opt$annot_file, header=F, row.names=NULL, sep="\t")
  reference_table <- annot  
} else {
  reference_table <- countdata_file
}

biomaRt_organism_info <- get_specific_dataframe_names(biomaRt_query_info, rownames(biomaRt_query_info), opt$model_organism)


if (opt$biomaRt_filter == "E"){
  opt$biomaRt_filter <- as.character(biomaRt_organism_info[1,7])
}

if (opt$biomaRt_filter == "R"){
  opt$biomaRt_filter <- as.character(biomaRt_organism_info[1,8])
}


###### Input control ######
if(opt$List_organisms == TRUE){
  print(as.character(rownames(biomaRt_query_info)))
  stop('program ends')
  } else if (is.null(opt$model_organism)){
    stop('No model organism indicated. Please indicate the
     model organism using parameter -m. Use -L to list all model organisms available')
}


############ CREATE FOLDERS #########3

paths <- list()
dir.create(opt$output_files)
paths$root <-opt$output_files

subfolders <- c('Results_topGO')
subfolders <- c(subfolders, 'Results_KEGG')


################# PREPROCESSING INPUT DATA FOR FUNCTIONAL ANALYSIS #############

############# Setting biomaRt query #############

#### Setting biomaRt query
organism_mart <- as.character(biomaRt_organism_info[,"Mart"])
organism_dataset <- as.character(biomaRt_organism_info[,"Dataset"])
organism_host <- as.character(biomaRt_organism_info[,"biomaRt_Host"])

#### Setting biomaRt attributes
attr <- c(opt$biomaRt_filter, as.character(biomaRt_organism_info[,"Attribute_GOs"]))
attr <- c(attr, as.character(biomaRt_organism_info[,"Attribute_entrez"]))


# Getting biomaRt information
# -----------------------------------------------

reference_table <- obtain_info_from_biomaRt(as.character(reference_table[,1]),opt$biomaRt_filter, organism_mart, organism_dataset, organism_host, attr)  

if (opt$save_query == TRUE){
    query <- reference_table
    saveRDS(query, file=file.path("query_results_temp"))
}



# Getting DEG subunits 
common_DEGs_df <- subset(DEG_annot_table, common_DEGs == TRUE)
common_DEGs <- unique(as.character(common_DEGs_df[,1]))

pos_logFC_common_DEGs <- get_specific_condition_dataframe_names(common_DEGs_df, common_DEGs_df["logFC_DESeq2"] > 0)
neg_logFC_common_DEGs <- get_specific_condition_dataframe_names(common_DEGs_df, common_DEGs_df["logFC_DESeq2"] < 0)

union_DEGs <- unique(as.character(DEG_annot_table$Row.names))

pos_logFC_union_DEGs <- get_specific_condition_dataframe_names(DEG_annot_table, DEG_annot_table["logFC_DESeq2"] > 0)
neg_logFC_union_DEGs <- get_specific_condition_dataframe_names(DEG_annot_table, DEG_annot_table["logFC_DESeq2"] < 0)


union_annot_DEGs_df <- get_specific_dataframe_names(reference_table, reference_table[,1], union_DEGs)
common_annot_DEGs_df <- get_specific_dataframe_names(reference_table, reference_table[,1], common_DEGs)



############ PERFORMING FUNCTIONAL ANALYSIS ######################

module_selected <- grepl("G", opt$functional_analysis)

if (module_selected == TRUE){

  module_selected <- grepl("M", opt$GO_graphs)
  if (module_selected==TRUE){
    perform_GSEA_analysis(common_DEGs, union_annot_DEGs_df, "MF", "topGO_common_MF.pdf")
    perform_GSEA_analysis(pos_logFC_common_DEGs, union_annot_DEGs_df, "MF", "topGO_pos_logFC_common_MF.pdf")
    perform_GSEA_analysis(neg_logFC_common_DEGs, union_annot_DEGs_df, "MF", "topGO_neg_logFC_common_MF.pdf")

    perform_GSEA_analysis(as.character(DEG_annot_table$Row.names), reference_table, "MF", "topGO_union_MF.pdf")
    perform_GSEA_analysis(pos_logFC_union_DEGs, reference_table, "MF", "topGO_pos_logFC_union_MF.pdf")
    perform_GSEA_analysis(neg_logFC_union_DEGs, reference_table, "MF", "topGO_neg_logFC_union_MF.pdf")    
  }

  module_selected <- grepl("B", opt$GO_graphs)
  if (module_selected==TRUE){
    perform_GSEA_analysis(common_DEGs, union_annot_DEGs_df, "BP", "topGO_common_BP.pdf")
    perform_GSEA_analysis(pos_logFC_common_DEGs, union_annot_DEGs_df, "BP", "topGO_pos_logFC_common_BP.pdf")
    perform_GSEA_analysis(neg_logFC_common_DEGs, union_annot_DEGs_df, "BP", "topGO_neg_logFC_common_BP.pdf")

    perform_GSEA_analysis(as.character(DEG_annot_table$Row.names), reference_table, "BP", "topGO_union_BP.pdf")
    perform_GSEA_analysis(pos_logFC_union_DEGs, reference_table, "BP", "topGO_pos_logFC_union_BP.pdf")
    perform_GSEA_analysis(neg_logFC_union_DEGs, reference_table, "BP", "topGO_neg_logFC_union_BP.pdf")    
  }

  module_selected <- grepl("C", opt$GO_graphs)
  if (module_selected==TRUE){
    perform_GSEA_analysis(common_DEGs, union_annot_DEGs_df, "CC", "topGO_common_CC.pdf")
    perform_GSEA_analysis(pos_logFC_common_DEGs, union_annot_DEGs_df, "CC", "topGO_pos_logFC_common_CC.pdf")
    perform_GSEA_analysis(neg_logFC_common_DEGs, union_annot_DEGs_df, "CC", "topGO_neg_logFC_common_CC.pdf")

    perform_GSEA_analysis(as.character(DEG_annot_table$Row.names), reference_table, "CC", "topGO_union_CC.pdf")
    perform_GSEA_analysis(pos_logFC_union_DEGs, reference_table, "CC", "topGO_pos_logFC_union_CC.pdf")
    perform_GSEA_analysis(neg_logFC_union_DEGs, reference_table, "CC", "topGO_neg_logFC_union_CC.pdf")    
  }
}


module_selected <- grepl("K", opt$functional_analysis)

if (module_selected == TRUE){
  
  common_unique_entrez <- unique(common_annot_DEGs_df[,"entrezgene"])
  pathway_common_id_EC_df <- find_interesting_pathways(biomaRt_organism_info, common_unique_entrez)     
  visualizing_KEGG_pathways("common_results_KEGG.html", pathway_common_id_EC_df) 

  union_unique_entrez <- unique(union_annot_DEGs_df[,"entrezgene"])
  pathway_union_id_EC_df <- find_interesting_pathways(biomaRt_organism_info, union_unique_entrez)
  visualizing_KEGG_pathways("union_results_KEGG.html", pathway_union_id_EC_df)
}


################## CREATE FINAL ANNOTATION TABLE ##########

write.table(DEG_annot_table, file=file.path(paths$root, "Annotated_table.txt"), quote=F, col.names=NA, sep="\t")

