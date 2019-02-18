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
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(biomaRt)) 
suppressPackageStartupMessages(library(topGO))
suppressPackageStartupMessages(library(KEGGREST))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(knitr))

# Load custom libraries
source(file.path(main_path_script, 'lib', 'general_functions.R'))
source(file.path(main_path_script, 'lib', 'functional_analysis_library.R'))

biomaRt_query_info <- read.table(file.path(main_path_script, "lib/biomaRt_organism_table.txt"), header=T, row.names=1, sep="\t")

#############################################
### MAIN 
#############################################

### Input/Output (I/O)
#############################################

# Parse command line
#------------------------------------------------

option_list <- list(
  make_option(c("-i", "--input_hunter_file"), type="character", default="hunter_DE_results/Common_results/hunter_results_table.txt",
    help="DEgenes Hunter's differential expression analysis output file"), 
  make_option(c("-c", "--countdata_file"), type="character", default="hunter_DE_results/filtered_count_data.txt",
    help="Filtered count data file"), 
  make_option(c("-a", "--annot_file"), type="character", default=NULL,
    help="Two column file with annotations for functional analysis of a non-model organism. First column must be a gene ensembl id or a refseq id from a model organism. This id must be a orthologue of the gene id of the second column, that is the custom id from a non-model organism (whose functional analysis is desired)"),
  make_option(c("-L", "--List_organisms"), action="store_true", type="logical", default=FALSE, 
    help="List all organisms provided"),
  make_option(c("-m", "--model_organism"), type="character", default=NULL,
    help="Ortologue Species. Default=%default"),
  make_option(c("-t", "--biomaRt_filter"), type="character", default="E",
    help="IDtype, 'E' for 'ensembl_gene_id' or 'R' for 'refseq_peptide', Default: E"),       	
  make_option(c("-f", "--functional_analysis"), type="character", default=c("GK"),
    help="Selection of type of functional analyses to be performed (G = GO, K = KEGG). By default the following modules Default=%default are performed"),
  make_option(c("-G", "--GO_graphs"), type="character", default=c("M"),
    help="Modules to able/disable (M = Molecular Function, B = Biological Process, C = Celular Components). By default Default=%default GO cathegory is performed"),
  make_option(c("-K", "--Kegg_organism"), type="character", default=NULL, 
    help="Indicate organism to look for in the Kegg database for doing the path enrichment"),
  make_option(c("-q", "--save_query"), type="logical", default=TRUE,
    help="Save biomaRt query. By default the biomaRt query is saved"), 
  make_option(c("-o", "--output_files"), type="character", default="results",
    help="Output path. Default=%default")
)
opt <- parse_args(OptionParser(option_list=option_list))


###### Parsing input data ######
if(opt$List_organisms == TRUE){
  print(as.character(rownames(biomaRt_query_info)))
  stop('program ends')
} else if (is.null(opt$model_organism)){
  stop('No model organism indicated. Please indicate the
     model organism using parameter -m. Use -L to list all model organisms available')
}

biomaRt_organism_info <- get_specific_dataframe_names(biomaRt_query_info, rownames(biomaRt_query_info), opt$model_organism)
KEGG_enrich_selected <- grepl("K", opt$functional_analysis)

if ((biomaRt_organism_info["KeggCode"] == "NotInKEGG") & (is.null(opt$Kegg_organism)) & (KEGG_enrich_selected==TRUE)){
  stop('This organism is not yet available in the KEGG database. You can select an alternative organism to perform
    the path enrichment analysis using parameter -K. Use -L to list all model organisms available.')
}

if (!is.null(opt$annot_file)){
  annot <- read.table(opt$annot_file, header=F, row.names=NULL, sep="\t")
  reference_table <- annot  
} else {
  if(!file.exists(opt$countdata_file)) stop('Count file not exists, check the PATH given to the -c flag')
  reference_table <- read.table(opt$countdata_file, header=T, row.names=NULL, sep="\t")
}

if (opt$biomaRt_filter == "E"){
  opt$biomaRt_filter <- 'ensembl_gene_id'
  #opt$biomaRt_filter <- as.character(biomaRt_organism_info[1,7])
}else if (opt$biomaRt_filter == "R"){
  opt$biomaRt_filter <- 'refseq_peptide'
  #opt$biomaRt_filter <- as.character(biomaRt_organism_info[1,8])
}

DEG_annot_table <- read.table(opt$input_hunter_file, header=T, row.names=1, sep="\t")

############ CREATE FOLDERS #########3
paths <- list()
dir.create(opt$output_files)
paths$root <-opt$output_files

############# BIOMART QUERY  #############
organism_mart <- as.character(biomaRt_organism_info[,"Mart"])
organism_dataset <- as.character(biomaRt_organism_info[,"Dataset"])
organism_host <- as.character(biomaRt_organism_info[,"biomaRt_Host"])
attr <- c(
  opt$biomaRt_filter, 
  as.character(biomaRt_organism_info[,"Attribute_GOs"]),
  as.character(biomaRt_organism_info[,"Attribute_entrez"])
)

reference_table <- obtain_info_from_biomaRt(as.character(reference_table[,1]),opt$biomaRt_filter, organism_mart, organism_dataset, organism_host, attr)  
if (opt$save_query == TRUE) saveRDS(reference_table, file=file.path("query_results_temp"))

################# PREPROCESSING INPUT DATA FOR FUNCTIONAL ANALYSIS #############
common_DEGs_df <- subset(DEG_annot_table, genes_tag=="PREVALENT_DEG")
if (nrow(common_DEGs_df) == 0) stop('Not detected genes tagged as PREVALENT_DEG. Likely some modules of degenes_Hunter.R failed and not generated output. Please, revise the input data')
common_DEGs <- rownames(common_DEGs_df)
common_annot_DEGs_df <- get_specific_dataframe_names(reference_table, reference_table[,1], common_DEGs)
common_unique_entrez <- unique(common_annot_DEGs_df[,"entrezgene"])

union_DEGs_df <- subset(DEG_annot_table, genes_tag=="POSSIBLE_DEG")
union_DEGs_df <- rbind(common_DEGs_df, union_DEGs_df)
union_DEGs <- rownames(union_DEGs_df)

union_annot_DEGs_df <- get_specific_dataframe_names(reference_table, reference_table[,1], union_DEGs)

############ PERFORMING FUNCTIONAL ANALYSIS ######################

module_selected <- grepl("G", opt$functional_analysis)

if (module_selected == TRUE){

    pos_logFC_common_DEGs <- rownames(subset(common_DEGs_df, common_DEGs_df["logFC_DESeq2"] > 0))
    neg_logFC_common_DEGs <- rownames(subset(common_DEGs_df, common_DEGs_df["logFC_DESeq2"] < 0))

    pos_logFC_union_DEGs <- rownames(subset(union_DEGs_df, union_DEGs_df["logFC_DESeq2"] > 0))
    neg_logFC_union_DEGs <- rownames(subset(union_DEGs_df, union_DEGs_df["logFC_DESeq2"] < 0))

  module_selected <- grepl("M", opt$GO_graphs)
  if (module_selected==TRUE){

    perform_GSEA_analysis(common_DEGs, union_annot_DEGs_df, "MF", "GOgraph_preval_MF.pdf")
    perform_GSEA_analysis(pos_logFC_common_DEGs, union_annot_DEGs_df, "MF", "GOgraph_preval_overex_MF.pdf")
    perform_GSEA_analysis(neg_logFC_common_DEGs, union_annot_DEGs_df, "MF", "GOgraph_preval_underex_MF.pdf")

    perform_GSEA_analysis(union_DEGs, reference_table, "MF", "GOgraph_allpos_MF.pdf")
    perform_GSEA_analysis(pos_logFC_union_DEGs, reference_table, "MF", "GOgraph_allpos_overex_MF.pdf")
    perform_GSEA_analysis(neg_logFC_union_DEGs, reference_table, "MF", "GOgraph_allpos_underex_MF.pdf")    
  }

  module_selected <- grepl("B", opt$GO_graphs)
  if (module_selected==TRUE){
    perform_GSEA_analysis(common_DEGs, union_annot_DEGs_df, "BP", "GOgraph_preval_BP.pdf")
    perform_GSEA_analysis(pos_logFC_common_DEGs, union_annot_DEGs_df, "BP", "GOgraph_preval_overex_BP.pdf")
    perform_GSEA_analysis(neg_logFC_common_DEGs, union_annot_DEGs_df, "BP", "GOgraph_preval_underex_BP.pdf")

    perform_GSEA_analysis(union_DEGs, reference_table, "BP", "GOgraph_allpos_BP.pdf")
    perform_GSEA_analysis(pos_logFC_union_DEGs, reference_table, "BP", "GOgraph_allpos_overex_BP.pdf")
    perform_GSEA_analysis(neg_logFC_union_DEGs, reference_table, "BP", "GOgraph_allpos_underex_BP.pdf")    
  }

  module_selected <- grepl("C", opt$GO_graphs)
  if (module_selected==TRUE){
    perform_GSEA_analysis(common_DEGs, union_annot_DEGs_df, "CC", "GOgraph_preval_CC.pdf")
    perform_GSEA_analysis(pos_logFC_common_DEGs, union_annot_DEGs_df, "CC", "GOgraph_preval_overex_CC.pdf")
    perform_GSEA_analysis(neg_logFC_common_DEGs, union_annot_DEGs_df, "CC", "GOgraph_preval_underex_CC.pdf")

    perform_GSEA_analysis(union_DEGs, reference_table, "CC", "GOgraph_allpos_CC.pdf")
    perform_GSEA_analysis(pos_logFC_union_DEGs, reference_table, "CC", "GOgraph_allpos_overex_CC.pdf")
    perform_GSEA_analysis(neg_logFC_union_DEGs, reference_table, "CC", "GOgraph_allpos_underex_CC.pdf")    
  }
}

module_selected <- grepl("K", opt$functional_analysis)

if (module_selected == TRUE){
  if(!(is.null(opt$Kegg_organism))){ 
    Kegg_organism_line <- get_specific_dataframe_names(biomaRt_query_info, rownames(biomaRt_query_info), opt$Kegg_organism)
    biomaRt_KEGG_organism_info["KeggCode"] <- Kegg_organism_line["KeggCode"]
    KEGG_table <- obtain_info_from_biomaRt(as.character(reference_table[,1]),opt$biomaRt_filter, organism_mart, organism_dataset, organism_host, "entrezgene")
  
    common_annot_DEGs_df <- get_specific_dataframe_names(KEGG_table, KEGG[,1], common_DEGs)
  
    common_unique_entrez <- unique(common_annot_DEGs_df[,"entrezgene"])
    pathway_common_id_EC_df <- find_interesting_pathways(biomaRt_KEGG_organism_info, common_unique_entrez)     
    visualizing_KEGG_pathways("KEGG_paths.html", pathway_common_id_EC_df) 
  } else { 
    pathway_common_id_EC_df <- find_interesting_pathways(biomaRt_organism_info, common_unique_entrez)     
    visualizing_KEGG_pathways("KEGG_paths.html", pathway_common_id_EC_df) 
  }
}

################## CREATE FINAL ANNOTATION TABLE ##########
write.table(DEG_annot_table, file=file.path(paths$root, "Annotated_table.txt"), quote=F, col.names=NA, sep="\t")
write.table(reference_table, file=file.path(paths$root, "entrez_Gos.txt"), quote=F, col.names=NA, sep="\t")

#### functional
generate_FA_report()
