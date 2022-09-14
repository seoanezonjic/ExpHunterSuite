#' Perform enrichment analysis on multiple lists of genes
#' @param input_file input file path. File must contain 1 cluster per line clusterA ..tab.. gene1,gene2
#' @param gene_attribute_file file to map genes to attributes e.g. fold change
#' @param org_db org_db object
#' @param gene_keytype ID used in the input file
#' @param temp_file enrichments will be read from this file if it exists, written to it if it doesnt
#' @param force if true, to run enrichments even if the temp_file exists
#' @param all_funsys list of GOs/dbs to use for enrichment analysis
#' @param task_size number of elements per packages used
#' @param workers cores to be used if parallel features are going to be used
#' @param current_organism_info organism info from text table,
#' @param pvalcutoff adjusted p-value threshold to use for clusterProfiler,
#' @param qvalcutoff q-value threshold to use for clusterProfiler,
#' @param all_custom_gmt list of custom annotation sets
#' @param kegg_data_file if KEGG included in enrich_dbs path needs
#' @param simplify simplify the merged enrichments
#' @param clean_parentals clean parental terms in merged enrichment
#' @export
main_clusters_to_enrichment <- function(
  input_file=NULL,
  gene_attribute_file = NULL,
  org_db = NULL,
  gene_keytype ="ENTREZID",
  temp_file = "enr_tmp.RData",
  force = FALSE,
  all_funsys = c("BP","MF","CC"),
  task_size = 1,
  current_organism_info = NULL,
  workers = 1,
  pvalcutoff = 0.1,
  qvalcutoff = 0.2,
  all_custom_gmt = NULL,
  kegg_data_file = NULL,
  simplify = FALSE,
  clean_parentals = FALSE){

  gene_attributes <- NULL
  gene_attribute_name <- NULL
  if (!is.null(gene_attribute_file)){
    gene_attributes_object <- obtain_gene_attributes(gene_attribute_file, org_db, gene_keytype)
    gene_attributes <- gene_attributes_object[["gene_attributes"]]
    gene_attribute_name <-  gene_attributes_object[["gene_attribute_name"]]
  }

  if (!file.exists(temp_file) || force) {
    cluster_genes <- read.table(input_file, header=FALSE, sep="\t")
    cluster_genes_list <- strsplit(cluster_genes[,2], ",")
    names(cluster_genes_list) <- cluster_genes[,1]
    if(gene_keytype != "ENTREZID") {
      cluster_genes_list <- lapply(cluster_genes_list, function(x){
          translate_ids_orgdb(ids=x, input_id = gene_keytype,
            org_db = org_db, just_output_ids=TRUE)}) 
    }
    names(cluster_genes_list) <- cluster_genes[,1]

    enrichments_ORA <- multienricher_ora(all_funsys =  all_funsys, 
                                  genes_list =  cluster_genes_list, 
                                  task_size = task_size, 
                                  organism_info= current_organism_info,
                                  workers = workers, 
                                  pvalueCutoff =  pvalcutoff, 
                                  qvalueCutoff = qvalcutoff,
                                  custom_sets = all_custom_gmt,
                                  kegg_file = kegg_data_file)

    enrichments_ORA_merged <- process_cp_list(enrichments_ORA, simplify_results = simplify, 
      clean_parentals = clean_parentals)
    save(enrichments_ORA,enrichments_ORA_merged, file = temp_file)

  } else {
    load(temp_file)
  }
  return(list(enrichments_ORA=enrichments_ORA, 
    enrichments_ORA_merged=enrichments_ORA_merged, 
    gene_attributes=gene_attributes, 
    gene_attribute_name=gene_attribute_name))
}