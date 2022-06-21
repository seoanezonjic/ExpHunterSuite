#' Perform Functional analysis with different 
#' enrichment corpora.
#' @param hunter_results list with DEG analysis results
#' @param model_organism organisms which genes are being studied
#' @param annot_table (OPTIONAL) annotation table to translate gene IDs
#' @param organisms_table (OPTIONAL) configuration table for given organism.
#'  Use see get_organism_table()
#' @param input_gene_id (OPTIONAL) type og gene IDs given. Allowed: ENSEMBL (E)
#' , entrezgene (e), TAIR/Arabidopsis (T), Gene Names (G).
#' @param func_annot_db (OPTIONAL) modules to be enriched
#' @param GO_subont (OPTIONAL) specific GO subontologies to be enriched
#' @param custom (OPTIONAL) list of dataframes with GMT corpus
#' @param analysis_type (OPTIONAL) enrichment type to be performed. Allowed: 
#' ORA (o) and GSEA (g)
#' @param remote (OPTIONAL) remote allowed actions. Are: Biomart query (b) and
#'  KEGG enrichment (k)
#' @param save_query (OPTIONAL) flag to indicate if query must be stored
#' @param pthreshold pvalue threshold
#' @param qthreshold qvalue threshold
#' @param cores cores to be used if parallel features are going to be used.
#'  Default: 1
#' @param task_size number of elements per packages used
#' @param output_files output folder
#' @param fc_colname main logFC colname (into hunter_results dataframe)
#' @return functional result object with enrichments performed
#' @keywords method
#' @export
#' @importFrom magrittr %>%
#' @importFrom clusterProfiler merge_result
#' @examples
#' # Load DE analysis result
#' data(degh_output)
#' annots <- "" # Include here wanted ontologies to be used
#' func_results <- functional_hunter(hunter_results = degh_output,
#' model_organism = "Mouse",func_annot_db = annots)
main_functional_hunter <- function(
    hunter_results,
    model_organism,
    annot_table = NULL,
    organisms_table = get_organism_table(),
    input_gene_id = "ENSEMBL", # Common examples: ENSEMBL, ENTREZID etc.
    #func_annot_db = "gKR",
    #GO_subont = "BMC",
    custom = NULL,
    #analysis_type = "o", #g
    enrich_dbs = c("MF", "Reactome"),
    kegg_data_file = "",
    enrich_methods = "ORA",
    # remote = "",
    annotation_source = "orgdb", # Other option Biomart, to be added
    save_query = FALSE,
    pthreshold = 0.1,
    qthreshold = 0.2,
    cores = 1,
    task_size = 1,
    output_files = "results",
    fc_colname = "mean_logFCs"
    ){

    ############################################################
    ##                                                        ##
    ##                INITIAL SETUP                           ## 
    ##                                                        ##
    ############################################################

    # Get initial parameters
    final_main_params <- as.list(environment(), all=TRUE)
    func_results <- list()
    func_results[["final_main_params"]] <- final_main_params

    # Check organism selected and obtain the organisms table
    if (any(is.null(model_organism), 
            !model_organism %in% rownames(organisms_table))) {
        stop(paste0('Model organism does not appear in organism table. 
            Currently available organisms are:', rownames(organisms_table)))
    } else { 
        current_organism_info <- subset(organisms_table, 
                rownames(organisms_table) %in% model_organism)  
    }

    # Get sample class info - JRP - to remove when updating expression hunter
    case_treat <- ifelse(hunter_results$sample_groups$class == "C", 
        "Control", "Treatment")
    sample_classes <- paste0("* [", case_treat, "] ", 
        hunter_results$sample_groups$name)

    # Load DEGenesHunter results
    DEGH_results <- hunter_results$DE_all_genes
    DEGH_results <- DEGH_results[DEGH_results$genes_tag != "FILTERED_OUT", ]

    # Chapuza to replicate previous 
    func_results$flags <- list(GO_cp    = any(c("MF","BP","CC") %in% enrich_dbs),
                  REACT = "Reactome" %in% enrich_dbs,
                  KEGG = "KEGG" %in% enrich_dbs,
                  GSEA  = "GSEA" %in% enrich_methods,
                  ORA   = "ORA" %in% enrich_methods)
    clusters_flag <- "Cluster_ID" %in% colnames(DEGH_results)
    # JRP to remove - needed for testing.
    if(clusters_flag) func_results$flags$WGCNA <- TRUE
    else func_results$flags$WGCNA <- FALSE
    ##############################################################
    ##                                                          ##
    ##  LOAD COMPLEMENTARY FILES AND ADD ENTREZ AND SYMBOL IDs  ## 
    ##                                                          ##
    ##############################################################

    if(!is.null(annot_table)){
        input_IDs <- translate_from_table(rownames(DEGH_results), annot_table)
    } else {
        input_IDs <- rownames(DEGH_results)
    }

    if(annotation_source == "biomart") {
        stop("biomart functionality for ID translation not implemented")
    } else if (annotation_source == "orgdb") {
        gene_translation_tables <- get_translation_tables_bm(
            input_gene_id=input_gene_id, 
            input_ids=input_IDs, 
            current_organism_info=current_organism_info)
    }

    DEGH_results <- add_translated_gene_ids(DEGH_results, input_IDs, 
        input_gene_id, gene_translation_tables)

    ############################################################
    ##                                                        ##
    ##              GET SIG GENES AND GENE LISTS              ## 
    ##                                                        ##
    ############################################################

    prev_genes <- get_sig_genes(DEGH_results)[["prev_genes"]]
    geneList <- get_gene_lists(DEGH_results, fc_colname)
    
    if(clusters_flag) {

        cl_genes <- get_sig_genes_cl(DEGH_results)
        cl_geneList <- get_gene_lists_cl(DEGH_results, fc_colname)
        cl_geneList <- cl_geneList[! names(cl_geneList) == 0]
    }

    ############################################################
    ##                                                        ##
    ##                     PERFORM ENRICHMENTS                ## 
    ##                                                        ##
    ############################################################

    if("ORA" %in% enrich_methods){
        deg_enr_ora  <- multienricher_ora(all_funsys=enrich_dbs, 
            genes_list=prev_genes, organism_info = current_organism_info, 
            pvalueCutoff = pthreshold, qvalueCutoff = qthreshold, 
            custom_sets=custom, kegg_file = kegg_data_file, symbols_in_plots=TRUE)
        
            deg_enr_ora <- add_term_sim_ora(deg_enr_ora)

        # JRP TO MAKE OUTPUT LIKE PREVIOUS
        func_results$ORA <- deg_enr_ora
        # names(func_results$ORA) <-  gsub("(MF|BP|CC)", "GO_\\1", 
        #                                  names(func_results$ORA)) 
        # names(func_results$ORA) <-  gsub("Reactome", "REACT", 
        #                                  names(func_results$ORA))

        if(clusters_flag){

            clusters_enr_ora <- multienricher_ora(all_funsys=enrich_dbs, 
                genes_list=cl_genes, organism_info = current_organism_info, 
                pvalueCutoff = pthreshold, qvalueCutoff = qthreshold, 
                custom_sets=custom, kegg_file = kegg_data_file, symbols_in_plots=TRUE) 
            
            clusters_enr_ora_compact <- merge_clusters(clusters_enr_ora)
            
            func_results$WGCNA_ORA <- add_term_sim_ora(clusters_enr_ora_compact)
            
            func_results$WGCNA_ORA_expanded <- add_term_sim_ora(clusters_enr_ora)


            # JRP TO MAKE OUTPUT LIKE PREVIOUS
            # names(func_results$WGCNA_ORA_expanded) <-  gsub("(MF|BP|CC)", 
            #   "GO_\\1", names(func_results$WGCNA_ORA_expanded))
            # names(func_results$WGCNA_ORA_expanded) <-  gsub("Reactome", 
            #   "REACT", names(func_results$WGCNA_ORA_expanded)) 
            # names(func_results$WGCNA_ORA) <-  gsub("(MF|BP|CC)", 
            #   "GO_\\1", names(func_results$WGCNA_ORA))
            # names(func_results$WGCNA_ORA) <-  gsub("Reactome", 
            #   "REACT", names(func_results$WGCNA_ORA)) 
        }
    }

    if("GSEA" %in% enrich_methods){
        deg_enr_gsea <- multienricher_gsea(all_funsys=enrich_dbs, 
            genes_list=geneList, 
            organism_info = current_organism_info, kegg_file = kegg_data_file,
            pvalueCutoff = pthreshold, symbols_in_plots=TRUE)

        deg_enr_gsea <- add_term_sim_ora(deg_enr_gsea)
        func_results$GSEA <- deg_enr_gsea
  
        # func_results$GSEA <- deg_enr_gsea[names(deg_enr_gsea) %in% enrich_dbs]        
        # names(func_results$GSEA) <-  gsub("(MF|BP|CC)", "GO_\\1", 
        #                                  names(func_results$GSEA)) 
        # names(func_results$GSEA) <-  gsub("Reactome", "REACT", 
        #                                  names(func_results$GSEA))

      if (clusters_flag) {
            message("GSEA not performed on clusters as non-sensical")
        }
    }

    ############################################################
    ##                                                        ##
    ##                     EXPORT DATA                        ## 
    ##                                                        ##
    ############################################################

    # JRP Silly thing to rename columns to use original report template
    names(DEGH_results)[names(DEGH_results) %in% c("ENTREZID", "SYMBOL")] <- c("entrezgene", "Symbol")

    # Reorder rows by combined FDR
    DEGH_results <- DEGH_results[order(DEGH_results[,"combined_FDR"]),] 
    func_results$DEGH_results_annot <- DEGH_results
    
    return(func_results)
}
