#' Perform Functional analysis with different 
#' enrichment corpora.
#' @param hunter_results output of main_degenes_hunter - list with DEG analysis results
#' @param model_organism organisms which genes are being studied
#' @param annot_table (OPTIONAL) annotation table to translate gene IDs
#' @param organisms_table (OPTIONAL) configuration table for given organism.
#'  Use see get_organism_table()
#' @param input_gene_id (OPTIONAL) type og gene IDs given. Allowed: ENSEMBL (E)
#' , entrezgene (e), TAIR/Arabidopsis (T), Gene Names (G).
#' @param custom (OPTIONAL) list of dataframes with GMT corpus
#' @param enrich_dbs annotation dbs to be used: BP,MF,CC,Reactome,KEGG
#' - MF and Reactome default
#' @param kegg_data_file if KEGG included in enrich_dbs path needs
#' @param enrich_methods enrichment methods to be used: ORA or GSEA
#' @param annotation_source annote with orgdbs or Biomart? Currently on orgdbs
#' @param pthreshold pvalue threshold
#' @param qthreshold qvalue threshold
#' @param cores cores to be used if parallel features are going to be used.
#'  Default: 1
#' @param task_size number of elements per packages used
#' @param output_files output folder
#' @param fc_colname main logFC colname (into hunter_results dataframe)
#' @param universe whether to use all genes as the background, or expressed only
#' @param simplify simplify the merged enrichments
#' @param clean_parentals clean parental terms in merged enrichment
#' @param top_categories numbers of categories from each cluster to use for merge
#' @param sim_thr value to use when combining similar categories in summary mode
#' @param summary_common_name 'significant' to use the most significant term to label each summarized group
#' @param clusters_flag execute clusters enrichments
#' @return functional result object with enrichments performed
#' @keywords method
#' @export
#' @importFrom magrittr %>%
#' @importFrom clusterProfiler merge_result
#' @examples
#' # Load DE analysis result
#' data(degh_output)
#' func_results <- main_functional_hunter(hunter_results = degh_output,
#' model_organism = "Mouse", enrich_dbs = "MF")
main_functional_hunter <- function(
    hunter_results,
    model_organism,
    annot_table = NULL,
    organisms_table = get_organism_table(),
    input_gene_id = "ENSEMBL", # Common examples: ENSEMBL, ENTREZID etc.
    custom = NULL,
    enrich_dbs = c("MF", "Reactome"),
    kegg_data_file = "",
    enrich_methods = "ORA",
    annotation_source = "orgdb", # Other option Biomart, to be added
    pthreshold = 0.1,
    qthreshold = 0.2,
    cores = 1,
    task_size = 1,
    output_files = "results",
    fc_colname = "mean_logFCs",
    universe = NULL,
    clean_parentals = FALSE,
    simplify = FALSE,
    top_categories = 50,
    sim_thr = NULL,
    summary_common_name = "ancestor",
    clusters_flag = FALSE
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

    clusters_flag <- "Cluster_ID" %in% colnames(DEGH_results) && clusters_flag
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
        gene_translation_tables <- get_translation_tables_orgdb(
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

    if(! is.null(universe)) {
        if(universe == "expressed") {
            universe <- DEGH_results$ENTREZID[DEGH_results$genes_tag != "FILTERED_OUT"]
        }
    }

    if (!is.null(hunter_results$pca_data)) {
    # prepare PCA data for topGO

        pca_eigenvectors <-  lapply(hunter_results$pca_data, 
                                    get_and_parse_pca_eigenvectors, 
                                    cor_pval = 1, 
                                    eig_abs_cutoff = NULL)
        pca_clusters <- lapply(hunter_results$pca_data, get_and_parse_clusters)

    }

    ############################################################
    ##                                                        ##
    ##                     PERFORM ENRICHMENTS                ## 
    ##                                                        ##
    ############################################################

    if("GSEA" %in% enrich_methods){
        deg_enr_gsea <- multienricher_gsea(all_funsys=enrich_dbs, 
            genes_list=geneList, 
            organism_info = current_organism_info, kegg_file = kegg_data_file,
            pvalueCutoff = pthreshold)

        deg_enr_gsea <- add_term_sim_ora(deg_enr_gsea)
        func_results$GSEA <- deg_enr_gsea
        if (clusters_flag) message("GSEA not performed on clusters as non-sensical")
    }

    if("ORA" %in% enrich_methods){

        deg_enr_ora  <- multienricher_ora(all_funsys = enrich_dbs, 
            genes_list = prev_genes, organism_info = current_organism_info, 
            pvalueCutoff = pthreshold, qvalueCutoff = qthreshold, 
            custom_sets = custom, kegg_file = kegg_data_file,
            return_all = TRUE, universe = universe)
        func_results$ORA  <- add_term_sim_ora(deg_enr_ora)

        if(clusters_flag){
            clusters_enr_ora <- multienricher_ora(all_funsys=enrich_dbs, 
                genes_list=cl_genes, organism_info = current_organism_info, 
                pvalueCutoff = pthreshold, qvalueCutoff = qthreshold, 
                custom_sets=custom, kegg_file = kegg_data_file, workers=cores, 
                task_size=task_size, return_all = TRUE, universe=universe)

            clusters_enr_ora_merged<- process_cp_list(clusters_enr_ora, simplify, clean_parentals)
            func_results$WGCNA_ORA <- clusters_enr_ora_merged
            func_results$WGCNA_ORA_expanded <- add_term_sim_ora(clusters_enr_ora)

            # TO improve and add to reports
            if(! is.null(sim_thr))
              func_results$summarized_ora <- summarize_merged_ora(clusters_enr_ora_merged, 
                  sim_thr, summary_common_name, pthreshold)
        }
    }

    if ("topGO" %in% enrich_methods){
        func_results$topGO <- multienricher_topGO(all_funsys=enrich_dbs, 
                                                   universe=universe,
                                                   genes_list=prev_genes,
                                                   organism_info=current_organism_info,
                                                   p_value_cutoff = pthreshold)
    }

    if (!is.null(hunter_results$pca_data)) {

        func_results$PCA_enrichments <- lapply(pca_eigenvectors, 
                                                multienricher_topGO, 
                                                all_funsys = enrich_dbs,
                                                organism_info = current_organism_info,
                                                p_value_cutoff = pthreshold,
                                                algorithm = "classic",
                                                statistic = "t",
                                                gene_id = "ensembl", 
                                                scoreOrder = "decreasing",
                                                workers = cores, 
                                                task_size = task_size)
        func_results$PCA_clusters_results <- lapply(pca_clusters, 
                                                multienricher_topGO, 
                                                all_funsys = enrich_dbs,
                                                organism_info = current_organism_info,
                                                p_value_cutoff = pthreshold,
                                                algorithm = "classic",
                                                statistic = "t",
                                                gene_id = "ensembl", 
                                                scoreOrder = "decreasing",
                                                workers = cores, 
                                                task_size = task_size)
    }
    ############################################################
    ##                                                        ##
    ##                     EXPORT DATA                        ## 
    ##                                                        ##
    ############################################################

    # Reorder rows by combined FDR
    DEGH_results <- DEGH_results[order(DEGH_results[,"combined_FDR"]),] 
    func_results$DEGH_results_annot <- DEGH_results
    func_results$pca_data <- hunter_results$pca_data
    return(func_results)
}
