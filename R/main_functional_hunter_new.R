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
    enrich_methods = "ORA",
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


    # Obtain parameters
    func_results <- list()
    func_results[["final_main_params"]] <- as.list(environment(), all=TRUE)

    # Check organism selected and obtain the organisms table
    if (any(is.null(model_organism), 
            !model_organism %in% rownames(organisms_table))) {
        stop(paste0('Model organism does not appear in organism table. Currently available organisms are:', rownames(organisms_table)))
    } else { 
        current_organism_info <- subset(organisms_table, 
                rownames(organisms_table) %in% model_organism)  
    }

    # JRP Sanitize later - should reuse original target
    experiments <- hunter_results$sample_groups
    experiments$class[experiments$class == "C"] <- "Control"
    experiments$class[experiments$class == "T"] <- "Treatment"
    sample_classes <- apply(experiments, 1, function(x){
        paste0("* [", x[1], "] ", x[2])
    })

    ## Load DEGenesHunter results
    ## DEGH_results => DEgenes Hunter results table and modifications
    DEGH_results <- hunter_results$DE_all_genes
    DEGH_results <- DEGH_results[DEGH_results$genes_tag != "FILTERED_OUT", ]

    # JRP only flag necessary
    clusters_flag <- "Cluster_ID" %in% colnames(DEGH_results)

    ##############################################################
    ##                                                          ##
    ##  LOAD COMPLEMENTARY FILES AND ADD ENTREZ AND SYMBOL IDs  ## 
    ##                                                          ##
    ##############################################################

    #Load annotation gene table for translation ID
    if(!is.null(annot_table)){
          DEGH_results$input_IDs <- translate_id(rownames(DEGH_results), 
                                                 annot_table)
    } else {
        DEGH_results$input_IDs <- rownames(DEGH_results)
    }

    if(annotation_source == "biomart") {

        stop("biomart functionality for ID translation remains to be implemented")

    } else if (annotation_source == "orgdb") {

        gene_translation_tables <- get_translation_tables_bm(input_gene_id=input_gene_id, DEGH_results=DEGH_results, current_organism_info=current_organism_info)
    
    }

    input_to_entrezgene <- gene_translation_tables[["input_to_entrezgene"]]    
    input_to_symbol <- gene_translation_tables[["input_to_symbol"]]    

    DEGH_results <- data.frame(ENTREZID = input_to_entrezgene[
      match(DEGH_results$input_IDs, input_to_entrezgene[[input_gene_id]]), "ENTREZID"], DEGH_results)
    
    if(!is.null(input_to_symbol)) {
      DEGH_results$SYMBOL <- input_to_symbol[
        match(DEGH_results$input_IDs, input_to_symbol[[input_gene_id]]), "SYMBOL"]
    }
    save(list = ls(all.names = TRUE), file = "~/environment_new.RData")

    ############################################################
    ##                                                        ##
    ##              PREPARE DATA FOR ENRICHMENTS              ## 
    ##                                                        ##
    ############################################################

    ## topGO & ORA
    likely_degs_entrez <- DEGH_results[DEGH_results$genes_tag == "PREVALENT_DEG" &
                                   !is.na(DEGH_results$ENTREZID), "ENTREZID"]

                                   print("entrez genes")
                                   print(likely_degs_entrez)
    # Verbose point
    message(paste("IDs used in ORA:",length(likely_degs_entrez)))
    ## TODO => ESTARIA BIEN REFLEJAR ESTA INFORMACION EN EL REPORT
    "%>%" <- magrittr::"%>%"
    union_DEGs_df <- subset(DEGH_results, genes_tag %in% c("POSSIBLE_DEG",
                             "PREVALENT_DEG"))
    union_DEGs <- union_DEGs_df[!is.na(union_DEGs_df$input_IDs), 
                                     "input_IDs"] %>% unique
    union_annot_DEGs_df <- subset(input_to_entrezgene, 
                        input_to_entrezgene[,1] %in% union_DEGs)

    ###############################
    ## GSEA
    # Calculate geneList
    #geneList: A NAMED VECTOR OF logFC; USED IN GSEA AND PLOTS
    geneList <- DEGH_results[!is.na(DEGH_results$ENTREZID),  fc_colname]
    names(geneList) <- DEGH_results[!is.na(DEGH_results$ENTREZID), 
                                    "ENTREZID"]
    # Sort FC
    geneList <- sort(geneList, decreasing = TRUE)

    #############################################
    ### WGCNA MODULES
    if(clusters_flag) {
        cls <- unique(DEGH_results$Cluster_ID)
        # DELETE GREY MODULE
        if (any(c(0,"grey") %in% cls)) {
            cls <- cls[!cls %in% c(0,"grey")]
        } else {
            warning("Module Zero/Grey not found")
        }

        clgenes <- lapply(cls,function(cl) { # Find
            unique(DEGH_results$entrezgene[which(DEGH_results$Cluster_ID == 
                                                  cl)])
        }) 
        names(clgenes) <- cls

        if ("GSEA" %in% enrich_methods) {
#            enrichments_GSEA_expanded <- list()
            cl_gene_fc <- lapply(cls, function(cl) {
                cl_results <- DEGH_results[!is.na(DEGH_results$entrezgene) & 
                                DEGH_results$Cluster_ID == cl,]            
                gene_fc <- as.vector(cl_results[,fc_colname])
                names(gene_fc) <- cl_results$entrezgene
                gene_fc <- sort(gene_fc, decreasing = TRUE)

                return(gene_fc)
            })
            names(cl_gene_fc) <- as.character(cls)
        }
    }

    message("Data prepared for enrichments")



    ############################################################
    ##                                                        ##
    ##                     PERFORM ENRICHMENTS                ## 
    ##                                                        ##
    ############################################################



    if("ORA" %in% enrich_methods){

        deg_enr_ora  <- multienricher_ora(all_funsys=enrich_dbs, genes=list(likely_degs_entrez), organism_info = current_organism_info, 
            pvalueCutoff = pthreshold, qvalueCutoff = qthreshold)
        # JRP Add this semsim stuff, important for some reason, to check later.
        for(funsys in names(deg_enr_ora)) {
            print(funsys)
            enrichment <- deg_enr_ora[[funsys]][[1]]
            if(!is.null(enrichment)) {
                if(nrow(enrichment) > 0) {
                    deg_enr_ora[[funsys]] <- catched_pairwise_termsim(enrichment)
                }
            }
        }

        func_results$ORA <- deg_enr_ora 

        if(clusters_flag){
            clusters_enr_ora <- multienricher_ora(all_funsys=enrich_dbs, genes=clgenes, organism_info = current_organism_info, 
            pvalueCutoff = pthreshold, qvalueCutoff = qthreshold)
            # JRP Add this semsim stuff, important for some reason, to check and clean up later.
            func_results$ORA$clusters <- clusters_enr_ora
            for(funsys in names(clusters_enr_ora)) {
                print(funsys)
                enrichment <- clusters_enr_ora[[funsys]]
                cat("length of enrichment is:", length(enrichment), "\n")
                for(i in 1:length(enrichment)) {
                    if(!is.null(enrichment[[i]])) {
                        #if(nrow(enrichment[[i]]) > 0) {
                            func_results$ORA$clusters[[funsys]][[i]] <- catched_pairwise_termsim(enrichment[[i]])
                        #}
                    }
                }
            }
        }
    }
    if("GSEA" %in% enrich_methods){
        deg_enr_gsea <- multienricher_gsea(all_funsys=enrich_dbs, genes_list=list(geneList), 
            organism_info = current_organism_info, 
            pvalueCutoff = pthreshold)
        
        print("check names")
        print(names(deg_enr_gsea))
        print("enrich_dbs")
        print(enrich_dbs)

        for(funsys in names(deg_enr_gsea)) {
            enrichment <- deg_enr_gsea[[funsys]][[1]]
            func_results[["GSEA"]][[funsys]] <- enrichment
        }   
      if (clusters_flag) {
            clusters_enr_gsea <- multienricher_gsea(all_funsys=enrich_dbs, genes_list = cl_gene_fc, 
                organism_info = current_organism_info, 
                pvalueCutoff = pthreshold)
            func_results[["GSEA"]][["clusters"]] <- clusters_enr_gsea
            # if(length(clusters_enr_gsea$WGCNA$GSEA)){
            #     aux <- clusters_enr_gsea$WGCNA
            #     func_results$WGCNA_GSEA <- aux$GSEA
            #     func_results$WGCNA_GSEA_expanded <- aux$GSEA_expanded
            # }
        }
    }

    #############################################
    ### CUSTOM ENRICHMENT
    if (!is.null(custom)) {
        message("Performing CUSTOM enrichments")
        ###########
        ### CUSTOM ENRICHMENTS
        custom_enrichments <- enrich_all_customs(custom_sets = custom, 
                                                 p_val_threshold = pthreshold, 
                                                 genes = likely_degs_entrez)
        func_results$CUSTOM <- custom_enrichments
        if (flags$WGCNA){
            custom_cls_ORA_expanded <- lapply(custom, function(gmt){
                enrich_clusters_with_gmt(custom_set = gmt, 
                                        genes_in_modules = clgenes, 
                                        p_val_threshold = pthreshold, 
                                        cores = cores,
                                        task_size= task_size)
            })
        }
        message("CUSTOM enrichments finished")
    } else {
        custom_enrichments <- NULL
    }

    ############################################################
    ##                                                        ##
    ##                     EXPORT DATA                        ## 
    ##                                                        ##
    ############################################################
    # Return here to sort out WGCNA stuff
    if (clusters_flag) {
        if (!is.null(custom)){
            custom_cls_ORA <- lapply(custom_cls_ORA_expanded, 
                                     clusterProfiler::merge_result)
            custom_cls_ORA <- lapply(custom_cls_ORA, function(res){
                if(nrow(res) > 0) res <- catched_pairwise_termsim(res)
                return(res)
            })
            func_results$WGCNA_CUSTOM <- custom_cls_ORA
            func_results$WGCNA_CUSTOM_expanded <-custom_cls_ORA_expanded
        }
    }


    # Reorder rows by combined FDR
    DEGH_results <- DEGH_results[order(DEGH_results[,"combined_FDR"]),] 
    func_results$DEGH_results_annot <- DEGH_results


    #func_results$flags <- flags
    return(func_results)
}
