#' Write Main DEgenes Hunter report
#'
#' This function allows you to report the DEG analysis.
#' @export
#' @examples
#' write_expression_report()
write_expression_report <- function(exp_results, output_files,template_folder = file.path(find.package('DEgenesHunter'), 'templates'), opt){
    DEG_pack_columns <- exp_results[['DEG_pack_columns']] 
    all_counts_for_plotting <- exp_results[['all_counts_for_plotting']] 
    all_FDR_names <- exp_results[['all_FDR_names']]
    all_LFC_names <- exp_results[['all_LFC_names']] 
    all_pvalue_names <- exp_results[['all_pvalue_names']]
    final_pvalue_names <- exp_results[['final_pvalue_names']]
    final_logFC_names <- exp_results[['final_logFC_names']]
    final_FDR_names <- exp_results[['final_FDR_names']]
    package_objects <- exp_results[['package_objects']]
    results_WGCNA <- exp_results[['WGCNA_all']]
    index_control_cols <- exp_results[['index_control_cols']] 
    index_treatmn_cols <- exp_results[['index_treatmn_cols']] 
    raw_filter <- exp_results[['raw_filter']] 
    design_vector <- exp_results[['design_vector']]
    all_data_normalized <- exp_results[['all_data_normalized']] 
    replicatesC <- exp_results[['replicatesC']] 
    replicatesT <- exp_results[['replicatesT']] 
    DE_all_genes <- exp_results[['DE_all_genes']] 

    outf <- file.path(normalizePath(output_files),"DEG_report.html")
    rmarkdown::render(file.path(template_folder, 'main_report.Rmd'), 
                      output_file = outf, intermediates_dir = output_files)	
}


#' Write Main DEgenes Hunter functional report
#'
#' This function allows you to report the Functional analysis.
#' @param
#' @param hunter_results DEG analysis results
#' @param func_results functional results
#' @param output_files output folder
#' @param template_folder (OPTIONAL) RMD templates folder
#' @param cores (OPTIONAL) cores for parallel features
#' @param model_organism model organism which genes are being studied
#' @param organisms_table configuration table for given organism. Use see get_organism_table()
#' @param fc_colname main logFC colname (into hunter_results dataframe)
#' @export
#' @examples
write_functional_report <- function(hunter_results, func_results, output_files, template_folder = file.path(find.package('DEgenesHunter'), 'templates'), cores = 1, model_organism, organisms_table, fc_colname){
    # TODO: update names into Rmd files instead this
    ############################################################
    ##               CREATE NECESSARY VARIABLES               ##
    ############################################################
    DEGH_results <- func_results$DEGH_results_annot
    # -
    experiments <- hunter_results$sample_groups
    sample_classes <- apply(experiments, 1, function(x) paste0("* [", x[1], "] ", x[2]))
    # -
    norm_counts <- hunter_results[["all_data_normalized"]][["DESeq2"]]
    scaled_counts <- scale_data_matrix(data_matrix = norm_counts, transpose = TRUE)
    scaled_counts_table <- as.data.frame(as.table(scaled_counts))
    colnames(scaled_counts_table) <- c("Gene","Sample","Count")
    # -
    flags <- func_results$flags
    if(flags$WGCNA){
        cls  <- unique(DEGH_results$Cluster_ID)
        cl_eigvalues <- as.matrix(hunter_results$WGCNA_all$plot_objects$trait_and_module[,grepl("^ME",colnames(hunter_results$WGCNA_all$plot_objects$trait_and_module))])
        cl_eigvalues <- as.data.frame(as.table(cl_eigvalues),stringsAsFactors = FALSE)
        colnames(cl_eigvalues) <- c("Sample","Cluster_ID","Count") 
        cl_eigvalues_gnorm <- cl_eigvalues
        cl_eigvalues_gnorm$Count <- (cl_eigvalues_gnorm$Count + 1) / 2 
        wgcna_pval_cl_trait <- as.matrix(hunter_results$WGCNA_all$package_objects$module_trait_cor_p)
        wgcna_corr_cl_trait <- as.matrix(hunter_results$WGCNA_all$package_objects$module_trait_cor)
        wgcna_count_sample_trait <- as.matrix(hunter_results$WGCNA_all$plot_objects$trait_and_module[,!grepl("^ME",colnames(hunter_results$WGCNA_all$plot_objects$trait_and_module))])
        wgcna_count_sample_trait <- scale_data_matrix(wgcna_count_sample_trait)        
    }
    #-
    if(any(grepl("WGCNA_ORA",names(func_results)))){
        enrichments_ORA <- func_results$WGCNA_ORA
        enrichments_ORA_expanded <- func_results$WGCNA_ORA_expanded
    }
    if(any(grepl("WGCNA_GSEA",names(func_results)))){
        enrichments_GSEA <- func_results$WGCNA_GSEA
        enrichments_GSEA_expanded <- func_results$WGCNA_GSEA_expanded
    }
    if(any(grepl("WGCNA_CUSTOM",names(func_results)))){
        custom_cls_ORA <- func_results$WGCNA_CUSTOM
        custom_cls_ORA_expanded <- func_results$WGCNA_CUSTOM_expanded
    }
    #-
    current_organism_info <- subset(organisms_table, rownames(organisms_table) %in% model_organism)  
    geneList <- func_results$DEGH_results_annot[!is.na(func_results$DEGH_results_annot$entrezgene), fc_colname]
    names(geneList) <- func_results$DEGH_results_annot[!is.na(func_results$DEGH_results_annot$entrezgene), "entrezgene"]
    geneList <- sort(geneList, decreasing = TRUE)
    # -
    custom_enrichments <- func_results$CUSTOM
    if("GO_ORA" %in% names(func_results)) enrich_go <- func_results$GO_ORA
    if("GO_GSEA" %in% names(func_results)) enrich_go_gsea <- func_results$GO_GSEA
    if("KEGG_ORA" %in% names(func_results)) enrich_ora <- func_results$KEGG_ORA
    if("KEGG_GSEA" %in% names(func_results)) enrich_gsea <- func_results$KEGG_GSEA
    if("REACT_ORA" %in% names(func_results)) enrich_react <- func_results$REACT_ORA
    if("REACT_GSEA" %in% names(func_results)) enrich_react_gsea <- func_results$REACT_GSEA
    # TODO: topGO is not bein loaded because files are not been search
    ############################################################
    ##                GENERATE CLUSTER REPORTS                ##
    ############################################################
    results_path <- normalizePath(output_files)

    if (any(grepl("WGCNA",names(func_results)))) { # Clustered
        message("Rendering specific cluster reports")
        invisible(parallel::mclapply(cls, function(cl) {
            # Take output name
            aux <- paste0("cl_func_",cl,".html")
            outf_cls_i <- file.path(results_path, aux)
            # Generate report
            rmarkdown::render(file.path(template_folder, 'cl_func_report.Rmd'), output_file = outf_cls_i, intermediates_dir = results_path)
        }, mc.cores = cores))

        message("\tRendering clustered report")
        outf_cls <- file.path(results_path, "clusters_func_report.html")
        rmarkdown::render(file.path(template_folder, 'clusters_main_report.Rmd'),output_file = outf_cls, intermediates_dir = results_path)
    }

    ############################################################
    ##              GENERATE DEG FUNCTIONAL REPORT            ##
    ############################################################
    message("\tRendering regular report")
    outf <- file.path(results_path, "functional_report.html")
    rmarkdown::render(file.path(template_folder, 'functional_report.Rmd'), output_file = outf, intermediates_dir = results_path)
}