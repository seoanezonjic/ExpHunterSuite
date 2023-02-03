#' Write Main DEgenes Hunter report
#'
#' This function allows you to report the DEG analysis.
#' @param exp_results DEG analysis results object
#' @param output_files output file path
#' @param template_folder RMD templates folder
#' @param opt option list
#' @return void
#' @export
#' @importFrom rmarkdown render
#' @examples
#' # Load DE analysis results
#' degh_output <- list() # data(degh_output)
#' write_expression_report(degh_output)
write_expression_report <- function(exp_results, 
    output_files=getwd(),
    template_folder = file.path(find.package('ExpHunterSuite'), 'templates'), 
    opt=NULL){
    if(length(exp_results) == 0){
        warning("Experiment results is not complete")
        return(NULL)
    }
    if(is.null(opt)){ opt <- exp_results[['final_main_params']]}
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
    final_results <- exp_results[['final_results']] 
    var_filter <-  exp_results[['var_filter']] 
    cpm_table <- exp_results[['cpm_table']]

    outf <- file.path(normalizePath(output_files),"DEG_report.html")
    rmarkdown::render(file.path(template_folder, 'main_report.Rmd'),
                      output_file = outf, intermediates_dir = output_files)
}


#' @importFrom heatmaply heatmaply
write_summarize_heatmaps <- function(summarized_ORA, output_path) {
  for (funsys in names(summarized_ORA)) {
   summ_ora_funsys <- summarized_ORA[[funsys]]
   heatmaply::heatmaply(summ_ora_funsys$summ_enr_table, 
                         grid_color = "gray50",
                         seriate = "mean",
                         grid_width = 0.00001,
                         dendrogram = "both",
                         scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
                         low = "#EE8291", 
                         high = "white", 
                         midpoint = 0.5, 
                         limits = c(0, 1)),
                         file = file.path(output_path, paste0("summ_",funsys,'_heatmap.html')))
    heatmaply::heatmaply(summ_ora_funsys$summ_enr_clean, 
                         grid_color = "gray50",
                         seriate = "mean",
                         grid_width = 0.00001,
                         fontsize_row = 11,
                          fontsize_col = 13,
                         dendrogram = "both",
                         scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
                         low = "#EE8291", 
                         high = "white", 
                         midpoint = 0.5, 
                         limits = c(0, 1)),
                         file = file.path(output_path, paste0("summ_rem_parent_",funsys,'_heatmap.html')))
    heatmaply::heatmaply(summ_ora_funsys$full_enr_table, 
                        grid_color = "gray50",
                        seriate = "mean",
                        dendrogram = "both",
                        grid_width = 0.00001,
                        scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
                        low = "#EE8291", 
                        high = "white", 
                        midpoint = 0.5, 
                        limits = c(0, 1)),
                        file = file.path(output_path, paste0("full_",funsys,'_heatmap.html')))
  }
}


write_merged_cluster_report <- function(enrichments_ORA, results_path, template_folder, 
    sample_classes=NULL, DEGH_results=NULL, showCategories, group_results) {
    message("\tRendering full cluster reports")
    if(is.null(enrichments_ORA)) {
        message("No WGCNA ORA results, not printing cluster report")
    } else {
        flags_cluster <- sapply(enrichments_ORA, function(x) nrow(x@compareClusterResult)) != 0
        names(flags_cluster) <- names(enrichments_ORA)
        outf_cls <- file.path(results_path, "clusters_func_report.html")
        rmarkdown::render(file.path(template_folder, 
          'clusters_main_report.Rmd'), output_file = outf_cls, 
          intermediates_dir = results_path)
    }
}


#' Write Main clusters to enrichment output
#' This function allows you to report the Functional analysis.
#' @param output_path output folder
#' @param output_file output file name for heatmaps
#' @param mode type of output to produce
#' @param enrichments_ORA list of enrich results for all clusters
#' @param enrichments_ORA_merged merged enrich results
#' @param task_size number of elements per packages used
#' @param workers (OPTIONAL) cores for parallel features
#' @param template_folder (OPTIONAL) RMD templates folder
#' @param top_categories numbers of categories from each cluster to use for merge
#' @param group_results experimental - whether to group results in the emap plot
#' @param n_category number of categories in the figures (per cluster)
#' @param sim_thr value to use when combining similar categories in summary mode
#' @param summary_common_name 'significant' to use the most significant term to label each summarized group
#' 'ancestor' to use the common ancestor of the group"
#' @param pvalcutoff used to select terms for summarizing
#' @param gene_attributes named list of attributes e.g. FCs for emap plot coloured nodes (genes)
#' @param gene_attribute_name name for the legend in the emap plot for the nodes (genes)
#' @param max_genes maximum number of genes to plot in cnet plot
#' @return void
#' @importFrom enrichplot emapplot
#' @importFrom enrichplot dotplot
#' @importFrom ggplot2 ggsave
#' @export
write_clusters_to_enrichment <- function( 
  output_path="results",
  output_file="results",
  mode="PR",
  enrichments_ORA=NULL,
  enrichments_ORA_merged=NULL,
  task_size = 1,
  workers = 1,
  template_folder = file.path(find.package('ExpHunterSuite'), 'templates'),
  top_categories = 50,
  group_results = FALSE,
  n_category = 30,
  sim_thr = 0.7, 
  summary_common_name = "ancestor", 
  pvalcutoff = 0.1, 
  gene_attributes=NULL,
  gene_attribute_name=NULL, 
  max_genes = 200) {

  if (grepl("R", mode)){
      enrichments_for_reports <- parse_results_for_report(enrichments_ORA)  
      write_enrich_clusters(enrichments_ORA, output_path)
      write_func_cluster_report(enrichments_for_reports, output_path, gene_attributes, 
        workers = workers, task_size = task_size, template_folder=template_folder, gene_attribute_name=gene_attribute_name)
  }

  if (grepl("P", mode) || grepl("S", mode)) 
    enrichments_ORA_merged <- filter_top_categories(enrichments_ORA_merged, top_categories)

  if (grepl("P", mode)) {
    for (funsys in names(enrichments_ORA_merged)){
      if (length(unique(enrichments_ORA_merged[[funsys]]@compareClusterResult$Description)) < 2 ) next

      if (group_results == TRUE){
        pp <- enrichplot::emapplot(enrichments_ORA_merged[[funsys]], showCategory= n_category, pie="Count", layout = "nicely", 
                    shadowtext = FALSE, node_label = "group", group_category = TRUE, 
                    nCluster = min(floor(nrow(enrichments_ORA_merged[[funsys]])/7), 20), nWords = 6, repel = TRUE)
      }else{
        pp <- enrichplot::emapplot(enrichments_ORA_merged[[funsys]], showCategory= n_category, pie="Count", layout = "nicely", 
                    shadowtext = FALSE, repel = TRUE)
      }

      ggplot2::ggsave(filename = file.path(output_path,paste0("emaplot_",funsys,".png")), pp, width = 30, height = 30, dpi = 300, units = "cm", device='png')

      pp <- enrichplot::dotplot(enrichments_ORA_merged[[funsys]], showCategory= n_category, label_format = 70)
      ggplot2::ggsave(filename = file.path(output_path,paste0("dotplot_",funsys,".png")), pp, width = 60, height = 40, dpi = 300, units = "cm", device='png')

    }
  }

  if (grepl("S", mode)){
    summarized_merged_ora <- summarize_merged_ora(enrichments_ORA_merged, sim_thr, summary_common_name, pvalcutoff)
    write_summarize_heatmaps(summarized_merged_ora, output_path)
  }
  if(grepl("R", mode)) {
    write_merged_cluster_report(enrichments_ORA_merged, results_path=output_path, template_folder, 
            showCategories=n_category, group_results=group_results)
  }
}



#' Write Main DEgenes Hunter functional report
#' This function allows you to report the Functional analysis.
#' @param hunter_results DEG analysis results
#' @param func_results functional results
#' @param output_files output folder.
#' @param template_folder (OPTIONAL) RMD templates folder
#' @param cores (OPTIONAL) cores for parallel features
#' @param organisms_table (OPTIONAL) configuration table for given organism.
#'  Use see get_organism_table()
#' @param fc_colname (OPTIONAL) main logFC colname (into hunter_results
#'  dataframe)
#' @param task_size number of elements per packages used
#' @param report string with reports to be written. Allowed: clusters (c)
#' @param showCategories number of categories in the figures (per cluster) 
#' @param group_results experimental - whether to group results in the emap plot
#'  and functional (f). Default = "fc"
#' @param max_genes maximum number of genes to plot in cnet plot
#' @return void
#' @importFrom rmarkdown render
#' @export
#' @examples
#' # Load func and DE results
#' data(degh_output)
#' func_results <- list() 
#' func_results <- main_functional_hunter(degh_output, "Mouse")
#' write_functional_report(degh_output, func_results)
write_functional_report <- function(hunter_results, 
                                    func_results, 
                                    output_files=getwd(), 
                                    fc_colname="mean_logFCs", 
                                    organisms_table=NULL, 
                                    template_folder = file.path(find.package('ExpHunterSuite'), 'templates'),
                                    cores = 2,
                                    task_size = 1, 
                                    report = "fc",
                                    showCategories = 30,
                                    group_results = FALSE,
                                    max_genes = 200
                                    ){
    # TO parallelize properly
    clean_tmpfiles_mod <- function() {
      message("Calling clean_tmpfiles_mod()")
    }
    assignInNamespace("clean_tmpfiles", clean_tmpfiles_mod, ns = "rmarkdown")

    if(!any(grepl("WGCNA", names(func_results))) && grepl("c|i", report)) {
        message("Cluster reports chosen but no cluster results available. Reports wont be plotted")
    }
    results_path <- normalizePath(output_files)
    model_organism <- func_results$final_main_params$model_organism

    if(is.null(organisms_table)){
        organisms_table <- get_organism_table()
    }
    current_organism_info <- subset(organisms_table, 
    rownames(organisms_table) %in% model_organism) 
    
    # Prepare the flag lists
    flags_ora <- sapply(func_results$ORA, nrow) != 0
    names(flags_ora) <- names(func_results$ORA)
    flags_gsea <- sapply(func_results$GSEA, nrow) != 0
    names(flags_gsea) <- names(func_results$GSEA)

    # JRP to clean up - should take target directly
    experiments <- hunter_results$sample_groups
    sample_classes <- apply(experiments, 1, function(x) paste0("* [", x[1],
                      "] ", x[2]))

    fc_vector <- func_results$DEGH_results_annot[
       !is.na(func_results$DEGH_results_annot$ENTREZID), fc_colname]
    names(fc_vector) <- func_results$DEGH_results_annot[
       !is.na(func_results$DEGH_results_annot$ENTREZID), "ENTREZID"]

    enrichments_ORA <- func_results$WGCNA_ORA
    DEGH_results <- func_results$DEGH_results_annot
    enrichments_ORA_expanded <- func_results$WGCNA_ORA_expanded

    # JRP This will get us one day
    norm_counts <- hunter_results[["all_data_normalized"]][["DESeq2"]]
    scaled_counts <- scale_data_matrix(data_matrix = as.matrix(norm_counts))
    scaled_counts_table <- as.data.frame(as.table(scaled_counts))
    colnames(scaled_counts_table) <- c("Gene","Sample","Count")

    if(grepl("f", report)){
        message("\tRendering regular report")
        outf <- file.path(results_path, "functional_report.html")
        rmarkdown::render(file.path(template_folder, 'functional_report.Rmd'), 
            output_file = outf, intermediates_dir = results_path)        
    }

    if(grepl("c", report)){
        write_merged_cluster_report(enrichments_ORA, results_path, template_folder, 
            sample_classes, DEGH_results, showCategories, group_results)
        write_summarize_heatmaps(func_results$summarized_ora, results_path)
    }

    if(grepl("a", report)){

         filter_compareCluster <- function(compareCluster, clusters_to_fil){
             fil_obj <- compareCluster
             fil_obj@compareClusterResult <- fil_obj@compareClusterResult[
                                     fil_obj@compareClusterResult$Cluster %in% clusters_to_fil,]
             return(fil_obj)
         }

         lapply(func_results$WGCNA_ORA, function(funsys){
             fil_funsys <- filter_compareCluster(funsys, fil_clusters)
         })
         write_merged_cluster_report(enrichments_ORA, results_path, template_folder,
             sample_classes, DEGH_results, showCategories, group_results)
    }
    if(grepl("i", report)) {
        message("\tRendering individual cluster reports")
        if(is.null(enrichments_ORA)) {
          message("No WGCNA ORA results, not printing individual cluster report")
        } else {
        cls  <- unique(DEGH_results$Cluster_ID)
        cls <- cls[cls != 0]
        trait_module <- hunter_results$WGCNA_all$plot_objects$trait_and_module
        cl_eigvalues <- as.matrix(trait_module[,grepl("^ME",colnames(trait_module))])
        cl_eigvalues <- as.data.frame(as.table(cl_eigvalues),
          stringsAsFactors = FALSE)
        colnames(cl_eigvalues) <- c("Sample","Cluster_ID","Count")
        cl_eigvalues_gnorm <- cl_eigvalues
        cl_eigvalues_gnorm$Count <- (cl_eigvalues_gnorm$Count + 1) / 2
        pack_obj <- hunter_results$WGCNA_all$package_objects
        wgcna_pval_cl_trait <- as.matrix(pack_obj$module_trait_cor_p)
        wgcna_corr_cl_trait <- as.matrix(pack_obj$module_trait_cor)
        wgcna_count_sample_trait <- as.matrix(trait_module[, !grepl("^ME",
          colnames(trait_module))])
        wgcna_count_sample_trait <- scale_data_matrix(wgcna_count_sample_trait, 
          norm_by_col = TRUE)

        message("\tRendering specific cluster reports")
        parallel_list(cls, function(cl) {
        #lapply(cls, function(cl) {
            cl_flags_ora <- lapply(enrichments_ORA_expanded, function(x) {
              nrow(x[[which(names(x) == cl)]]) != 0 
            })
            temp_path_cl <- file.path(results_path, paste0(cl,"_temp_cl_rep"))
            outf_cls_i <- file.path(results_path, paste0("cl_func_",cl,".html"))
            DEGH_results <- DEGH_results[which(DEGH_results$Cluster_ID == cl), ]
            rmarkdown::render(file.path(template_folder, 
                   'cl_func_report.Rmd'), output_file = outf_cls_i, 
                   clean=TRUE, intermediates_dir = temp_path_cl)
        #}) 
        }, workers=cores, task_size=task_size)
        unlink(list.files(results_path, pattern="_temp_cl_rep$", full.names=TRUE), recursive=TRUE) 
      }
    }
}

#' @importFrom rmarkdown render
write_miRNA_cor_report <- function(
report_name,
template_folder, 
output_files , 
p_val_cutoff,
corr_cutoff, 
strategies, 
unsig_strategies,
cont_tables,
filters_summary,
score_comp,
all_pairs,
mirna_names,
gene_id_translation,
sample_proportion,
selected_predicted_databases,
all_cor_dist,
miRNAseq, 
miRNA_cont_tables,
eval_method,
miRNA_cont_tables_adj,
RNAseq,
sig_pairs,
raw_databases_scores,
p_fisher,
mapping_output,
output_pairs           
){
 "%>%" <- magrittr::"%>%"

 rmarkdown::render(
               file.path(template_folder, 'miRNA_RNA.Rmd'), 
               output_file = file.path(output_files, report_name), 
               intermediates_dir = file.path(output_files))
    
}

parse_strat_text <- function(strategies){
  o_text <- c()
  default_strats <-  c(
    "Eigengene_0_RNA_vs_miRNA_normalized_counts", 
    "normalized_counts_RNA_vs_miRNA_Eigengene_0", 
    "DEGs_DEMs_permutated")
  strategies <- strategies[!strategies %in% default_strats]
  dictionary <- list(
    "Eigengene" = "Eigengene profiles for coexpression modules",
    "hub_1" = "hub gene profile for coexpression modules",
    "normalized_counts" = "normalized counts"
    )
  for (strategy_name in strategies){
    strategy <- unlist(strsplit(strategy_name, "_RNA_vs_miRNA_"))
    o_text <- c(o_text,
      paste0("\t+ **", 
             strategy_name, 
             ":** correlates ", 
             dictionary[[strategy[1]]], 
             " of RNAseq data with ", 
             dictionary[[strategy[1]]],
             " of miRNAseq data."))
  }
  return(paste(o_text, collapse = "\n"))
} 


write_func_cluster_report <- function(enrichments_for_reports, output_path, 
  gene_attributes, workers, task_size, template_folder, gene_attribute_name="fold change", max_genes = 200){
  clean_tmpfiles_mod <- function() {
    message("Calling clean_tmpfiles_mod()")
  }
  assignInNamespace("clean_tmpfiles", clean_tmpfiles_mod, ns = "rmarkdown")

  parallel_list(names(enrichments_for_reports), function(cluster) {
    temp_path_cl <- file.path(output_path, paste0(cluster,"_temp"))
    func_results <- enrichments_for_reports[[cluster]]
    cl_flags_ora <- sapply(func_results, nrow) > 0
    fc_vector <- gene_attributes[[cluster]]
    outfile <- file.path(output_path, paste0(cluster, "_func_report.html"))
    test_env <- list2env(list(func_results=func_results, 
      cl_flags_ora=cl_flags_ora))
    rmarkdown::render(file.path(template_folder, 
                   'clusters_to_enrichment.Rmd'), output_file = outfile, 
               clean=TRUE, intermediates_dir = temp_path_cl, envir=test_env)
  }, workers=workers, task_size=task_size)
  # temp files not deleted properly in parallel 
  unlink(list.files(output_path, pattern="_temp$", full.names=TRUE), 
    recursive=TRUE)
}
