#' Write Main DEgenes Hunter report
#'
#' This function allows you to report the DEG analysis.
#' @param exp_results DEG analysis results object
#' @param output_files output file path
#' @param template_folder htmlreportR templates folder
#' @param opt option list
#' @param source_folder Folder where js and css templates to load in htmlreportR
#' object are located.
#' @return void
#' @export
#' @importFrom rmarkdown render
#' @examples
#' \dontrun{
#'      # Load DE analysis results
#'      degh_output <- list() # data(degh_output)
#'      write_expression_report(degh_output)
#' }
write_expression_report <- function(exp_results, 
    output_files = getwd(),
    template_folder = NULL, 
    opt = NULL,
    source_folder = NULL){
    if (is.null(template_folder))
      template_folder <- file.path(find.package('ExpHunterSuite'), 'templates')
    
    if(length(exp_results) == 0){
      warning("Experiment results is not complete")
      return(NULL)
    }
    if(is.null(opt)){
      opt <- exp_results[['final_main_params']]}
    if(is.null(source_folder)){
      source_folder <- find.package("htmlreportR")
    }
    template <- file.path(template_folder, "main_report.txt")
    tmp_folder <- file.path(output_files, "tmp_lib")
    container <- list(opt = opt, template_folder = template_folder,
    final_main_params = exp_results[['final_main_params']],
    DEG_pack_columns = exp_results[['DEG_pack_columns']],
    all_counts_for_plotting = exp_results[['all_counts_for_plotting']],
    all_FDR_names = exp_results[['all_FDR_names']],
    all_LFC_names = exp_results[['all_LFC_names']],
    all_pvalue_names = exp_results[['all_pvalue_names']],
    final_pvalue_names = exp_results[['final_pvalue_names']],
    final_logFC_names = exp_results[['final_logFC_names']],
    final_FDR_names = exp_results[['final_FDR_names']],
    package_objects = exp_results[['package_objects']],
    results_WGCNA = exp_results[['WGCNA_all']],
    index_control_cols = exp_results[['index_control_cols']],
    index_treatmn_cols = exp_results[['index_treatmn_cols']],
    raw_filter = exp_results[['raw_filter']],
    design_vector = exp_results[['design_vector']],
    all_data_normalized = exp_results[['all_data_normalized']],
    replicatesC = exp_results[['replicatesC']],
    replicatesT = exp_results[['replicatesT']],
    DE_all_genes = exp_results[['DE_all_genes']],
    final_results = exp_results[['final_results']],
    var_filter =  exp_results[['var_filter']],
    cpm_table = exp_results[['cpm_table']],
    coverage_df = exp_results[['coverage_df']],
    mean_counts_df = exp_results[['mean_counts_df']],
    exp_genes_df = exp_results[['exp_genes_df']],
    numeric_factors = exp_results[["numeric_factors"]],
    string_factors = exp_results[["string_factors"]],
    PCA_res = exp_results[["PCA_res"]],
    library_sizes = exp_results[["library_sizes"]],
    target = exp_results[["target"]],
    WGCNA_results = exp_results[["WGCNA_results"]])
    
    outf <- file.path(normalizePath(output_files), "DEG_report.html")
    plotter <- htmlreportR::htmlReport$new(title_doc = "DEG_report", 
                                            container = container,
                                            tmp_folder = tmp_folder,
                                            src = source_folder,
                                            compress_obj = TRUE,
                                            type_index = "contents_list")
    plotter$build(template)
    plotter$write_report(outf)
    message(paste0("Report written in ", outf))
}


#' Function to create the expression hunter folder
#' Used by the degenes_Hunter script
#' @param final_results object containing expression results,
#' as created by main_degenes_Hunter
#' @param output_files Output path
#' @return used for side effect of creating the folder and files
#' @export
#' @examples
#' # Can only be run after running main_degenes_Hunter 
#' # eg as part of the degenes_Hunter.R script
#' \dontrun{
#'      write_expression_data(final_results, opt$output_files)
#' }
write_expression_data <- function(final_results, output_files){

  write.table(final_results[['raw_filter']], 
    file=file.path(output_files, "filtered_count_data.txt"), quote=FALSE, 
    col.names=NA, sep="\t")
  write.table(final_results[['sample_groups']], file=file.path(output_files, 
    "control_treatment.txt"), row.names=FALSE, quote=FALSE, sep="\t")
  write_df_list_as_tables(final_results[['all_data_normalized']], 
    prefix = 'Normalized_counts_', root = output_files)
  write_df_list_as_tables(final_results[['all_counts_for_plotting']], 
    prefix = 'allgenes_', root = output_files)
  dir.create(file.path(output_files, "Common_results"))
  write.table(final_results[['DE_all_genes']], file=file.path(output_files, 
    "Common_results", "hunter_results_table.txt"), quote=FALSE, 
  row.names=TRUE, sep="\t")

  
  write_pca_data(final_results[['PCA_res']], output_files)
}

write_pca_data <- function(PCA_res, output_files){

    pca_output <- file.path(output_files, "PCA_results")
    dir.create(pca_output)
    all_genes_pca <- PCA_res$all_genes
    
    write_general_pca(all_genes_pca, pca_output, "all_genes_")


    prevalent_pca <- PCA_res$DEGs
    if (!is.null(prevalent_pca)){
       write_general_pca(prevalent_pca, pca_output, "prevalent_")
    }

    saveRDS(PCA_res, file = file.path(pca_output, "pca_data.rds"))
}

write_general_pca <- function(pca_data, output_files, tag = ""){
  pca_res <- pca_data$dim_data_merged
  
  merged_metrics <- merge_dim_table_metrics(pca_res)
  if (!is.null(merged_metrics))
      write.table(merged_metrics, file = file.path(output_files, paste0(tag, "dim_metrics.txt")),sep = "\t", quote = FALSE, row.names=FALSE)
  
  dimnames(pca_data$pca_data$svd$V) <- dimnames(pca_data$pca_data$var$cor)
  write.table(pca_data$pca_data$svd$V, quote = FALSE, file = file.path(output_files, paste0(tag, "eigenvectors.txt")), sep = "\t")
  hcpc_table <- pca_data$res.hcpc$call$X
  hcpc_table$samples <- rownames(hcpc_table)
  write.table(hcpc_table, file = file.path(output_files, paste0(tag, "hcpc.txt")),sep = "\t", quote = FALSE, row.names=FALSE)
  if (!is.null(pca_data$pca_data$var$cor)){
    write.table(pca_data$pca_data$var$cor, file = file.path(output_files, paste0(tag, "pca_vars.txt")), quote = FALSE, sep = "\t")
  } else {
    write.table(pca_data$pca_data$var$eta2, file = file.path(output_files, paste0(tag, "pca_vars.txt")), quote = FALSE, sep = "\t")
  }

}

merge_dim_table_metrics <- function(merged_dim_table){
        if (is.null(merged_dim_table))
            return(NULL)
        
        if(nrow(merged_dim_table$qualitative) > 0) {
            names(merged_dim_table$qualitative)[names(merged_dim_table$qualitative) == "R2"] <- "metric"
            merged_dim_table$qualitative$metric_type <- "R2"
        } else { merged_dim_table$qualitative <- NULL}
        
        if(nrow(merged_dim_table$quantitative) > 0) {
            names(merged_dim_table$quantitative)[names(merged_dim_table$quantitative) == "correlation"] <- "metric"
            merged_dim_table$quantitative$metric_type <- "correlation"
        } else { merged_dim_table$quantitative <- NULL}

        if(nrow(merged_dim_table$qual_category) > 0) {
            names(merged_dim_table$qual_category)[names(merged_dim_table$qual_category) == "Estimate"] <- "metric"
            merged_dim_table$qual_category$metric_type <- "coord_var_barycentre"
         } else { merged_dim_table$qual_category <- NULL}

        dim_data_merged <-data.table::rbindlist(merged_dim_table, use.names = TRUE,idcol = "var_type")
        dim_data_merged <- as.data.frame(dim_data_merged)
        dim_data_merged <- dim_data_merged[,c("factor","var_type" ,"metric_type", "dimension","metric","p.value")]
        return(dim_data_merged)
} 

#' Write global coRmiT report
#'
#' This function writes a report of coRmiT results.
#' @param cormit_res coRmiT results object
#' @param output_files output file path
#' @param template_folder htmlreportR templates folder
#' @param source_folder Folder where js and css templates to load in htmlreportR
#' object are located.
#' @return void
#' @export
#' @importFrom rmarkdown render
#' @examples
#' # Load DE analysis results
#' \dontrun{
#'      degh_output <- list() # data(degh_output)
#'      write_global_cormit(degh_output)
#' }

write_global_cormit <- function(cormit_res, 
    output_files=getwd(),
    template_folder = NULL,
    source_folder = NULL){
    integrated_pairs <- as.data.frame(cormit_res$integrated_pairs)
    miRNA_cont_tables <- as.data.frame(cormit_res$miRNA_cont_tables)
    mirna_names <- cormit_res$mirna_names
    integrated_pairs <- integrated_pairs[integrated_pairs$miRNAseq %in% unique(miRNA_cont_tables[miRNA_cont_tables$db_group == "multimir", "miRNA"]),]
    miRNA_cont_tables$miRNA <- mirna_names[match(miRNA_cont_tables$miRNA, mirna_names$ACCESSION), "NAME"]
    cormit_res$integrated_pairs <- integrated_pairs
    cormit_res$miRNA_cont_tables <- miRNA_cont_tables
    if(is.null(template_folder)){
      template_folder <- file.path(find.package('ExpHunterSuite'), 'templates')
    }
    if(length(cormit_res) == 0){
        warning("coRmiT object is incomplete")
        return(NULL)
    }
    if(is.null(opt)){
        opt <- list(report_name = cormit_res$opt$report,
        output_files = cormit_res$opt$output_files,
        # p_val_cutoff = cormit_res$opt$p_val_cutoff,
        # corr_cutoff = cormit_res$opt$corr_cutoff,
        # eval_method = cormit_res$opt$eval_method,
        # sample_proportion = cormit_res$opt$sample_proportion,
        genomic_ranges = cormit_res$opt$genomic_ranges,
        genome_ref = cormit_res$opt$genome_ref,
        mapping_output = cormit_res$opt$mapping_output,
        output_pairs = cormit_res$opt$output_pairs)
    }
    if(is.null(source_folder)){
      source_folder <- find.package("htmlreportR")
    }
    template <- file.path(template_folder, "global_cormit.txt")
    tmp_folder <- file.path(output_files, "tmp_lib")
    container <- list(miRNA_cor_results = cormit_res$miRNA_cor_results,
                      strategies = cormit_res$strategies,
                      all_pairs = cormit_res$all_pairs,
                      miRNA_cont_tables = cormit_res$miRNA_cont_tables,
                      cont_tables = cormit_res$cont_tables,
                      mirna_names = cormit_res$mirna_names,
                      int_miRNA_cont_tables = cormit_res$int_miRNA_cont_tables,
                      int_cont_tables = cormit_res$int_cont_tables,
                      integrated_pairs = cormit_res$integrated_pairs,
                      gene_id_translation = cormit_res$gene_id_translation,
                      miRNAseq = cormit_res$miRNAseq,
                      RNAseq = cormit_res$RNAseq,
                      sel_pred_data = cormit_res$selected_predicted_databases,
                      all_cor_dist = cormit_res$all_cor_dist,
                      multimir_stats = cormit_res$multimir_stats,
                      # p_val_cutoff = opt$p_val_cutoff,
                      # corr_cutoff = opt$corr_cutoff,
                      # eval_method = opt$eval_method,
                      # sample_proportion = opt$sample_proportion,
                      genomic_ranges = opt$genomic_ranges,
                      genome_ref = opt$genome_ref,
                      mapping_output = opt$mapping_output,
                      output_pairs = opt$output_pairs)
    outf <- file.path(normalizePath(opt$output_files), "coRmiT_report.html")
    plotter <- htmlreportR::htmlReport$new(title_doc = "coRmiT report", 
                                            container = container,
                                            tmp_folder = tmp_folder,
                                            src = source_folder,
                                            compress_obj = TRUE,
                                            type_index = "contents_list")
    plotter$build(template)
    plotter$write_report(outf)
    message(paste0("Report written in ", outf))
    if (!is.null(opt$genomic_ranges)){
        g_ranges <- read.table(opt$genomic_ranges, header = FALSE)
        g_ranges <- g_ranges[,c(1,2,3,4,6)]
        colnames(g_ranges) <- c("chromosome", "start","end","miRNA","strand")
        annotated <- annotate_genomic_ranges(g_ranges, opt$genome_ref)
        annotated <- annotated[!grepl("MIR",annotated$annot.symbol),]
        miRNA_annot <- aggregate(annot.symbol ~ miRNA, annotated, unique)
        miRNA_annot$annot.symbol <- unlist(lapply(miRNA_annot$annot.symbol,
                                                  paste, collapse = ","))
        miRNA_loci <- miRNA_annot[match(integrated_pairs$miRNAseq,
                                        miRNA_annot$miRNA),"annot.symbol"]
        integrated_pairs$miRNA_loci <- miRNA_loci
        integrated_pairs[is.na(integrated_pairs$miRNA_loci),"miRNA_loci"] <- ""
    }    

    integrated_pairs$miRNA <- mirna_names[match(cormit_res$integrated_pairs$miRNAseq,
                                                mirna_names$ACCESSION), "NAME"]
    output_pairs_all <- add_attrib_to_pairs(integrated_pairs, cormit_res$RNAseq, cormit_res$miRNAseq)
    colnames(output_pairs_all)[match(c("normalized_counts_RNA_vs_miRNA_normalized_counts_correlation","normalized_counts_RNA_vs_miRNA_normalized_counts_pval"), colnames(output_pairs_all))] <- c("Raw_correlation", "Raw_cor_Pcalue")
    gene_id_translation <- as.data.frame(cormit_res$gene_id_translation)
    output_pairs_all$Target_SYMBOL <- gene_id_translation[match(output_pairs_all$Target_ID, gene_id_translation$ensembl_gene_id), "Symbol"]
   
    if (!cormit_res$output_pairs %in% c("predicted", "validated", "multimir"))   
        cormit_res$output_pairs <- "multimir"
  
    integrated_pairs$db_type <- ifelse(integrated_pairs[,cormit_res$output_pairs], "DB","ND")
   
    out_pairs <- data.frame()
    genes_attr <- data.frame()
    attr <- cormit_res$all_pairs[,c("miRNAseq", "RNAseq", "validated_c", "predicted_c")]
    for (miRNA in unique(integrated_pairs$miRNA)){ 

        DB <- data.frame(miRNA = paste0(miRNA, "_DB"), 
                         genes= paste(integrated_pairs$RNAseq[integrated_pairs$db_type== "DB" & integrated_pairs$miRNA == miRNA], collapse = ","))
        
        if (cormit_res$output_pairs == "All") {

            ALL <- data.frame(miRNA = paste0(miRNA, "_ALL"), 
            genes= paste(integrated_pairs$RNAseq,collapse = ","))
            out_pairs <- rbind(out_pairs, DB, ALL)
        } else {
            out_pairs <- rbind(out_pairs, DB)
        }
    }

best_strats <- select_best_strategy(cormit_res$int_miRNA_cont_tables[cormit_res$int_miRNA_cont_tables$db_group == "multimir",])[,c("miRNA", "Odds_ratio")]
output_pairs_all <- merge(output_pairs_all, best_strats, by.x ="miRNA_ID", by.y = "miRNA", all.x = TRUE)


    write.table(out_pairs, col.names = FALSE, sep = "\t",file = file.path(output_files,"integrated_miRNA.txt"), quote = FALSE, row.names = FALSE)
    write.table(output_pairs_all, col.names = TRUE, sep = "\t",file = file.path(output_files,"target_results_table.txt"), quote = FALSE, row.names = FALSE)

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
                        file = file.path(output_path, paste0("full_",
                                         funsys,'_heatmap.html')))
  }
}

#' Write clustering report for functional Hunter
#' This function will usually be called by write_functional_report.
#' @inheritParams write_functional_report
#' @param enrichments_ORA ORA enrichment resunts.
#' @param results_path Path where report will be written.
#' @param sample_classes Sample classes.
#' @param DEGH_results DEGenes Hunter results object.
#' @returns invisible(NULL)
write_merged_cluster_report <- function(enrichments_ORA, results_path,
    template_folder, sample_classes=NULL, DEGH_results=NULL, showCategories,
    group_results, func_results = NULL, source_folder = NULL) {
    message("\tRendering full cluster reports")
    if(is.null(enrichments_ORA)) {
        message("No WGCNA ORA results, not printing cluster report")
    } else {
        if(is.null(source_folder)) {
            source_folder <- find.package("htmlreportR")
        }
        flags_cluster <- sapply(enrichments_ORA,
                                function(x) nrow(x@compareClusterResult)) != 0
        names(flags_cluster) <- names(enrichments_ORA)
        outf_cls <- file.path(results_path, "clusters_func_report.html")
        tmp_folder <- file.path(results_path, "tmp_lib")
        template <- file.path(template_folder, "clusters_main_report.txt")
        container <- list(flags_cluster = flags_cluster,
            DEGH_results = DEGH_results, sample_classes = sample_classes,
            func_results = func_results, enrichments_ORA = enrichments_ORA,
            group_results = group_results, showCategories = showCategories)
        plotter <- htmlreportR::htmlReport$new(container = container,
                                 title_doc = "clusters functional report",
                                 tmp_folder = tmp_folder, src = source_folder,
                                 compress_obj = TRUE)
        plotter$build(template)
        plotter$write_report(outf_cls)
        message(paste0("Report written in ", outf_cls))
    }
}


#' Write Main clusters to enrichment output
#' This function allows you to report the Functional analysis.
#' @param output_path output folder
#' @param output_file output file name for heatmaps
#' @param mode type of output to produce - P for plots, R for reports, and S for summarized heatmaps
#' @param enrichments_ORA list of enrich results for all clusters
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
#' @param simplify Activate for process ClusterProfiler results with simplify function
#' @param clean_parentals Activate to reduce significant terms in GO by removin the parentals
#' @param source_folder htmlreportR source folder to load js and css libraries
#' @return void
#' @importFrom enrichplot emapplot
#' @importFrom enrichplot dotplot
#' @importFrom ggplot2 ggsave
#' @examples
#' \dontrun{
#'    write_clusters_to_enrichment(enrichments_ORA = ORA_results)
#'}
#' @export
write_clusters_to_enrichment <- function(output_path="results",
      output_file="results", mode="PR", enrichments_ORA=NULL, task_size = 1,
      workers = 1, top_categories = NULL, group_results = FALSE,
      n_category = 30, sim_thr = 0.7, pvalcutoff = 0.1, gene_attributes=NULL,
      max_genes = 200, summary_common_name = "ancestor", simplify = FALSE,
      gene_attribute_name=NULL, clean_parentals = FALSE, source_folder = NULL,
      template_folder = file.path(find.package('ExpHunterSuite'), 'templates')){
      enrichments_ORA_merged <- process_cp_list(enrichments_ORA,
                simplify_results = simplify, clean_parentals = clean_parentals)
      if(!is.null(top_categories)) {
        enrichments_ORA_merged <- filter_top_categories(enrichments_ORA_merged,
                                                        top_categories)
      }
      if(is.null(source_folder)) {
        source_folder <- find.package("htmlreportR")
      }
      if (grepl("R", mode)){
          enrichments_for_reports <- parse_results_for_report(enrichments_ORA)  
          write_enrich_clusters(enrichments_ORA, output_path)
          write_func_cluster_report(enrichments_for_reports, output_path,
            gene_attributes, workers = workers, task_size = task_size,
            template_folder = template_folder, source_folder = source_folder,
            gene_attribute_name = gene_attribute_name)
      }
      if (grepl("P", mode)) {
        for (funsys in names(enrichments_ORA_merged)){
            enrich <- enrichments_ORA_merged[[funsys]]
            desc <- enrich@compareClusterResult$Description
          if (length(unique(desc)) < 2 ) next
          if (group_results == TRUE){
            n_Cluster <- min(floor(nrow(enrich)/7), 20)
            pp <- enrichplot::emapplot(enrich, showCategory= n_category,
                                      pie="Count", layout = "nicely",
                                      shadowtext = FALSE, node_label = "group",
                                      group_category = TRUE, 
                                      nCluster = min(floor(nrow(enrich)/7), 20),
                                                     nWords = 6, repel = TRUE)
          }else{
            pp <- enrichplot::emapplot(enrich, showCategory= n_category,
                                       pie="Count", layout = "nicely", 
                                       shadowtext = FALSE, repel = TRUE)
          }
          ggplot2::ggsave(filename = file.path(output_path,
                            paste0("emaplot_", funsys, ".png")), pp, width = 30,
                            height = 30, dpi = 300, units = "cm", device='png')
          pp <- enrichplot::dotplot(enrich, showCategory= n_category,
                                    label_format = 70)
          ggplot2::ggsave(filename = file.path(output_path,
                            paste0("dotplot_", funsys, ".png")), pp, width = 60,
                            height = 40, dpi = 300, units = "cm", device='png')
        }
      }
      if (grepl("S", mode)) {
        summarized_merged_ora <- summarize_merged_ora(enrichments_ORA_merged,
                                       sim_thr, summary_common_name, pvalcutoff)
        write_summarize_heatmaps(summarized_merged_ora, output_path)
      }
      if(grepl("R", mode)) {
        write_merged_cluster_report(enrichments_ORA_merged,
                         results_path=output_path, template_folder, 
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
#' @param corr_threshold minimum module eigengene-trait vector absolute Pearson
#' R value  
#' @param pvalcutoff maximum module eigengene-trait vector correlation P value
#' @param source_folder htmlreportR source folder to load js and css libraries
#' @return void
#' @importFrom rmarkdown render
#' @export
#' @examples
#' \dontrun{
#' # Load func and DE results
#' data(degh_output)
#' func_results <- list() 
#' func_results <- main_functional_hunter(degh_output, "Mouse")
#' write_functional_report(degh_output, func_results)
#'}
write_functional_report <- function(hunter_results, func_results, cores = 2,
                                    output_files = getwd(),
                                    organisms_table = NULL,
                                    fc_colname = "mean_logFCs", task_size = 1, 
                                    report = "fc", source_folder = NULL,
                                    group_results = FALSE, showCategories = 30,
                                    corr_threshold = 0.8, pvalcutoff = 0.05,
                                    template_folder = NULL, max_genes = 200) {
    if(is.null(template_folder)){
      template_folder <- file.path(find.package('ExpHunterSuite'), 'templates')
    }
    if(length(hunter_results) == 0){
      warning("Incomplete hunter results")
      return(NULL)
    }
    if(length(func_results) == 0){
      warning("Incomplete functional results")
      return(NULL)
    }
    if(is.null(source_folder)){
      source_folder <- find.package("htmlreportR")
    }
    tmp_folder <- file.path(output_files, "tmp_lib")
    # # TO parallelize properly
    # clean_tmpfiles_mod <- function() {
    #   message("Calling clean_tmpfiles_mod()")
    # }
    # assignInNamespace("clean_tmpfiles", clean_tmpfiles_mod, ns = "rmarkdown")
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
    attr_vector <- func_results$DEGH_results_annot[
       !is.na(func_results$DEGH_results_annot$ENTREZID), fc_colname]
    names(attr_vector) <- func_results$DEGH_results_annot[
       !is.na(func_results$DEGH_results_annot$ENTREZID), "ENTREZID"]
    gene_attribute_name <- "Log2FC"
    enrichments_ORA <- func_results$WGCNA_ORA
    DEGH_results <- func_results$DEGH_results_annot
    enrichments_ORA_expanded <- func_results$WGCNA_ORA_expanded
    container <- list(hunter_results = hunter_results, max_genes = max_genes,
                  func_results = func_results,
                  model_organism = model_organism,
                  current_organism_info = current_organism_info,
                  flags_ora = flags_ora, flags_gsea = flags_gsea,
                  sample_classes = sample_classes,
                  attr_vector = attr_vector,
                  gene_attribute_name = gene_attribute_name,
                  enrichments_ORA = enrichments_ORA,
                  ennrichments_ORA_expanded = enrichments_ORA_expanded,
                  DEGH_results = DEGH_results, pvalcutoff = pvalcutoff)
    if(grepl("f", report)){
        message("\tRendering regular report")
        template <- file.path(template_folder, "functional_report.txt")
        outf <- file.path(results_path, "functional_report.html")
        plotter <- htmlreportR::htmlReport$new(title_doc = "functional report",
                                                container = container,
                                                tmp_folder = tmp_folder,
                                                src = source_folder,
                                                compress_obj = TRUE,
                                                type_index = "contents_list")
        plotter$build(template)
        plotter$write_report(outf)
        message(paste0("Report written in ", outf))
    }
    if(!any(grepl("WGCNA", names(func_results))) && grepl("c|i", report)) {
        message("Cluster reports chosen but no cluster results available. 
                Reports wont be plotted")
    } else {
        res <- hunter_results$WGCNA_all$package_objects$module_trait_cor_p
        mod_t_cor_p <- res
        mod_t_cor <- hunter_results$WGCNA_all$package_objects$module_trait_cor
        corr_cl <- mod_t_cor[abs(mod_t_cor[,"treat_Ctrl"]) > corr_threshold 
                                & mod_t_cor_p[,"treat_Ctrl"] < 0.05,]
        corr_cl <- rownames(corr_cl)
        corr_cl <- gsub("Cluster_","",corr_cl)
        if (length(corr_cl) > 0) {
            enrichments_ORA <- lapply(enrichments_ORA, 
                                      filter_cluster_enrichment,
                                      filter_list = corr_cl)
        } else {
            warning(paste0(c("There are not clusters with higher absolute ",
                             "correlation with treat/control hinger than ",
                             corr_threshold, ". Modify corr_threshold option. ",
                             "Reporting enrichments of all clusters...")))
        }
    }
    if(grepl("c", report)){
            message("Launching cluster report")
            write_merged_cluster_report(enrichments_ORA = enrichments_ORA,
                results_path = results_path, template_folder = template_folder,
                sample_classes = sample_classes, DEGH_results = DEGH_results, 
                showCategories = showCategories, group_results = group_results,
                func_results = func_results, source_folder = source_folder)
            write_summarize_heatmaps(func_results$summarized_ora, results_path)
    }

   
    if(grepl("i", report)) {
        #PCC adaptar todo esto para usar write_func_cluster_report
        # JRP This will get us one day
        norm_counts <- hunter_results[["all_data_normalized"]][["DESeq2"]]
        scaled_counts <- scale_data_matrix(data_matrix = as.matrix(norm_counts))
        scaled_counts_table <- as.data.frame(as.table(scaled_counts))
        colnames(scaled_counts_table) <- c("Gene","Sample","Count")

        message("\tRendering individual cluster reports")
        if(is.null(enrichments_ORA)) {
          message("No WGCNA ORA results, not printing individual cluster report")
        } else {
        cls  <- unique(DEGH_results$Cluster_ID)
        cls <- cls[cls != 0]
        trait_module <- hunter_results$WGCNA_all$plot_objects$trait_and_module
        cl_eigvalues <- as.matrix(trait_module[,grepl("^ME",
                                            colnames(trait_module))])
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
        c_update <- list(norm_counts = norm_counts,
                         scaled_counts_table = scaled_counts_table, cls = cls,
                         cl_eigvalues = cl_eigvalues,
                         cl_eigvalues_gnorm = cl_eigvalues_gnorm,
                         wgcna_pval_cl_trait = wgcna_pval_cl_trait,
                         wgcna_corr_cl_trait = wgcna_corr_cl_trait,
                         wgcna_count_sample_trait = wgcna_count_sample_trait)
        container <- c(container, c_update)
        template <- file.path(template_folder, "cl_func_report.txt")
        parallel_list(cls, function(cl) {
        #lapply(cls, function(cl) {
            cl_flags_ora <- lapply(enrichments_ORA_expanded, function(x) {
              nrow(x[[which(names(x) == cl)]]) != 0 
            })
            temp_path_cl <- file.path(results_path, paste0(cl,"_temp_cl_rep"))
            outf_cls_i <- file.path(results_path, paste0("cl_func_",cl,".html"))
            DEGH_subset <- DEGH_results[which(DEGH_results$Cluster_ID == cl), ]
            container$DEGH_results <- DEGH_subset
            container$cl <- cl
            plotter <- htmlreportR::htmlReport$new(title_doc = "functional report",
                                                container = container,
                                                tmp_folder = temp_path_cl,
                                                src = source_folder,
                                                compress_obj = TRUE,
                                                type_index = "contents_list")
            plotter$hash_vars$cl_flags_ora <- cl_flags_ora
            plotter$build(template)
            plotter$write_report(outf_cls_i)
            message(paste0("Report written in ", outf_cls_i))
        #}) 
        }, workers=cores, task_size=task_size)
        unlink(list.files(results_path, pattern="_temp_cl_rep$", 
                        full.names=TRUE), recursive=TRUE) 
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
  gene_attributes, workers, task_size, template_folder, max_genes = 200,
  gene_attribute_name="fold change", source_folder = NULL){
  clean_tmpfiles_mod <- function() {
    message("Calling clean_tmpfiles_mod()")
  }
  if(is.null(source_folder)) {
    source_folder <- find.package("htmlreportR")
  }
  #assignInNamespace("clean_tmpfiles", clean_tmpfiles_mod, ns = "rmarkdown")
  parallel_list(names(enrichments_for_reports), function(cluster) {
    temp_path_cl <- file.path(output_path, paste0(cluster, "_temp"))
    func_results <- enrichments_for_reports[[cluster]]
    cl_flags_ora <- sapply(func_results, nrow) > 0
    attr_vector <- gene_attributes[[cluster]]
    outfile <- file.path(output_path, paste0(cluster, "_to_enrichment",
                                             "_report.html"))
    container <- list(func_results = func_results, cl_flags_ora = cl_flags_ora,
                      max_genes = max_genes)
    message("\tRendering regular report")
    template <- file.path(template_folder, "clusters_to_enrichment.txt")
    plotter <- htmlreportR::htmlReport$new(title_doc = "func cluster",
                                            container = container,
                                            tmp_folder = temp_path_cl,
                                            src = source_folder,
                                            compress_obj = TRUE,
                                            type_index = "contents_list")
    plotter$build(template)
    plotter$write_report(outfile)
    message(paste0("Report written in ", outfile))
  }, workers=workers, task_size=task_size)
  # temp files not deleted properly in parallel 
  unlink(list.files(output_path, pattern="_temp$", full.names=TRUE), 
    recursive=TRUE)
}


render_multivar_report <- function(multivar_res, output_files, template_folder, string_factors = NULL, numeric_factors = NULL,
                                   opt) {
  source_folder <- find.package("htmlreportR")
  if( Sys.getenv('HTMLREPORTER_MODE') == 'DEVELOPMENT' )
    source_folder <- file.path(source_folder, "inst")

  plotter <- htmlreportR::htmlReport$new(title_doc = "PCA report", 
                          container = multivar_res, 
                          tmp_folder = file.path(normalizePath(output_files), "tmp"),
                          src = source_folder,
                          compress_obj = FALSE,
                          type_index = "menu")
  
  plotter$build(file.path(template_folder, 'multivar_main.txt'))
  plotter$write_report(file.path(opt$output_files, "PCA_report.html"))
}