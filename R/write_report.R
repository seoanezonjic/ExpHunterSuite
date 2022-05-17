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
#'  and functional (f). Default = "fc"
#' @return void
#' @importFrom rmarkdown render
#' @export
#' @examples
#' # Load func and DE results
#' data(degh_output)
#' func_results <- list() 
#' # func_results <- functional_hunter(degh_output,"Mouse")
#' write_functional_report(degh_output, func_results)
write_functional_report <- function(hunter_results, 
                                    func_results, 
                                    output_files=getwd(), 
                                    fc_colname="mean_logFCs", 
       organisms_table=NULL, 
       template_folder = file.path(find.package('ExpHunterSuite'), 'templates'),
                                    cores = 2,
                                    task_size = 1, 
                                    report = "fc"){
    if(is.null(organisms_table)){
        organisms_table <- get_organism_table()
    }
    if(length(hunter_results) == 0 || length(func_results) == 0){
        warning("Results objects are not complete")
        return(NULL)
    }
    model_organism <- func_results$final_main_params$model_organism
    # TODO: update names into Rmd files instead this
    ############################################################
    ##               CREATE NECESSARY VARIABLES               ##
    ############################################################
    degh_exp_threshold <- hunter_results$final_main_params$p_val_cutoff
    DEGH_results <- func_results$DEGH_results_annot
    # -
    experiments <- hunter_results$sample_groups
    sample_classes <- apply(experiments, 1, function(x) paste0("* [", x[1],
                      "] ", x[2]))
    # -
    if(! "externalDEA" %in% names(hunter_results[["all_data_normalized"]])) {
        norm_counts <- hunter_results[["all_data_normalized"]][["DESeq2"]]
        scaled_counts <- scale_data_matrix(data_matrix = as.matrix(norm_counts))
        scaled_counts_table <- as.data.frame(as.table(scaled_counts))
        colnames(scaled_counts_table) <- c("Gene","Sample","Count")
    }

    # -
    flags <- func_results$flags
    if(flags$WGCNA){
        cls  <- unique(DEGH_results$Cluster_ID)
        aux <- hunter_results$WGCNA_all$plot_objects$trait_and_module
        cl_eigvalues <- as.matrix(aux[,grepl("^ME",colnames(aux))])
        cl_eigvalues <- as.data.frame(as.table(cl_eigvalues),
            stringsAsFactors = FALSE)
        colnames(cl_eigvalues) <- c("Sample","Cluster_ID","Count") 
        cl_eigvalues_gnorm <- cl_eigvalues
        cl_eigvalues_gnorm$Count <- (cl_eigvalues_gnorm$Count + 1) / 2
        aux2 <- hunter_results$WGCNA_all$package_objects
        wgcna_pval_cl_trait <- as.matrix(aux2$module_trait_cor_p)
        wgcna_corr_cl_trait <- as.matrix(aux2$module_trait_cor)
        wgcna_count_sample_trait <- as.matrix(aux[,!grepl("^ME",
            colnames(aux))])
        wgcna_count_sample_trait <- scale_data_matrix(wgcna_count_sample_trait, 
            norm_by_col = TRUE)
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
    current_organism_info <- subset(organisms_table, 
        rownames(organisms_table) %in% model_organism)  
    geneList <- func_results$DEGH_results_annot[
       !is.na(func_results$DEGH_results_annot$entrezgene), fc_colname]
    names(geneList) <- func_results$DEGH_results_annot[
       !is.na(func_results$DEGH_results_annot$entrezgene), "entrezgene"]
    geneList <- sort(geneList, decreasing = TRUE)
    # -
    custom_enrichments <- func_results$CUSTOM
    if("ORA" %in% names(func_results)){
        aux <- grepl("GO", names(func_results$ORA))
        if(any(aux)) enrich_go <- func_results$ORA[aux]
        if("KEGG" %in% names(func_results$ORA)) 
            enrich_ora <- func_results$ORA$KEGG
        if("REACT" %in% names(func_results$ORA)) 
            enrich_react <- func_results$ORA$REACT
    }
    if("GSEA" %in% names(func_results)){
        aux <- grepl("GO", names(func_results$GSEA))
        if(any(aux)) enrich_go_gsea <- func_results$GSEA[aux]
        if("KEGG" %in% names(func_results$GSEA)) 
            enrich_gsea <- func_results$GSEA$KEGG
        if("REACT" %in% names(func_results$GSEA)) 
            enrich_react_gsea <- func_results$GSEA$REACT   
    }
    # TODO: topGO is not bein loaded because files are not been search
    ############################################################
    ##                GENERATE CLUSTER REPORTS                ##
    ############################################################
    clean_tmpfiles_mod <- function() {
        message("Calling clean_tmpfiles_mod()")
    }

#https://github.com/rstudio/rmarkdown/issues/1632#issuecomment-545824711
    utils::assignInNamespace("clean_tmpfiles", 
                             clean_tmpfiles_mod, ns = "rmarkdown") 

    results_path <- normalizePath(output_files)
    results_temp <- file.path(paste0(results_path, "_tmp"))
    check_and_create_dir(results_temp)
    if(grepl("c", report)){
        if (any(grepl("WGCNA",names(func_results)))) { # Clustered
            message("Rendering specific cluster reports")
            invisible(parallel_list(cls, function(cl) {
                # Take output name
                aux <- paste0("cl_func_",cl,".html")
                outf_cls_i <- file.path(results_path, aux)
                # Generate report
                rmarkdown::render(file.path(template_folder, 
                    'cl_func_report.Rmd'), output_file = outf_cls_i, 
                intermediates_dir = file.path(results_temp, cl), quiet=TRUE)
            }, workers = cores, task_size= task_size))

            message("\tRendering clustered report")
            outf_cls <- file.path(results_path, "clusters_func_report.html")
            rmarkdown::render(file.path(template_folder, 
                'clusters_main_report.Rmd'),output_file = outf_cls, 
            intermediates_dir = results_path)
        }        
    }
    unlink(results_temp, recursive = TRUE)

    ############################################################
    ##              GENERATE DEG FUNCTIONAL REPORT            ##
    ############################################################
    if(grepl("f", report)){
        message("\tRendering regular report")
        outf <- file.path(results_path, "functional_report.html")
        rmarkdown::render(file.path(template_folder, 'functional_report.Rmd'), 
            output_file = outf, intermediates_dir = results_path)        
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
RNAseq,
sig_pairs,
raw_databases_scores,
p_fisher           
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
    "DEGs_RNA_vs_miRNA_DEMs",
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

write_functional_targets <- function(
    enrich_GO,
    enrich_KEGG,
    enrich_react,
    launch_default,
    launch_expanded,
    output_path,
    strategy,
    enrichments_ORA_expanded,
    enrichments_ORA,
    RNAseq_folder,
    miRNAseq_folder,
    templates_folder,
    nomenclatures,
    current_organism_info,
    geneList,
    enrich_custom = NULL,
    strat,
    unique_miRNAs,
    raw_data,
    cores,
    task_size

){
    message("\tRendering regular report")
    results_temp <- file.path(paste0(output_path, "_tmp"))
    check_and_create_dir(output_path)
    check_and_create_dir(results_temp)
    if (launch_default) {
        outf <- file.path(output_path, paste0(strategy, 
            "_targets_functional.html"))

        rmarkdown::render(file.path(templates_folder, 
            'targets_functional.Rmd'), output_file = outf, 
            intermediates_dir = results_temp)
    }
   
    if (launch_expanded) {
        message("Rendering specific miRNA reports")
     
        RNAseq <- load_DEGH_information(RNAseq_folder)
        miRNAseq <- load_DEGH_information(miRNAseq_folder)
        miRNA_strat <- unlist(strsplit(strat, ""))[2]
        if (miRNA_strat == "h") {
            hub_miRNAs <- get_hub_genes_by_MM(miRNAseq[["normalized_counts"]], 
            miRNAseq[["DH_results"]], top = 1)
        } else if (miRNA_strat == "E") {
            miRNAseq$Eigengene <- as.data.frame(as.table(miRNAseq$Eigengene), 
                               stringsAsFactors = FALSE)
            colnames(miRNAseq$Eigengene) <- c("Sample","Cluster_ID","Count") 
            tgt_eigvalues_gnorm <- miRNAseq$Eigengene
            tgt_eigvalues_gnorm$Count <- (tgt_eigvalues_gnorm$Count + 1) / 2 
        }
        unique_miRNAs <- unique_miRNAs[!is.na(unique_miRNAs)]

        invisible(parallel_list(unique_miRNAs, function(miRNA) {
            # Take output name
            target_outf <- file.path(output_path, paste0("targets_", 
                miRNA,".html"))
            # Generate report
            rmarkdown::render(file.path(templates_folder, 
                'miRNA_target_func.Rmd'), output_file = target_outf, 
            intermediates_dir = file.path(results_temp, miRNA), quiet=TRUE)
            
        }, workers = cores, task_size= task_size))
        message("\tRendering merged miRNA report")
        outf_cls <- file.path(output_path, "expanded_targets_func.html")
        rmarkdown::render(file.path(templates_folder, 
            'targets_global_report.Rmd'),output_file = outf_cls,
             intermediates_dir = results_temp)
    }
    unlink(results_temp, recursive = TRUE)
}