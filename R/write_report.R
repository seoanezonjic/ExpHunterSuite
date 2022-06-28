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
                                    report = "fc"){
    # report <- "i"
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
       !is.na(func_results$DEGH_results_annot$entrezgene), fc_colname]
    names(fc_vector) <- func_results$DEGH_results_annot[
       !is.na(func_results$DEGH_results_annot$entrezgene), "entrezgene"]

    enrichments_ORA <- func_results$WGCNA_ORA
    DEGH_results <- func_results$DEGH_results_annot
    enrichments_ORA_expanded <- func_results$WGCNA_ORA_expanded

    # JRP This will get us one day
    norm_counts <- hunter_results[["all_data_normalized"]][["DESeq2"]]
    scaled_counts <- scale_data_matrix(data_matrix = as.matrix(norm_counts))
    scaled_counts_table <- as.data.frame(as.table(scaled_counts))
    colnames(scaled_counts_table) <- c("Gene","Sample","Count")

    rm("func_results")

    if(grepl("f", report)){
        message("\tRendering regular report")
        outf <- file.path(results_path, "functional_report.html")
        rmarkdown::render(file.path(template_folder, 'functional_report.Rmd'), 
            output_file = outf, intermediates_dir = results_path)        
    }

    if(grepl("c", report)){
        message("\tRendering full cluster reports")
        if(is.null(enrichments_ORA)) {
            message("No WGCNA ORA results, not printing cluster report")
        } else {
            flags_cluster <- sapply(enrichments_ORA, function(x) nrow(x@compareClusterResult)) != 0
            names(flags_cluster) <- names(enrichments_ORA)
            outf_cls <- file.path(results_path, "clusters_func_report.html")
            rmarkdown::render(file.path(template_folder, 
              'clusters_main_report.Rmd'),output_file = outf_cls, 
              intermediates_dir = results_path)
        }
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

        results_temp <- file.path(paste0(results_path, "_tmp"))
        check_and_create_dir(results_temp)

        message("\tRendering specific cluster reports")
        parallel_list(cls, function(cl) {
        #lapply(cls, function(cl) {
            cl_flags_ora <- lapply(enrichments_ORA_expanded, function(x) {
              nrow(x[[which(names(x) == cl)]]) != 0 
            })
            temp_path_cl <- file.path(results_path, paste0(cl,"_temp_cl_rep"))
            print("temp path:")
            print(temp_path_cl)
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