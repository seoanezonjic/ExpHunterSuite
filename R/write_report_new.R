

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
write_functional_report_new <- function(hunter_results, 
                                    func_results, 
                                    output_files=getwd(), 
                                    fc_colname="mean_logFCs", 
                                    organisms_table=NULL, 
                                    template_folder = file.path(find.package('ExpHunterSuite'), 'templates'),
                                    cores = 2,
                                    task_size = 1, 
                                    report = "fc"){

    if(!any(grepl("WGCNA", names(func_results_new))) && grepl("c|i", report)) {
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

    if(grepl("f", report)){
        message("\tRendering regular report")

        DEGH_res_report <- DEGH_results
        outf <- file.path(results_path, "functional_report_new.html")
        rmarkdown::render(file.path(template_folder, 'functional_report_new.Rmd'), 
            output_file = outf, intermediates_dir = results_path)        
    }


    if(grepl("c", report)){

        message("\tRendering full cluster reports")

        if(is.null(func_results$WGCNA_ORA)) {
            message("No WGCNA ORA results, not printing cluster report")
        } else {
            flags_cluster <- sapply(func_results$WGCNA_ORA, function(x) nrow(x@compareClusterResult)) != 0
            names(flags_cluster) <- names(func_results$WGCNA_ORA)
            enrichments_ORA_expanded <- func_results$WGCNA_ORA_expanded
            norm_counts <- hunter_results[["all_data_normalized"]][["DESeq2"]]
            scaled_counts <- scale_data_matrix(data_matrix = as.matrix(norm_counts))
            scaled_counts_table <- as.data.frame(as.table(scaled_counts))
            colnames(scaled_counts_table) <- c("Gene","Sample","Count")

            outf_cls <- file.path(results_path, "clusters_func_report_new.html")
            rmarkdown::render(file.path(template_folder, 
              'clusters_main_report_new.Rmd'),output_file = outf_cls, 
              intermediates_dir = results_path, clean=FALSE)
        }
    }

    if(grepl("i", report)) {
      message("\tRendering individual cluster reports")
      if(is.null(func_results$WGCNA_ORA)) {
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

      #if(any(grepl("WGCNA",names(func_results)))) { # Clustered
      message("Rendering specific cluster reports")
      for(cl in cls) {
        cl_flags_ora <- parallel_list(enrichments_ORA_expanded, function(x) {
          nrow(x[[which(names(x) == cl)]]) != 0 
        }, workers = cores, task_size = task_size)
        outf_cls_i <- file.path(results_path, paste0("cl_func_",cl,".html"))
        DEGH_res_report <- DEGH_results[which(DEGH_results$Cluster_ID == cl), ]
        rmarkdown::render(file.path(template_folder, 
          'cl_func_report_new.Rmd'), output_file = outf_cls_i, 
          intermediates_dir = file.path(results_temp, cl))
      }
      }
    }

}