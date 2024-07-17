#' Function to generate targets file
#' Load target file if it exists, otherwise use the -C and -T flags. 
#' Note target takes precedence over target.
#' @param from_file input targets file
#' @param ctrl_samples control samples
#' @param treat_samples target samples
#' @return targets info loaded
#' @keywords design
#' @export
#' @importFrom utils read.table
#' @examples
#' # To generate target or load from file. 
#' # Example shows generation based on names of treament/control samples.
#' # These should match the column names in the table of counts.
#' case_sample_names <- paste("case", seq(5), sep="_")
#' control_sample_names <- paste("control", seq(5), sep="_")
#' target <- target_generation(ctrl_samples=control_sample_names, 
#'                             treat_samples=case_sample_names)
target_generation <- function(from_file=NULL, 
                              ctrl_samples=NULL, 
                              treat_samples=NULL){
    target <- NULL
    if(! is.null(from_file)) {
      target <- utils::read.table(from_file, 
                                    header=TRUE, 
                                    sep="\t", 
                                    check.names = FALSE)
      # Check there is a column named treat
      if(! "treat" %in% colnames(target)) {
        stop(cat("No column named treat in the target file.\nPlease resubmit"))
      }
    } else {
      index_control_cols <- unlist(strsplit(ctrl_samples, ",")) 
      index_treatmn_cols <- unlist(strsplit(treat_samples, ","))
      # Create the target data frame needed in calls to the DE detection methods
      target <- data.frame(sample=c(index_control_cols, index_treatmn_cols), 
        treat=c(rep("Ctrl", length(index_control_cols)), 
                rep("Treat", length(index_treatmn_cols))))
    }    
    return(target)
}    


#' Function to write expression package results.
#' Write to disk the results of each diff expression package
#' @param df_list dataframes list 
#' @param prefix prefix concatenated to output file name
#' @param root path where write
#' @return void
#' @keywords output
#' @export
#' @importFrom utils write.table
#' @examples
#' # Typically used with the @all_counts_for_plotting slot of the main output
#' deg_res <- list(aa=data.frame(a=seq(5),b=seq(5)), 
#'                 bb=data.frame(a=seq(6,10),b=seq(6,10)))
#' write_df_list_as_tables(deg_res, prefix="test", root=tempdir())
write_df_list_as_tables <- function(df_list, prefix, root=getwd()) {
  for(pack in names(df_list)) {
      folder <- file.path(root, paste0("Results_", pack))
      if(!dir.exists(folder)) dir.create(folder)
      utils::write.table(df_list[[pack]],
      file=file.path(folder, paste0(prefix, pack, '.txt')), 
      quote=FALSE, col.names = NA, sep="\t")
  }
}

#' @importFrom utils tail
load_hunter_folder <- function(path = NULL){
    packages_results <- list.dirs(path,recursive=FALSE)
    packages_results <- packages_results[grepl("Results_",packages_results)]
    packages_results <- packages_results[!grepl("WGCNA$",packages_results)]

    dgh_exp_results <- list()
    dgh_exp_results[["DE_all_genes"]] <- utils::read.table(
        file.path(path,"Common_results","hunter_results_table.txt"))

    dgh_exp_results[["sample_groups"]] <- utils::read.table(
        file.path(path,"control_treatment.txt"))

    WGCNA_res <- load_WGCNA_results(path,dgh_exp_results[["DE_all_genes"]])
    if(length(WGCNA_res) > 0) dgh_exp_results[["WGCNA_all"]] <- WGCNA_res

    dgh_exp_results[["all_data_normalized"]] <- lapply(packages_results, function(pack_path) {
        #pack_name <- utils::tail(unlist(strsplit(pack_path,"_")),1)
        exp_matrix_path <- grep("Norm", dir(file.path(pack_path)), value=TRUE)
        exp_matrix <- utils::read.table(file.path(pack_path, exp_matrix_path))
        return(exp_matrix)
    })
    names(dgh_exp_results[["all_data_normalized"]]) <- sapply(strsplit(packages_results,"_"), utils::tail, 1)

    dgh_exp_results$pca_all_genes <- read.table(file.path(path, "PCA_results/all_genes_eigenvectors.txt"), 
                                                row.names = 1, header = TRUE)
    dgh_exp_results$pca_degs<- read.table(file.path(path, "PCA_results/prevalent_eigenvectors.txt"), 
                                                row.names = 1, header = TRUE)

    return(dgh_exp_results)
}

#' Loads stored WGCNA results in DEgenes Hunter expression analysis package 
#' folder and returns in compact object
#' @param path of expression analysis results (Hunter folder)
#' @param main_deg_table dataframe with main DEG analysis results
#' @return WGCNA results loaded
#' @keywords input
#' @importFrom utils read.table
load_WGCNA_results <- function(path, main_deg_table){
    info <- list()
    if(dir.exists(file.path(path,"Results_WGCNA"))){
        # Package objects
        package_objects <- list()
        package_objects[["gene_module_cor"]] <- utils::read.table(
            file.path(path,"Results_WGCNA","gene_MM.txt"))
        package_objects[["gene_module_cor_p"]] <- utils::read.table(
            file.path(path,"Results_WGCNA","gene_MM_p_val.txt"))
        package_objects[["gene_trait_cor"]] <- utils::read.table(
            file.path(path,"Results_WGCNA","gene_trait.txt"))
        package_objects[["gene_trait_cor_p"]] <- utils::read.table(
            file.path(path,"Results_WGCNA","gene_trait_p_val.txt"))
        package_objects[["module_trait_cor"]] <- utils::read.table(
            file.path(path,"Results_WGCNA","module_trait.txt"))
        package_objects[["module_trait_cor_p"]] <- utils::read.table(
           file.path(path,"Results_WGCNA","module_trait_p_val.txt"))
        # gene_cluster_info
        gene_cluster_info <- cbind(rownames(main_deg_table), 
            main_deg_table[,c("Cluster_ID", "Cluster_MM", "Cluster_MM_pval")])
        colnames(gene_cluster_info) <- c("ENSEMBL_ID", "Cluster_ID", 
                                         "Cluster_MM", "Cluster_MM_pval")
        gene_cluster_info <- gene_cluster_info[rownames(
                                        package_objects[["gene_module_cor"]]),]
        # plot_objects
        plot_objects <- list()

        sample_module <- utils::read.table(
            file.path(path,"Results_WGCNA","eigen_values_per_samples.txt"))
        sample_trait <- utils::read.table(
            file.path(path,"Results_WGCNA","sample_trait.txt"))
        plot_objects[["trait_and_module"]] <- cbind(sample_module,sample_trait)

        # Store into info
        info$package_objects <- package_objects
        info$gene_cluster_info <- gene_cluster_info
        info$plot_objects <- plot_objects
    }

    return(info)
}

#' get path for KEGG db for downloading/usage
#' @param kegg_data_file file path. if null, inst/kegg_data_files will be used
#' @param current_organism_info information for the organism, including KEGG code
get_kegg_db_path <- function(kegg_data_file, current_organism_info){
  if(is.null(kegg_data_file)) {
    kegg_code <- current_organism_info$KeggCode[1]
    kegg_data_file <- paste0(kegg_code, "_KEGG.rds")

    if( Sys.getenv('DEGHUNTER_MODE') == 'DEVELOPMENT' ) {
      # Root path must be defined outside function i.e. in script that calls it
      kegg_data_file <- file.path(root_path, 'inst', 'kegg_data_files', kegg_data_file)
    }
    else { 
      kegg_data_file <- system.file("kegg_data_files", kegg_data_file, package="ExpHunterSuite")
    }
  }
  return(kegg_data_file)
}

#' download kegg db for a given organism
#' @param current_organism_info organism info for which to download the file
#' @param file where to save the file
download_latest_kegg_db <- function(current_organism_info, file) {
  organism <- current_organism_info$KeggCode[1]
  prepare_KEGG <- get("prepare_KEGG", envir=asNamespace("clusterProfiler"), inherits = FALSE)
  ENRICH_DATA <- prepare_KEGG(organism, "KEGG", "kegg")
  saveRDS(ENRICH_DATA, file=file)
}

write_table_ehs <- function(x, file) {
    utils::write.table(x=x, file=file, quote=FALSE, col.names=TRUE, 
        row.names = FALSE, sep="\t")
}

write_enrich_tables <- function(func_res_tables, method_type, output_files){
    for(res in names(func_res_tables)) {
        filename <- file.path(output_files, 
                paste(res, method_type, "results", sep="_"))
        res_to_print <- func_res_tables[[res]]
        write_table_ehs(res_to_print, file=filename)
    }
}
write_enrich_clusters <- function(func_clusters, output_path) {
    lapply(names(func_clusters), function(funsys) {
      func_res <- func_clusters[[funsys]]
      merged_func_res <- clusterProfiler::merge_result(func_res)
      write_table_ehs(merged_func_res,
        file=file.path(output_path, paste0("enr_cls_", funsys, ".csv")))
      # write_enrich_tables(func_res, 
      #   paste0(funsys, "_cluster"), output_path)
    })
}

merge_and_write_tables <- function(tables_list, output_path, idcol = "group") {
    merged_table <- data.table::rbindlist(tables_list, fill=FALSE, idcol= idcol)
    write.table(merged_table, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)
}

write_pca_enrichments <- function(pca_enrichments, output_path) {
    for (pca_type in names(pca_enrichments)) {
        enrichments <- pca_enrichments[[pca_type]]
        for (funsys in names(enrichments)) {
            file_name <- file.path(output_path, paste(pca_type, funsys, sep = "_"))
            merge_and_write_tables(enrichments[[funsys]], paste0(file_name, ".txt"), idcol = "PC")
        }
    } 
}

#' Write enrichment files related to functional_hunter results list
#' @param func_results functional enrichment results
#' @param output_path output folder path
#' @return void
#' @keywords export
#' @export
#' @importFrom utils write.table
#' @examples
#' # Use a @functional_hunter or @multienricher result object to 
#' # create real files
#' write_enrich_files(list(),"./")
write_enrich_files <- function(func_results, output_path=getwd()) {
    if(!dir.exists(output_path)) dir.create(output_path)
    final_params <- func_results$final_main_params[! names(func_results$final_main_params) %in%
        c("hunter_results", "organisms_table", "annot_table", "custom")]
    write_table_ehs(data.frame(A = names(final_params), 
                               B = sapply(final_params, paste, collapse=" ")), 
                       file = file.path(output_path,"functional_opt.txt"))
    write_pca_enrichments(func_results$PCA_enrichments, output_path)
    if("ORA" %in% names(func_results)) {
        write_enrich_tables(func_results$ORA, "ORA", output_path)
    }
    if("GSEA" %in% names(func_results)) {
        write_enrich_tables(func_results$GSEA, "GSEA", output_path)
    }
    if("DEGH_results_annot" %in% names(func_results)) 
        write_table_ehs(func_results$DEGH_results_annot, 
                           file=file.path(output_path, 
                                          "hunter_results_table_annotated.txt"))
    if("WGCNA_ORA" %in% names(func_results)) {
        write_enrich_clusters(func_results$WGCNA_ORA_expanded, output_path)
        # func_res_cls <- lapply(func_results$WGCNA_ORA, as.data.frame)
        # write_enrich_tables(func_res_cls, "cls_ORA", output_path)

        # for(cl in unique(func_results$DEGH_results_annot$Cluster_ID)) {
        #     DEGH_res_cl <- func_results$DEGH_results_annot[func_results$DEGH_results_annot$Cluster_ID == cl,]              
        #     write_table_ehs(DEGH_res_cl[DEGH_res_cl$genes_tag == "PREVALENT_DEG", 
        #                         c("Symbol", "entrezgene", "mean_logFCs", "combined_FDR")],
        #                     file=file.path(output_files, paste0("cluster_", cl, "_DEGs.txt")))
        # }
    }
}



#' Loads Hunter DEG analysis table and simulated TRUE DEG vector. Concat in
#' training/testing ML data frame format and returns it.
#' @param folder where files are placed
#' @param hunter_table_file DEG analys result table file name
#' @param prediction_file_regex regex to find TRUE DEG file generated
#' @return a data frame with DEG result and real predictions column
#' @export
#' @examples
#' # Used to load a synth + DE analysis result folder
#' ht <- load_synth_dataset() 
#' # Will return an empty object, you need a correct folder
load_synth_dataset <- function(
    folder = NULL,
    hunter_table_file = "hunter_results_table.txt",
    prediction_file_regex = "_predv"){
    # Check
    if(is.null(folder)){
      warning("Folder has not been specified. Returning NULL")
      return(NULL)
    }
    # Load hunter table
    ht <- read.table(file = file.path(folder,hunter_table_file), sep = "\t",
                     header = TRUE)
    # Load prediction
    prediction_file <- list.files(folder)
    prediction_file <- prediction_file[which(grepl(prediction_file_regex,
                                                   prediction_file))]
    if(length(prediction_file) > 1) prediction_file <- prediction_file[1]
    true_preds <- read.table(file = file.path(folder,prediction_file),
                             sep = "\t", header = TRUE)
    # Concat
    ht$Prediction <- unlist(lapply(seq(nrow(ht)),function(i){
        return(true_preds[which(true_preds[,1] == rownames(ht)[i])[1],2])
    }))
    # Return
    return(ht)
}

