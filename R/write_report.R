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
