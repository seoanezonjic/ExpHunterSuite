#!/usr/bin/env Rscript

options(warn=1)
if( Sys.getenv('DEGHUNTER_MODE') == 'DEVELOPMENT' ){
  # Obtain this script directory
  full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile), 
                 error=function(e) # works when using R CMD
                normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', 
                  commandArgs())], '='))[2]))
  main_path_script <- dirname(full.fpath)
  root_path <- file.path(main_path_script, '..', '..')
  # Load custom libraries
  custom_libraries <- c('general_functions.R', "io_handling.R",
    'functional_analysis_library.R', 'functional_analysis_library_new.R', 'plotting_functions.R', "clusters_to_enrichments_functions.R")
  for (lib in custom_libraries){
    source(file.path(root_path, 'R', lib))
  }
  template_folder <- file.path(root_path, 'inst', 'templates')
  organisms_table_file <- file.path(root_path, "inst", "external_data", 
        "organism_table.txt")
} else {
  require('ExpHunterSuite')
  root_path <- find.package('ExpHunterSuite')
  template_folder <- file.path(root_path, 'templates')
  organisms_table_file <- file.path(root_path, "external_data", 
        "organism_table.txt")
}

col <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
         "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
         "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
         "#8A7C64", "#599861")



########################## OPTIONS
option_list <- list(
  optparse::make_option(c("-i", "--input_file"), type="character", default=NULL,
                        help="2 columns - cluster and comma separated gene ids"),
  optparse::make_option(c("-w", "--workers"), type="integer", default=1,
                        help="number of processes for parallel execution. Default=%default"),
  optparse::make_option(c("-p", "--pvalcutoff"), type="double", default=0.05,
                        help="Cutoff for P value and adjusted P value for enrichments. Default=%default"),
  optparse::make_option(c("-q", "--qvalcutoff"), type="double", default=0.02,
                        help="Cutoff for Q value for enrichments. Default=%default"),
  optparse::make_option(c("-t", "--task_size"), type="integer", default=1,
                        help="number of clusters per task. Default=%default"),
  optparse::make_option(c("-F", "--force"), type="logical", default=FALSE, 
                        action = "store_true", help="Ignore temporal files"),
  optparse::make_option(c("-f", "--funsys"), type="character", default="BP,MF,CC", 
                        help="Funsys to execute: MF => GO Molecular Function, BP => GO Biological Process, CC => GO Celular Coponent, KEGG => KEGG. Default=%default"),
  optparse::make_option(c("-K", "--kegg_data_file"), ,type = "character", default=NULL,
                        help=paste0("KEGG database file. Can download with download_latest_kegg_db().",
                          "If required but not provided, it will be downloaded to the current working directory")), 
  optparse::make_option(c("--custom"), type="character", default=NULL,
                        help="Comma separated path of custom enrichment sets."),
  optparse::make_option(c("--showCategories"), type="integer", default=30, 
                        help="Number of top categories to show on clusterProfiler dorplot and emaplot"),
  optparse::make_option(c("-c", "--clean_parentals"), type="logical", default=FALSE, 
                        action = "store_true", help="Clean parentals GO terms that appears on the same clusters than child."),
  optparse::make_option(c("-s", "--simplify"), type="logical", default=FALSE, 
                        action = "store_true", help="Apply simplify function from cluster profiler to enrichment."),
  optparse::make_option(c("-O", "--model_organism"), type="character", default="Human", 
                        help="Model organism. Human or Mouse"),
  optparse::make_option(c("-M", "--mode"), type="character", default="PR", 
                        help="String indicating report modes to execute. R: Generate report for each cluster and all clusters combined. P: Plots (dotplot and emaplot) combining all cluster enrichments. S: Summary mode (sumamrized categories heatmap)"),
  optparse::make_option(c("-S", "--sim_thr"), type="double", default=0.7,
                        help="Similarity cutoff for grouping categories in Summary mode. Default=%default"),
  optparse::make_option(c("-C", "--common_name"), type="character", default="ancestor", 
                        help="Name of the term groups. 'significant' to use the most significant term of each group. 'ancestor' to use the common ancestor of the group"),
  optparse::make_option(c("-T", "--top_categories"), type="integer", default=50,
                        help="Number of top categories for each cluster. Default=%default"),
  optparse::make_option(c("-d", "--description_file"), type="character", default=NULL,
                        help="Markdown file describing of the enriched clusters."),
  optparse::make_option(c("-k", "--gene_keytype"), type="character", default="ENTREZID",
                        help="What identifier is being used for the genes in the clusters?. Default=%default"),
  optparse::make_option(c("--gmt_id"), type="character", default="ENTREZID",
                        help="What identifier is being used for the genes in the custom gmt file. Default=%default"),
  optparse::make_option(c("-g", "--gene_mapping"), type="character", default=NULL,
                        help="3 columns tabular file- Cluster - InputGeneID - NumericGeneMapping. Header must be indicated as cluster - geneid - [numeric_mapping]"),
  optparse::make_option(c("-G", "--group_results"), type="logical", default=FALSE, 
                        action = "store_true", help="Functions are gropuped in most frequent words in emaplots."),
  optparse::make_option(c("-o", "--output_file"), type="character", default="results",
                        help="Define the output path.")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

##################################### INITIALIZE ##
library(clusterProfiler)
output_path <- paste0(opt$output_file, "_functional_enrichment")
dir.create(output_path)
output_path <- normalizePath(output_path)
temp_file <- file.path(output_path, "enr_tmp.RData")
all_funsys <- unlist(strsplit(opt$funsys, ","))
n_category <- opt$showCategories
organisms_table <- get_organism_table(organisms_table_file)
current_organism_info <- organisms_table[rownames(organisms_table) %in% opt$model_organism,]
gene_mapping <- NULL
org_db <- get_org_db(current_organism_info)
# Load customs
all_custom_gmt <- NULL
if (!is.null(opt$custom)) {
    custom_files <- unlist(strsplit(opt$custom, ","))
    all_custom_gmt <- lapply(custom_files, load_and_parse_gmt)
    names(all_custom_gmt) <- custom_files
  names(all_custom_gmt) <- basename(names(all_custom_gmt))

  if(opt$gmt_id != "ENTREZID") {
    all_custom_gmt <- lapply(all_custom_gmt, function(gmt){
        tr_gmt <- translate_gmt(gmt, opt$gmt_id, org_db)
        return(tr_gmt)
      })
    
  }
}

if("KEGG" %in% all_funsys) {
  kegg_data_file <- get_kegg_db_path(opt$kegg_data_file, current_organism_info=current_organism_info)
  if(! file.exists(kegg_data_file)) stop(paste("KEGG file:", kegg_data_file, "not found"))
}

################################### MAIN ##

if (!is.null(opt$gene_mapping)){
  parsed_mappings <- parse_mappings(opt$gene_mapping, org_db, opt$gene_keytype)
  gene_mapping <- parsed_mappings[["gene_mapping"]]
  gane_mapping_name <-  parsed_mappings[["gane_mapping_name"]]
}

if (!file.exists(temp_file) || opt$force) {

  cluster_genes <- read.table(opt$input_file, header=FALSE, sep="\t")
  cluster_genes_list <- strsplit(cluster_genes[,2], ",")
  if(opt$gene_keytype != "ENTREZID") {
    cluster_genes_list <- lapply(cluster_genes_list, function(x){
                                      convert_ids_to_entrez(ids=x, 
                                                            gene_keytype=opt$gene_keytype,
                                                            org_db = org_db)}) 
  }

  names(cluster_genes_list) <- cluster_genes[,1]
  enrichments_ORA <- multienricher_ora(all_funsys =  all_funsys, 
                                  genes_list =  cluster_genes_list, 
                                  task_size = opt$task_size, 
                                  organism_info= current_organism_info,
                                  workers = opt$workers, 
                                  pvalueCutoff =  opt$pvalcutoff, 
                                  qvalueCutoff = opt$qvalcutoff,
                                  custom = all_custom_gmt,
                                  kegg_file = kegg_data_file)


  enrichments_ORA_merged <- parse_cluster_results(enrichments_ORA, simplify_results = opt$simplify, clean_parentals = opt$clean_parentals)
  save(enrichments_ORA,enrichments_ORA_merged, file = temp_file)

} else {
  load(temp_file)
}

if (!is.null(opt$custom))
  all_funsys <- c(all_funsys, names(all_custom_gmt))


if (grepl("R", opt$mode)){
    enrichments_for_reports <- parse_results_for_report(enrichments_ORA)
    write_fun_enrichments(enrichments_ORA, output_path, all_funsys)
    write_func_cluster_report(enrichments_for_reports,output_path,gene_mapping, workers = opt$workers, task_size = opt$task_size)
}

if (grepl("P", opt$mode) || grepl("S", opt$mode)) 
  enrichments_ORA_merged <- filter_top_categories(enrichments_ORA_merged, opt$top_categories)

if (grepl("P", opt$mode)) {
  for (funsys in names(enrichments_ORA_merged)){
    if (length(unique(enrichments_ORA_merged[[funsys]]@compareClusterResult$Description)) < 2 ) next

    if (opt$group_results == TRUE){
      pp <- enrichplot::emapplot(enrichments_ORA_merged[[funsys]], showCategory= n_category, pie="Count", layout = "nicely", 
                  shadowtext = FALSE, node_label = "group", group_category = TRUE, 
                  nCluster = min(floor(nrow(enrichments_ORA_merged[[funsys]])/7), 20), nWords = 6, repel = TRUE)
    }else{
      pp <- enrichplot::emapplot(enrichments_ORA_merged[[funsys]], showCategory= n_category, pie="Count", layout = "nicely", 
                  shadowtext = FALSE, repel = TRUE)
    }

    ggplot2::ggsave(filename = file.path(output_path,paste0("emaplot_",funsys,"_",opt$output_file,".png")), pp, width = 30, height = 30, dpi = 300, units = "cm", device='png')

    pp <- enrichplot::dotplot(enrichments_ORA_merged[[funsys]], showCategory= n_category, label_format = 70)
    ggplot2::ggsave(filename = file.path(output_path,paste0("dotplot_",funsys,"_",opt$output_file,".png")), pp, width = 60, height = 40, dpi = 300, units = "cm", device='png')

  }
}




if (grepl("S", opt$mode)){

# So that it doesn't break with non GO categories
enrichments_ORA_merged <- enrichments_ORA_merged[names(enrichments_ORA_merged) %in% c("BP","MF","CC")]

  print("time summarize categories")

  print(system.time(
  sum_enrichments <- summarize_categories(enrichments_ORA_merged, sim_thr = opt$sim_thr,common_name = opt$common_name)
 ))
        save(list = ls(all.names = TRUE), file = "environment_to_summarize.RData")

  for (funsys in names(sum_enrichments)) {

    sum_table <- sum_enrichments[[funsys]]
    sum_table <- (sum_table > opt$pvalcutoff) + 0
    sum_table <- sum_table[rownames(sum_table) != "to_remove",]
    sum_table_to_plot <- sum_table
    rownames(sum_table_to_plot) <- get_GOid_term(rownames(sum_table_to_plot))
    heatmaply::heatmaply(sum_table_to_plot, 
                         grid_color = "gray50",
                         seriate = "mean",
                         grid_width = 0.00001,
                         dendrogram = "both",
                         scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
                         low = "#EE8291", 
                         high = "white", 
                         midpoint = 0.5, 
                         limits = c(0, 1)),
                         file = file.path(output_path, paste0("sum_",funsys,'_heatmap.html')))
 
 print("time clean parentals")
print(system.time(
    sum_table_clean <-  clean_parentals_in_matrix(sum_table, funsys)
    ))
print("time get goid")
print(system.time(

    rownames(sum_table_clean) <- get_GOid_term(rownames(sum_table_clean))
))
    sum_table_clean <- sum_table_clean[rowSums(sum_table_clean < opt$pvalcutoff) != 0,]

    heatmaply::heatmaply(sum_table_clean, 
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
                         file = file.path(output_path, paste0("sum_clean_",funsys,'_heatmap.html')))


print("cluster enr to matrix")
print(system.time(
    all_enrichments <- cluster_enr_to_matrix(enrichments_ORA_merged[[funsys]]@compareClusterResult)
))
    all_enrichments <- (all_enrichments > opt$pvalcutoff) + 0 

    heatmaply::heatmaply(all_enrichments, 
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

