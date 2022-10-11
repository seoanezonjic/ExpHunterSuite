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
    'functional_analysis_library.R', 'plotting_functions.R', "main_clusters_to_enrichment.R",
    "write_report.R")
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

########################## OPTIONS
option_list <- list(
  optparse::make_option(c("-i", "--input_file"), type="character", default=NULL,
                        help="2 columns - cluster and comma separated gene ids"),
  optparse::make_option(c("-w", "--workers"), type="integer", default=1,
                        help="number of processes for parallel execution. Default=%default"),
  optparse::make_option(c("-p", "--pvalcutoff"), type="double", default=0.1,
                        help="Cutoff for P value and adjusted P value for enrichments. Default=%default"),
  optparse::make_option(c("-q", "--qvalcutoff"), type="double", default=0.2,
                        help="Cutoff for Q value for enrichments. Default=%default"),
  optparse::make_option(c("-t", "--task_size"), type="integer", default=1,
                        help="number of clusters per task. Default=%default"),
  optparse::make_option(c("-F", "--force"), type="logical", default=FALSE, 
                        action = "store_true", help="Ignore temporal files"),
  optparse::make_option(c("-f", "--funsys"), type="character", default="BP,MF,CC,KEGG,Reactome", 
                        help="Funsys to execute: MF => GO Molecular Function, BP => GO Biological Process, CC => GO Celular Coponent, 
                        KEGG => KEGG, Reactome => Reactome, DOSE => DOSE, DGN => DisGeNET. Default=%default"),
  optparse::make_option(c("-K", "--kegg_data_file"), ,type = "character", default=NULL,
                        help=paste0("KEGG database file. Can download with download_latest_kegg_db().",
                          "If required but not provided, it will be downloaded to the current working directory")), 
  optparse::make_option(c("--custom"), type="character", default=NULL,
                        help="Comma separated path of custom enrichment sets."),
  optparse::make_option(c("--showCategories"), type="integer", default=30, 
                        help="Number of top categories to show on clusterProfiler dotplot and emaplot"),
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
  optparse::make_option(c("-C", "--summary_common_name"), type="character", default="ancestor", 
                        help="Name of the term groups. 'significant' to use the most significant term of each group. 'ancestor' to use the common ancestor of the group"),
  optparse::make_option(c("-T", "--top_categories"), type="integer", default=50,
                        help="Number of top categories for each cluster. Default=%default"),
  optparse::make_option(c("-d", "--description_file"), type="character", default=NULL,
                        help="Markdown file describing of the enriched clusters."),
  optparse::make_option(c("-k", "--gene_keytype"), type="character", default="ENTREZID",
                        help="What identifier is being used for the genes in the clusters?. Default=%default"),
  optparse::make_option(c("--gmt_id"), type="character", default="ENTREZID",
                        help="What identifier is being used for the genes in the custom gmt file. Default=%default"),
  optparse::make_option(c("-g", "--gene_attribute_file"), type="character", default=NULL,
                        help="3 columns tabular file- Cluster - InputGeneID - NumericAttribute. Header must be indicated as cluster - geneid - [numeric_atribute]"),
  optparse::make_option(c("-G", "--group_results"), type="logical", default=FALSE, 
                        action = "store_true", help="Functions are grouped in most frequent words in emaplots."),
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
organisms_table <- get_organism_table(organisms_table_file)
current_organism_info <- organisms_table[rownames(organisms_table) %in% opt$model_organism,]
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
ce_list <- main_clusters_to_enrichment(
  input_file = opt$input_file,
  gene_attribute_file = opt$gene_attribute_file,
  org_db = org_db,
  gene_keytype = opt$gene_keytype,
  temp_file = temp_file,
  force = opt$force,
  all_funsys = all_funsys,
  task_size = opt$task_size,
  current_organism_info = current_organism_info,
  workers = opt$workers,
  pvalcutoff = opt$pvalcutoff,
  qvalcutoff = opt$qvalcutoff,
  all_custom_gmt = all_custom_gmt,
  kegg_data_file = kegg_data_file,
  simplify = opt$simplify,
  clean_parentals = opt$clean_parentals
)
enrichments_ORA <- ce_list[["enrichments_ORA"]]
enrichments_ORA_merged <- ce_list[["enrichments_ORA_merged"]]

write_clusters_to_enrichment(
  output_path=output_path,
  output_file=opt$output_file,
  mode=opt$mode,
  enrichments_ORA=enrichments_ORA,
  enrichments_ORA_merged=enrichments_ORA_merged,
  task_size = opt$task_size,
  workers = opt$workers,
  template_folder = template_folder,
  top_categories = opt$top_categories,
  group_results = opt$group_results,
  n_category = opt$showCategories,
  sim_thr = opt$sim_thr, 
  summary_common_name = opt$summary_common_name, 
  pvalcutoff = opt$pvalcutoff,
  gene_attributes = ce_list[["gene_attributes"]],
  gene_attribute_name = ce_list[["gene_attribute_name"]]
)