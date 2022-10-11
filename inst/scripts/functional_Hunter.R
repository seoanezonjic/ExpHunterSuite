#! /usr/bin/env Rscript
#############################################
############## FUNCTIONAL HUNTER ###########
#############################################
if( Sys.getenv('DEGHUNTER_MODE') == 'DEVELOPMENT' ){
    full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile),  
               error=function(e) # works when using R CMD
              normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', 
                commandArgs())], '='))[2]))
    main_path_script <- dirname(full.fpath)
    root_path <- file.path(main_path_script, '..', '..')

    custom_libraries <- c('general_functions.R', 
        'functional_analysis_library.R', 'plotting_functions.R', 
        'main_functional_hunter.R', "io_handling.R", "plotting_functions.R", 
        "write_report.R")
    for (lib in custom_libraries){
        source(file.path(root_path, 'R', lib))
    }
    organisms_table_file <- file.path(root_path, "inst", "external_data", 
        "organism_table.txt")
    template_folder <- file.path(root_path, 'inst', 'templates')

}else{
    require('ExpHunterSuite')
    root_path <- find.package('ExpHunterSuite')
    template_folder <- file.path(root_path, 'templates')
    organisms_table_file <- file.path(root_path, "external_data", 
        "organism_table.txt")
}
 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
##                                                                      ##
##                            INITIALIZE                                ##                                                     
##                                                                      ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

option_list <- list(
  optparse::make_option(c("-i", "--input_hunter_folder"), type="character",
    help="DEgenes Hunter differential expression analysis output folder"), 
  optparse::make_option(c("-m", "--model_organism"), type="character",
    help=paste0("Species to use for functional enrichment. You can see all",
        " available species running functional_Hunter.R -L")),
  optparse::make_option(c("-L", "--List_organisms"), action="store_true", 
    type="logical", default=FALSE, 
    help="Print all organisms available and ends the program"),
  optparse::make_option(c("-a", "--annot_file"), type="character", 
    default = NULL,
    help=paste0("If the species used does not exist in organism_table.txt,",
        " please add a two-column file mapping orthologue ENSEMBL gene ids",
        " from a model organism (first column) to the original IDs",
        " (second column).")),
  optparse::make_option(c("-t", "--input_gene_id"), type="character", 
    default="E",
    help=paste0("Input gene IDs. Available IDs are: ENSEMBL (E), entrezgene",
        " (e), TAIR/Arabidopsis (T), Gene Names (G). [Default:%default]")),          
  optparse::make_option(c("-f", "--func_annot_db"), type="character", 
    default="gKR",
    help=paste0("Functional annotation database and enrichment method(s) to",
        " use (topGO: G = GO | clusterProfiler: K = KEGG, g = GO, R = ",
        "Reactome). D = DOSE, d = DGN. [Default=%default]")),
  optparse::make_option(c("-k", "--kegg_data_file"), ,type = "character", default=NULL,
    help=paste0("KEGG database file. Can download with download_latest_kegg_db().",
        "If not required but not provided, it will be downloaded to working directory")), 
  optparse::make_option(c("-G", "--GO_subont"), type="character",
    default=c("BMC"),
    help=paste0("GO sub-ontologies to use for functional analysis ",
        "(M = Molecular Function, B = Biological Process, C = Celular",
        " Component). Default=%default")), # Not Checked
  optparse::make_option(c("-C", "--custom"), ,type = "character", default=NULL,
    help=paste0("Files with custom functional annotation database ",
        "(in GMT format) separated by commas (,)")),
  optparse::make_option(c("--gmt_id"), type="character", default="ENTREZID",
    help="What identifier is being used for the genes in the custom gmt file. Default=%default"),
  optparse::make_option(c("-A", "--analysis_type"), type="character", 
    default="o",
    help=paste0("Analysis performance (g = Gene Set Enrichment Analysis,",
    " o = Over Representation Analysis). Default=%default")), # Not Checked
   optparse::make_option(c("-r", "--remote"), ,type = "character", default="",
    help=paste0("Flags to activate remote query from enrichments and Genes",
        " translation. Use (b) to launch biomaRt translation; (k) to use Kegg",
        " remote data base")),
  optparse::make_option(c("-q", "--save_query"), type="logical", 
    action = "store_true", default=FALSE,
    help="Flag to save biomaRt query."), 
  optparse::make_option(c("-P", "--pthreshold"), type="double", default=0.1,
    help="Enrichment p-value threshold. [Default = %default]"),
  optparse::make_option(c("-Q", "--qthreshold"), type="double", default=0.2,
    help="Enrichment q-value threshold. [Default = %default]"),
  optparse::make_option(c("-c", "--cores"), ,type = "numeric", default=1,
    help=paste0("Cores to be used to parallelize clusters enrichments.",
        " Default : %default")),
  optparse::make_option(c("-s", "--task_size"), ,type = "numeric", default=10,
    help=paste0("Number of items to be processed in each parallel task.",
        " Default : %default")),  
  optparse::make_option(c("-o", "--output_files"), type="character", 
    default="results",
    help="Output path. Default=%default"),
  optparse::make_option(c("-R", "--report_modes"), type="character", 
    default="fci",
    help="HTML report modes. 'f' for functional_report, 'c' for cluster_main_report and 'i' for individual module report. Default=%default"),
  optparse::make_option(c("-u", "--universe"), type="character", default=NULL, 
    help="Background genes for enrichment. Default all. Alternative = expressed"),
  optparse::make_option("--clean_parentals", type="logical", default=FALSE, 
    action = "store_true", help="Clean parentals GO terms that appears on the same clusters as child."),
  optparse::make_option("--simplify", type="logical", default=FALSE, 
    action = "store_true", help="Apply simplify function from cluster profiler to enrichment."),
  optparse::make_option("--showCategories", type="integer", default=30, 
    help="Number of top categories to show on clusterProfiler dotplot and emaplot"),
  optparse::make_option(c("-T", "--top_categories"), type="integer", default=50,
    help="Number of top categories for each cluster. Default=%default"),
  optparse::make_option(c("-S", "--sim_thr"), type="double", default=NULL,
    help="Similarity cutoff for grouping categories in Summary mode. Default=%default"),
  optparse::make_option("--group_results", type="logical", default=FALSE, 
    action = "store_true", help="Functions are grouped in most frequent words in emaplots."),
   optparse::make_option("--summary_common_name", type="character", default="ancestor", 
                        help="Name of the term groups. 'significant' to use the most significant term of each group. 'ancestor' to use the common ancestor of the group")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

# NB remote and save_query don't do anything - need to fix remote so it can toggle Biomart

# Special IDs
fc_colname <- "mean_logFCs"

organisms_table <- get_organism_table(organisms_table_file)
if(opt$List_organisms || ! opt$model_organism %in% rownames(organisms_table)){
    print(as.character(rownames(organisms_table)))
    stop('Check this list and choose one model species.')
}

hunter_results <- load_hunter_folder(opt$input_hunter_folder)
if(is.null(opt$annot_file)) {
    annot_table <- NULL
} else {
    annot_table <- read.table(opt$annot_file, header=FALSE, row.names=NULL, 
        sep="\t", stringsAsFactors = FALSE, quote = "")
}

current_organism_info <- organisms_table[rownames(organisms_table) %in% opt$model_organism,]
org_db <- get_org_db(current_organism_info)
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

if(opt$input_gene_id == "e") input_gene_id <- "ENTREZID"
if(opt$input_gene_id == "E") input_gene_id <- "ENSEMBL"
if(opt$input_gene_id == "T") input_gene_id <- "TAIR"
if(opt$input_gene_id == "G") input_gene_id <- "GENENAME"

# Simplest option just to grow the vectors, given complexity of input arguments & interplay
enrich_dbs <- vector()
if(grepl("R", opt$func_annot_db)) enrich_dbs = c(enrich_dbs, "Reactome")
if(grepl("K", opt$func_annot_db)) enrich_dbs = c(enrich_dbs, "KEGG")
if(grepl("G", opt$func_annot_db) || grepl("g", opt$func_annot_db)) {
    if(grepl("B", opt$GO_subont)) enrich_dbs = c(enrich_dbs, "BP")
    if(grepl("C", opt$GO_subont)) enrich_dbs = c(enrich_dbs, "CC")
    if(grepl("M", opt$GO_subont)) enrich_dbs = c(enrich_dbs, "MF")
}
if(grepl("D", opt$func_annot_db)) enrich_dbs = c(enrich_dbs, "DOSE")
if(grepl("d", opt$func_annot_db)) enrich_dbs = c(enrich_dbs, "DGN")
enrich_methods <- vector()
if(grepl("G", opt$func_annot_db)) enrich_methods <- c(enrich_methods, "topGO")
if(grepl("o", opt$analysis_type)) enrich_methods <- c(enrich_methods, "ORA")
if(grepl("g", opt$analysis_type)) enrich_methods <- c(enrich_methods, "GSEA")

kegg_data_file <- opt$kegg_data_file
if("KEGG" %in% enrich_dbs) {
    kegg_data_file <- get_kegg_db_path(opt$kegg_data_file, current_organism_info=current_organism_info)
    if(! file.exists(kegg_data_file)) stop(paste("KEGG file:", kegg_data_file, "not found"))
}

clusters_flag <- grepl("c", opt$report_modes) || grepl("i", opt$report_modes)


func_results <- main_functional_hunter(
       hunter_results = hunter_results,
       model_organism = opt$model_organism,
       annot_table = annot_table,
       input_gene_id = input_gene_id,
       custom = all_custom_gmt,
       enrich_dbs = enrich_dbs,
       kegg_data_file = kegg_data_file,
       enrich_methods = enrich_methods,
       annotation_source = "orgdb", # Other option Biomart, to be added
       pthreshold = opt$pthreshold,
       qthreshold = opt$qthreshold,
       cores = opt$cores,
       task_size = opt$task_size,
       output_files = opt$output_files,
       organisms_table = organisms_table,
       fc_colname = fc_colname,
       universe = opt$universe,
       clean_parentals = opt$clean_parentals,
       simplify = opt$simplify,
       top_categories = opt$top_categories,
       sim_thr = opt$sim_thr,
       summary_common_name = opt$summary_common_name,
       clusters_flag = clusters_flag
)

write_enrich_files(func_results, opt$output_files)
write_functional_report(hunter_results = hunter_results, 
                            func_results = func_results, 
                            output_files = opt$output_files,
                            organisms_table = organisms_table,
                            template_folder = template_folder,
                            cores =  opt$cores,
                            task_size = opt$task_size,
                            report = opt$report_modes,
                            showCategories = opt$showCategories,
                            group_results = opt$group_results
)


