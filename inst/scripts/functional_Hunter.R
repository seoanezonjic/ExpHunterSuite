#! /usr/bin/env Rscript
#############################################
############## FUNCTIONAL HUNTER ###########
#############################################

if( Sys.getenv('DEGHUNTER_MODE') == 'DEVELOPMENT' ){
	full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile),  # works when using source
               error=function(e) # works when using R CMD
              normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', commandArgs())], '='))[2]))
	main_path_script <- dirname(full.fpath)
	root_path <- file.path(main_path_script, '..', '..')

	custom_libraries <- c('general_functions.R', 'functional_analysis_library.R', 'plotting_functions.R', 'main_functional_hunter.R', "io_handling.R", "plotting_functions.R", "write_report.R")
	for (lib in custom_libraries){
		source(file.path(root_path, 'R', lib))
	}
	organisms_table_file <- file.path(root_path, "inst", "external_data", "organism_table.txt")
	template_folder <- file.path(root_path, 'inst', 'templates')

	#############################################
	### LOAD LIBRARIES
	#############################################
	# suppressPackageStartupMessages(require(optparse))
	# suppressPackageStartupMessages(require(biomaRt)) 
	# suppressPackageStartupMessages(require(topGO))
	# suppressPackageStartupMessages(require(KEGGREST))
	# suppressPackageStartupMessages(require(stringr))
	# suppressPackageStartupMessages(require(plyr))
	# suppressPackageStartupMessages(require(knitr))
	# suppressPackageStartupMessages(require(clusterProfiler))
	# suppressPackageStartupMessages(require(dplyr))
}else{
	require('DEgenesHunter')
	root_path <- find.package('DEgenesHunter')
	template_folder <- file.path(root_path, 'templates')
	organisms_table_file <- file.path(root_path, "external_data", "organism_table.txt")
}


 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
##                                                                                                                   ##
##                                                     INITIALIZE                                                    ##                                                     
##                                                                                                                   ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

option_list <- list(
  optparse::make_option(c("-i", "--input_hunter_folder"), type="character",
    help="DEgenes Hunter differential expression analysis output folder"), 
  optparse::make_option(c("-m", "--model_organism"), type="character",
    help="Species to use for functional enrichment. You can see all available species running functional_Hunter.R -L"),
  optparse::make_option(c("-L", "--List_organisms"), action="store_true", type="logical", default=FALSE, 
    help="Print all organisms available and ends the program"),
  optparse::make_option(c("-a", "--annot_file"), type="character", default = NULL,
  	help="If the species used does not exist in organism_table.txt, please add a two-column file mapping orthologue ENSEMBL gene ids from a model organism (first column) to the original IDs (second column)."),
  optparse::make_option(c("-t", "--input_gene_id"), type="character", default="E",
    help="Input gene IDs. Available IDs are: ENSEMBL (E), entrezgene (e), TAIR/Arabidopsis (T), Gene Names (G). [Default:%default]"),      	
  optparse::make_option(c("-f", "--func_annot_db"), type="character", default="gKR",
    help="Functional annotation database and enrichment method(s) to use (topGO: G = GO | clusterProfiler: K = KEGG, g = GO, R = Reactome). [Default=%default]"),
  optparse::make_option(c("-G", "--GO_subont"), type="character", default=c("BMC"),
    help="GO sub-ontologies to use for functional analysis (M = Molecular Function, B = Biological Process, C = Celular Component). Default=%default"), # Not Checked
  optparse::make_option(c("-C", "--custom"), ,type = "character", default=NULL,
    help="Files with custom functional annotation database (in GMT format) separated by commas (,)"),
  optparse::make_option(c("-A", "--analysis_type"), type="character", default=c("go"),
    help="Analysis performance (g = Gene Set Enrichment Analysis, o = Over Representation Analysis). Default=%default"), # Not Checked
   optparse::make_option(c("-r", "--remote"), ,type = "character", default="",
    help="Flags to activate remote query from enrichments and Genes translation. Use (b) to launch biomaRt translation; (k) to use Kegg remote data base"),
  optparse::make_option(c("-q", "--save_query"), type="logical", action = "store_true", default=FALSE,
    help="Flag to save biomaRt query."), 
  optparse::make_option(c("-P", "--pthreshold"), type="double", default=0.1,
    help="Enrichment p-value threshold. [Default = %default]"),
  optparse::make_option(c("-Q", "--qthreshold"), type="double", default=0.2,
    help="Enrichment q-value threshold. [Default = %default]"),
  optparse::make_option(c("-c", "--cores"), ,type = "numeric", default=1,
    help="Cores to be used to parallelize clusters enrichments. Default : %default"),
  optparse::make_option(c("-o", "--output_files"), type="character", default="results",
    help="Output path. Default=%default")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))



# # Special IDs
fc_colname <- "mean_logFCs"

############ CREATE OUTPUT FOLDERS #########
paths <- list()
dir.create(opt$output_files)
paths$root <- opt$output_files


organisms_table <- get_organism_table(organisms_table_file)


# Load DEGenesHunter config
DEGenesHunter_expression_opt <- read.table(file.path(opt$input_hunter_folder, "opt_input_values.txt"), header = FALSE, stringsAsFactors = FALSE, sep = "\t")
degh_exp_threshold <- as.numeric(DEGenesHunter_expression_opt[which(DEGenesHunter_expression_opt[,1] == "p_val_cutoff"), 2])

# Load Hunter Folder
hunter_results <- load_hunter_folder(opt$input_hunter_folder)
if(is.null(opt$annot_file)){
	annot_table <- NULL
}else{
	annot_table <- read.table(opt$annot_file, header=FALSE, row.names=NULL, sep="\t", stringsAsFactors = FALSE, quote = "")
}

# Load customs
if (!is.null(opt$custom)) {
	custom_files <- unlist(strsplit(opt$custom, ","))
	all_custom_gmt <- lapply(custom_files, load_and_parse_gmt)
	names(all_custom_gmt) <- custom_files
}else{
	all_custom_gmt <- NULL
}


if(opt$List_organisms){
	# print(as.character(rownames(all_organisms_info)))
	message('Check this list and choose one model species.')
}else{
	func_results <- functional_hunter(
		hunter_results = hunter_results,
		model_organism = opt$model_organism,
		annot_table = annot_table,
		input_gene_id = opt$input_gene_id,
		func_annot_db = opt$func_annot_db,
		GO_subont = opt$GO_subont,
		custom = all_custom_gmt,
		analysis_type = opt$analysis_type,
		remote = opt$remote,
		save_query = opt$save_query,
		pthreshold = opt$pthreshold,
		qthreshold = opt$qthreshold,
		cores = opt$cores,
		output_files = opt$output_files,
		organisms_table = organisms_table,
		fc_colname = fc_colname)

	# Write outputs
	if("GO_ORA" %in% names(func_results)) write.table(as.data.frame(do.call(rbind,lapply(func_results$GO_ORA,function(res) {as.data.frame(res)}))), file=file.path(opt$output_files, "GO_CL_ora"), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")
	if("GO_GSEA" %in% names(func_results)) write.table(as.data.frame(do.call(rbind,lapply(func_results$GO_GSEA,function(res) {as.data.frame(res)}))), file=file.path(opt$output_files, "GO_CL_gsea"), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")	
	if("KEGG_ORA" %in% names(func_results)) write.table(func_results$KEGG_ORA, file=file.path(opt$output_files, "KEGG_results"), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")
	if("KEGG_GSEA" %in% names(func_results)) write.table(func_results$KEGG_GSEA, file=file.path(opt$output_files, "KEGG_GSEA_results"), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")
	if("REACT_ORA" %in% names(func_results)) write.table(func_results$REACT_ORA, file=file.path(opt$output_files, "REACT_results"), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")
	if("REACT_GSEA" %in% names(func_results)) write.table(func_results$REACT_GSEA, file=file.path(opt$output_files, "REACT_GSEA_results"), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")
	if("DEGH_results_annot" %in% names(func_results)) write.table(func_results$DEGH_results_annot, file=file.path(opt$output_files, "hunter_results_table_annotated.txt"), quote=FALSE, col.names=NA, sep="\t")
	if("CUSTOM" %in% names(func_results)){
		invisible(lapply(func_results$CUSTOM,function(res){
			write.table(res$Result, file=file.path(opt$output_files, paste0(basename(res$File),"_ora_results")), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")
		}))
	}
	if("WGCNA_ORA" %in% names(func_results)){
		for(enrichment_i in 1:length(func_results$WGCNA_ORA)) {
			df <- enrichplot:::fortify.compareClusterResult(func_results$WGCNA_ORA[[enrichment_i]])
			write.table(df, file=file.path(opt$output_files, paste0(names(func_results$WGCNA_ORA[enrichment_i]),"_cls_ora")), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")
		}
	}
	if("WGCNA_GSEA" %in% names(func_results)){
		for (enrichment_i in 1:length(func_results$WGCNA_GSEA)) {
			df <- func_results$WGCNA_GSEA[[enrichment_i]]@compareClusterResult
			write.table(df, file=file.path(opt$output_files, paste0(names(func_results$WGCNA_GSEA[enrichment_i]),"_cls_gsea")), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")
		}
	}
	if("WGCNA_CUSTOM" %in% names(func_results)){
		for(enrichment_i in 1:length(func_results$WGCNA_CUSTOM)) {
			df <- enrichplot:::fortify.compareClusterResult(func_results$WGCNA_CUSTOM[[enrichment_i]])
			write.table(df, file=file.path(opt$output_files, paste0(basename(names(func_results$WGCNA_CUSTOM[enrichment_i])),"_cls_ORA")), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")
		}
	}

	write_functional_report(hunter_results, func_results, opt$output_files, template_folder, opt$cores, opt$model_organism, organisms_table, fc_colname)
}
