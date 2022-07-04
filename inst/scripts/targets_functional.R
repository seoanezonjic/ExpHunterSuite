#!/usr/bin/env Rscript

#' @import dplyr


 
 #### NOTAS PARA EL SCRIPT

full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile), 
               error=function(e) # works when using R CMD
              normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', 
                commandArgs())], '='))[2]))
main_path_script <- dirname(full.fpath)

root_path <- file.path(main_path_script, '..', '..')
custom_libraries <- c('functional_analysis_library.R', 'plotting_functions.R',
                      'miRNA_RNA_functions.R','general_functions.R',
                      'main_target_functional.R', 'write_report.R', 'io_handling.R')

for (lib in custom_libraries){
    source(file.path(root_path, 'R', lib))
}
templates_folder <- file.path(main_path_script, '..', 'templates')

option_list <- list(
    optparse::make_option(c("-t", "--target_miRNA_results"), type="character", 
        default=NULL,
        help="Set miRNA-target results folder"),
    optparse::make_option(c("-s", "--strat"), type="character", 
        default="EE,Eh,Ed,hd,hE", help = "dd = All possible correlations between miRNA counts and RNA counts, EE = Anticorrelation between RNAseq modules Eigengene and miRNAseq modules Eigengene, Eh = Anticorrelation between RNAseq modules Eigengene and miRNAseq hub genes, Ed = Anticorrelation between RNAseq modules Eigengene and miRNAseq DEG expression profiles, hd = Anticorrelation between RNAseq hub genes and miRNAseq DEG expression profiles, hE = Anticorrelation between RNAseq hub genes and miRNAseq modules Eigengene. Default : %default"),
    optparse::make_option(c("--organism"), type ="character", default="Human",
        help= "Reference organism to use. Available 'Human' for human and 'Mouse' for mouse."),
    optparse::make_option(c("-r", "--RNAseq_folder"), type="character", default=NULL,
        help="DEgenesHunter RNAseq execution folder"),
    optparse::make_option(c("-m", "--miRNAseq_folder"), type="character", default=NULL,
        help="DEgenesHunter miRNAseq execution folder"),
    optparse::make_option(c("--mode"), type ="character", default="d",
        help= "Set execution mode. 'd' = default mode, enrichments are performed with all targets. 'e' = expanded mode, targets of each miRNA are enriched separatedly."),
    optparse::make_option(c("-f", "--func_annot_db"), type="character", default="g,K,R",
        help="Comma separated nomenclature to perform analysis (clusterProfiler: K = KEGG, g = GO, R = Reactome). [Default=%default]"),
    optparse::make_option(c("-G", "--GO_subont"), type="character", 
        default=c("BP,MF,CC"), help="GO sub-ontologies to use for functional analysis (MF = Molecular Function, BP = Biological Process, CC = Celular Component). Default=%default"), 
    optparse::make_option(c("-C", "--custom"), ,type = "character", default="",
        help="Files with custom functional annotation database (in GMT format) separated by commas (,)"),
    optparse::make_option(c("-P", "--pthreshold"), type="double", default=0.1,
        help="Enrichment p-value threshold. [Default = %default]"),
    optparse::make_option(c("-Q", "--qthreshold"), type="double", default=0.2,
        help="Enrichment q-value threshold. [Default = %default]"),
    optparse::make_option(c("-c", "--cores"), type="double", 
        help="Cores available for parallelization", default=1),
    optparse::make_option(c("-T", "--task_size"), ,type = "numeric", default=10,
        help="Number of items to be processed in each parallel task. Default : %default"),
    optparse::make_option(c("-o", "--output_files"), type="character", 
        default="results", help="Output path. Default=%default"))

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
opt$help <- NULL
opt$output_files <- file.path(opt$output_files, "targets_functional")


target_func_results <- do.call("main_targets_functional", opt)



for (funsys in c("enrich_GO", "enrich_react", "enrich_KEGG")){
    if (is.null(target_func_results[[funsys]]))
        next
    if (funsys == "enrich_GO"){
        for (subont in c("MF", "BP", "CC")) {
          utils::write.table(target_func_results[[funsys]][[subont]], 
                file=file.path(opt$output_files, paste0("targets_", unlist(strsplit(funsys, "_"))[2], "_", subont, "_enrichment")), 
                quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")
        }
    } else {
      utils::write.table(target_func_results[[funsys]], 
        file=file.path(opt$output_files, paste0("targets_", unlist(strsplit(funsys, "_"))[2], "_enrichment")), 
        quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")
    
    }
} 
   


for (funsys in names(target_func_results$enrichments_ORA)){
    utils::write.table(target_func_results$enrichments_ORA[[funsys]], 
                file=file.path(opt$output_files, paste0("miRNAs_targets_", funsys, "_enrichments")), 
                quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")
}


target_func_results$enrich_custom <- NULL
target_func_results <- c(target_func_results,
    list(output_path = opt$output_files,
         RNAseq_folder = opt$RNAseq_folder,
         miRNAseq_folder = opt$miRNAseq_folder,
         templates_folder = templates_folder,
         strat = opt$strat,
         task_size = opt$task_size,
         cores = opt$cores
))

do.call("write_functional_targets", target_func_results)