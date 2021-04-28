#!/usr/bin/env Rscript

#' @import dplyr

no_pkg_messages <- suppressPackageStartupMessages
suppressPackageStartupMessages(library(optparse)) 

 
 #### NOTAS PARA EL SCRIPT

full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile), 
               error=function(e) # works when using R CMD
              normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', 
                commandArgs())], '='))[2]))
main_path_script <- dirname(full.fpath)


templates_folder <- file.path(main_path_script, '../templates')
option_list <- list(
    optparse::make_option(c("-t", "--target_miRNA_results"), type="character", 
        default=NULL,
        help="Set miRNA-target results folder"),
    optparse::make_option(c("-s", "--strat"), type="character", 
        default="EE,Eh,Ed,hd,hE",
        help = paste0("zero = All possible correlations between miRNA counts",
            " and RNA counts, EE = Anticorrelation between RNAseq modules ",
            "Eigengene and miRNAseq modules Eigengene, Eh = Anticorrelation ",
            "between RNAseq modules Eigengene and miRNAseq hub genes, Ed = ",
            "Anticorrelation between RNAseq modules Eigengene and miRNAseq DEG",
            " expression profiles, hd = Anticorrelation between RNAseq hub ",
            "genes and miRNAseq DEG expression profiles, hE = Anticorrelation",
            " between RNAseq hub genes and miRNAseq modules Eigengene. ",
            "Default : %default")),
    optparse::make_option(c("--organism"), type ="character", default="Human",
        help= paste0("Reference organism to use. Available 'Human' for human",
            " and 'Mouse' for mouse.")),
    optparse::make_option(c("-r", "--RNAseq_folder"), type="character", default=NULL,
        help="DEgenesHunter RNAseq execution folder"),
    optparse::make_option(c("-m", "--miRNAseq_folder"), type="character", default=NULL,
        help="DEgenesHunter miRNAseq execution folder"),
    optparse::make_option(c("--mode"), type ="character", default="d",
        help= paste0("Set execution mode. 'd' = default mode, enrichments are",
            " performed with all targets. 'e' = expanded mode, targets of each",
            " miRNA are enriched separatedly.")),
    optparse::make_option(c("-f", "--func_annot_db"), type="character", default="g,K,R",
        help=paste0("Comma separated nomenclature to perform analysis ",
            "(clusterProfiler: K = KEGG, g = GO, R = Reactome). ",
            "[Default=%default]")),
    optparse::make_option(c("-G", "--GO_subont"), type="character", 
        default=c("BP,MF,CC"),
        help=paste0("GO sub-ontologies to use for functional analysis ",
            "(MF = Molecular Function, BP = Biological Process, CC = Celular ",
            "Component). Default=%default")), # Not Checked
    optparse::make_option(c("-C", "--custom"), ,type = "character", default="",
        help=paste0("Files with custom functional annotation database ",
            "(in GMT format) separated by commas (,)")),
    optparse::make_option(c("-P", "--pthreshold"), type="double", default=0.1,
        help="Enrichment p-value threshold. [Default = %default]"),
    optparse::make_option(c("-Q", "--qthreshold"), type="double", default=0.2,
        help="Enrichment q-value threshold. [Default = %default]"),
      optparse::make_option(c("-c", "--cores"), type="double", 
        help=paste0("Cores available for parallelization"),
        default=1),
      optparse::make_option(c("-o", "--output_files"), type="character", 
        default="results",
        help="Output path. Default=%default"))

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
opt$help <- NULL

source(file.path(main_path_script, '../../R', 'functional_analysis_library.R'))
source(file.path(main_path_script, '../../R', 'plotting_functions.R'))
source(file.path(main_path_script, '../../R', 'miRNA_RNA_functions.R'))
source(file.path(main_path_script, '../../R', 'general_functions.R'))
source(file.path(main_path_script, '../../R', 'main_target_functional.R'))

suppressPackageStartupMessages(require(clusterProfiler))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(enrichplot))
suppressPackageStartupMessages(require(stringr))  
suppressPackageStartupMessages(require(dplyr))  
main_targets_functional <- function(
target_miRNA_results, 
strat,
organism,
RNAseq_folder,
miRNAseq_folder,
func_annot_db,
GO_subont,
custom,
pthreshold,
qthreshold,
output_files,
mode

    ){

str_names <- list("dd" = "normalized_counts_RNA_vs_miRNA_normalized_counts",  
                    "EE" = "Eigengene_RNA_vs_miRNA_Eigengene" ,
                    "Eh" = "Eigengene_RNA_vs_miRNA_hub_1", 
                    "Ed" = "Eigengene_RNA_vs_miRNA_normalized_counts", 
                    "hd" = "hub_1_RNA_vs_miRNA_normalized_counts", 
                    "hE" = "hub_1_RNA_vs_miRNA_Eigengene") 

    #Configure execution 

    launch_default <- grepl("d", mode)
    launch_expanded <- grepl("e", mode)


    #Load biomart table
    all_organisms_info <- read.table(file.path(main_path_script, "../external_data", 
        "organism_table.txt"), header = TRUE, row.names=1, sep="\t", 
    stringsAsFactors = FALSE, fill = NA)
    current_organism_info <- subset(all_organisms_info, rownames(
        all_organisms_info) %in% organism)  
    nomenclatures <- unlist(strsplit(func_annot_db, ","))

    strategy <- str_names[[strat]]

    output_files <- file.path(output_files, strategy)
    dir.create(output_files)
    #Enrichments
    message(paste0("\tPerforming functional analysis for ", strategy))
    strategy_folder <- file.path(target_miRNA_results, strategy)

    if (!"target_results_table.txt" %in% list.files(strategy_folder)) {
        warning(paste0("No files have been found in results folder for ", strategy, 
            " strategy. Aborting..."))
        next
    }
    raw_data <- read.table(file.path(strategy_folder, "target_results_table.txt"), 
        header=TRUE, row.names=NULL, sep="\t")
        unique_genes <- unique(raw_data[!raw_data$entrezgene %in% c("", NA), 
            c("entrezgene", "mean_logFCs")])
        enrich_genes <- unique_genes$entrezgene
        geneList <- unique_genes$mean_logFCs
        names(geneList) <- unique_genes$entrezgene

    if (launch_expanded) {
        unique_miRNAs <- unique(raw_data$miRNA_names)
        expanded_targets <- lapply(unique_miRNAs, function(miRNA){
                mirna_targets <- raw_data[raw_data$miRNA_names == miRNA, 
                "entrezgene"]
                return(unique(mirna_targets))
            })
        names(expanded_targets) <- unique_miRNAs
        enrichments_ORA_expanded <- list()
    }

    if ("g" %in% nomenclatures){
       
        GO_subontologies <- unlist(strsplit(GO_subont, ","))
        if (launch_default) {
            enrich_GO <- lapply(GO_subontologies, function(mod) {
                message(paste0("Performing ", mod, "_subonthology of GO"))
                enrich <- no_pkg_messages(enrichment_ORA(genes = enrich_genes, 
                                                    organism = current_organism_info$Bioconductor_DB[1], 
                                                    keyType = "ENTREZID",
                                                    pvalueCutoff = pthreshold, 
                                                    pAdjustMethod = "BH", 
                                                    ont = paste0("GO_", mod), 
                                                    qvalueCutoff = qthreshold))
                return(enrich)
            })
            names(enrich_GO) <- GO_subontologies

        }
        if (launch_expanded) {
            for (subont in GO_subontologies) {
                subont_name <- paste0("GO_", subont)
                enrichments_ORA_expanded[[subont_name]] <- enrichment_clusters_ORA(
                                    genes = expanded_targets,
                                    organism = current_organism_info$Bioconductor_DB[1],
                                    keyType = "ENTREZID",
                                    pvalueCutoff = pthreshold,
                                    pAdjustMethod = "BH",
                                    ont = subont_name,
                                    qvalueCutoff = qthreshold,
                                    useInternal = FALSE)
            }
        }


    }

    if ("K" %in% nomenclatures) {
        message("\tPerforming KEGG enrichments")
        if (launch_default) {
            enrich_KEGG <- enrichment_ORA(genes = enrich_genes, 
                                    organism = current_organism_info$KeggCode[1], 
                                    keyType = "kegg",
                                    pvalueCutoff = pthreshold,
                                    pAdjustMethod = "BH",
                                    ont = "KEGG", 
                                    useInternal = TRUE, 
                                    qvalueCutoff = qthreshold)
            
                
        }
        if (launch_expanded) {
            enrichments_ORA_expanded[["KEGG"]] <- enrichment_clusters_ORA(
                                     genes = expanded_targets,
                                     organism = current_organism_info$KeggCode[1],
                                     keyType = "kegg",
                                     pvalueCutoff = pthreshold,
                                     pAdjustMethod = "BH",
                                     ont = "KEGG",
                                     qvalueCutoff = qthreshold,
                                     useInternal = TRUE)
        }

    }

    if ("R" %in% nomenclatures) {
        message("\tPerforming Reactome enrichments")
        if (launch_default) {
            enrich_react <- enrichment_ORA(genes = enrich_genes,
                                    organism = current_organism_info$Reactome_ID[1],
                                    keyType = "ENTREZID",
                                    pvalueCutoff = pthreshold, 
                                    pAdjustMethod = "BH",
                                    ont = "REACT", 
                                    qvalueCutoff = qthreshold)        
                
        }
        if (launch_expanded) {
            enrichments_ORA_expanded[["REACT"]] <- enrichment_clusters_ORA(
                                genes = expanded_targets,
                                organism = current_organism_info$Reactome_ID[1],
                                keyType = "ENTREZID",
                                pvalueCutoff = pthreshold,
                                pAdjustMethod = "BH",
                                ont = "REACT",
                                qvalueCutoff = qthreshold,
                                useInternal = FALSE)
        }

    }

    enrich_custom <- NULL
    if (custom != "") {
        message("\tPerforming custom enrichments")
        if (launch_default) {

            custom_files <- unlist(strsplit(custom,","))

            enrich_custom <- enrich_all_customs(custom_files = custom_files, 
                                                p_val_threshold = pthreshold, 
                                                genes = enrich_genes)
        }
        if (launch_expanded) {
            custom_cls_ORA_expanded <- lapply(custom_files, function(gmt){ 
    # custom_cls_ORA_expanded name is given for make compatible 
    # ora_customEnrichResult.Rmd 
                custom_set <- load_and_parse_gmt(gmt)
                enrich_clusters_with_gmt(custom_set = custom_set, 
                                        genes_in_modules = expanded_targets, 
                                        p_val_threshold = pthreshold)
            })
            custom_targets_ORA <- lapply(custom_cls_ORA_expanded, merge_result)
        }
    }
    save(list = ls(all.names = TRUE), file = "/home/bio267lab/proyectos/target_miRNA_2020/test_func.RData")
    q()
    return(
        target_func_results <- 
        list(enrich_GO = enrich_GO,
        enrich_KEGG = enrich_KEGG,
        enrich_react = enrich_react,
        # enrich_custom = enrich_custom,
        launch_default = launch_default,
        launch_expanded = launch_expanded,
        strategy = strategy,
        enrichments_ORA_expanded = enrichments_ORA_expanded,
        nomenclatures = nomenclatures,
        current_organism_info = current_organism_info,
        geneList = geneList,
        unique_miRNAs = unique_miRNAs,
raw_data=raw_data

    ))

}

write_functional_targets <- function(
    enrich_GO,
    enrich_KEGG,
    enrich_react,
    launch_default,
    launch_expanded,
    output_files,
    strategy,
    enrichments_ORA_expanded,
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

    if (launch_default) {

        outf <- file.path(output_files, paste0(strategy, 
            "_targets_functional.html"))

        rmarkdown::render(file.path(templates_folder, 
            'targets_functional.Rmd'), output_file = outf, 
            intermediates_dir = output_files)
    }
    if (launch_expanded) {
        message("Rendering specific miRNA reports")
        enrichments_ORA <- lapply(enrichments_ORA_expanded, merge_result)
        enrichments_ORA <- lapply(enrichments_ORA, function(res){
                                if(nrow(res@compareClusterResult) > 0)
                                   res <- catched_pairwise_termsim(res, 200)
                             return(res)
                           })
        RNAseq <- load_DEGH_information(RNAseq_folder)
        miRNAseq <- load_DEGH_information(miRNAseq_folder)
        RNAseq[['normalized_counts']] <- as.data.frame(as.table(
            scale_data_matrix(data_matrix = as.matrix(RNAseq[['normalized_counts']]))))
        colnames(RNAseq[['normalized_counts']]) <- c("Sample","Gene","Count")
        miRNAseq[['normalized_counts']] <- scale_data_matrix(
            data_matrix = as.matrix(miRNAseq[['normalized_counts']]))
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
        results_temp <- file.path(paste0(output_files, "_tmp"))
        check_and_create_dir(results_temp)

        invisible(parallel_list(unique_miRNAs, function(miRNA) {
            # Take output name
            target_outf <- file.path(output_files, paste0("targets_", 
                miRNA,".html"))
            # Generate report
            rmarkdown::render(file.path(templates_folder, 
                'miRNA_target_func.Rmd'), output_file = target_outf, 
            intermediates_dir = file.path(results_temp, miRNA), quiet=TRUE)
            
        }, workers = cores, task_size= task_size))
       
        message("\tRendering merged miRNA report")
        outf_cls <- file.path(output_files, "expanded_targets_func.html")
        rmarkdown::render(file.path(templates_folder, 
            'targets_global_report.Rmd'),output_file = outf_cls,
             intermediates_dir = output_files)
    }
}


target_func_results <- do.call("main_targets_functional", opt)
#  save(target_func_results, file = "/home/bio267lab/proyectos/target_miRNA_2020/test_func_2.RData")
# q()
# load("/home/bio267lab/proyectos/target_miRNA_2020/test_func_2.RData")
task_size <- 4
target_func_results$enrich_custom <- NULL
target_func_results <- c(target_func_results,
    list(output_files = opt$output_files,


    RNAseq_folder = opt$RNAseq_folder,

    miRNAseq_folder = opt$miRNAseq_folder,

    templates_folder = templates_folder,
    strat = opt$strat,
    task_size = task_size,
    cores = opt$cores

))

do.call("write_functional_targets", target_func_results)