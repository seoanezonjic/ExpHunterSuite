#!/usr/bin/env Rscript

#' @import dplyr

suppressPackageStartupMessages(library(optparse)) 

 
 #### NOTAS PARA EL SCRIPT

full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile), 
               error=function(e) # works when using R CMD
              normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', 
                commandArgs())], '='))[2]))
main_path_script <- dirname(full.fpath)

str_names <- list("dd" = "counts_RNA_vs_counts_miRNA",  
                    "EE" = "Eigen_RNA_v_Eigen_miRNA" ,
                    "Eh" = "Eigen_RNA_v_Hub1_miRNA", 
                    "Ed" = "Eigen_RNA_v_counts_miRNA", 
                    "hd" = "Hub1_RNA_v_counts_miRNA", 
                    "hE" = "Hub1_RNA_v_Eigen_miRNA"
            ) 


option_list <- list(
    make_option(c("-t", "--target_miRNA_results"), type="character", 
        default=NULL,
        help="Set miRNA-target results folder"),
    make_option(c("-a", "--aproaches"), type="character", 
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
    make_option(c("--organism"), type ="character", default="hsa",
        help= paste0("Reference organism to use. Available 'Human' for human",
            " and 'Mouse' for mouse.")),
    make_option(c("-r", "--RNAseq_folder"), type="character", default=NULL,
        help="DEgenesHunter RNAseq execution folder"),
    make_option(c("-m", "--miRNAseq_folder"), type="character", default=NULL,
        help="DEgenesHunter miRNAseq execution folder"),
    make_option(c("--mode"), type ="character", default="d",
        help= paste0("Set execution mode. 'd' = default mode, enrichments are",
            " performed with all targets. 'e' = expanded mode, targets of each",
            " miRNA are enriched separatedly.")),
    make_option(c("-f", "--func_annot_db"), type="character", default="g,K,R",
        help=paste0("Comma separated nomenclature to perform analysis ",
            "(clusterProfiler: K = KEGG, g = GO, R = Reactome). ",
            "[Default=%default]")),
    make_option(c("-G", "--GO_subont"), type="character", 
        default=c("BP,MF,CC"),
        help=paste0("GO sub-ontologies to use for functional analysis ",
            "(MF = Molecular Function, BP = Biological Process, CC = Celular ",
            "Component). Default=%default"), # Not Checked
    make_option(c("-C", "--custom"), ,type = "character", default=NULL,
        help=paste0("Files with custom functional annotation database ",
            "(in GMT format) separated by commas (,)")),
    make_option(c("-P", "--pthreshold"), type="double", default=0.1,
        help="Enrichment p-value threshold. [Default = %default]"),
    make_option(c("-Q", "--qthreshold"), type="double", default=0.2,
        help="Enrichment q-value threshold. [Default = %default]"),
      make_option(c("-o", "--output_files"), type="character", 
        default="results",
        help="Output path. Default=%default")

    )
opt <- parse_args(OptionParser(option_list=option_list))

print(opt)
source(file.path(main_path_script, 'lib', 'functional_analysis_library.R'))
source(file.path(main_path_script, 'lib', 'plotting_functions.R'))
source(file.path(main_path_script, 'lib', 'miRNA_RNA_functions.R'))
suppressPackageStartupMessages(require(clusterProfiler))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(enrichplot))
suppressPackageStartupMessages(require(stringr))  
suppressPackageStartupMessages(require(dplyr))  




#Configure execution 

launch_default <- grepl("d", opt$mode)
launch_expanded <- grepl("e", opt$mode)


#Load biomart table
all_organisms_info <- read.table(file.path(main_path_script, "lib", 
    "organism_table.txt"), header = TRUE, row.names=1, sep="\t", 
stringsAsFactors = FALSE, fill = NA)
current_organism_info <- subset(all_organisms_info, rownames(
    all_organisms_info) %in% opt$organism)  
nomenclatures <- unlist(strsplit(opt$func_annot_db, ","))

strategy <- str_names[[opt$aproach]]

output_files <- file.path(opt$output_files, strategy)
dir.create(output_files)
#Enrichments
message(paste0("\tPerforming functional analysis for ", strategy))
strategy_folder <- file.path(opt$target_miRNA_results, strategy)

if (!"target_results_table.txt" %in% list.files(strategy_folder)) {
    warning(paste0("No files have been fond in results folder for ", strategy, 
        " strategy. Aborting..."))
    next
}
raw_data <- read.table(file.path(strategy_folder, "target_results_table.txt"), 
    header=TRUE, row.names=NULL, sep="\t")

    unique_genes <- unique(raw_data[!raw_data$entrezgene %in% c("", NA), 
        c("entrezgene", "mean_logFCs")])
    enrich_genes <- unique_genes$entrezgene
    geneList <- unique_genes$mean_logFCs
    names(geneList) <-unique_genes$entrezgene

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
    GO_subontologies <- unlist(strsplit(opt$GO_subont, ","))
    if (launch_default) {
    enrich_GO <- lapply(GO_subontologies, function(mod) {
        message(paste0("Performing ", mod, "_subonthology of GO"))
        enrich <- enrichment_ORA(genes = enrich_genes, 
                            organism = current_organism_info$Bioconductor_DB[1], 
                            keyType = "ENTREZID",
                            pvalueCutoff = opt$pthreshold, 
                            pAdjustMethod = "BH", 
                            ont = paste0("GO_", mod), 
                            qvalueCutoff = opt$qthreshold)
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
                            pvalueCutoff = opt$pthreshold,
                            pAdjustMethod = "BH",
                            ont = subont_name,
                            qvalueCutoff = opt$qthreshold,
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
                                pvalueCutoff = opt$pthreshold,
                                pAdjustMethod = "BH",
                                ont = "KEGG", 
                                useInternal = TRUE, 
                                qvalueCutoff = opt$qthreshold)
        
            
    }
    if (launch_expanded) {
        enrichments_ORA_expanded[["KEGG"]] <- enrichment_clusters_ORA(
                                 genes = expanded_targets,
                                 organism = current_organism_info$KeggCode[1],
                                 keyType = "kegg",
                                 pvalueCutoff = opt$pthreshold,
                                 pAdjustMethod = "BH",
                                 ont = "KEGG",
                                 qvalueCutoff = opt$qthreshold,
                                 useInternal = TRUE)
    }

}

if ("R" %in% nomenclatures) {
    message("\tPerforming Reactome enrichments")
    if (launch_default) {
        enrich_react <- enrichment_ORA(genes = enrich_genes,
                                organism = current_organism_info$Reactome_ID[1],
                                keyType = "ENTREZID",
                                pvalueCutoff = opt$pthreshold, 
                                pAdjustMethod = "BH",
                                ont = "REACT", 
                                qvalueCutoff = opt$qthreshold)        
            
    }
    if (launch_expanded) {
        enrichments_ORA_expanded[["REACT"]] <- enrichment_clusters_ORA(
                            genes = expanded_targets,
                            organism = current_organism_info$Reactome_ID[1],
                            keyType = "ENTREZID",
                            pvalueCutoff = opt$pthreshold,
                            pAdjustMethod = "BH",
                            ont = "REACT",
                            qvalueCutoff = opt$qthreshold,
                            useInternal = FALSE)
    }

}

enrich_custom <- NULL
if (!is.null(opt$custom)) {
    message("\tPerforming custom enrichments")
    if (launch_default) {

        custom_files <- unlist(strsplit(opt$custom,","))

        enrich_custom <- enrich_all_customs(custom_files = custom_files, 
                                            p_val_threshold = opt$pthreshold, 
                                            genes = enrich_genes)
    }
    if (launch_expanded) {
        custom_cls_ORA_expanded <- lapply(custom_files, function(gmt){ 
# custom_cls_ORA_expanded name is given for make compatible 
# ora_customEnrichResult.Rmd 
            custom_set <- load_and_parse_gmt(gmt)
            enrich_clusters_with_gmt(custom_set = custom_set, 
                                    genes_in_modules = expanded_targets, 
                                    p_val_threshold = opt$pthreshold)
        })
        custom_targets_ORA <- lapply(custom_cls_ORA_expanded, merge_result)

    }
}

message("\tRendering regular report")
debug_folder <- file.path(output_files, "debug_files")
dir.create(debug_folder)

save.image(file = file.path(debug_folder, "debug.RData"))

if (launch_default) {

    outf <- file.path(output_files, paste0(strategy, 
        "_targets_functional.html"))

    rmarkdown::render(file.path(main_path_script, 'templates', 
        'targets_functional.Rmd'), output_file = outf, 
        intermediates_dir = output_files)
}

if (launch_expanded) {
    message("Rendering specific miRNA reports")
    enrichments_ORA <- lapply(enrichments_ORA_expanded, merge_result)

    if (!is.null(opt$RNAseq_folder) && !is.null(opt$RNAseq_folder)) {
        RNAseq <- load_DEGH_information(opt$RNAseq_folder)
        miRNAseq <- load_DEGH_information(opt$miRNAseq_folder)


        RNAseq[['normalized_counts']] <- as.data.frame(as.table(
            scale_data_matrix(data_matrix = as.matrix(RNAseq[['normalized_counts']]))))
        colnames(RNAseq[['normalized_counts']]) <- c("Sample","Gene","Count")
        

        miRNAseq[['normalized_counts']] <- scale_data_matrix(
            data_matrix = as.matrix(miRNAseq[['normalized_counts']]))

        miRNA_strat <- unlist(strsplit(opt$aproaches, ""))[2]
        
        if (miRNA_strat == "h") {
            print(miRNA_strat)
            hub_miRNAs <- get_hub_genes_by_MM(miRNAseq[["normalized_counts"]], 
                miRNAseq[["DH_results"]], top = 1)
            
        } else if (miRNA_strat == "E") {
            print(miRNA_strat)

            miRNAseq$Eigengene <- as.data.frame(as.table(miRNAseq$Eigengene), 
                stringsAsFactors = FALSE)
            colnames(miRNAseq$Eigengene) <- c("Sample","Cluster_ID","Count") 
            tgt_eigvalues_gnorm <- miRNAseq$Eigengene
            tgt_eigvalues_gnorm$Count <- (tgt_eigvalues_gnorm$Count + 1) / 2 

        }
    }

    for(miRNA in unique_miRNAs) {
        
        target_outf <- file.path(output_files, paste0("targets_", 
            miRNA,".html"))
        # Generate report
        rmarkdown::render(file.path(main_path_script, 'templates', 
            'miRNA_target_func.Rmd'), output_file = target_outf, 
        intermediates_dir = output_files)
    }

    message("\tRendering merged miRNA report")
    outf_cls <- file.path(output_files, "expanded_targets_func.html")
    rmarkdown::render(file.path(main_path_script, 'templates', 
        'targets_global_report.Rmd'),output_file = outf_cls,
         intermediates_dir = output_files)
}








