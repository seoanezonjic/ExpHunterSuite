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
mode,
cores,
task_size
){
no_pkg_messages <- suppressPackageStartupMessages

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
    all_organisms_info <- read.table(file.path(main_path_script, 
        "../external_data", "organism_table.txt"), header = TRUE, row.names=1, 
        sep="\t", stringsAsFactors = FALSE, fill = NA)
    curr_org_info <- subset(all_organisms_info, rownames(
        all_organisms_info) %in% organism)  
    nomenclatures <- unlist(strsplit(func_annot_db, ","))

    strategy <- str_names[[strat]]

    output_files <- file.path(output_files, strategy)
    check_and_create_dir(output_files)
    #Enrichments
    message(paste0("\tPerforming functional analysis for ", strategy))
    strategy_folder <- file.path(target_miRNA_results, strategy)

    if (!"target_results_table.txt" %in% list.files(strategy_folder)) {
        stop(paste0("No files have been found in results folder for ", 
            strategy, " strategy. Aborting..."))
    }
    raw_data <- read.table(file.path(strategy_folder, 
                                    "target_results_table.txt"), 
        header=TRUE, row.names=NULL, sep="\t")
     str(raw_data)
        unique_genes <- 
        unique(raw_data[!raw_data$ENTREZGENE_target %in% c("", NA), 
            c("ENTREZGENE_target", "Target_log2FC")])
        
        enrich_genes <- unique_genes$ENTREZGENE_target
        geneList <- unique_genes$Target_log2FC
        names(geneList) <- unique_genes$ENTREZGENE_target

    if (launch_expanded) {
        all_miRNA <- raw_data$miRNA
        all_miRNA[is.na(all_miRNA)] <- raw_data[is.na(all_miRNA), "miRNA_ID"]
        unique_miRNAs <-unique(all_miRNA)
        expanded_targets <- lapply(unique_miRNAs, function(miRNA){
                # miRNA_check <- all_miRNA == miRNA
                miRNA_check <- all_miRNA == miRNA & !is.na(raw_data$miRNA) 
          
                if(any(miRNA_check)){
                    mirna_targets <- raw_data[raw_data$miRNA == miRNA, 
                "ENTREZGENE_target"]
                } else {
                    mirna_targets <- raw_data[raw_data$miRNA_ID == miRNA, 
                "ENTREZGENE_target"]
                }
                return(unique(mirna_targets))
            })
        names(expanded_targets) <- unique_miRNAs
        enr_ORA_expanded <- list()

    }

    if ("g" %in% nomenclatures){
       
        GO_subontologies <- unlist(strsplit(GO_subont, ","))
        if (launch_default) {
            enrich_GO <- lapply(GO_subontologies, function(mod) {
                message(paste0("Performing ", mod, "_subonthology of GO"))
                enrich <- no_pkg_messages(
                    enrichment_ORA(genes = enrich_genes, 
                                   organism = curr_org_info$Bioconductor_DB[1], 
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
                enr_ORA_expanded[[subont_name]] <- enrichment_clusters_ORA(
                                    genes = expanded_targets,
                                    organism = curr_org_info$Bioconductor_DB[1],
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
                                    organism = curr_org_info$KeggCode[1], 
                                    keyType = "kegg",
                                    pvalueCutoff = pthreshold,
                                    pAdjustMethod = "BH",
                                    ont = "KEGG", 
                                    useInternal = TRUE, 
                                    qvalueCutoff = qthreshold)
            
                
        }
        if (launch_expanded) {
            enr_ORA_expanded[["KEGG"]] <- enrichment_clusters_ORA(
                                     genes = expanded_targets,
                                     organism = curr_org_info$KeggCode[1],
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
                                    organism = curr_org_info$Reactome_ID[1],
                                    keyType = "ENTREZID",
                                    pvalueCutoff = pthreshold, 
                                    pAdjustMethod = "BH",
                                    ont = "REACT", 
                                    qvalueCutoff = qthreshold)        
                
        }
        if (launch_expanded) {
            enr_ORA_expanded[["REACT"]] <- enrichment_clusters_ORA(
                                genes = expanded_targets,
                                organism = curr_org_info$Reactome_ID[1],
                                keyType = "ENTREZID",
                                pvalueCutoff = pthreshold,
                                pAdjustMethod = "BH",
                                ont = "REACT",
                                qvalueCutoff = qthreshold,
                                useInternal = FALSE)
        }

    }
    # save(enr_ORA_expanded, file = "/mnt/scratch/users/bio_267_uma/josecordoba/NGS_projects/LaforaRNAseq/target_miRNA/test_fun.RData")
    # q()
    enrichments_ORA <- lapply(enr_ORA_expanded, clusterProfiler::merge_result)
    enrichments_ORA <- lapply(enrichments_ORA, function(res){
                             if(nrow(res@compareClusterResult) > 0)
                                res <- catched_pairwise_termsim(res, 200)
                             return(res)
                             })
   
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
                custom_set <- load_and_parse_gmt(gmt)
                enrich_clusters_with_gmt(custom_set = custom_set, 
                                        genes_in_modules = expanded_targets, 
                                        p_val_threshold = pthreshold)
            })
            custom_targets_ORA <- lapply(custom_cls_ORA_expanded, merge_result)
        }
    }

    return(
        target_func_results <- 
        list(enrich_GO = enrich_GO,
        enrich_KEGG = enrich_KEGG,
        enrich_react = enrich_react,
        # enrich_custom = enrich_custom,
        enrichments_ORA_expanded = enr_ORA_expanded,
        launch_default = launch_default,
        launch_expanded = launch_expanded,
        strategy = strategy,
        enrichments_ORA = enrichments_ORA,
        nomenclatures = nomenclatures,
        current_organism_info = curr_org_info,
        geneList = geneList,
        unique_miRNAs = unique_miRNAs,
        raw_data=raw_data
    ))
}