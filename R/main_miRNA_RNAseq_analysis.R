#' @importFrom utils read.table 
#' @importFrom rmarkdown render
coRmiT <- function(
RNAseq_folder,
miRNAseq_folder,
output_files,
strat_names,
sample_proportion, 
organism,
multimir_db,
p_val_cutoff,
corr_cutoff,
permutations, 
MM_cutoff,
report, #this option can be set to default value
translation_file,
databases,
translate_ensembl = FALSE,
filter_db_theshold = 0,
database_to_filter = NULL,
mc_cores = 1,
template_folder = file.path(find.package('ExpHunterSuite'), "templates"),
organism_table_path = file.path(find.package('ExpHunterSuite'), "inst", 
    "external_data", "organism_table.txt") 
){
 #create output folder
 "%>%" <- magrittr::"%>%"
  pred_dbs <- c("diana_microt", "elmmo", "microcosm",
   "miranda","mirdb", "pictar", "pita", "targetscan")
  
 #parse strategies and add default strategies
 strat_names <- c(
    "normalized_counts_RNA_vs_miRNA_normalized_counts", 
    "Eigengene_0_RNA_vs_miRNA_normalized_counts", 
    "normalized_counts_RNA_vs_miRNA_Eigengene_0", 
    # # "DEGs_RNA_vs_miRNA_DEMs",
    "DEGs_DEMs_permutated",
     parse_strategies(strat_names))
  
 # Prepare for RNA ID translation
 organism_info <- utils::read.table(organism_table_path, 
    header = TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE, 
    fill = NA) 
 organism_info <- organism_info[organism_info$KeggCode %in% organism,]
 #Prepare multiMiR
 selected_predicted_databases <- unlist(strsplit(databases, ","))
 selected_predicted_databases <- pred_dbs[pred_dbs %in% unlist(
                                                selected_predicted_databases)]
 if (!is.null(database_to_filter)){
    database_to_filter <- unlist(strsplit(database_to_filter, ","))
 }
 multimir_info <- load_and_parse_multimir(multimir_path = multimir_db, 
    selected_predicted_databases = selected_predicted_databases, 
    filter_db_theshold= filter_db_theshold, #filter_db_theshold
     database_to_filter = database_to_filter  #database_to_filter
    )
 multimir <- multimir_info[["multimir_table"]]
 raw_databases_scores <- multimir_info[["raw_databases_scores"]]
 message("multiMiR database has been parsed and summarized")
 
 #Load and prepare DGH data
 RNAseq   <- load_DEGH_information(RNAseq_folder)
 miRNAseq <- load_DEGH_information(miRNAseq_folder)
 message("RNAseq and miRNAseq data has been loaded")
 
 all_pairs <- prepare_all_pairs(RNAseq = RNAseq, miRNAseq = miRNAseq, 
    multimir = multimir, organism = organism)
 multimir <- NULL
 # Perform strategies
 miRNA_cor_results <- perform_all_strategies(strat_names = strat_names, 
    RNAseq = RNAseq, miRNAseq = miRNAseq, corr_cutoff=  corr_cutoff, 
    p_val_cutoff = p_val_cutoff, MM_cutoff = MM_cutoff,
    permutations = permutations, all_pairs = all_pairs, 
    selected_predicted_databases = selected_predicted_databases, 
    sample_proportion = sample_proportion, 
    raw_databases_scores=raw_databases_scores)
 
 miRNA_cor_results$cont_tables <- v_get_stats(miRNA_cor_results$cont_tables)
  
 miRNA_cor_results$filters_summary$type <- factor(
    miRNA_cor_results$filters_summary$type, 
    levels=c("novel_miRNAs","known_miRNAs","multimir", "multimir_random", 
        "predicted", "predicted_random", "validated","validated_random", 
        "pred_and_val", "pred_and_val_random"))
 miRNA_cor_results$score_comp$log.p.value <- -log10( 
    miRNA_cor_results$score_comp$p.value)
 miRNA_cor_results$score_comp$log.boot.p.value <- -log10(
    miRNA_cor_results$score_comp$boot.p.value)
 
 miRNA_cor_results$p_fisher$fisher.log.p.value <- -log10( 
    miRNA_cor_results$p_fisher$fisher.p.value)

 ################# GET IDS TRANSLATIONS
 translated_ids <- translate_all_id(
    miRNA_IDs = miRNA_cor_results$all_pairs$miRNAseq, 
    RNA_IDs = miRNA_cor_results$all_pairs$RNAseq, 
    organism_info = organism_info, translate_ensembl = translate_ensembl)
 write.table(miRNA_cor_results$cont_tables, 
    file.path(output_files, "strategies_stats.txt"), 
    quote=FALSE, col.names=TRUE, sep="\t")
 
 miRNA_cor_results <- c(
    miRNA_cor_results,
    translated_ids,
    list(
     miRNAseq = miRNAseq,
     RNAseq = RNAseq,
     selected_predicted_databases = selected_predicted_databases
 ))
 return(miRNA_cor_results)
}