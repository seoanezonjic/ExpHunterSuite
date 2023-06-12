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
corr_cutoffs,
permutations, 
MM_cutoff,
report, #this option can be set to default value
translation_file,
databases,
translate_ensembl = FALSE,
filter_db_theshold = 0,
database_to_filter = NULL,
mc_cores = 1,
tag_filter,
corr_type,
corr_coef,
f_p_val,
selected_targets_file,
eval_method = "aggregated",  
template_folder = file.path(find.package('ExpHunterSuite'), "templates"),
organism_table_path = file.path(find.package('ExpHunterSuite'), "inst", 
    "external_data", "organism_table.txt"),
compare_pred_scores = FALSE
){


#corr_cutoffs <- gsub("'","", corr_cutoffs)
corr_cutoffs<- eval(parse(text=paste('c(', gsub("'","", corr_cutoffs), ')')))
#corr_cutoffs <- as.numeric(unlist(strsplit(corr_cutoffs, ",")))



 #create output folder
 "%>%" <- magrittr::"%>%"
 pred_dbs <- c("diana_microt", "elmmo", "microcosm",
   "miranda","mirdb", "pictar", "pita", "targetscan")
  
 #parse strategies and add default strategies
 strat_names <- c(
    "DEGs_RNA_vs_miRNA_DEMs_opp",
    "DEGs_RNA_vs_miRNA_DEMs_sim",
    "normalized_counts_RNA_vs_miRNA_normalized_counts", 
    "Eigengene_0_RNA_vs_miRNA_normalized_counts", 
    "normalized_counts_RNA_vs_miRNA_Eigengene_0", 
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
    filter_db_theshold= filter_db_theshold, 
     database_to_filter = database_to_filter  
    )
 multimir <- multimir_info[["multimir_table"]]
 raw_databases_scores <- multimir_info[["raw_databases_scores"]]
 
 #Load and prepare DGH data
 RNAseq <- list()

 if (RNAseq_folder != "") {
   RNAseq <- load_DEGH_information(RNAseq_folder) 
 } else {
   RNAseq <- NULL
 } 

 selected_targets <- load_selected_targets(selected_targets_file)
 
 miRNAseq <- load_DEGH_information(miRNAseq_folder) 
 
 all_pairs <- prepare_all_pairs(RNAseq = RNAseq, miRNAseq = miRNAseq, 
    multimir = multimir, selected_targets = selected_targets, 
    organism = organism)
 multimir <- NULL


 # Perform strategies
 miRNA_cor_results <- perform_all_strategies(strat_names = strat_names, 
    RNAseq = RNAseq, miRNAseq = miRNAseq, MM_cutoff = MM_cutoff,
    permutations = permutations, all_pairs = all_pairs, 
    selected_predicted_databases = selected_predicted_databases, 
    tag_filter = tag_filter, sample_proportion = sample_proportion, 
    raw_databases_scores = raw_databases_scores, corr_type = corr_type, 
    corr_coef = corr_coef, selected_targets = selected_targets, 
    compare_pred_scores = compare_pred_scores, eval_method = eval_method)


#remove eval_method
miRNA_cont_tables <- data.frame()
cont_tables <- data.frame()
all_pairs <- as.data.frame(miRNA_cor_results$all_pairs)
message("Preparing data for stats computing")
 for (corr_cutoff in corr_cutoffs){
 #mkdir folder
  message(paste0("Computing ", corr_cutoff, " threshold"))
   for (strategy in strat_names){
       
        contingency_tables <- prepare_for_stats(all_pairs, strategy, corr_cutoff, p_val_cutoff, corr_type = corr_type)
        if (is.null(contingency_tables)) next
          miRNA_cont_tables <- rbind(miRNA_cont_tables, contingency_tables$strat_miRNA_ct)
          cont_tables <- rbind(cont_tables, contingency_tables$strat_ct)
   }
   message("test")
 
 } 
 contingency_tables <- list(miRNA_cont_tables = miRNA_cont_tables,
                              cont_tables = cont_tables) 
 miRNA_cor_results <- c(miRNA_cor_results, 
   contingency_tables)


 miRNA_cor_results$cont_tables <- v_get_stats(miRNA_cor_results$cont_tables)
 miRNA_cor_results$miRNA_cont_tables <- 
                              v_get_stats(miRNA_cor_results$miRNA_cont_tables, 
                              selected_stats = c("v.fisher.test","odds_ratio"))
 miRNA_cont_tables <- miRNA_cor_results$miRNA_cont_tables
 integrated_stats <- filter_and_integrate_OR(miRNA_cont_tables,
                                             p.adjust.method = "BH")
   
 integrated_OR <- integrated_stats[["integrated_stats"]]
 miRNA_cont_tables <- integrated_stats[["miRNA_ct"]]
 

 miRNA_cont_tables <- miRNA_cont_tables[!grepl("Eigengene_0", miRNA_cont_tables$strategy) & 
                                          !grepl("sim", miRNA_cont_tables$strategy) &
                                          !grepl("opp", miRNA_cont_tables$strategy) &
                                          miRNA_cont_tables$Odds_ratio >0.0001 &
                                          miRNA_cont_tables$Odds_ratio < 10000 &
                                          miRNA_cont_tables$Pvalue <= f_p_val,]


 if(nrow(miRNA_cont_tables) == 0 ){
    stop(paste0("ERROR: Any miRNA has significant overlapping with databases at any ",
                "given strategy/correlation threshold. Please try to modify parametres,"))
 }
 miRNA_cor_results$miRNA_cont_tables <-  miRNA_cont_tables                                
 miRNA_cor_results$cont_tables <- merge(miRNA_cor_results$cont_tables, 
                                          integrated_OR, 
                                          by=c("strategy","db_group", "corr_cutoff"), 
                                          all.x = TRUE)


  integrated_strat <- get_optimal_pairs(miRNA_cor_results$all_pairs, 
                                         miRNA_cor_results$miRNA_cont_tables,
                                             p_val_cutoff,
                                             corr_type = corr_type
                                             )
 all_pairs$integrated_strat <- integrated_strat
 

 integrated_ct <- prepare_for_stats(all_pairs, "integrated_strat", 0, p_val_cutoff, 
                                                sig_pairs = TRUE)

integrated_ct$cont_tables <- v_get_stats(integrated_ct$strat_ct)

 integrated_ct$miRNA_cont_tables <- 
                              v_get_stats(integrated_ct$strat_miRNA_ct, 
                              selected_stats = c("v.fisher.test","odds_ratio"))
 #integrated_ct <- integrated_ct$miRNA_cont_tables
 integrated_stats <- filter_and_integrate_OR(integrated_ct$miRNA_cont_tables,
                                            p.adjust.method = "BH")
   
 integrated_OR <- integrated_stats[["integrated_stats"]]
 int_miRNA_cont_tables <- integrated_stats[["miRNA_ct"]]

 miRNA_cor_results$int_miRNA_cont_tables <- rbind(miRNA_cor_results$miRNA_cont_tables,
                                                   integrated_stats[["miRNA_ct"]])


 integrated_ct$cont_tables <- merge(integrated_ct$cont_tables, 
                                          integrated_OR, 
                                          by=c("strategy","db_group", "corr_cutoff"), 
                                          all.x = TRUE)
 miRNA_cor_results$int_cont_tables <- rbind(miRNA_cor_results$cont_tables,
                                                   integrated_ct$cont_tables)


 ## add crossvalidation

 miRNA_cont_tables$strategy <- set_strats_readable(miRNA_cont_tables$strategy)
 all_cor_dist <- get_cor_dist(all_pairs, strat_names)
 ################# GET IDS TRANSLATIONS
 RNA_IDs <- selected_targets
 if (!is.null(RNAseq_folder)){
   RNA_IDs <- all_pairs$RNAseq
 }

 translated_ids <- translate_all_id(
    miRNA_IDs = all_pairs$miRNAseq, 
    RNA_IDs = RNA_IDs, 
    organism_info = organism_info, translate_ensembl = translate_ensembl)
    
 miRNA_cor_results$integrated_pairs <-  all_pairs[
  all_pairs$integrated_strat &
  !is.na(all_pairs$integrated_strat),]  


 miRNA_cor_results$all_pairs <- all_pairs

 miRNA_cor_results <- c(
    miRNA_cor_results,
    translated_ids,
    list(
     miRNAseq = miRNAseq,
     RNAseq = RNAseq,
     selected_predicted_databases = selected_predicted_databases,
     all_cor_dist = all_cor_dist
 ))
 
 return(miRNA_cor_results)
}