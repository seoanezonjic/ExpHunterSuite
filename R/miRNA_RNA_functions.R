################### FUNCTIONS

#' @importFrom utils read.table
#' @importFrom stringr str_remove
load_DEGH_information <- function(execution_path){ 
  DEGH_exec <- list()
  DEGH_exec[['DH_results']] <- utils::read.table(
                                     file.path(execution_path, 
                                     "Common_results/hunter_results_table.txt"),
                                     header=TRUE, row.names=1, sep="\t")
  DEGH_exec[['DH_results']]$gene_name <- rownames(DEGH_exec[['DH_results']])
  rownames(DEGH_exec[['DH_results']]) <- NULL
  DEGH_exec$DH_results <- DEGH_exec$DH_results[
                                 DEGH_exec$DH_results$gene_name != "id",]
  DEGH_exec[['Eigengene']] <- as.matrix(utils::read.table(
                                  file.path(execution_path, 
                                  "Results_WGCNA/eigen_values_per_samples.txt"),
                                  header=TRUE, row.names=1, sep="\t"))
  colnames(DEGH_exec[["Eigengene"]]) <- stringr::str_remove(
                                             colnames(DEGH_exec[["Eigengene"]]),
                                             "ME")
  DEGH_exec[['normalized_counts']] <- as.matrix(utils::read.table(
                                 file.path(execution_path, 
                                 "Results_DESeq2/Normalized_counts_DESeq2.txt"),
                                 header=TRUE, row.names=1, sep="\t"))
  DEGH_exec[['normalized_counts']] <- t(DEGH_exec[['normalized_counts']]) 
  DEGH_exec[['hub_1']] <- get_hub_genes_by_MM(DEGH_exec[['normalized_counts']], 
                                              DEGH_exec[['DH_results']])
  return(DEGH_exec)
}


filter_DEGH_data <- function(DGH_data, MM_cutoff, tag_filter){

    DGH_results <- DGH_data$DH_results 
    if (tag_filter == "prevalent") {
        tags_to_filter <- c("PREVALENT_DEG")
    } else if (tag_filter %in% c("putative", "possible")){
        tags_to_filter <- c("PREVALENT_DEG", "POSSIBLE_DEG")
    }
    
    all_degs <- DGH_results$genes_tag %in% tags_to_filter

    DGH_data$Eigengene_0 <- DGH_data$Eigengene[,"0", drop = FALSE]
    if (tag_filter == "putative") {
        modules_with_DEGS <- unique(DGH_results[all_degs & 
            DGH_results$Cluster_MM >= MM_cutoff, "Cluster_ID"])
        #Removing module 0
        modules_with_DEGS <- modules_with_DEGS[modules_with_DEGS != 0]
        DGH_data$Eigengene <- DGH_data$Eigengene[,
               colnames(DGH_data$Eigengene) %in% modules_with_DEGS]
        DGH_data$hub_1 <- DGH_data$hub_1[,
               colnames(DGH_data$hub_1) %in% modules_with_DEGS]
        
        candidate_not_deg <- DGH_results$Cluster_ID %in% modules_with_DEGS &
                                           DGH_results$genes_tag == "NOT_DEG" &
                                           DGH_results$Cluster_MM >= MM_cutoff
        candidate_not_deg <- DGH_results[candidate_not_deg, "gene_name"]
        relevant_genes <- c(candidate_not_deg, DGH_results[all_degs, 
                                                            "gene_name"]) 

    }else {
        relevant_genes <- DGH_results[all_degs, "gene_name"] 
    }

    DGH_data$DH_results$relevant_genes <- DGH_results$gene_name %in% 
                                             relevant_genes
    DGH_data$normalized_counts <- DGH_data$normalized_counts[, 
                    colnames(DGH_data$normalized_counts) %in% relevant_genes ]
                    
    return(DGH_data)
}

#' @importFrom data.table as.data.table 
load_and_parse_multimir <- function(
    multimir_table = NULL, 
    multimir_path = NULL,
    selected_predicted_databases = c("diana_microt", "elmmo", "microcosm", 
                        "miranda","mirdb", "pictar", "pita", "targetscan"),
    database_to_filter = NULL,
    filter_db_theshold = 0
    ) {
    multimir_summary <- NULL
    if (is.null(multimir_table) && file.exists(multimir_path)) {
        load(multimir_path)
        multimir_table <- multimir_summary 
    } 
    if (is.null(multimir_table)){
        stop("You have to give a parsed multiMiR data frame or an RData path")
    }
    multimir_table <- as.data.frame(multimir_table)
    dbs <- list(
    "validated_c" = c("mirecords", "mirtarbase", "tarbase"),
    "predicted_c" = selected_predicted_databases
    )
    multimir_table[["validated_c"]] <- rowSums(multimir_table[,
        dbs[["validated_c"]]], na.rm = TRUE)
    multimir_table[["predicted_c"]] <- rowSums(multimir_table[,
        dbs[["predicted_c"]]] > 0, na.rm = TRUE)
    gc()
    if ( ! is.null(database_to_filter) && filter_db_theshold > 0) {
      multimir_table <-  multimir_table[
      rowSums(!is.na(multimir_table[,database_to_filter])
      ) > filter_db_theshold | multimir_table$validated_c > 0,]
    }
    raw_databases_scores <- lapply(selected_predicted_databases, 
      function(database){
          database_scores <- multimir_table[,database]
          database_scores <- database_scores[!is.na(database_scores)]
      })
    names(raw_databases_scores) <- selected_predicted_databases
    return(list(raw_databases_scores = raw_databases_scores,
    multimir_table = data.table::as.data.table(multimir_table)))
}

parse_strategies <- function(strategies){
    strat_equivalence <- list(
    "d" = "normalized_counts",
    "E" = "Eigengene",
    "h" = "hub_1")
    parsed_strategies <- c()
    for (strategy in strategies) {
      splitted_strategy <- unlist(strsplit(strategy, ""))
      RNA_profile <- strat_equivalence[[splitted_strategy[1]]]
      miRNA_profile <- strat_equivalence[[splitted_strategy[2]]]
      parsed_strategies <- c(parsed_strategies, paste0(RNA_profile, 
          "_RNA_vs_miRNA_", miRNA_profile))
    }
    return(parsed_strategies)
}


#' @importFrom utils read.table
load_selected_targets <- function(selected_targets_file){
    selected_targets <- NULL
    if (!is.null(selected_targets_file))
        selected_targets <- utils::read.table(selected_targets_file)[,1] 
    return(selected_targets)
}

#' @importFrom data.table as.data.table
#' @importFrom tidyr unite
#' @importFrom plyr rbind.fill
perform_all_strategies <- function(
strat_names, 
RNAseq, 
miRNAseq, 
all_pairs,
tag_filter,
selected_predicted_databases,
raw_databases_scores,
MM_cutoff = 0.7,
permutations = 10, 
sample_proportion = 0.01,
selected_targets = c(),
corr_type,
corr_coef,
compare_pred_scores = FALSE,
eval_method = "aggregated"
){
 ######### PREPARE PAIRS SCAFOLD
 if (is.null(RNAseq)) {
    if(!is.null(selected_targets)){
        strat_names <- c("selected_targets_RNA_vs_miRNA_DEMs")
        message(paste0("RNAseq folder has been not set, only ",
 "selected_targets_RNA_vs_miRNA_DEMs strategy will be performed"))
    } else {
        strat_names <- c("all_posible_targets_RNA_vs_miRNA_DEMs")
        message(paste0("RNAseq folder has been not set, only ",
        " strategy will be performed"))
        selected_targets <- unique(all_pairs$RNAseq)
    }
} 
 strategies <- list()
 tag_filter <- unlist(strsplit(tag_filter, ","))
 if (length(tag_filter) == 1)
 tag_filter <- c(tag_filter,tag_filter)
 ###### FILTER DGH DATA
 RNAseq <- filter_DEGH_data(RNAseq, MM_cutoff, tag_filter[1])
 miRNAseq <- filter_DEGH_data(miRNAseq, MM_cutoff, tag_filter[2])
 message("Data has been filtered")
 gc()
 std_positions <- tidyr::unite(all_pairs, "pairs", RNAseq:miRNAseq, sep = "_")
 std_positions <- std_positions$pairs
 

 ###### BUILD EMPTY TABLES
 #Vector with strategies without significant results
 unsig_strategies <- c() 
 #Dataframe with all contingency tables
 cont_tables <- data.frame()  
 #Datafarme with stats of strategy pairs and randoms 
 miRNA_cont_tables <- data.frame()
 filters_summary <- data.frame()
 #Dataframe with stats from prediction scores comparison
 score_comp<- data.frame()
 #Dataframe with weighted fisher pvalues
 p_fisher <- data.frame()

 ####PERFORM STRATEGIES
 for(strategy in strat_names){ 
   message(paste0("\n", strategy, " strategy is been performed"))
    message("Computing stats by group")
    strategy_data <- perform_correlations(strategy = strategy, 
                    RNAseq = RNAseq, miRNAseq = miRNAseq, 
                    std_positions = std_positions, corr_coef = corr_coef,
                    corr_type  = corr_type, selected_targets = selected_targets)



   strategy_data <- data.table::as.data.table(strategy_data)
   
   # if (sum(strategy_data$correlated_pairs) == 0){
   #     unsig_strategies <- c(strategy, unsig_strategies)
   #     next
   # }
   strategies[[strategy]] <- strategy_data
  
   # all_pairs[,strategy] <- FALSE
   # all_pairs[strategy_data$pair_n, strategy] <- strategy_data$correlated_pairs
   
   all_pairs[strategy_data$pair_n, 
              paste0(strategy,"_correlation")] <- strategy_data$correlation
   all_pairs[strategy_data$pair_n, 
              paste0(strategy,"_pval")] <- strategy_data$pval
  
 }

 all_pairs <- data.table::as.data.table(all_pairs)
 # if (all(names(strategies) %in% unsig_strategies)){
 #  stop("There is not significant pairs for any strategy selected.")
 # } 

  # all_cor_dist <- parse_correlations(strategies)
  # sig_pairs <- get_sig_pairs(strategies, all_pairs)
 miRNA_cor_results <- list(strategies = strategies,
    # unsig_strategies =  unsig_strategies, 
   all_pairs = all_pairs 
   #all_cor_dist = all_cor_dist,   sig_pairs = sig_pairs
   )
 return(miRNA_cor_results)
}





get_degs_sets <- function(DH_results_table){

     degs <- list()
     degs[["pos"]] <- DH_results_table[DH_results_table$genes_tag %in%
               c("PREVALENT_DEG", "POSSIBLE_DEG") & 
                DH_results_table$mean_logFCs > 0, "gene_name"]
     degs[["neg"]] <- DH_results_table[DH_results_table$genes_tag %in%
               c("PREVALENT_DEG", "POSSIBLE_DEG") &
                DH_results_table$mean_logFCs < 0, "gene_name"]  
    return(degs)
}


prepare_for_stats <- function(
all_pairs, 
strategy, 
corr_cutoff, 
p_val_cutoff,
db_groups = c("predicted", "validated", "multimir"),
sig_pairs = FALSE,
crossval = FALSE, 
test_sample = 0.25,
corr_type = "lower"
){
 cont_tables <- data.frame()
 miRNA_cont_tables <- data.frame()
 # message("Computing contingency tables")

    strat_data <- all_pairs[,c("miRNAseq", "predicted", "validated", "multimir")]
    
    if (!sig_pairs){
        if (corr_type == "lower"){

        strat_data$sig_pairs <- all_pairs[,paste0(strategy,
                                            "_correlation")] <= corr_cutoff &
                           all_pairs[,paste0(strategy,
                                            "_pval")] <= p_val_cutoff &
                           all_pairs$known_miRNA
        } else if (corr_type == "higher"){
            strat_data$sig_pairs <- all_pairs[,paste0(strategy,
                                            "_correlation")] >= corr_cutoff &
                           all_pairs[,paste0(strategy,
                                            "_pval")] <= p_val_cutoff &
                           all_pairs$known_miRNA
        }

    } else {
        strat_data$sig_pairs <- all_pairs[,strategy]
    }

    strat_data$sig_pairs[is.na(strat_data$sig_pairs)] <- FALSE

    if (sum(strat_data$sig_pairs) == 0) return(NULL)
 
    strat_data <- data.table::as.data.table(strat_data)

    all_cts <- get_ct_by_group(strat_data = strat_data, 
                              db_groups = db_groups,
                              crossval = crossval, test_sample = test_sample)
    strat_miRNA_ct <-  all_cts$strat_miRNA_cts 
    strat_ct <-  all_cts$strat_cts
       
       if (nrow(strat_miRNA_ct) > 0) {
          strat_miRNA_ct$strategy <- strategy
          strat_miRNA_ct$corr_cutoff <- corr_cutoff
       }

       if (nrow(strat_miRNA_ct) > 0) {
          strat_ct$strategy <- strategy
          strat_ct$corr_cutoff <- corr_cutoff
      }

      return(list( strat_miRNA_ct = strat_miRNA_ct,
        strat_ct = strat_ct
        ))
}




get_cor_dist <- function(all_pairs, strategies){
    corr_dist <- data.frame()
    for (strategy in strategies){
        strat_corr <- data.frame(
        correlation = all_pairs[,paste0(strategy, "_correlation")],
        pval = all_pairs[,paste0(strategy, "_pval")],
        strategy = strategy)

        strat_corr <- strat_corr[!is.na(strat_corr$correlation),]   
        corr_dist <- rbind(corr_dist, strat_corr)
    }
    return(corr_dist)
}

#' @importFrom data.table as.data.table
perform_correlations <- function(
      strategy = "normalized_counts_RNA_vs_miRNA_normalized_counts", 
      RNAseq, miRNAseq, std_positions = NULL, corr_coef,
      corr_type, selected_targets = c()){ #correct_positions 
    #is a mirna_RNA pairs vector
    # the parsed strategy name text is used to subset RNAseq/miRNAseq obects
    
    if (strategy == "DEGs_RNA_vs_miRNA_DEMs_opp"){
      #   q()   
      
       DEGs <- get_degs_sets(RNAseq$DH_results)
       DEMs <- get_degs_sets(miRNAseq$DH_results)
        
       all_pairs <- permute_pairs(DEGs[["pos"]], DEMs[["neg"]], 
                                    corr_type = corr_type)
       all_pairs <- rbind(all_pairs, permute_pairs(DEGs[["neg"]], DEMs[["pos"]], 
                                                    corr_type = corr_type))

    } else if ( strategy == "DEGs_RNA_vs_miRNA_DEMs_sim" ){
       DEGs <- get_degs_sets(RNAseq$DH_results)
       DEMs <- get_degs_sets(miRNAseq$DH_results)
       
       all_pairs <- permute_pairs(DEGs[["pos"]], DEMs[["pos"]],
                                     corr_type = corr_type)
       all_pairs <- rbind(all_pairs, permute_pairs(DEGs[["neg"]], DEMs[["neg"]], 
                                                        corr_type = corr_type))
    
    } else if ( strategy %in% c("selected_targets_RNA_vs_miRNA_DEMs",
                                "all_posible_targets_RNA_vs_miRNA_DEMs")) {
        DEMS <- miRNAseq$DH_results[miRNAseq$DH_results$genes_tag %in% 
               c("PREVALENT_DEG", "POSSIBLE_DEG"), "gene_name"]    
        all_pairs <- permute_pairs(selected_targets, DEMS, corr_type = corr_type)
    } else {
        strat_description <- unlist(strsplit(strategy, "_RNA_vs_miRNA_"))
        RNA_profiles <- as.matrix(RNAseq[[strat_description[1]]])
        miRNA_profiles <- as.matrix(miRNAseq[[strat_description[2]]])
        if (ncol(RNA_profiles) == 0 || ncol(miRNA_profiles) == 0) {
          return(data.frame(NULL))
        }
        all_pairs <- correlate_profiles(RNA_profiles, miRNA_profiles, 
                                        corr_coef= corr_coef)
        if (strat_description[1] != "normalized_counts"){
            if (strat_description[1] == "Eigengene_0") {
                RNAseq$DH_results$relevant_genes <- TRUE
            }
            all_pairs <- expand_module(all_pairs = all_pairs, tag = "RNAseq", 
                DH_results = 
                RNAseq$DH_results[RNAseq$DH_results$relevant_genes, ])
        }
        if (strat_description[2] != "normalized_counts"){
            if (strat_description[2] == "Eigengene_0") {
                miRNAseq$DH_results$relevant_genes <- TRUE
            }
            all_pairs <- expand_module(all_pairs = all_pairs, tag = "miRNAseq", 
                DH_results = 
                miRNAseq$DH_results[miRNAseq$DH_results$relevant_genes, ])
        }
    }
    all_pairs <- data.table::as.data.table(all_pairs)

    # if(corr_type == "higher"){
    #         all_pairs$correlated_pairs <- all_pairs$correlation >= cor_cutoff & 
    #                                  all_pairs$pval < pval_cutoff
    # } else if ( corr_type == "lower"){
    #         all_pairs$correlated_pairs <- all_pairs$correlation <= cor_cutoff & 
    #                                 all_pairs$pval < pval_cutoff
    # } else {
    #     stop("Non valid --corr_type arguments")
    # }
    
     # Add extra column to indicate number of pair in std_positions
    all_pairs_str <- all_pairs[,c("RNAseq", "miRNAseq")]
    all_pairs_str <- tidyr::unite(all_pairs_str, "pairs", RNAseq:miRNAseq, 
                               sep = "_")
    all_pairs_str <- all_pairs_str$pairs  
    all_pairs$pair_n <- match(all_pairs_str, std_positions) 
    all_pairs[,c("RNAseq", "miRNAseq")] <- NULL

    return(all_pairs)
}


permute_selected_targets <- function(miRNAseq, selected_targets,corr_type){
    DEMS <- miRNAseq[miRNAseq$genes_tag %in% 
     c("PREVALENT_DEG", "POSSIBLE_DEG"), "gene_name"]
     all_pairs <- expand.grid(selected_targets, DEMS, 
                               stringsAsFactors= FALSE, 
                               KEEP.OUT.ATTRS = FALSE)
    colnames(all_pairs) <- c("RNAseq", "miRNAseq")
    if (corr_type == "higher"){
       all_pairs$correlation <- 1
    } else if (corr_type == "lower"){
       all_pairs$correlation <- -1
    }
    all_pairs$pval <- 0
    return(all_pairs)
}


permute_pairs <- function(RNAseq, miRNAseq, corr_type){
    all_pairs <- expand.grid(RNAseq, miRNAseq, 
                               stringsAsFactors= FALSE, 
                               KEEP.OUT.ATTRS = FALSE)
    colnames(all_pairs) <- c("RNAseq", "miRNAseq")
    if (corr_type == "higher"){
       all_pairs$correlation <- 1
    } else if (corr_type == "lower"){
       all_pairs$correlation <- -1
    }
    all_pairs$pval <- 0

    return(all_pairs)
}

#' @importFrom WGCNA corPvalueStudent
#' @importFrom stats cor
correlate_profiles <- function(RNA_profiles, miRNA_profiles,corr_coef) {

    if (corr_coef == "bicor") {
        correlations <- WGCNA::bicor(RNA_profiles, miRNA_profiles, 
                                maxPOutliers = round(ncol(RNA_profiles)/5, 0))    
    } else if (corr_coef %in% c("pearson", "spearman", "kendall")){

        correlations <- WGCNA::cor(RNA_profiles, miRNA_profiles, 
            method = corr_coef)
    }
    
    nSamples <- nrow(RNA_profiles)
    cor_pval <- as.data.frame(WGCNA::corPvalueStudent(
        as.matrix(correlations), nSamples))

    names_matrix <- as.matrix(correlations)  
    all_pairs <- as.data.frame(as.table(correlations)) 
    colnames(all_pairs) <- c("RNAseq", "miRNAseq", "correlation")
    
    all_pairs$correlation <- round(all_pairs$correlation, 2)

    all_pairs$pval <- as.data.frame(as.table(as.matrix(cor_pval)))$Freq

    all_pairs$miRNAseq <- as.character(all_pairs$miRNAseq)
    all_pairs$RNAseq <- as.character(all_pairs$RNAseq)
    return(all_pairs)
}



#' @importFrom dplyr desc between row_number filter arrange group_by
#' @importFrom rlang .data
get_hub_genes_by_MM <- function(normalized_counts, hunter_results, top = 1){
    #maybe we can change this function by  splitting all DGH by modules
    "%>%" <- magrittr::"%>%"
    hub_genes <- hunter_results %>% 
    dplyr::filter(.data$Cluster_MM_pval <= 0.05) %>% 
    dplyr::arrange(dplyr::desc(.data$Cluster_MM)) %>% 
    dplyr::group_by(.data$Cluster_ID) %>% 
    dplyr::filter(dplyr::between(dplyr::row_number(), 1, top))
    hub_genes_profile <- as.matrix(
      normalized_counts[,colnames(normalized_counts) %in% hub_genes$gene_name])
    
    colnames(hub_genes_profile) <- as.data.frame(hub_genes)[match(
                                            colnames(hub_genes_profile), 
                                            hub_genes$gene_name), "Cluster_ID"]
    return(hub_genes_profile)
}


#' @importFrom data.table as.data.table merge.data.table  
#' @importFrom rlang .data
expand_module <- function(all_pairs, tag, DH_results){
    mod_tag <- paste0(tag, "_mod")
    names(all_pairs)[names(all_pairs)== tag] <- mod_tag

    partial_expanded <- DH_results[,c("Cluster_ID", "gene_name")]
    partial_expanded <- data.table::as.data.table(partial_expanded)
    colnames(partial_expanded) <- c("module", tag)
    partial_expanded$module <- as.character(partial_expanded$module)
    all_pairs <- data.table::as.data.table(all_pairs)
    all_pairs <- data.table::merge.data.table(x = all_pairs, 
        y = partial_expanded, by.x = mod_tag, 
        by.y = "module", allow.cartesian  = TRUE)
    all_pairs[,mod_tag] <- NULL
    return(all_pairs)
}

#' @importFrom data.table as.data.table merge.data.table
add_multimir_info <- function(all_pairs, multiMiR = NULL) {
    all_pairs <- data.table::as.data.table(all_pairs)
    
    if (!is.null(multiMiR)){
     multiMiR <- data.table::as.data.table(multiMiR)
     all_pairs <- data.table::merge.data.table(all_pairs, 
        multiMiR, by.x = c("RNAseq","miRNAseq"), 
        by.y =c("target_ensembl","mature_mirna_acc"), all.x = TRUE, 
        suffixes = c("", ""), no.dups = TRUE)
    }
    all_pairs <- as.data.frame(all_pairs)
    for (pred_db in c("predicted_c", "validated_c")){

       all_pairs[is.na(all_pairs[,pred_db]), pred_db] <- 0

    }
    for (pred_db in c("mirecords", "mirtarbase", "tarbase")) {
       all_pairs[is.na(all_pairs[,pred_db]), pred_db] <- FALSE
    }


    all_pairs$predicted <- all_pairs$predicted_c > 0
    all_pairs$validated <- all_pairs$validated_c > 0
    all_pairs$pred_and_val <- all_pairs$predicted & all_pairs$validated
    all_pairs$multimir <- all_pairs$predicted | all_pairs$validated
    gc()
    all_pairs <- data.table::as.data.table(all_pairs)
 
    return(all_pairs)
}


#' @importFrom data.table as.data.table merge.data.table
get_entrez_symbol_translations <- function(ensembl_ids, organism_info){
    org_db <- get_org_db(organism_info)
    input_to_entrezgene <- translate_ids_orgdb(ids = ensembl_ids, 
                    input_id="ENSEMBL", output_id = "ENTREZID", org_db=org_db)
    # Fix names
    colnames(input_to_entrezgene) <- c("ensembl_gene_id", "entrezgene")     
    
    input_to_symbol <- translate_ids_orgdb(input_id="ENTREZID", 
                    output_id = "SYMBOL", org_db=org_db,
                    ids = unique(input_to_entrezgene$entrezgene))
    colnames(input_to_symbol) <- c("entrezgene", "Symbol")

    input_to_entrezgene <- data.table::as.data.table(input_to_entrezgene)
    input_to_symbol <- data.table::as.data.table(input_to_symbol)
    gene_id_translation <- data.table::merge.data.table(input_to_entrezgene, 
                             input_to_symbol, by = "entrezgene", all.x = TRUE)

    return(gene_id_translation)
}

score_comparison <- function(databases, all_pairs, sample_size, strategy){
  db_comp <- data.frame()
  all_pairs <- data.table::as.data.table(all_pairs)
  for (database in databases) {
    bckg_scores <- all_pairs[, get(database)]
    bckg_scores <- bckg_scores[!is.na(bckg_scores)]
    strat_scores <- all_pairs[get(strategy), get(database)] #Revise that
    strat_scores <- strat_scores[!is.na(strat_scores)]

    # sampled_bckg <- sample(bckg_scores, 
    #         size = 100, replace = FALSE)
    res <- data.frame(p.value = 1, boot.p.value = 1)
    if (length(strat_scores) > 1) {
      # res <- stats::t.test(x = strat_scores, y = bckg_scores, 
        res <- boot.t.test.2(x = strat_scores, y = bckg_scores, 
          alternative = "greater", 
          vectorial = FALSE, 
          R = 1000,
          chunk_size = 100)
    } 
    actual_db_comp <- data.frame(
                        database = database,
                        strategy = strategy,
                        boot.p.value = res$boot.p.value,
                        p.value = res$p.value)
    if (actual_db_comp$boot.p.value == 0){
      if(actual_db_comp$p.value == 0){
        actual_db_comp[,c("boot.p.value", "p.value")] <- 1e-16
      } else {
        actual_db_comp$boot.p.value <- actual_db_comp$p.value
      }
    }
    db_comp <- rbind(db_comp,actual_db_comp) 
  }
  db_comp <- rbind(
              db_comp, 
              data.frame(
                database = "comb_pval", strategy = strategy,
                boot.p.value = vectorial_fisher_method(db_comp$boot.p.value),
                p.value = vectorial_fisher_method(db_comp$p.value)
                ))
  return(db_comp)
}

get_miRNA_ct <- function(strat_data, db_group, crossval = FALSE, test_sample = 0.25, bootstrap = 40){
    db_gs <- strat_data[,get(db_group)]
    experiment_pairs <- strat_data[,sig_pairs]

    if (sum(experiment_pairs) == 0) return(NULL)

    if (crossval) {
        strat_group_ct <- data.frame()
        seed <- 1234
        for (bt_it in seq(1, bootstrap)){
                set.seed(seed + bt_it)
                db_gs_trues <- which(db_gs)
                db_gs_to_falsify <- sample(db_gs_trues, ceiling(test_sample*length(db_gs_trues)))

                db_gs[db_gs_to_falsify] <- FALSE
                              
                train_ct <- calc_ct(experiment = experiment_pairs, gold_standard = db_gs)
                train_ct$crossval_it <- bt_it
                train_ct$crossval_status <- "train"

                strat_group_ct <- rbind(strat_group_ct, train_ct)
        }

    } else {
     strat_group_ct <- calc_ct(experiment = experiment_pairs, 
                              gold_standard = db_gs)
    
    }




    strat_group_ct$db_group <- db_group
    return(strat_group_ct)
}


#' @importFrom plyr rbind.fill
#' @importFrom data.table setnames
get_ct_by_group <- function(
strat_data, 
db_groups = c("predicted", "validated", "multimir"),
crossval = FALSE, test_sample = 0.25
){
 strat_cts <- data.frame() #contingency table
 strat_miRNA_cts <- data.frame() #contingency table for every miRNA

 for (db_group in db_groups){

    strat_group_ct <- get_miRNA_ct(strat_data, db_group)
    strat_cts <- rbind(strat_cts, strat_group_ct)
   all_miRNA_ct <- lapply(split(strat_data, strat_data$miRNAseq), 
                            get_miRNA_ct, db_group = db_group, crossval = crossval,
                            test_sample = test_sample)
    strat_miRNA_cts <- rbind(strat_miRNA_cts, 
                              data.table::rbindlist(all_miRNA_ct, 
                                                    use.names = TRUE, 
                                                    idcol = "miRNA"))
 }
 return(list(strat_cts= strat_cts, 
    strat_miRNA_cts = strat_miRNA_cts))
}


prepare_all_pairs <- function(
RNAseq, 
miRNAseq, 
multimir,
selected_targets,
organism){
 
 

 if (is.null(RNAseq)){
    if (!is.null(selected_targets)) {
        putative_targets <- selected_targets
    } else {
        putative_targets <- unique(multimir[multimir$mature_mirna_acc  %in% 
                    miRNAseq$DH_results$gene_name, "target_ensembl"])[[1]]
    }
 } else {
    putative_targets <- colnames(RNAseq$normalized_counts)
 }
 all_pairs <- expand.grid(
                putative_targets, 
                colnames(miRNAseq$normalized_counts), 
                stringsAsFactors = TRUE)
 names(all_pairs) <- c("RNAseq", "miRNAseq")
 all_pairs <- add_multimir_info(all_pairs = all_pairs, 
                                multiMiR = multimir)

 
 all_pairs$known_miRNA <- grepl("MIMAT", all_pairs$miRNAseq)
 
 #possible_positives: Pairs with RNA and miRNA included in multiMiR. 
 possible_positives <- all_pairs$known_miRNA & 
                                 all_pairs$miRNAseq %in% 
                                  multimir$mature_mirna_acc &
                                 all_pairs$RNAseq %in% 
                                   multimir$target_ensembl 

 all_pairs$predicted <- all_pairs$predicted & possible_positives
 all_pairs$validated <- all_pairs$validated & possible_positives
 all_pairs$multimir <- all_pairs$predicted | all_pairs$validated

 return(all_pairs)
}

#' @importFrom miRBaseVersions.db miRBaseVersions.db
#' @importFrom AnnotationDbi select
translate_all_id <- function(
miRNA_IDs, 
RNA_IDs, 
organism_info, 
translate_ensembl){
 mirna_names <- unique(miRNA_IDs)
 mirna_names <- mirna_names[grepl("MIMAT", mirna_names)]
 mirna_names <-  AnnotationDbi::select(
                            miRBaseVersions.db::miRBaseVersions.db, 
                            keys = mirna_names, keytype =  "VW-MIMAT-22.0",
                            columns = c("ACCESSION", "NAME"))
 gene_id_translation <- NULL
 if (! organism_info$Bioconductor_VarName_SYMBOL[1] == "" && 
     translate_ensembl) {
     gene_id_translation <- get_entrez_symbol_translations(
         ensembl_ids = unique(RNA_IDs), 
         organism_info = organism_info)
 } else {
    message("RNA IDs were not translated")
 }
 return(list(mirna_names = mirna_names, 
    gene_id_translation= gene_id_translation))
}

#' @importFrom miRBaseVersions.db miRBaseVersions.db
#' @importFrom AnnotationDbi select
translate_miRNA_ids <- function(miRNA_IDs, miRBaseVersion = "VW-MIMAT-22.0") {
    miRNA_tr_table <- AnnotationDbi::select(
                            miRBaseVersions.db::miRBaseVersions.db, 
                            keys = miRNA_IDs, keytype = miRBaseVersion,
                            columns = c("ACCESSION", "NAME"))
    translated_miRNA <- miRNA_tr_table[match(miRNA_IDs, 
                                        miRNA_tr_table$ACCESSION) ,"NAME"]
    return(translated_miRNA)
}

#' @importFrom data.table rbindlist
parse_correlations <- function(strategies){
 all_cor_dist <- strategies
 all_cor_dist["DEGs_DEMs_permutated"] <- NULL
 all_cor_dist <- lapply(all_cor_dist, function(strategy_data){
    dist <- strategy_data[,c("correlation", "pval", "correlated_pairs")]
    return(dist)
 })
 all_cor_dist <- as.data.frame(data.table::rbindlist(all_cor_dist, 
                                                    use.names=TRUE, 
                                                    idcol= "strategy"))
 return(all_cor_dist)
}
 
get_sig_pairs <- function(strategies, all_pairs){
  groups_db <- c("all","multimir", "predicted", "validated", "pred_and_val")
  all_overlap <- list()
  strat_names <- names(strategies)
  strat_names <- strat_names[strat_names %in% colnames(all_pairs)]
  for (group_db in groups_db){
    group_matrix <- matrix(
      data = NA, nrow = length(strat_names), ncol = length(strat_names),
      dimnames = list(strat_names, strat_names))
    for(xstrategy in strat_names){
      for(ystrategy in strat_names){
        group <- rep(TRUE, nrow(all_pairs))
        if (group_db != "all"){
            group <- all_pairs[, get(group_db)]
        }
        group_matrix[xstrategy, ystrategy] <- sum(
          all_pairs[,get(xstrategy)] & all_pairs[,get(ystrategy)] & group)
      }
    }
    all_overlap[[group_db]] <- group_matrix
  }
  names(all_overlap) <- groups_db
  return(all_overlap)
}


calc_pred_ctables <- function(
 all_pairs,
 strategy,
 selected_predicted_databases)
{
 cont_tables <- data.frame()
 for (database in selected_predicted_databases){
   c_table <- calc_ct(experiment = all_pairs[,get(strategy)],
                       gold_standard = !is.na(all_pairs[,get(database)]))
   c_table$database <- database
   cont_tables <- rbind(cont_tables, c_table)
 }
 return(cont_tables)
}





select_best_strategy <- function(miRNA_cont_tables){
    #taken from https://stackoverflow.com/questions/24558328/select-the-row-with-the-maximum-value-in-each-group
    max_mIRNA_ID <- with(miRNA_cont_tables, ave(Odds_ratio, miRNA, FUN=function(x) seq_along(x) == which.max(x))) == 1
    miRNA_cont_tables_max <- miRNA_cont_tables[max_mIRNA_ID,]
    return(miRNA_cont_tables_max)
}


set_strats_readable <- function(strategies){
dictionary <- list(
    "Eigengene" = "E",
    "Eigengene_0" = "E0",
    "hub_1" = "H",
    "normalized_counts" = "C"
    )
  translated_strategies <- c()
  # str(dictionary)
  for (strategy_name in strategies){

    strategies <- unlist(strsplit(strategy_name, "_RNA_vs_miRNA_"))
    summ_name <- c()
    for (strategy in strategies) {
        if (!is.null(dictionary[[strategy]])){
            summ_name <- c(summ_name, dictionary[[strategy]])
        } else {
            if(stringr::str_detect(strategy, "_opp")) strategy <- gsub("_opp", " (opp)", strategy)
            if(stringr::str_detect(strategy, "_sim")) strategy <- gsub("_sim", " (sim)", strategy)
            summ_name <- c(summ_name, strategy)
        }
    }
    translated_str <- paste0(summ_name, c("g ","m"), collapse = "")
    translated_strategies <- c(translated_strategies, translated_str)
  }
  return(translated_strategies)
}



 get_optimal_pairs <- function(all_pairs, miRNA_cont_tables, pval_thr, corr_type){
    all_pairs <- as.data.frame(all_pairs)
   selected_miRNA_strat <- select_best_strategy(miRNA_cont_tables)
   miRNA_strat_pairs <- rep(FALSE, nrow(all_pairs))
   for(selected_miRNA in selected_miRNA_strat$miRNA) {
         best_miRNA_strat <- selected_miRNA_strat[selected_miRNA_strat$miRNA == selected_miRNA,]
      strat_corr <- paste0(best_miRNA_strat$strategy, "_correlation")
      strat_pval <- paste0(best_miRNA_strat$strategy, "_pval")
      if (corr_type == "lower"){
             miRNA_best_pairs  <- all_pairs$miRNAseq == selected_miRNA &
                        all_pairs[,strat_corr] <= best_miRNA_strat$corr_cutoff &
                        all_pairs[,strat_pval] <= pval_thr
      } else if (corr_type == "higher"){
                 miRNA_best_pairs  <- all_pairs$miRNAseq == selected_miRNA &
                        all_pairs[,strat_corr] >= best_miRNA_strat$corr_cutoff &
                        all_pairs[,strat_pval] <= pval_thr
      }
      miRNA_strat_pairs <- miRNA_strat_pairs | miRNA_best_pairs 
   }
   return(miRNA_strat_pairs)
 }



 add_attrib_to_pairs <- function(pairs, RNAseq = NULL, miRNAseq = NULL, translate_targets = FALSE){
    fields_to_extract <- c("miRNA", "miRNAseq", "RNAseq", "validated_c","predicted_c")
    if ("miRNA_loci" %in% colnames(pairs))
        fields_to_extract <- c(fields_to_extract, "miRNA_loci")  
 
       output_pairs <- pairs[,fields_to_extract]

    if (!is.null(RNAseq)) {
        RNA_data <- RNAseq$DH_results[,c("gene_name","mean_logFCs","genes_tag")]
        colnames(RNA_data) <- c("RNAseq", "Target_log2FC", "DEG_tag_target")
        output_pairs <- merge(output_pairs, RNA_data, by = "RNAseq", all.x = TRUE)
    }
    if (!is.null(miRNAseq)) {
        miRNA_data <- miRNAseq$DH_results[,c("gene_name","mean_logFCs","genes_tag")]
        colnames(miRNA_data) <- c("miRNAseq", "miRNA_log2FC", "DEM_tag_miRNA")
        output_pairs <- merge(output_pairs, miRNA_data, by = "miRNAseq", all.x = TRUE)
    }
    colnames(output_pairs)[match(c("miRNAseq", "RNAseq", "validated_c","predicted_c"),
                                   colnames(output_pairs))] <- 
            c("miRNA_ID", "Target_ID","Validated_DB_count", "Predicted_DB_count")

    output_pairs$All_DB_count <- output_pairs$Validated_DB_count + output_pairs$Predicted_DB_count


    if ("miRNA_loci" %in% colnames(output_pairs))
    colnames(output_pairs)[colnames(output_pairs) == "miRNA_loci"] <- "miRNA_source_gene"
    
    return(output_pairs)
}
