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


filter_DEGH_data <- function(DGH_data, MM_cutoff){
    DGH_results <- DGH_data$DH_results 
    all_degs <- DGH_results$genes_tag %in% c("PREVALENT_DEG", "POSSIBLE_DEG")
    modules_with_DEGS <- unique(DGH_results[all_degs & 
        DGH_results$Cluster_MM >= MM_cutoff, "Cluster_ID"])

    candidate_not_deg <- DGH_results$Cluster_ID %in% modules_with_DEGS &
                                       DGH_results$genes_tag == "NOT_DEG" &
                                       DGH_results$Cluster_MM >= MM_cutoff

    candidate_not_deg <- DGH_results[candidate_not_deg, "gene_name"]

    relevant_genes <- c(candidate_not_deg, DGH_results[all_degs, "gene_name"]) 

    DGH_data$DH_results$relevant_genes <- DGH_results$gene_name %in% 
                                             relevant_genes
    
    DGH_data$normalized_counts <- DGH_data$normalized_counts[, 
                    colnames(DGH_data$normalized_counts) %in% relevant_genes ]


    DGH_data$Eigengene_0 <- DGH_data$Eigengene[,"0", drop = FALSE]

    DGH_data$Eigengene <- DGH_data$Eigengene[,
               colnames(DGH_data$Eigengene) %in% modules_with_DEGS]
    DGH_data$hub_1 <- DGH_data$hub_1[,
               colnames(DGH_data$hub_1) %in% modules_with_DEGS]

    return(DGH_data)
}


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



#' @importFrom data.table as.data.table
#' @importFrom tidyr unite
#' @importFrom plyr rbind.fill
perform_all_strategies <- function(
strat_names, 
RNAseq, 
miRNAseq, 
all_pairs,
selected_predicted_databases,
raw_databases_scores,
corr_cutoff = -0.8, 
p_val_cutoff = 0.05,
MM_cutoff = 0.7,
permutations = 10, 
sample_proportion = 0.01
){
 ######### PREPARE PAIRS SCAFOLD
 strategies <- list()
 
 ###### FILTER DGH DATA
 RNAseq <- filter_DEGH_data(RNAseq, MM_cutoff)
 miRNAseq <- filter_DEGH_data(miRNAseq, MM_cutoff)
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
 filters_summary <- data.frame()
 #Dataframe with stats from prediction scores comparison
 score_comp<- data.frame()
 #Dataframe with weighted fisher pvalues
 p_fisher <- data.frame()

 ####PERFORM STRATEGIES
 for(strategy in strat_names){ 
   message(paste0("\n", strategy, " strategy is been performed"))
    strategy_data <- perform_correlations(strategy = strategy, 
                        RNAseq = RNAseq, miRNAseq = miRNAseq, 
                        std_positions = std_positions, 
                        cor_cutoff = corr_cutoff, pval_cutoff = p_val_cutoff)
   strategy_data <- data.table::as.data.table(strategy_data)
   strategies[[strategy]] <- strategy_data
   if (sum(strategy_data$correlated_pairs) == 0){
       unsig_strategies <- c(strategy, unsig_strategies)
       next
   }
   all_pairs[,strategy] <- FALSE
   all_pairs[strategy_data$pair_n, strategy] <- strategy_data$correlated_pairs
   all_pairs[strategy_data$pair_n, 
              paste0(strategy,"_correlation")] <- strategy_data$correlation

   all_stats <- get_stats_by_group(all_pairs = all_pairs, strategy = strategy,
                                   permutations = permutations)
   filters_summary <- plyr::rbind.fill(filters_summary,all_stats$strat_summary)
   cont_tables <- rbind(cont_tables, all_stats$strat_cts)
   all_stats <- NULL
   pred_cont_tables <- calc_pred_ctables(
                 all_pairs =all_pairs, 
                 strategy= strategy,
                 selected_predicted_databases = selected_predicted_databases)
   # , 
   #                          raw_databases_scores = raw_databases_scores)
 
   pred_cont_tables$fisher.p.value <- v.fisher.test(pred_cont_tables)
   pred_fisher <- pred_cont_tables[,c("database", "fisher.p.value")]
   pred_comb_pval <- vectorial_fisher_method(pred_fisher$fisher.p.value)
   pred_fisher <- rbind(pred_fisher, 
                        data.frame(database = "comb_pval",
                                   fisher.p.value = pred_comb_pval
                              ))
   pred_fisher$strategy <- strategy
   p_fisher <- rbind(p_fisher, pred_fisher)
   score_comp <- rbind(score_comp,
                       score_comparison(
                        databases = selected_predicted_databases,  
                        all_pairs = all_pairs,
                        sample_size = sample_proportion, strategy = strategy))
   gc()
 }
 if (all(names(strategies) %in% unsig_strategies)){
  stop("There is not significant pairs for any strategy selected.")
 } 
  all_cor_dist <- parse_correlations(strategies)
  sig_pairs <- get_sig_pairs(strategies, all_pairs)
 miRNA_cor_results <- list(strategies = strategies,
     unsig_strategies =  unsig_strategies, cont_tables = cont_tables,
     filters_summary = filters_summary, score_comp = score_comp,
     all_pairs = all_pairs, all_cor_dist = all_cor_dist, sig_pairs = sig_pairs,
     p_fisher = p_fisher, raw_databases_scores = raw_databases_scores)
}

#' @importFrom data.table as.data.table
perform_correlations <- function(
      strategy = "normalized_counts_RNA_vs_miRNA_normalized_counts", 
      RNAseq, miRNAseq, std_positions = NULL, cor_cutoff = 0, 
      pval_cutoff = 0.05){ #correct_positions is a mirna_RNA pairs vector
    # the parsed strategy name text is used to subset RNAseq/miRNAseq obects
    if (strategy == "DEGs_DEMs_permutated"){
        all_pairs <- perform_deg_dem_strat(RNAseq$DH_results,
                                          miRNAseq$DH_results)
    } else {
        strat_description <- unlist(strsplit(strategy, "_RNA_vs_miRNA_"))
        RNA_profiles <- as.matrix(RNAseq[[strat_description[1]]])
        miRNA_profiles <- as.matrix(miRNAseq[[strat_description[2]]])
        if (ncol(RNA_profiles) == 0 || ncol(miRNA_profiles) == 0) {
          return(data.frame(NULL))
        }
        all_pairs <- correlate_profiles(RNA_profiles, miRNA_profiles)
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
    all_pairs$correlated_pairs <- all_pairs$correlation <= cor_cutoff & 
                                     all_pairs$pval < pval_cutoff
    if (!is.null(std_positions)){
     # Add extra column to indicate number of pair in std_positions
     all_pairs_str <- all_pairs[,c("RNAseq", "miRNAseq")]
     all_pairs_str <- tidyr::unite(all_pairs_str, "pairs", RNAseq:miRNAseq, 
                                   sep = "_")
     all_pairs_str <- all_pairs_str$pairs  
     all_pairs$pair_n <- match(all_pairs_str, std_positions) 
     all_pairs[,c("RNAseq", "miRNAseq")] <- NULL
    }
    return(all_pairs)
}

perform_deg_dem_strat <- function(RNAseq, miRNAseq){
    DEGS <- RNAseq[RNAseq$genes_tag %in%
     c("PREVALENT_DEG", "POSSIBLE_DEG"), "gene_name"]
    DEMS <- miRNAseq[miRNAseq$genes_tag %in% 
     c("PREVALENT_DEG", "POSSIBLE_DEG"), "gene_name"]
    all_pairs <- expand.grid(DEGS, DEMS, 
                               stringsAsFactors= FALSE, 
                               KEEP.OUT.ATTRS = FALSE)
    colnames(all_pairs) <- c("RNAseq", "miRNAseq")
    all_pairs$correlation <- -1
    all_pairs$pval <- 0
    return(all_pairs)
}

#' @importFrom WGCNA corPvalueStudent
#' @importFrom stats cor
correlate_profiles <- function(RNA_profiles, miRNA_profiles) {

    correlations <- WGCNA::cor(RNA_profiles, miRNA_profiles)
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
#' @importFrom dplyr select
#' @importFrom rlang .data
expand_module <- function(all_pairs, tag, DH_results){
    "%>%" <- magrittr::"%>%"

    mod_tag <- paste0(tag, "_mod")
    names(all_pairs)[names(all_pairs)== tag] <- mod_tag

    partial_expanded <- DH_results %>% dplyr::select(.data$Cluster_ID,
    .data$gene_name)
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

    input_to_entrezgene <- id_translation_orgdb(input_ids = ensembl_ids,
     organism_db = organism_info$Bioconductor_DB[1],
     org_var_name = organism_info$Bioconductor_VarName[1])

    # Fix names
    colnames(input_to_entrezgene) <- c("ensembl_gene_id", "entrezgene") 

    input_to_symbol <- id_translation_orgdb(
        input_ids = input_to_entrezgene$entrezgene,
        organism_db = organism_info$Bioconductor_DB[1],
    org_var_name = organism_info$Bioconductor_VarName_SYMBOL[1])
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



#' @importFrom plyr rbind.fill
#' @importFrom data.table setnames
get_stats_by_group <- function(
all_pairs, 
strategy, 
permutations = 50, 
db_groups = c("predicted", "validated", "pred_and_val")
){
 strat_cts <- data.frame() #contingency table
 known_miRNAs_c <- sum(all_pairs[,get(strategy) & known_miRNA])
 novel_miRNAs_c <- sum(all_pairs[,get(strategy) & !known_miRNA])
 background <- all_pairs[,possible_positives]
 strat_summary <- data.frame(
              pairs = c(known_miRNAs_c, novel_miRNAs_c),
              type = c("known_miRNAs", "novel_miRNAs"))
 for (db_group in db_groups){
    db_gs <- all_pairs[,get(db_group)]
    strat_group_ct <- calc_ct(
                        experiment = all_pairs[,get(strategy) & known_miRNA], 
                        gold_standard = db_gs)
    strat_group_ct$db_group <- db_group
   
    strat_cts <- rbind(strat_cts, strat_group_ct)
    perm_stats <- permutations_stats(permutations = permutations, 
        db_gs = db_gs, background = background, 
        sample_size = known_miRNAs_c)
    
    data.table::setnames(perm_stats, "mean", "pairs")
    perm_stats$type <- paste0(db_group, "_random")
    group_summary <- data.frame(type = db_group,
                                pairs = sum(all_pairs[,get(strategy) & 
                                                    get(db_group)]))
    strat_summary <- plyr::rbind.fill(strat_summary, group_summary, perm_stats)
 }
 strat_summary$strategy <- strategy
 strat_cts$strategy <- strategy
 return(list(strat_summary = strat_summary, strat_cts= strat_cts))
}


#' @importFrom miRBaseConverter getAllMiRNAs
prepare_all_pairs <- function(
RNAseq, 
miRNAseq, 
multimir, 
organism){

 all_pairs <- expand.grid(
                colnames(RNAseq$normalized_counts), 
                colnames(miRNAseq$normalized_counts), 
                stringsAsFactors = TRUE)
 names(all_pairs) <- c("RNAseq", "miRNAseq")
 all_pairs <- add_multimir_info(all_pairs = all_pairs, 
                                multiMiR = multimir)

 
 all_known_miRNAs <- miRBaseConverter::getAllMiRNAs(version = "v22", 
                                                    type = "all", 
                                                    species = organism)
 all_known_miRNAs <- all_known_miRNAs$Accession
 all_pairs$known_miRNA <- all_pairs$miRNAseq %in% all_known_miRNAs
 
 #possible_positives: Pairs with RNA and miRNA included in multiMiR. 
 all_pairs$possible_positives <- all_pairs$known_miRNA & 
                                 all_pairs$miRNAseq %in% 
                                  multimir$mature_mirna_acc &
                                 all_pairs$RNAseq %in% 
                                   multimir$target_ensembl 
 return(all_pairs)
}

#' @importFrom miRBaseConverter miRNA_AccessionToName
translate_all_id <- function(
miRNA_IDs, 
RNA_IDs, 
organism_info, 
translate_ensembl){
 mirna_names <- unique(miRNA_IDs)
 mirna_names <- mirna_names[grepl("MIMAT", mirna_names)]
 mirna_names <- miRBaseConverter::miRNA_AccessionToName(mirna_names, 
     targetVersion = "v22")
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
   # background_weights <- raw_databases_scores[[database]]
   # background_weights <- background_weights[!is.na(background_weights)]
   c_table <- calc_ct(experiment = all_pairs[,get(strategy)],
                       gold_standard = !is.na(all_pairs[,get(database)]))
                       # weights = all_pairs[,get(database)],
                       # background_weights = background_weights)
   c_table$database <- database
   cont_tables <- rbind(cont_tables, c_table)
 }
 return(cont_tables)
}
