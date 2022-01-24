convert_ids_to_entrez <- function(ids, gene_keytype, org_db){
  possible_ids <- AnnotationDbi::columns(org_db)
  if(! gene_keytype %in% possible_ids) 
    stop(paste(c("gene keytype must be one of the following:", possible_ids), collapse=" "))
  ids <- tryCatch(
    ids <- AnnotationDbi::mapIds(org_db, keys=ids, column="ENTREZID", keytype=gene_keytype),
    error=function(cond){
      ids <- NULL
    }
  )
  return(ids[!is.na(ids)])
}


get_organism_id <- function(organism_info, funsys){
  if(funsys %in% c("BP","CC","MF")){
    org_id <- organism_info$Bioconductor_DB[1]
  }else if(funsys == "KEGG"){
    org_id <- organism_info$KeggCode[1]
  }else if(funsys == "REACT"){
    org_id <- organism_info$Reactome_ID[1]
  }
  return(org_id)
}



get_org_db <- function(current_organism_info){
  org_db <- current_organism_info$Bioconductor_DB[1]
  org_db <- eval(parse(text = paste0(org_db,"::",org_db)))
  return(org_db)
}


multienricher_2 <- function(funsys, cluster_genes_list, organism_info,org_db = NULL,task_size, workers, pvalcutoff = 0.05, qvalcutoff = 0.02, readable = TRUE){
  enrichments_ORA <- list()
  if (is.null(org_db)){
    org_db <- get_org_db(current_organism_info)
  }

  for(funsys in all_funsys) {
    if (funsys %in% c("CC","BP","MF")){
      enrf <- clusterProfiler::enrichGO
      get_enr_data <- get("get_GO_data", envir = asNamespace("clusterProfiler"),inherits = FALSE)
      pattern_to_remove <- "GO_DATA *<-"
      ENRICH_DATA <- get_enr_data(org_db, funsys, "ENTREZID")
      ltorem <- grep(pattern_to_remove, body(enrf))
      body(enrf)[[ltorem]] <- substitute(GO_DATA <- ENRICH_DATA)

    } else  if (funsys == "REACT"){
    
    }else if (funsys == "KEGG"){

    }
    params <- list(OrgDb = org_db,
                   pAdjustMethod = "BH",minGSSize = 5, ont = funsys,
                   pvalueCutoff  = pvalcutoff, qvalueCutoff = qvalcutoff, readable = readable) #Now Param is configured Only for GO
    enriched_cats <- parallel_list(cluster_genes_list, function(cl_genes){

       params_cl <- c(params, list(gene = cl_genes) )
       enr <- do.call("enrf", params_cl)
      }, 
      workers= workers, task_size = task_size
    )
    enrichments_ORA[[funsys]] <- enriched_cats
  }
  return(enrichments_ORA)
}

parse_cluster_results <- function(enrichments_ORA, simplify_results, clean_parentals){
  enrichments_ORA_tr <- list()
  for (funsys in names(enrichments_ORA)){
    enr_obj <- clusterProfiler::merge_result(enrichments_ORA[[funsys]])
    if(nrow(enr_obj@compareClusterResult) > 0){
      if (funsys %in% c("MF", "CC", "BP") && clean_parentals){
        enr_obj@fun <- "enrichGO"
        enr_obj <- clean_all_parentals(enr_obj, subont = funsys) 
      } 
      if (funsys %in% c("MF", "CC", "BP") && simplify_results){
        enr_obj@fun <- "enrichGO"
        enr_obj <- clusterProfiler::simplify(enr_obj) 
      } 
      enr_obj <- catched_pairwise_termsim(enr_obj, 200)
    }                              
    enrichments_ORA_tr[[funsys]] <- enr_obj 
  }
  return(enrichments_ORA_tr)
}


clean_all_parentals <- function(enr_obj, subont){
   ##ADD control for enrichresults or comparecluster
  if (subont=="BP"){
    GO_ancestors <- GO.db::GOBPANCESTOR
  } else if (subont=="MF"){
    GO_ancestors <-  GO.db::GOMFANCESTOR
  } else if (subont=="CC"){
    GO_ancestors <- GO.db::GOCCANCESTOR
  }
  GO_ancestors <- as.list(GO_ancestors)
  GO_ancestors <- GO_ancestors[!is.na(GO_ancestors)]

  enrich_obj <- enr_obj@compareClusterResult
  pre_hamming_m <- matrix(0, nrow = length(unique(enrich_obj$Cluster)), 
                              ncol = length(unique(enrich_obj$ID)), 
                            dimnames= list(unique(enrich_obj$Cluster), 
                                          unique(enrich_obj$ID)))
  for(pair in seq(nrow(enrich_obj))){
    pair_text <- enrich_obj[pair, c("Cluster","ID")]
    pair_text$Cluster <- as.character(pair_text$Cluster)
    pre_hamming_m[pair_text[1,1], pair_text[1,2]] <- 1
  }
  hamming_matrix <- hamming_binary(pre_hamming_m)
  hamming_0 <- as.data.frame(as.table(hamming_matrix)) 
  hamming_0 <- hamming_0[hamming_0[,3] == 0,c(1,2)]
  hamming_0 <- hamming_0[hamming_0[,1] != hamming_0[,2], ]
  hamming_0[] <- lapply(hamming_0, "as.character")
  terms_to_discard <- c()
  for(pair in seq(nrow(hamming_0))){
    possible_parentals <- unlist(hamming_0[pair, c(1,2)])
    parental <- are_parentals(possible_parentals[1], possible_parentals[2], GO_ancestors)
    if (!is.na(parental)){
      terms_to_discard  <- c(terms_to_discard,parental)
    }
  }
  enr_obj@compareClusterResult <- enrich_obj[!enrich_obj$ID %in% terms_to_discard,]
  return(enr_obj)
}

are_parentals <- function(term_A, term_B, GO_ancestors){
  parental <- NA
   if (term_A %in% GO_ancestors[[term_B]]){
     parental <- term_A
  } else if (term_B %in% GO_ancestors[[term_A]]) {
    parental <- term_B
  }
  return(parental)
}

hamming_binary <- function(X, Y = NULL) {
    if (is.null(Y)) {
        D <- t(1 - X) %*% X
        D + t(D)
    } else {
        t(1 - X) %*% Y + t(X) %*% (1 - Y)
    }
} #taken from https://johanndejong.wordpress.com/2015/10/02/faster-hamming-distance-in-r-2/

write_fun_enrichments <- function(enrichments_ORA, output_path, all_funsys){
  for(funsys in all_funsys) {
    enriched_cats <- enrichments_ORA[[funsys]]
    enriched_cats_dfs <- lapply(enriched_cats, data.frame)
    enriched_cats_bound <- data.table::rbindlist(enriched_cats_dfs, use.names= TRUE, idcol= "Cluster_ID" )
    if (nrow(enriched_cats_bound) == 0) next 
    utils::write.table(enriched_cats_bound, 
                       file=file.path(output_path, paste0("enrichment_",funsys,".csv")),
                       quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")
  }
}

parse_results_for_report <- function(enrichments, simplify_results = FALSE){
  enrichments_for_reports <- list()
  for (funsys in names(enrichments)){
    for (cluster in names(enrichments[[funsys]])){
      if (is.null(enrichments_for_reports[[cluster]]))
      enrichments_for_reports[[cluster]] <- list()
      enr <- enrichments[[funsys]][[cluster]]
      if (simplify_results)
        enr <- clusterProfiler::simplify(enr) 
      if (length(enr$Description) < 3 ) next
      enr <- catched_pairwise_termsim(enr, 200)
      enrichments_for_reports[[cluster]][[funsys]] <- enr 
    }
  }
  return(enrichments_for_reports)
}


parse_mappings <- function(gene_mapping, org_db, keytype){
  cl_genes_mapping <- read.table(gene_mapping, header=TRUE)
  gane_mapping_name <- colnames(cl_genes_mapping)[3]
  ids <- AnnotationDbi::select(org_db , keys=cl_genes_mapping[,2], column="ENTREZID", keytype=keytype)
  ids <- unique(ids)
  cl_genes_mapping <- merge(cl_genes_mapping, ids, by.x = "geneid", by.y = keytype)
  cl_genes_mapping <- cl_genes_mapping[!is.na(cl_genes_mapping$ENTREZID),]
  gene_mapping <- cl_genes_mapping[,2:4]
  gene_mapping <- split(gene_mapping, gene_mapping$cluster)
  gene_mapping <-lapply(gene_mapping, function(cluster_mapping){
    gene_map <- cluster_mapping[,gane_mapping_name]
    names(gene_map) <- cluster_mapping$ENTREZID
    return(gene_map)
  })
  return(list(gane_mapping_name = gane_mapping_name, gene_mapping = gene_mapping))
}

write_func_cluster_report <- function(enrichments_for_reports, output_path, gene_mapping){
  temp_path <- file.path(output_path, "temp")
  dir.create(temp_path) #perform parallel
  for (cluster in names(enrichments_for_reports)){
      geneList <- NULL
      if (!is.null(gene_mapping)){ 
        geneList <- gene_mapping[[cluster]]
      }
      outfile <- file.path(output_path, paste0(cluster, "_func_report.html"))
      temp_path_cl <- file.path(temp_path, paste0(cluster,"_temp"))
      rmarkdown::render(file.path(template_folder, 
                  'clusters_to_enrichment.Rmd'),output_file = outfile, 
              intermediates_dir = temp_path_cl)
  }
  unlink(temp_path, recursive = TRUE)
}


summarize_categories <- function(all_enrichments, sim_thr = 0.7, common_name = "significant"){
   enrichment_summary <- list()

  clusterized_terms <- clusterize_terms(all_enrichments, threshold = sim_thr, common_name = common_name)
  for (funsys in names(clusterized_terms)){
      enrichments_cl <- all_enrichments[[funsys]]
    combined_enrichments <- combine_terms_by_cluster(enrichments_cl@compareClusterResult, clusterized_terms[[funsys]])
    combined_enrichments_tmp <-  reshape2::acast(combined_enrichments, term_cluster~cluster, value.var="p.adjust", fill = 1)
    combined_enrichments <- apply(combined_enrichments_tmp, 2, FUN=as.numeric)
    dimnames(combined_enrichments) <- dimnames(combined_enrichments_tmp)
    enrichment_summary[[funsys]] <- combined_enrichments
  }
  return(enrichment_summary)
}

clusterize_terms <- function(all_enrichments, threshold = 0.7, common_name = "significant"){
  clusterized_terms <- list()
  for (funsys in names(all_enrichments)){
    if (funsys=="BP"){
      GO_parents <- GO.db::GOBPPARENTS
      GO_ancestor <- GO.db::GOBPANCESTOR
    } else if (funsys=="MF"){
      GO_parents <-  GO.db::GOMFPARENTS
      GO_ancestor <-  GO.db::GOMFANCESTOR
    } else if (funsys=="CC"){
      GO_parents <- GO.db::GOCCPARENTS
      GO_ancestor <- GO.db::GOCCANCESTOR
    }

    GO_ancestor <- as.list(GO_ancestor)
    GO_parents <- as.list(GO_parents)
    for (ancestor in names(GO_ancestor)){
      GO_ancestor[[ancestor]] <- c(GO_ancestor[[ancestor]], ancestor)
    }
    paths <- calc_all_paths(ids = unique(names(GO_parents)), GO_parents = GO_parents)
    levels <- unlist(lapply(paths, function(id){min(lengths(id))}))

    enrichments_cl <- all_enrichments[[funsys]]
    term_sim <- enrichments_cl@termsim
    term_sim <- Matrix::forceSymmetric(term_sim, uplo="U")
    term_dis <- as.dist(1 - term_sim)
    term_clust <- fastcluster::hclust(term_dis, method = "average")
    threshold = 1 - threshold
    clust_term <- cutree(term_clust, h = threshold)
    cluster_names <- lapply(clust_term, function(cluster) {
       terms_in_cl <- names(clust_term)[clust_term == cluster]
       cluster_enrichments <- enrichments_cl@compareClusterResult
       cluster_enrichments <- cluster_enrichments[cluster_enrichments$Description %in% terms_in_cl,]
       if (common_name == "ancestor") {
          main_name <- get_common_ancestor(unique(cluster_enrichments$ID), GO_ancestor = GO_ancestor, levels = levels)
       } else if (common_name == "significant"){
          main_name <- cluster_enrichments[which.min(cluster_enrichments$p.adjust), "Description"]
       } 
       return(main_name)
    })
    clust_term <- unlist(cluster_names)

    clusterized_terms[[funsys]] <- clust_term
  }
  return(clusterized_terms)
}




combine_terms_by_cluster <- function(clusters_enr, term_clustering) {
  combined_terms <- data.frame()
  for (cluster in unique(clusters_enr$Cluster)) {
    clusters_enr$Description <- as.character(clusters_enr$Description)
    cluster_terms <- clusters_enr[clusters_enr$Cluster == cluster, "Description"]
    term_clusters <- unique(term_clustering[names(term_clustering) %in% cluster_terms])
    for (term_cluster in term_clusters) {
      clustered_terms <- names(term_clustering)[term_clustering == term_cluster]
      combined_terms <- rbind(combined_terms, 
        c(cluster, term_cluster, min(clusters_enr[clusters_enr$Description %in% clustered_terms, "p.adjust"])))
    }
  }
  colnames(combined_terms) <- c("cluster", "term_cluster", "p.adjust")
  return(combined_terms)
}


get_common_ancestor <- function(terms, GO_ancestor, levels, common_ancestor_method = "lowest"){
  terms_ancestors <- GO_ancestor[terms]
  common_ancestors <- Reduce(intersect, terms_ancestors)
  common_ancestors <- levels[names(levels) %in% common_ancestors]
  common_ancestor <- NULL
  if (common_ancestor_method == "lowest") {
      common_ancestors <- sort(common_ancestors,  decreasing = TRUE)
  } else if (common_ancestor_method == "highest") {
      common_ancestors <- sort(common_ancestors,  decreasing = FALSE)
  }

  common_ancestor <- common_ancestors[1]
  if (common_ancestor > 2) {
      
      common_ancestor <- names(common_ancestor) #get_GOid_term(names(common_ancestor))

  } else {
      common_ancestor <- "to_remove"
  }
  return(common_ancestor)
}

get_all_parentals <- function(term, GO_parents, relative_level = 0){
  all_parents <- GO_parents[[term]]
  isa_parents <- all_parents[names(all_parents) == "isa"]
  term_level <- data.frame(term = term, 
                               relative_level = relative_level)
  if (any(isa_parents %in% "all")) {
    return(term_level)
  } else {
    all_branches <- term_level
    for (parent in isa_parents) {
      parental_branches <- get_all_parentals(term = parent, 
                                             GO_parents = GO_parents, 
                                             relative_level = (relative_level + 1))
      all_branches <- rbind(all_branches, parental_branches)
    }
    return(all_branches)
  }
}



get_GOid_term <- function(GOid, output = "term"){
 term <- AnnotationDbi::Term(GO.db::GOTERM[GOid])
 names(term) <- NULL
 return(term)
}


cluster_enr_to_matrix <- function(compareClusterResult){
   compareClusterResult <- compareClusterResult[,c("Cluster", "Description", "p.adjust")]
    compareClusterResult_tmp <-  reshape2::acast(compareClusterResult, Description~Cluster, value.var="p.adjust", fill = 1)
    compareClusterResult <- apply(compareClusterResult_tmp, 2, FUN=as.numeric)
    dimnames(compareClusterResult) <- dimnames(compareClusterResult_tmp)
    return(compareClusterResult_tmp)
}

calc_all_paths <- function(ids, GO_parents){
    env <- environment()
    paths <- list()
    calc_path(ids = ids, GO_parents = GO_parents, env = env)
    return(paths)
}

calc_path <- function(ids, GO_parents, env){
  for (id in ids) {
    if (is.null(env$paths[[id]])) {
      env$paths[[id]] <- list()
      direct_ancestors <- GO_parents[[id]]
      direct_ancestors <- direct_ancestors[names(direct_ancestors) == "isa"]
      if (length(direct_ancestors) == 0) {
        env$paths[[id]] <- c(env$paths[[id]], list(c(id)))
      } else {
        for (anc in direct_ancestors) {
          calc_path(ids = anc, GO_parents = GO_parents,env = env)
          env$paths[[id]] <- c(env$paths[[id]], 
                               lapply(env$paths[[anc]], function(path){
                                c(id,unlist(path))}))
        }
      }
    }
  }
}

vectdist <- function(vectA, vectB){
  # VectA and B must have same length. Exception not handled
  return(sqrt(sum((vectA - vectB)^2)))
}

toDistances <- function(vectors_matrix, rows = TRUE){
  if(!rows){
    vectors_matrix = t(vectors_matrix)
  }
  # Calc similitudes of rows
  numItems = nrow(vectors_matrix)
  Mdist = matrix(Inf,nrow = numItems, ncol = numItems)
  invisible(lapply(seq(numItems), function(i){
    if(i != numItems){
      invisible(lapply(seq(i+1, numItems), function(j){
        v = vectdist(vectors_matrix[i,],vectors_matrix[j,])
        Mdist[i,j] <<- v
        Mdist[j,i] <<- v
      }))
    }
    Mdist[i,i] <<- 0
  }))
  return(Mdist)
}


clean_parentals_in_matrix <- function(enrichment_mx, subont){
  if (subont=="BP"){
    GO_ancestors <- GO.db::GOBPANCESTOR
  } else if (subont=="MF"){
    GO_ancestors <-  GO.db::GOMFANCESTOR
  } else if (subont=="CC"){
    GO_ancestors <- GO.db::GOCCANCESTOR
  }
  for (cluster in colnames(enrichment_mx)) {

    enriched_cats <- rownames(enrichment_mx)
    for (category_ref in enriched_cats){
       parental_terms <- enriched_cats %in% GO_ancestors[[category_ref]]

       if (any(parental_terms)) {

          enrichment_mx[parental_terms, cluster] <- 1
       }

    }
  }

return(enrichment_mx)
}
