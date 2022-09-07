are_parentals <- function(term_A, term_B, GO_ancestors){
  parental <- NA
   if (term_A %in% GO_ancestors[[term_B]]){
     parental <- term_A
  } else if (term_B %in% GO_ancestors[[term_A]]) {
    parental <- term_B
  }
  return(parental)
}

parse_results_for_report <- function(enrichments, simplify_results = FALSE){
  enrichments_for_reports <- list()
  for (funsys in names(enrichments)){
    for (cluster in names(enrichments[[funsys]])){
      if (is.null(enrichments_for_reports[[cluster]]))
      enrichments_for_reports[[cluster]] <- list()
      enr <- enrichments[[funsys]][[cluster]]
      if (funsys %in% c("CC","MF","BP")){
        if (simplify_results) 
          enr <- clusterProfiler::simplify(enr) 
      }   
      if (length(enr$Description) > 2 ) 
      enr <- catched_pairwise_termsim(enr, 200)
      enrichments_for_reports[[cluster]][[funsys]] <- enr 
    }
  }
  return(enrichments_for_reports)
}


parse_mappings <- function(gene_mapping, org_db, keytype){
  cl_genes_mapping <- read.table(gene_mapping, header=TRUE)
  gane_mapping_name <- colnames(cl_genes_mapping)[3]
  ids <- AnnotationDbi::select(org_db , keys=cl_genes_mapping[,2], 
    column="ENTREZID", keytype=keytype)
  ids <- unique(ids)
  cl_genes_mapping <- merge(cl_genes_mapping, ids, by.x = "geneid", 
    by.y = keytype)
  cl_genes_mapping <- cl_genes_mapping[!is.na(cl_genes_mapping$ENTREZID),]
  gene_mapping <- cl_genes_mapping[,2:4]
  gene_mapping <- split(gene_mapping, gene_mapping$cluster)
  gene_mapping <-lapply(gene_mapping, function(cluster_mapping){
    gene_map <- cluster_mapping[,gane_mapping_name]
    names(gene_map) <- cluster_mapping$ENTREZID
    return(gene_map)
  })
  return(list(gane_mapping_name = gane_mapping_name, 
    gene_mapping = gene_mapping))
}


summarize_categories <- function(all_enrichments, sim_thr = 0.7, 
  common_name = "significant"){
   enrichment_summary <- list()
       # save(list = ls(all.names = TRUE), file = "summarize_environment.RData")

     print("time clusterize_terms")
     print(system.time(
  clusterized_terms <- clusterize_terms(all_enrichments, threshold = sim_thr, 
    common_name = common_name)
))

  for (funsys in names(clusterized_terms)){
    enrichments_cl <- all_enrichments[[funsys]]
    combined_enrichments <- combine_terms_by_cluster(
      enrichments_cl@compareClusterResult, clusterized_terms[[funsys]])
    combined_enrichments_tmp <-  reshape2::acast(combined_enrichments, 
      term_cluster~cluster, value.var="p.adjust", fill = 1)
    combined_enrichments <- apply(combined_enrichments_tmp, 2, FUN=as.numeric)
    dimnames(combined_enrichments) <- dimnames(combined_enrichments_tmp)
    enrichment_summary[[funsys]] <- combined_enrichments
  }
  return(enrichment_summary)
}

#' @importFrom stats as.dist
#' @importFrom Matrix forceSymmetric
#' @importFrom fastcluster hclust
#' @importFrom GO.db GOBPPARENTS GOBPANCESTOR GOMFPARENTS GOMFANCESTOR GOCCPARENTS GOCCANCESTOR
clusterize_terms <- function(all_enrichments, threshold = 0.7, 
  common_name = "significant"){
    # save(list = ls(all.names = TRUE), file = "clusterize_environment.RData")
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
    print("time calc_all_paths")
     print(system.time(
    paths <- calc_all_paths(ids = unique(names(GO_parents)), 
      GO_parents = GO_parents)
    ))
    levels <- unlist(lapply(paths, function(id){min(lengths(id))}))

    enrichments_cl <- all_enrichments[[funsys]]
    term_sim <- enrichments_cl@termsim
    term_sim <- Matrix::forceSymmetric(term_sim, uplo="U")
    term_sim[is.na(term_sim)] <- 0
    term_dis <- stats::as.dist(1 - term_sim)
    term_clust <- fastcluster::hclust(term_dis, method = "average")
    threshold = 1 - threshold
    clust_term <- stats::cutree(term_clust, h = threshold)

    cluster_names <- lapply(clust_term, function(cluster) {
       terms_in_cl <- names(clust_term)[clust_term == cluster]
       cluster_enrichments <- enrichments_cl@compareClusterResult
       cluster_enrichments <- cluster_enrichments[
         cluster_enrichments$Description %in% terms_in_cl,]
       if (common_name == "ancestor") {
          main_name <- get_common_ancestor(unique(cluster_enrichments$ID), 
            GO_ancestor = GO_ancestor, levels = levels)
       } else if (common_name == "significant"){
          main_name <- cluster_enrichments[which.min(
            cluster_enrichments$p.adjust), "Description"]
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
    cluster_terms <- clusters_enr[clusters_enr$Cluster == cluster, 
                                   "Description"]
    term_clusters <- unique(term_clustering[names(term_clustering) 
      %in% cluster_terms])
    for (term_cluster in term_clusters) {
      clustered_terms <- names(term_clustering)[term_clustering == term_cluster]
      combined_terms <- rbind(combined_terms, 
        c(cluster, term_cluster, min(clusters_enr[clusters_enr$Description %in% 
                                    clustered_terms, "p.adjust"])))
    }
  }
  colnames(combined_terms) <- c("cluster", "term_cluster", "p.adjust")
  return(combined_terms)
}


get_common_ancestor <- function(terms, GO_ancestor, levels, 
  common_ancestor_method = "lowest"){
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
      #get_GOid_term(names(common_ancestor))
      common_ancestor <- names(common_ancestor) 

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
                                             relative_level = (
                                              relative_level + 1))
      all_branches <- rbind(all_branches, parental_branches)
    }
    return(all_branches)
  }
}


#' @importFrom GO.db GOTERM
get_GOid_term <- function(GOid, output = "term"){
 term <- AnnotationDbi::Term(GO.db::GOTERM[GOid])
 names(term) <- NULL
 return(term)
}


cluster_enr_to_matrix <- function(compareClusterResult){
   compareClusterResult <- compareClusterResult[,c("Cluster", "Description", 
    "p.adjust")]
    compareClusterResult_tmp <-  reshape2::acast(compareClusterResult, 
      Description~Cluster, value.var="p.adjust", fill = 1)
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

#' @importFrom GO.db GOBPANCESTOR GOMFANCESTOR GOCCANCESTOR
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

