#!/usr/bin/env Rscript

options(warn=1)
if( Sys.getenv('DEGHUNTER_MODE') == 'DEVELOPMENT' ){
  # Obtain this script directory
  full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile), 
                 error=function(e) # works when using R CMD
                normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', 
                  commandArgs())], '='))[2]))
  main_path_script <- dirname(full.fpath)
  root_path <- file.path(main_path_script, '..', '..')
  # Load custom libraries
  custom_libraries <- c('general_functions.R', 
    'functional_analysis_library.R','plotting_functions.R')
  for (lib in custom_libraries){
    source(file.path(root_path, 'R', lib))
  }
  template_folder <- file.path(root_path, 'inst', 'templates')
  organisms_table_file <- file.path(root_path, "inst", "external_data", 
        "organism_table.txt")
}else{
  require('ExpHunterSuite')
  root_path <- find.package('ExpHunterSuite')
  template_folder <- file.path(root_path, 'templates')
  organisms_table_file <- file.path(root_path, "external_data", 
        "organism_table.txt")
}

col <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
         "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
         "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
         "#8A7C64", "#599861")


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
       #   main_name <- get_common_ancestor(unique(cluster_enrichments$ID), GO_offspring)
        #  main_name <- get_common_ancestor_path(unique(cluster_enrichments$ID), GO_parents)
          main_name <- get_common_ancestor_final(unique(cluster_enrichments$ID), GO_ancestor = GO_ancestor, levels = levels)
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


get_common_ancestor_final <- function(terms, GO_ancestor, levels, common_ancestor_method = "lowest"){
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
      
      common_ancestor <- get_GOid_term(names(common_ancestor))

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



get_GOid_term <- function(GOid){
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

########################## OPTIONS
option_list <- list(
  optparse::make_option(c("-i", "--input_file"), type="character", default=NULL,
                        help="2 columns - cluster and comma separated gene ids"),
  optparse::make_option(c("-w", "--workers"), type="integer", default=1,
                        help="number of processes for parallel execution. Default=%default"),
  optparse::make_option(c("-p", "--pvalcutoff"), type="double", default=0.05,
                        help="Cutoff for P value and adjusted P value for enrichments. Default=%default"),
  optparse::make_option(c("-q", "--qvalcutoff"), type="double", default=0.02,
                        help="Cutoff for Q value for enrichments. Default=%default"),
  optparse::make_option(c("-t", "--task_size"), type="integer", default=1,
                        help="number of clusters per task. Default=%default"),
  optparse::make_option(c("-F", "--force"), type="logical", default=FALSE, 
                        action = "store_true", help="Ignore temporal files"),
  optparse::make_option(c("-f", "--funsys"), type="character", default="BP,MF,CC", 
                        help="Funsys to execute: MF => GO Molecular Function, BP => GO Biological Process, CC => GO Celular Coponent. Default=%default"),
  optparse::make_option(c("-c", "--clean_parentals"), type="logical", default=FALSE, 
                        action = "store_true", help="Clean parentals GO terms that appears on the same clusters than child."),
  optparse::make_option(c("-s", "--simplify"), type="logical", default=FALSE, 
                        action = "store_true", help="Apply simplify function from cluster profiler to enrichment."),
  optparse::make_option(c("-O", "--model_organism"), type="character", default="Human", 
                        help="Model organism. Human or Mouse"),
  optparse::make_option(c("-M", "--mode"), type="character", default="D", 
                        help="String indicating report modes to execute. R: Generate report for each cluster and all clusters combined. P: Plots (dotplot and emaplot) combining all cluster enrichments. S: Summary mode (sumamrized categories heatmap)"),
  optparse::make_option(c("-S", "--sim_thr"), type="double", default=0.7,
                        help="Similarity cutoff for grouping categories in Summary mode. Default=%default"),
  optparse::make_option(c("-C", "--common_name"), type="character", default="significant", 
                        help="Name of the term groups. 'significant' to use the most significant term of each group. 'ancestor' to use the common ancestor of the group"),
  optparse::make_option(c("-d", "--description_file"), type="character", default=NULL,
                        help="Markdown file describing of the enriched clusters."),
  optparse::make_option(c("-k", "--gene_keytype"), type="character", default="ENTREZID",
                        help="What identifier is being used for the genes in the clusters?. Default=%default"),
  optparse::make_option(c("-g", "--gene_mapping"), type="character", default=NULL,
                        help="3 columns tabular file- Cluster - InputGeneID - NumericGeneMapping. Header must be indicated as cluster - geneid - [numeric_mapping]"),
  optparse::make_option(c("-G", "--group_results"), type="logical", default=FALSE, 
                        action = "store_true", help="Functions are gropuped in most frequent words in emaplots."),
  optparse::make_option(c("-o", "--output_file"), type="character", default="results",
                        help="Define the output path.")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

##################################### INITIALIZE ##
library(clusterProfiler)
output_path <- paste0(opt$output_file, "_functional_enrichment")
dir.create(output_path)
output_path <- normalizePath(output_path)
temp_file <- file.path(output_path, "enr_tmp.RData")
all_funsys <- c("MF", "CC", "BP") 
n_category <- 30
organisms_table <- get_organism_table(organisms_table_file)
current_organism_info <- organisms_table[rownames(organisms_table) %in% opt$model_organism,]
gene_mapping <- NULL
org_db <- get_org_db(current_organism_info)

#################################### MAIN ##


if (!is.null(opt$gene_mapping)){
  parsed_mappings <- parse_mappings(opt$gene_mapping, org_db, opt$gene_keytype)
  gene_mapping <- parsed_mappings[["gene_mapping"]]
  gane_mapping_name <-  parsed_mappings[["gane_mapping_name"]]
}

if (!file.exists(temp_file) || opt$force) {

  cluster_genes <- read.table(opt$input_file, header=FALSE)
  cluster_genes_list <- strsplit(cluster_genes[,2], ",")
  if(opt$gene_keytype != "ENTREZID") {
    cluster_genes_list <- lapply(cluster_genes_list, function(x){
                                      convert_ids_to_entrez(ids=x, 
                                                            gene_keytype=opt$gene_keytype,
                                                            org_db = org_db)}) 
  }

  names(cluster_genes_list) <- cluster_genes[,1]

  enrichments_ORA <- multienricher_2(funsys =  all_funsys, 
                                  cluster_genes_list =  cluster_genes_list, 
                                  task_size = opt$task_size, 
                                  org_db= org_db,
                                  workers = opt$workers, 
                                  pvalcutoff =  opt$pvalcutoff, 
                                  qvalcutoff = opt$qvalcutoff)
  save(enrichments_ORA, file = temp_file)

} else {
  load(temp_file)
}

enrichments_ORA_merged <- parse_cluster_results(enrichments_ORA, simplify_results = opt$simplify, clean_parentals = opt$clean_parentals)

if (grepl("R", opt$mode)){
    enrichments_for_reports <- parse_results_for_report(enrichments_ORA)
    write_fun_enrichments(enrichments_ORA, output_path, all_funsys)
    write_func_cluster_report(enrichments_for_reports,output_path,gene_mapping)
}

if (grepl("P", opt$mode)) {
  for (funsys in names(enrichments_ORA_merged)){
    if (length(unique(enrichments_ORA_merged[[funsys]]@compareClusterResult$Description)) < 2 ) next

    if (opt$group_results == TRUE){
      pp <- enrichplot::emapplot(enrichments_ORA_merged[[funsys]], showCategory= n_category, pie="Count", layout = "nicely", 
                  shadowtext = FALSE, node_label = "group", group_category = TRUE, 
                  nCluster = min(floor(nrow(enrichments_ORA_merged[[funsys]])/7), 20), nWords = 6, repel = TRUE)
    }else{
      pp <- enrichplot::emapplot(enrichments_ORA_merged[[funsys]], showCategory= n_category, pie="Count", layout = "nicely", 
                  shadowtext = FALSE, repel = TRUE)
    }

    ggplot2::ggsave(filename = file.path(output_path,paste0("emaplot_",funsys,"_",opt$output_file,".png")), pp, width = 30, height = 30, dpi = 300, units = "cm", device='png')

    pp <- enrichplot::dotplot(enrichments_ORA_merged[[funsys]], showCategory= n_category)
    ggplot2::ggsave(filename = file.path(output_path,paste0("dotplot_",funsys,"_",opt$output_file,".png")), pp, width = 60, height = 40, dpi = 300, units = "cm", device='png')

  }
}




if (grepl("S", opt$mode)){
  sum_enrichments <- summarize_categories(enrichments_ORA_merged, sim_thr = opt$sim_thr,common_name = opt$common_name)
 
  for (funsys in names(sum_enrichments)) {
    sum_table <- sum_enrichments[[funsys]]
    sum_table <- (sum_table > opt$pvalcutoff) + 0
    save(sum_table, file = "test.RData")
    sum_table <- sum_table[rownames(sum_table) != "to_remove",]
    heatmaply::heatmaply(sum_table, grid_color = "gray50",
                                    grid_width = 0.00001,
                                    dendrogram = "column",
                                    scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
                                    low = "#EE8291", 
                                    high = "white", 
                                    midpoint = 0.5, 
                                    limits = c(0, 1)),
                                    file = file.path(output_path, paste0("sum_",funsys,'_heatmap.html')))
 
    all_enrichments <- cluster_enr_to_matrix(enrichments_ORA_merged[[funsys]]@compareClusterResult)
    all_enrichments <- (all_enrichments > opt$pvalcutoff) + 0 
    heatmaply::heatmaply(all_enrichments, grid_color = "gray50",
                                    grid_width = 0.00001,
                                    dendrogram = "column",
                                    scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
                                    low = "#EE8291", 
                                    high = "white", 
                                    midpoint = 0.5, 
                                    limits = c(0, 1)),
                                    file = file.path(output_path, paste0("full_",funsys,'_heatmap.html')))
  
  }

}

