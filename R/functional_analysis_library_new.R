get_org_db <- function(current_organism_info){
  org_db <- current_organism_info$Bioconductor_DB[1]
  org_db <- eval(parse(text = paste0(org_db,"::",org_db)))
  return(org_db)
}

download_latest_kegg_db <- function(organism=hsa, file) {
  prepare_KEGG <- get("prepare_KEGG", envir=asNamespace("clusterProfiler"), inherits = FALSE)
  ENRICH_DATA <- prepare_KEGG(organism, "KEGG", "kegg")
  saveRDS(ENRICH_DATA, file=file)
}

#' Perform topGO enrichment analysis of a list of genes
#' @export
multienricher_topGO <- function(all_funsys, genes_list, universe=NULL, organism_info, gene_id="entrez", ...){

 org_db <- current_organism_info$Bioconductor_DB[1]

 for(funsys in all_funsys) {
  if (funsys %in% c("CC","BP","MF")){
    if(is.null(universe)) {
      go_to_genes <- annFUN.org(funsys, mapping = org_db, ID = gene_id)
      universe <- unique(unlist(go_to_genes))
    }

      de_genes <- sample(universe, 500)
      geneList <- factor(as.integer(universe %in% de_genes))
      names(geneList) <- universe
     
      save(list = ls(all.names = TRUE), file = "~/environment_topgo.RData")

     GOdata <- new("topGOdata",
                   ontology = funsys,
                   allGenes = geneList,
                   nodeSize = 5,
                   annot = go_to_genes, 
                   mapping = "org.Hs.eg.db",
                   ID = gene_id)
      save(list = ls(all.names = TRUE), file = "~/environment_topgo.RData")
  }
}
  #GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

}

#' Perform gsea enrichment analysis of a list of genes
#' @export
multienricher_gsea <- function(all_funsys, genes_list, organism_info, org_db = NULL, task_size=1, 
  workers=1, pvalueCutoff = 0.05, pAdjustMethod = "BH", kegg_file=NULL, 
  all_custom_sets=NULL, readable=TRUE, ...){

  common_params <- list(pvalueCutoff = pvalueCutoff, 
    pAdjustMethod = pAdjustMethod, ...)

  # if (is.null(org_db)){
  #   org_db <- get_org_db(organism_info)
  # }

  enrichments_gsea <- list()
  for(funsys in all_funsys) {
  
    if (funsys %in% c("CC","BP","MF")){
      org_db <- get_org_db(organism_info)
      enrf <- prepare_enrichment_GO(enrichment_type="gsea", subont = funsys, org_db = org_db)
      specific_params <- list(OrgDb = org_db, ont = funsys)
    } else  if (funsys == "Reactome"){
      enrf <- prepare_enrichment_Reactome(enrichment_type="gsea", 
                                          reactome_id = organism_info$Reactome_ID[1])
            specific_params <- list(organism = organism_info$Reactome_ID[1])

    } else  if (funsys == "KEGG"){
      enrf <- prepare_enrichment_KEGG(enrichment_type="gsea", kegg_file = kegg_file)
      specific_params <- list(organism = organism_info$KeggCode[1])

    } else if (funsys %in% names(all_custom_sets)) {

    } else {
      stop("funsys", funsys, "not recognized")
    }

    save(list = ls(all.names = TRUE), file = "~/environment_gsea.RData")

    enriched_cats <- parallel_list(genes_list, function(l_genes){
       params_genes <- c(specific_params, common_params, list(gene = l_genes))
       enr <- do.call("enrf", params_genes)
      }, 
      workers= workers, task_size = task_size
    )

    enrichments_gsea[[funsys]] <- enriched_cats
  }
  return(enrichments_gsea)
}

prepare_enrichment_GO <- function(enrichment_type, subont, org_db) {
 if(enrichment_type == "ora") enrf <- clusterProfiler::enrichGO
 if(enrichment_type == "gsea") enrf <- clusterProfiler::gseGO

 get_enr_data <- get("get_GO_data", envir = asNamespace("clusterProfiler"), inherits = FALSE)      
 pattern_to_remove  <- "GO_DATA *<-"
 ENRICH_DATA <- get_enr_data(org_db, subont, "ENTREZID")
 ltorem <- grep(pattern_to_remove, body(enrf))
 body(enrf)[[ltorem]] <- substitute(GO_DATA <- ENRICH_DATA)
 return(enrf)
}

prepare_enrichment_Reactome <- function(enrichment_type, reactome_id) {
  if(enrichment_type == "ora") enrf <- ReactomePA::enrichPathway
  if(enrichment_type == "gsea") enrf <- ReactomePA::gsePathway

  get_enr_data <- get("get_Reactome_DATA", envir = asNamespace("ReactomePA"), inherits = FALSE)
  ENRICH_DATA <- get_enr_data(reactome_id) 
  pattern_to_remove <- "Reactome_DATA *<-"
  ltorem <- grep(pattern_to_remove, body(enrf))
  body(enrf)[[ltorem]] <- substitute(Reactome_DATA <- ENRICH_DATA)
  return(enrf)
}

prepare_enrichment_KEGG <- function(enrichment_type, kegg_file) {
  if(! file.exists(kegg_file)) stop("kegg_file not found. Please download using download_latest_kegg_db or if using script ensure not using remote mode")

  if(enrichment_type == "ora") enrf <- clusterProfiler::enrichKEGG
  if(enrichment_type == "gsea") enrf <- clusterProfiler::gseKEGG

  ENRICH_DATA <- readRDS(kegg_file)
  pattern_to_remove <- "KEGG_DATA *<-"
  ltorem <- grep(pattern_to_remove, body(enrf))
  body(enrf)[[ltorem]] <- substitute(KEGG_DATA <- ENRICH_DATA)
  return(enrf)
}

#' Perform ORA enrichment analysis of a list of genes
#' @export
multienricher_ora <- function(all_funsys=NULL, genes_list, universe=NULL, organism_info, org_db = NULL, task_size=1, 
  workers=1, pvalueCutoff = 0.05, qvalueCutoff = 0.2, pAdjustMethod = "BH", kegg_file=NULL, 
  all_custom_sets=NULL, readable=TRUE, ...){

  common_params <- list(universe = universe, pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, 
    pAdjustMethod = pAdjustMethod, ...)

  if(! is.null(all_custom_sets)) {
    if(is.null(names(all_custom_sets))) stop("Custom sets enrichment object must be a named list")
    all_funsys <- c(all_funsys, names(all_custom_sets))
  }

  enrichments_ORA <- list()
  for(funsys in all_funsys) {

    if (funsys %in% c("CC","BP","MF")){
      org_db <- get_org_db(organism_info)
      enrf <- prepare_enrichment_GO(enrichment_type="ora", subont = funsys, org_db = org_db)
      specific_params <- list(OrgDb = org_db, ont = funsys, readable = readable)

    } else  if (funsys == "Reactome"){

      enrf <- prepare_enrichment_Reactome(enrichment_type="ora", reactome_id = organism_info$Reactome_ID[1])
      specific_params <- list(organism = organism_info$Reactome_ID[1])

    } else if (funsys == "KEGG"){

      enrf <- prepare_enrichment_KEGG(enrichment_type="ora", kegg_file = kegg_file)
      specific_params <- list(organism = organism_info$KeggCode[1])
                       
    } else if (funsys %in% names(all_custom_sets)) {

      enrf <- clusterProfiler::enricher
      specific_params <- list(TERM2GENE = all_custom_sets[[funsys]])

    } else {
      stop("funsys", funsys, "not recognized")
    }
    save(list = ls(all.names = TRUE), file = "~/environment_ora.RData")

    enriched_cats <- parallel_list(genes_list, function(l_genes){
       params_genes <- c(specific_params, common_params, list(gene = l_genes))
       enr <- do.call("enrf", params_genes)
      }, 
      workers= workers, task_size = task_size
    )
    enrichments_ORA[[funsys]] <- enriched_cats
  }
  return(enrichments_ORA)
}
