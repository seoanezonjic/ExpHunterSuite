get_org_db <- function(current_organism_info) {
  org_db <- current_organism_info$Bioconductor_DB[1]
  org_db <- eval(parse(text = paste0(org_db,"::",org_db)))
  return(org_db)
}


  check_id_valid_orgdb <- function(gene_id, id_type="input", organism_info, outcome_action="stop") {
          print("we")
    org_db <- get_org_db(organism_info)
         print("wa")
    if(id_type == "input") possible_ids <- AnnotationDbi::keytypes(org_db)
    else possible_ids <- AnnotationDbi::columns(org_db)

    if(! gene_id %in% possible_ids) {
      if(outcome_action=="stop") {
        stop(paste(c("gene id must be one of the following:", possible_ids), collapse=" "))
      } else if(outcome_action=="warn") {
        warning(paste(c("gene id must be one of the following:", possible_ids), collapse=" "))
        return(FALSE)
      }
    }
    return(TRUE)
  }

get_translation_tables_bm <- function(input_gene_id, DEGH_results, current_organism_info) {
  if(input_gene_id == "ENTREZID") {
    input_to_entrezgene <- data.frame(input=row.names(DEGH_results), 
                                      ENTREZID=row.names(DEGH_results))
  } else {
    # Check input gene ID valid
    check_id_valid_orgdb(gene_id=input_gene_id, id_type="input", organism_info=current_organism_info, outcome_action="stop")
        input_to_entrezgene <- translate_ids_orgdb(input_genes=row.names(DEGH_results), 
        input_id=input_gene_id, organism_info = current_organism_info) 
  }
  symbol_output_available <- check_id_valid_orgdb(gene_id=input_gene_id, id_type="output", organism_info=current_organism_info, outcome_action="warning")
            
  if(symbol_output_available == TRUE) {
    input_to_symbol <- translate_ids_orgdb(input_genes=row.names(DEGH_results), 
    input_id=input_gene_id, output_id="SYMBOL", organism_info = current_organism_info)
  } else {
    input_to_symbol <- NULL
  }
  return(list(input_to_entrezgene = input_to_entrezgene, input_to_symbol = input_to_symbol))
}

    translate_ids_orgdb <- function(input_genes, input_id, output_id="ENTREZID", organism_info){
         print("we")
         org_db <- get_org_db(organism_info)
         print("wa")
#         input_genes <- c("ENSMUSG00000051951", "ENSMUSG00000103377", "ENSMUSG00000103025", 
# "ENSMUSG00000103201", "ENSMUSG00000055493", "ENSMUSG00000024164", 
# "ENSMUSG00000026822", "ENSMUSG00000097971", "ENSMUSG00000069516", 
# "ENSMUSG00000073418")
#         input_id <- "ENSEMBL"
#         output_id <- "ENTREZID"
        possible_ids <- AnnotationDbi::columns(org_db)
        if(! input_id %in% possible_ids) 
        stop(paste(c("gene keytype must be one of the following:", possible_ids), collapse=" "))
        ids <- tryCatch(
        ids <- AnnotationDbi::select(org_db, keys=input_genes, column=output_id, keytype=input_id),
        error=function(cond){
            ids <- NULL
        }
        )
    #    return(ids)
    return(ids[!is.na(ids[,2]),])
    }

download_latest_kegg_db <- function(organism=hsa, file) {
  prepare_KEGG <- get("prepare_KEGG", envir=asNamespace("clusterProfiler"), inherits = FALSE)
  ENRICH_DATA <- prepare_KEGG(organism, "KEGG", "kegg")
  saveRDS(ENRICH_DATA, file=file)
}

  
#' Perform topGO enrichment analysis of a list of genes
#' @export
#' @importClassesFrom topGO topGOdata
#' @importFrom topGO annFUN
multienricher_topGO <- function(all_funsys, genes_list, universe=NULL, organism_info, gene_id="entrez",  
  algorithm = "classic", statistic = "fisher", nodeSize = 5, task_size=1, 
  workers=1, ...){

 org_db <- organism_info$Bioconductor_DB[1]
 enrichments_topGO <- list()

 for(funsys in all_funsys) {
  if (funsys %in% c("CC","BP","MF")){
    if(is.null(universe)) {
      go_to_genes <- topGO::annFUN.org(funsys, mapping = org_db, ID = gene_id)
      universe <- unique(unlist(go_to_genes))
    }

    geneList <- factor(as.integer(universe %in% genes_list[[1]]))
    names(geneList) <- universe
     
    save(list = ls(all.names = TRUE), file = "~/environment_topgo.RData")

    # Necessary to create environment variables used in the initialization of the topGOdata object.
    topGO::groupGOTerms()
    GOdata <- new("topGOdata",
                 ontology = funsys,
                 allGenes = geneList,
                 nodeSize = nodeSize,
                 mapping = org_db,
                 annotationFun = topGO::annFUN.org,
                 ID = gene_id)
    }
    save(list = ls(all.names = TRUE), file = "~/environment_topgo.RData")

    enriched_cats <- parallel_list(genes_list, function(l_genes) {
      geneList <- factor(as.integer(universe %in% l_genes))
      names(geneList) <- universe
      l_GOdata <- topGO::updateGenes(object = GOdata, geneList = geneList)
      resultFis <- topGO::runTest(l_GOdata, algorithm = algorithm, statistic = statistic)
    }, workers = workers, task_size = task_size )

    enrichments_topGO[[funsys]] <- enriched_cats
  }
  return(enrichments_topGO)
}

#' Perform gsea enrichment analysis of a list of genes
#' @export
multienricher_gsea <- function(all_funsys=NULL, genes_list, organism_info, org_db = NULL, task_size=1, 
  workers=1, pvalueCutoff = 0.05, pAdjustMethod = "BH", kegg_file=NULL, 
  all_custom_sets=NULL, readable=TRUE, ...){

  common_params <- list(pvalueCutoff = pvalueCutoff, 
    pAdjustMethod = pAdjustMethod, ...)

  if(! is.null(all_custom_sets)) {
    if(is.null(names(all_custom_sets))) stop("Custom sets enrichment object must be a named list")
    all_funsys <- c(all_funsys, names(all_custom_sets))
  }

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
      enrf <- clusterProfiler::GSEA
      specific_params <- list(TERM2GENE = all_custom_sets[[funsys]])
    } else {
      stop("funsys", funsys, "not recognized")
    }

    save(list = ls(all.names = TRUE), file = "~/environment_gsea.RData")

    #l_genes <- genes_list[[1]]
    print("parallel list for")
    print(funsys)
    print(names(genes_list))
    enriched_cats <- parallel_list(genes_list, function(l_genes){
      print("ltest")
       params_genes <- c(specific_params, common_params, list(gene = l_genes))
       enriched_cats <- do.call("enrf", params_genes)
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
  all_custom_sets=NULL, readable=FALSE, ...){
cat("first")
  common_params <- list(universe = universe, pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, 
    pAdjustMethod = pAdjustMethod, ...)

  if(! is.null(all_custom_sets)) {
    if(is.null(names(all_custom_sets))) stop("Custom sets enrichment object must be a named list")
    all_funsys <- c(all_funsys, names(all_custom_sets))
  }

  enrichments_ORA <- list()
cat("out")
  for(funsys in all_funsys) {
cat(funsys)
    if (funsys %in% c("CC","BP","MF")){
      cat("here\n")
      org_db <- get_org_db(organism_info)
      enrf <- prepare_enrichment_GO(enrichment_type="ora", subont = funsys, org_db = org_db)
      specific_params <- list(OrgDb = org_db, ont = funsys, readable = readable)

    } else  if (funsys == "Reactome"){

      enrf <- prepare_enrichment_Reactome(enrichment_type="ora", reactome_id = organism_info$Reactome_ID[1])
      specific_params <- list(organism = organism_info$Reactome_ID[1], readable = readable)

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
