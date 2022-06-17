add_translated_gene_ids <- function(DEGH_results, input_ids, input_gene_id, gene_translation_tables) {
        input_to_entrezgene <- gene_translation_tables[["input_to_entrezgene"]]    
        input_to_symbol <- gene_translation_tables[["input_to_symbol"]]    

        if(!is.null(input_to_symbol)) {
          DEGH_results <- data.frame(SYMBOL = input_to_symbol[
            match(input_ids, input_to_symbol[[input_gene_id]]), "SYMBOL"], 
            DEGH_results)
         }
         DEGH_results <- data.frame(ENTREZID = input_to_entrezgene[
           match(input_ids, input_to_entrezgene[[input_gene_id]]), "ENTREZID"], 
           DEGH_results)
         DEGH_results <- data.frame(input_IDs=input_ids, DEGH_results)
}

get_sig_genes <- function(DEGH_results) {
    prev_genes <- DEGH_results[DEGH_results$genes_tag == "PREVALENT_DEG" &
                               !is.na(DEGH_results$ENTREZID), "ENTREZID"]
                               "%>%" <- magrittr::"%>%"
    ## TODO => ESTARIA BIEN REFLEJAR ESTA INFORMACION EN EL REPORT
    union_DEGs_df <- subset(DEGH_results, genes_tag %in% c("POSSIBLE_DEG",
                         "PREVALENT_DEG"))
    union_DEGs <- union_DEGs_df[!is.na(union_DEGs_df$input_IDs), 
                                 "input_IDs"] %>% unique
    return(list(prev_genes=prev_genes, union_DEGs_df=union_DEGs_df,
        union_DEGs=union_DEGs))
}

get_gene_lists <- function(DEGH_results, fc_colname) {
    geneList <- DEGH_results[!is.na(DEGH_results$ENTREZID),  fc_colname]
    names(geneList) <- DEGH_results[!is.na(DEGH_results$ENTREZID), "ENTREZID"]
    geneList <- sort(geneList, decreasing = TRUE)
    return(geneList)
}

get_sig_genes_cl <- function(DEGH_results) {
    cls <- unique(DEGH_results$Cluster_ID)
    # DELETE GREY MODULE
    if (any(c(0,"grey") %in% cls)) {
        cls <- cls[!cls %in% c(0,"grey")]
    } else {
        warning("Module Zero/Grey not found")
    }

    clgenes <- lapply(cls,function(cl) { # Find
        unique(DEGH_results$ENTREZID[which(DEGH_results$Cluster_ID == 
                                              cl)])
    }) 
    names(clgenes) <- cls
    return(clgenes)
}

get_gene_lists_cl <- function(DEGH_results, fc_colname) {
    # JRP - note we don't delete grey module for GSEA, any reason?
    DEGH_res_list <- split(DEGH_results, DEGH_results$Cluster_ID)
    lapply(DEGH_res_list, function(x) get_gene_lists(x, fc_colname))
}

get_org_db <- function(current_organism_info) {
  org_db <- current_organism_info$Bioconductor_DB[1]
  org_db <- eval(parse(text = paste0(org_db,"::",org_db)))
  return(org_db)
}

check_id_valid_orgdb <- function(gene_id, id_type="input", organism_info, outcome_action="stop") {
  org_db <- get_org_db(organism_info)
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

#' Translates a given gene ID using a dictionary. Note: one unknown ID can
#' corresponds to many known ids. 
#' @param ids_to_translate set of IDs to be translated
#' @param annot_table dictionary to translate IDs
#' @keywords translate
#' @return translated IDs or NA if it's not possible to translate
translate_from_table <- function(ids_to_translate, annot_table){ 
  translated_ids <- unlist(lapply(ids_to_translate, function(id){
    indx <- which(annot_table[,2] == id)
    if(length(indx) == 0){
      return(NA)
    } else{
      return(annot_table[indx[1],1])
    }
  }))
  return(translated_ids)
}

get_translation_tables_bm <- function(input_gene_id, input_ids, current_organism_info) {
  if(input_gene_id == "ENTREZID") {
    input_to_entrezgene <- data.frame(input=input_ids, 
                                      ENTREZID=input_ids)
  } else {
    # Check input gene ID valid
    check_id_valid_orgdb(gene_id=input_gene_id, id_type="input", organism_info=current_organism_info, outcome_action="stop")

        input_to_entrezgene <- translate_ids_orgdb(input_genes=input_ids, 
        input_id=input_gene_id, organism_info = current_organism_info) 
  }
  symbol_output_available <- check_id_valid_orgdb(gene_id=input_gene_id, id_type="output", 
                                                  organism_info=current_organism_info, outcome_action="warning")
            
  if(symbol_output_available == TRUE) {
    input_to_symbol <- translate_ids_orgdb(input_genes=input_ids, 
    input_id=input_gene_id, output_id="SYMBOL", organism_info = current_organism_info)
  } else {
    input_to_symbol <- NULL
  }
  return(list(input_to_entrezgene = input_to_entrezgene, input_to_symbol = input_to_symbol))
}

translate_ids_orgdb <- function(input_genes, input_id, output_id="ENTREZID", organism_info){
  org_db <- get_org_db(organism_info)

  possible_ids <- AnnotationDbi::columns(org_db)
  if(! input_id %in% possible_ids) 
    stop(paste(c("gene keytype must be one of the following:", possible_ids), collapse=" "))

    ids <- tryCatch(
      ids <- AnnotationDbi::select(org_db, keys=input_genes, column=output_id, keytype=input_id),
      error=function(cond){
            ids <- NULL
        }
    )
    return(ids[!is.na(ids[,2]),])
    }
  
#' Perform topGO enrichment analysis of a list of genes
#' @export
#' @importClassesFrom topGO topGOdata
#' @importFrom topGO annFUN
multienricher_topGO <- function(all_funsys, genes_list, universe=NULL, organism_info, gene_id="entrez",  
  algorithm = "classic", statistic = "fisher", nodeSize = 5, task_size=1, 
  workers=1, ...){

  unlisted_input_flag <- FALSE
  if(! is.list(genes_list)) {
    unlisted_input_flag <- TRUE
    genes_list <- list(genes_list)
  }

 org_db <- organism_info$Bioconductor_DB[1]
 enrichments_topGO <- vector("list", length(all_funsys))
 names(enrichments_topGO) <- all_funsys
 for(funsys in all_funsys) {
  if (funsys %in% c("CC","BP","MF")){
    if(is.null(universe)) {
      go_to_genes <- topGO::annFUN.org(funsys, mapping = org_db, ID = gene_id)
      universe <- unique(unlist(go_to_genes))
    }

    geneList <- factor(as.integer(universe %in% genes_list[[1]]))
    names(geneList) <- universe
     

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

    enriched_cats <- parallel_list(genes_list, function(l_genes) {
      geneList <- factor(as.integer(universe %in% l_genes))
      names(geneList) <- universe
      l_GOdata <- topGO::updateGenes(object = GOdata, geneList = geneList)
      resultFis <- topGO::runTest(l_GOdata, algorithm = algorithm, statistic = statistic)
    }, workers = workers, task_size = task_size )
    
    if(unlisted_input_flag) enriched_cats <- enriched_cats[[1]]
    enrichments_topGO[[funsys]] <- enriched_cats
  }
  return(enrichments_topGO)
}

#' Perform gsea enrichment analysis of a list of genes
#' @export
multienricher_gsea <- function(all_funsys=NULL, genes_list, organism_info, org_db = NULL, task_size=1, 
  workers=1, pvalueCutoff = 0.05, pAdjustMethod = "BH", kegg_file=NULL, 
  custom_sets=NULL, readable=FALSE, symbols_in_plots=TRUE, ...){

  unlisted_input_flag <- FALSE
  if(! is.list(genes_list)) {
    unlisted_input_flag <- TRUE
    genes_list <- list(genes_list)
  }


  common_params <- list(pvalueCutoff = pvalueCutoff, 
    pAdjustMethod = pAdjustMethod, ...)

  if(! is.null(custom_sets)) {
    if(is.null(names(custom_sets))) stop("Custom sets enrichment object must be a named list")
    all_funsys <- c(all_funsys, names(custom_sets))
  }

  enrichments_gsea <- vector("list", length(all_funsys))
  names(enrichments_gsea) <- all_funsys
  org_db <- get_org_db(organism_info)
  for(funsys in all_funsys) {
    if (funsys %in% c("CC","BP","MF")){

      enrf <- prepare_enrichment_GO(enrichment_type="gsea", subont = funsys, org_db = org_db)
      specific_params <- list(OrgDb = org_db, ont = funsys)
    } else  if (funsys == "Reactome"){
      enrf <- prepare_enrichment_Reactome(enrichment_type="gsea", 
                                          reactome_id = organism_info$Reactome_ID[1])
            specific_params <- list(organism = organism_info$Reactome_ID[1])

    } else  if (funsys == "KEGG"){
      enrf <- prepare_enrichment_KEGG(enrichment_type="gsea", kegg_file = kegg_file)
      specific_params <- list(organism = organism_info$KeggCode[1])

    } else if (funsys %in% names(custom_sets)) {
      enrf <- clusterProfiler::GSEA
      specific_params <- list(TERM2GENE = custom_sets[[funsys]])
    } else {
      stop("funsys", funsys, "not recognized")
    }

    library("ReactomePA")
    enriched_cats <- parallel_list(genes_list, function(l_genes){
       params_genes <- c(specific_params, common_params, list(gene = l_genes))
       enriched_cats <- do.call("enrf", params_genes)
      }, 
      workers= workers, task_size = task_size
    )
    # enriched_cats[sapply(enriched_cats, is.null)] <- NULL
    # enriched_cats <- lapply(enriched_cats, function(x) { 
    #    DOSE::setReadable(x, OrgDb = org_db, 
    #    keyType="ENTREZID")
    # })

    if(unlisted_input_flag) enriched_cats <- enriched_cats[[1]]
    enrichments_gsea[[funsys]] <- enriched_cats
  }
  save(enrichments_gsea, file=paste0("~/", funsys, "_enrichments_gsea.RData"))

  return(enrichments_gsea)
}

prepare_enrichment_GO <- function(enrichment_type, subont, org_db) {
  if(enrichment_type == "ora")   enrf <- clusterProfiler::enrichGO
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
  if(is.null(kegg_file) || ! file.exists(kegg_file) ) stop("kegg_file not found or not provided. 
  It can be downloaded using download_latest_kegg_db()")

  if(enrichment_type == "ora") enrf <- clusterProfiler::enrichKEGG
  if(enrichment_type == "gsea") enrf <- clusterProfiler::gseKEGG

  ENRICH_DATA <- readRDS(kegg_file)
  pattern_to_remove <- "KEGG_DATA *<-"
  ltorem <- grep(pattern_to_remove, body(enrf))
  body(enrf)[[ltorem]] <- substitute(KEGG_DATA <- ENRICH_DATA)
  return(enrf)
}


check_multienricher_input <- function(genes_list, custom_sets) {
  unlisted_input_flag <- FALSE
  if(! is.list(genes_list)) {
    unlisted_input_flag <- TRUE
    genes_list <- list(genes_list)
  }

  if(! is.null(custom_sets)) {
  if(is.null(names(custom_sets))) stop("Custom sets enrichment object must be a named list")
    all_funsys <- c(all_funsys, names(custom_sets))
  }
}

#' Perform ORA enrichment analysis of a list of genes
#' @export
multienricher_ora <- function(all_funsys=NULL, genes_list, universe=NULL, organism_info, org_db = NULL, task_size=1, 
  workers=1, pvalueCutoff = 0.05, qvalueCutoff = 0.2, pAdjustMethod = "BH", kegg_file=NULL, 
  custom_sets=NULL, readable=FALSE, symbols_in_plots=TRUE, ...){
print("multi time")
  unlisted_input_flag <- FALSE
  if(! is.list(genes_list)) {
    unlisted_input_flag <- TRUE
    genes_list <- list(genes_list)
  }

  common_params <- list(universe = universe, pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, 
    pAdjustMethod = pAdjustMethod, ...)

  if(! is.null(custom_sets)) {
    if(is.null(names(custom_sets))) stop("Custom sets enrichment object must be a named list")
    all_funsys <- c(all_funsys, names(custom_sets))
  }

  org_db <- get_org_db(organism_info)

  enrichments_ORA <- vector("list", length(all_funsys))
  names(enrichments_ORA) <- all_funsys
  for(funsys in all_funsys) {
    if (funsys %in% c("CC","BP","MF")){
      enrf <- prepare_enrichment_GO(enrichment_type="ora", subont = funsys, org_db = org_db)
      specific_params <- list(OrgDb = org_db, ont = funsys, readable = readable)

    } else  if (funsys == "Reactome"){
      enrf <- prepare_enrichment_Reactome(enrichment_type="ora", reactome_id = organism_info$Reactome_ID[1])
      specific_params <- list(organism = organism_info$Reactome_ID[1], readable = readable)

    } else if (funsys == "KEGG"){
      enrf <- prepare_enrichment_KEGG(enrichment_type="ora", kegg_file = kegg_file)
      specific_params <- list(organism = organism_info$KeggCode[1])
                       
    } else if (funsys %in% names(custom_sets)) {
      enrf <- clusterProfiler::enricher
      specific_params <- list(TERM2GENE = custom_sets[[funsys]])

    } else {
      stop("funsys", funsys, "not recognized")
    }
    # print()
    #save(genes_list, enrf, specific_params, common_params, file = "/mnt/scratch/users/bio_267_uma/josecordoba/claudia_prot/unique_prots/get_cdf_table.rb_0000/test.Rdata")
    enriched_cats <- parallel_list(genes_list, function(l_genes){
       params_genes <- c(specific_params, common_params, list(gene = l_genes))
       enr <- do.call("enrf", params_genes)
      }, 
      workers= workers, task_size = task_size
    )
    enriched_cats[sapply(enriched_cats, is.null)] <- NULL
    enriched_cats <- lapply(enriched_cats, function(x) { 
      DOSE::setReadable(x, OrgDb = org_db, 
        keyType="ENTREZID")
    })

    if(unlisted_input_flag) enriched_cats <- enriched_cats[[1]]
    enrichments_ORA[[funsys]] <- enriched_cats
  }

  return(enrichments_ORA)
}

merge_clusters <- function(results_list) {
  merged_clusters <- lapply(results_list, clusterProfiler::merge_result)
}

add_term_sim_ora <- function(deg_enr_ora) {
  enr_with_termsim <- list()
  for(funsys in names(deg_enr_ora)) {
    enrichment <- deg_enr_ora[[funsys]]
    if(is.list(enrichment)) {
      enr_with_termsim[[funsys]] <- sapply(enrichment, function(enr) {
        if(nrow(enr) > 0) {
#         return(catched_pairwise_termsim(enr))
          return(trycatch_pairwise_termsim(enr))
        } else {
          return(enr)
        }
      })
    } else {
      if(nrow(enrichment) > 0) {
        enr_with_termsim[[funsys]] <- trycatch_pairwise_termsim(enrichment)
        #enr_with_termsim[[funsys]] <- catched_pairwise_termsim(enrichment)
      } else {
        enr_with_termsim[[funsys]] <- enrichment
      }
    }
  }
  return(enr_with_termsim)
}

#' if pairwise_termsim throws an error, remove dup cat descriptions
#' Only seems to be a problem with Reactome - to study further
#' Further investigation re: number of cats also required
#' @param enr enrichment object to be studied
#' @param num_cats number of categories to be shown
#' @return enrichment object after add termsim info
#' @importFrom enrichplot pairwise_termsim
trycatch_pairwise_termsim <- function(enr, num_cats = 200){
  enr <- tryCatch(
  {
    enr <-enrichplot::pairwise_termsim(enr)
  },
    error = function(cond){
    message("ERROR ADDING TERM SIMILARITY TO ENRICHMENT OBJECT")
    message(cond)
    message("ATTEMPTING TO FIX BY REMOVING DUPLICATED ID DESCRIPTIONS")
    ccr <- enr@compareClusterResult
    unique_desc_id <- unique(ccr[c("ID","Description")])
    IDs_with_dupl_desc <- unique_desc_id$ID[duplicated(unique_desc_id$Description)]
    enr@compareClusterResult <- ccr[! ccr$ID %in% IDs_with_dupl_desc, ]

    enrichplot::pairwise_termsim(enr)
  })
  return(enr)
}

