###############################################################################
########################## FUNCTIONAL ANALYSIS LIBRARY ########################
###############################################################################


#' Scale a matrix using minimum-maximum method
#' @param data_matrix to be scaled
#' @param norm_by_col boolean flag: if true scaling will be performed 
#' by columns intead of by rows. Default: FALSE
#' @keywords method
#' @return scaled matrix
scale_data_matrix <- function(data_matrix, norm_by_col = FALSE) {
  if ( norm_by_col == TRUE) {
    data_matrix <- t(data_matrix)
  } 
  dm_min_max <- matrixStats::rowRanges(data_matrix, na.rm = TRUE)
  dm_diffs <- as.vector(matrixStats::rowDiffs(dm_min_max))
  dm_mins <- dm_min_max[,1]
  #when all values are the same in the row, all values are turned to 0.
  dm_diffs[dm_diffs == 0] <- 1 
  #main scaling function
  scaled_counts <- (data_matrix - dm_mins) / dm_diffs
  if (norm_by_col == TRUE) {
    scaled_counts <- t(scaled_counts)
  } 
  return(scaled_counts)
}




#' Load a GMT format file and return a dataframe in correct format
#' @param gmt_file file to be loaded
#' @return GMT loaded info
#' @keywords file
#' @export
#' @examples
#' gmt_file <- system.file("extData", 
#' "toy_categories_1.gmt", package = "ExpHunterSuite")
#' load_and_parse_gmt(gmt_file)
load_and_parse_gmt <- function(gmt_file) {
    # Load file
    gmt <- readLines(con = gmt_file)
    gmt_list <- strsplit(gmt, "\t")
    parsed_gmt <- do.call(rbind, lapply(gmt_list, function(category) {
          category_name <- category[1]
          genes <- category[3:length(category)]
          parsedTerms <- data.frame(Term = category_name, 
                                    Gene= genes, 
                                    stringsAsFactors = FALSE)
          return(parsedTerms)
    }))
    return(parsed_gmt)
}

#' Table with information abaut all organism available
#' @param file to be loaded. If none given, internal organism table loaded
#' @return organism table
#' @keywords method
#' @export
#' @importFrom utils read.table
#' @examples
#' ot <- get_organism_table()
get_organism_table <- function(file = NULL){
  if(is.null(file)) {
    file <- system.file("external_data", "organism_table.txt", 
                          package="ExpHunterSuite")
  }
  return(utils::read.table(file, 
                           header = TRUE, 
                           row.names=1, 
                           sep="\t", 
                           stringsAsFactors = FALSE, 
                           fill = NA))
}

add_translated_gene_ids <- function(DEGH_results, 
                                    input_ids, 
                                    input_gene_id, 
                                    gene_translation_tables) {
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
    genes_tag <- NULL
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
    DEGH_res_list <- split(DEGH_results, DEGH_results$Cluster_ID)
    lapply(DEGH_res_list, function(x) get_gene_lists(x, fc_colname))
}

#' Obtain an org_db object from a current_organism_info, obtained from
#' subsetting the results of get_organism_table()
#' @param current_organism_info row from output of get_organism_table() 
#' @export
#' @keywords org_db
#' @return org_db 
get_org_db <- function(current_organism_info) {
  org_db <- current_organism_info$Bioconductor_DB[1]
  org_db <- eval(parse(text = paste0(org_db,"::",org_db)))
  return(org_db)
}

check_id_valid_orgdb <- function(gene_id, 
                                 id_type="input", 
                                 organism_info, 
                                 outcome_action="stop") {
  org_db <- get_org_db(organism_info)
  if(id_type == "input") possible_ids <- AnnotationDbi::keytypes(org_db)
  else possible_ids <- AnnotationDbi::columns(org_db)

  if(! gene_id %in% possible_ids) {
    if(outcome_action=="stop") {
      stop(paste(c("gene id must be one of the following:", possible_ids), 
                                                              collapse=" "))
    } else if(outcome_action=="warn") {
      warning(paste(c("gene id must be one of the following:", possible_ids), 
                                                                collapse=" "))
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

get_translation_tables_orgdb <- function(
  input_gene_id, 
  input_ids, 
  current_organism_info) {
  org_db <- get_org_db(current_organism_info)
  if(input_gene_id == "ENTREZID") {
    input_to_entrezgene <- data.frame(input=input_ids, 
                                      ENTREZID=input_ids)
  } else {
    # Check input gene ID valid
    check_id_valid_orgdb(gene_id=input_gene_id, 
                         id_type="input", 
                         organism_info=current_organism_info, 
                         outcome_action="stop")

        input_to_entrezgene <- translate_ids_orgdb(ids=input_ids, 
        input_id=input_gene_id, org_db=org_db) 
  }
  symbol_output_available <- check_id_valid_orgdb(
                                        gene_id=input_gene_id, 
                                        id_type="output", 
                                        organism_info=current_organism_info, 
                                        outcome_action="warning")
            
  if(symbol_output_available == TRUE) {
    if(input_gene_id == "SYMBOL") {
      input_to_symbol <- data.frame(input=input_ids, 
                                      SYMBOL=input_ids)
    } else {
      input_to_symbol <- translate_ids_orgdb(ids=input_ids, 
      input_id=input_gene_id, output_id="SYMBOL", org_db=org_db)
    }
  } else {
    input_to_symbol <- NULL
  }
  return(list(input_to_entrezgene = input_to_entrezgene, 
              input_to_symbol = input_to_symbol))
}

translate_ids_orgdb <- function(
ids, 
input_id, 
output_id="ENTREZID", 
org_db=org_db, 
just_output_ids=FALSE){
  
  possible_ids <- AnnotationDbi::columns(org_db)
  if(! input_id %in% possible_ids) 
    stop(paste(c("gene keytype must be one of the following:", possible_ids),
                                                               collapse=" "))
    ids <- tryCatch(
      ids <- AnnotationDbi::select(org_db, keys=ids, column=output_id, 
                                  keytype=input_id),
      error=function(cond){
            ids <- NULL
        }
    )
    ids <- ids[!is.na(ids[,output_id]),]
    if(just_output_ids == TRUE) {
      return(unique(ids[,output_id]))
    } else {
      return(ids)
    }
}

translate_gmt <- function(gmt, gene_keytype, org_db){
  splitted_gmt <- split(gmt$Gene, gmt$Term)
  tr_splitted_gmt  <- lapply(splitted_gmt, function(x) {
                      tr_table <- translate_ids_orgdb(ids=x, 
                                            input_id=gene_keytype,
                                            org_db = org_db)
                      return(unique(tr_table[,2])) })
  translated_gmt <- lapply(tr_splitted_gmt, as.data.frame)
  translated_gmt <- as.data.frame(data.table::rbindlist(translated_gmt , 
    use.names = TRUE, idcol = TRUE))
  names(translated_gmt) <- c("Term","Gene")
  return(translated_gmt)
}

#' Perform topGO enrichment analysis of a list of genes
#' @param all_funsys vector of funsys to use (e.g. MF, Reatcome)
#' @param genes_list vector of genes or list of vectors, to be enriched
#' @param universe background of genes to use as universe
#' @param organism_info from the annotation table: infor on db names etc
#' @param gene_id what identifier do the genes use - must be entrez for 
#' @param algorithm alrgorithm for decting enriched functions - see topGO docs
#' @param statistic method for detecting enriched functions
#' @param nodeSize related to the creation of the GOdata object
#' @param workers for parallelization
#' @param task_size for parallelization
#' @export
#' @importClassesFrom topGO topGOdata
#' @importFrom topGO annFUN
multienricher_topGO <- function(all_funsys, genes_list, universe=NULL, 
  organism_info, gene_id="entrez", p_value_cutoff = 0.05, algorithm = "classic", statistic = "fisher", 
  nodeSize = 5, task_size=1, workers=1, scoreOrder = "increasing", clean_parentals = FALSE) {


  #checking input format
  if (!is.list(genes_list)) 
    genes_list <- list(genes_list)
  

  if (!any(all_funsys %in% c("CC","BP","MF"))) {
    message("None GO subontology has been selected")
    return(NULL)
  }

  org_db <- organism_info$Bioconductor_DB[1]

  all_funsys <- all_funsys[all_funsys %in% c("CC","BP","MF")]
  enrichments_topGO <- as.list(all_funsys)
  names(enrichments_topGO) <- all_funsys
  enrichments_topGO <- lapply(enrichments_topGO, multi_topGOTest, genes_list = genes_list,
                            nodeSize = nodeSize, org_db = org_db, p_value_cutoff = p_value_cutoff, universe = universe, 
                            gene_id = gene_id, algorithm = algorithm, statistic = statistic,
                            scoreOrder = scoreOrder, clean_parentals = clean_parentals,workers = workers, task_size = task_size)
  return(enrichments_topGO)
}

multi_topGOTest <- function(funsys, genes_list, nodeSize = 5, org_db, 
                            universe = NULL, gene_id = "entrez",p_value_cutoff = 0.05, algorithm = "classic", statistic = "fisher",
                            scoreOrder = "increasing", workers = 1, task_size = 1, geneSel = NULL, clean_parentals = FALSE) {
  if(is.null(universe) && statistic == "fisher") {
    go_to_genes <- topGO::annFUN.org(funsys, mapping = org_db, ID = gene_id)
    universe <- unique(unlist(go_to_genes))
  }

  if (is.null(geneSel)) 
    geneSel <- get_default_geneSel(statistic)

  genes_list <- lapply(genes_list, function(l_genes){
    if(is.numeric(l_genes)) { 
      if(!statistic %in%  c("ks", "t"))
        stop("ERROR: numeric scores can be only computed using 'ks' and 't' statistic")
      
      return(l_genes)      
    } else {
      if (statistic %in% c("ks", "t"))
        stop("ERROR: 'ks' and 't' can be only computed using numeric scores") 
      l_genes <- factor(as.integer(universe %in% l_genes))
      names(l_genes) <- universe
      return(l_genes)
    }   
  })
  all_names <- unique(unlist(lapply(genes_list, names)))
  if (statistic == "fisher") {
    rand <- factor(c(1,2))
  } else {
    rand <- seq(0,1,0.1)
  }

  all_genes <- sample(rand, length(all_names), replace = TRUE)
  names(all_genes) <-  all_names 
  # Create environment variables used for initialization of the topGOdata obj
  topGO::groupGOTerms()
  enriched_cats <- parallel_list(genes_list, function(l_genes) {
    if(length(l_genes) == 0)
      return(data.frame())
    
    GOdata <- new("topGOdata",
                     ontology = funsys,
                     allGenes = l_genes,
                     nodeSize = nodeSize,
                     mapping = org_db,
                     geneSel = geneSel,
                     annotationFun = topGO::annFUN.org,
                     ID = gene_id)

    if (length(GOdata@graph@nodes) == 0)
      return(data.frame())

      result <- topGO::runTest(GOdata, algorithm = algorithm, 
      statistic = statistic, scoreOrder = scoreOrder)
    res_table <- topGO::GenTable(GOdata, 
                                 p.value = result,
                                 topNodes = length(result@score), 
                                 format.FUN = function(x, dig, eps){x})
    if(clean_parentals)
      res_table <- clean_topGO_parentals(res_table, funsys, p_value_cutoff)
    if("p.value" %in% colnames(res_table)) {
      res_table$fdr <- stats::p.adjust(res_table$p.value, method = "BH", n = length(res_table$p.value))  
    } else {
      res_table$fdr <- NA_real_
    }
    res_table
  
  }, workers = workers, task_size = task_size )

  if (length(enriched_cats) == 1)
    enriched_cats <- enriched_cats[[1]]
  return(enriched_cats)

}

get_default_geneSel <-function(statistic){
  if (statistic == "fisher"){
    return(NULL) 
  } else {
    return(function(x){x})
  } 
}





#' Perform gsea enrichment analysis of a list of genes
#' @export
#' @param all_funsys vector of funsys to use (e.g. MF, Reatcome)
#' @param genes_list vector of genes or list of vectors, to be enriched
#' @param organism_info from the annotation table: infor on db names etc
#' @param org_db org db file (optional) - normally obtained from organism_info
#' @param workers for parallelization
#' @param task_size for parallelization
#' @param pvalueCutoff p-value finding enriched categories
#' @param pAdjustMethod adjustment method for p-values
#' @param kegg_file must be provided if kegg annotation required
#' @param custom_sets custom set object (processed by load_and_parse_gmt)
#' @param ... other arguments passed to the enrichment function
multienricher_gsea <- function(all_funsys=NULL, genes_list, organism_info, 
  org_db = NULL, task_size=1, workers=1, pvalueCutoff = 0.05, 
  pAdjustMethod = "BH", kegg_file=NULL, custom_sets=NULL, ...){

  unlisted_input_flag <- FALSE
  if(! is.list(genes_list)) {
    unlisted_input_flag <- TRUE
    genes_list <- list(genes_list)
  }


  common_params <- list(pvalueCutoff = pvalueCutoff, 
    pAdjustMethod = pAdjustMethod, ...)

  if(! is.null(custom_sets)) {
    if(is.null(names(custom_sets))) 
      stop("Custom sets enrichment object must be a named list")
    all_funsys <- c(all_funsys, names(custom_sets))
  }

  enrichments_gsea <- vector("list", length(all_funsys))
  names(enrichments_gsea) <- all_funsys
  org_db <- get_org_db(organism_info)
  for(funsys in all_funsys) {
    if (funsys %in% c("CC","BP","MF")){

      enrf <- prepare_enrichment_GO(enrichment_type="gsea", subont = funsys, 
        org_db = org_db)
      specific_params <- list(OrgDb = org_db, ont = funsys)
    } else  if (funsys == "Reactome"){
      enrf <- prepare_enrichment_Reactome(enrichment_type="gsea", 
                                          reactome_id = 
                                          organism_info$Reactome_ID[1])
            specific_params <- list(organism = organism_info$Reactome_ID[1])

    } else  if (funsys == "KEGG"){
      enrf <- prepare_enrichment_KEGG(enrichment_type="gsea", 
        kegg_file = kegg_file)
      specific_params <- list(organism = organism_info$KeggCode[1])

    } else if (funsys %in% names(custom_sets)) {
      enrf <- clusterProfiler::GSEA
      specific_params <- list(TERM2GENE = custom_sets[[funsys]])
    } else {
      stop("funsys", funsys, "not recognized")
    }

    #enriched_cats <- parallel_list(genes_list, function(l_genes){
    enriched_cats <- lapply(genes_list, function(l_genes){
       params_genes <- c(specific_params, common_params, list(gene = l_genes))
       enriched_cats <- do.call("enrf", params_genes)
    })
    #},workers= 1, task_size = 1 )

    # enriched_cats[sapply(enriched_cats, is.null)] <- data.frame()
    # enriched_cats <- lapply(enriched_cats, function(x) { 
    #    DOSE::setReadable(x, OrgDb = org_db, 
    #    keyType="ENTREZID")
    # })

    if(unlisted_input_flag) enriched_cats <- enriched_cats[[1]]
    enrichments_gsea[[funsys]] <- enriched_cats
  }

  return(enrichments_gsea)
}

#' @importFrom clusterProfiler enrichGO gseGO
prepare_enrichment_GO <- function(enrichment_type, subont, org_db) {
  if(enrichment_type == "ora")   enrf <- clusterProfiler::enrichGO
  if(enrichment_type == "gsea") enrf <- clusterProfiler::gseGO

  get_enr_data <- get("get_GO_data", envir = asNamespace("clusterProfiler"), 
    inherits = FALSE)      
  pattern_to_remove  <- "GO_DATA *<-"
  ENRICH_DATA <- get_enr_data(org_db, subont, "ENTREZID")
  ltorem <- grep(pattern_to_remove, body(enrf))
  body(enrf)[[ltorem]] <- substitute(GO_DATA <- ENRICH_DATA)
  return(enrf)
}

#' @importFrom ReactomePA enrichPathway gsePathway
prepare_enrichment_Reactome <- function(enrichment_type, reactome_id) {
  if(enrichment_type == "ora") enrf <- ReactomePA::enrichPathway
  if(enrichment_type == "gsea") enrf <- ReactomePA::gsePathway

  get_enr_data <- get("get_Reactome_DATA", envir = asNamespace("ReactomePA"), 
    inherits = FALSE)
  ENRICH_DATA <- get_enr_data(reactome_id) 
  pattern_to_remove <- "Reactome_DATA *<-"
  ltorem <- grep(pattern_to_remove, body(enrf))
  body(enrf)[[ltorem]] <- substitute(Reactome_DATA <- ENRICH_DATA)
  return(enrf)
}

#' @importFrom clusterProfiler gseKEGG enrichKEGG
prepare_enrichment_KEGG <- function(enrichment_type, kegg_file) {
############# DEPRECATED FUNCTION

  if(is.null(kegg_file) || ! file.exists(kegg_file) ) 
                stop("kegg_file not found or not provided. 
  It can be downloaded using download_latest_kegg_db()")

  if(enrichment_type == "ora") enrf <- clusterProfiler::enrichKEGG

  if(enrichment_type == "gsea") enrf <- clusterProfiler::gseKEGG

  ENRICH_DATA <- readRDS(kegg_file)
  pattern_to_remove <- "KEGG_DATA *<-"
  ltorem <- grep(pattern_to_remove, body(enrf))
  body(enrf)[[ltorem]] <- substitute(KEGG_DATA <- ENRICH_DATA)
  return(enrf)
}




enrichKEGG_user_data <- function (
  gene, 
  organism = "hsa", 
  keyType = "kegg", 
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH", 
  universe, minGSSize = 10, 
  maxGSSize = 500,
  qvalueCutoff = 0.2, 
  user_data, ...)
{
  enr_int <- utils::getFromNamespace("enricher_internal", "clusterProfiler")
    res <- enr_int(gene, 
      pvalueCutoff = pvalueCutoff,
      pAdjustMethod = pAdjustMethod, 
      universe = universe, 
      minGSSize = minGSSize,
      maxGSSize = maxGSSize, 
      qvalueCutoff = qvalueCutoff, 
      USER_DATA = user_data,
      ...)
    if (is.null(res))
        return(res)
    res@ontology <- "KEGG"
    res@organism <- organism
    res@keytype <- keyType
    return(res)
}


#' Perform ORA enrichment analysis of a list of genes
#' @param all_funsys vector of funsys to use (e.g. MF, Reatcome)
#' @param genes_list vector of genes or list of vectors, to be enriched
#' @param universe background of genes to use as universe
#' @param organism_info from the annotation table: infor on db names etc
#' @param org_db org db file (optional) - normally obtained from organism_info
#' @param workers for parallelization
#' @param task_size for parallelization
#' @param pvalueCutoff p-value finding enriched categories
#' @param qvalueCutoff q-value finding enriched categories
#' @param pAdjustMethod adjustment method for p-values
#' @param kegg_file must be provided if kegg annotation required
#' @param custom_sets custom set object (processed by load_and_parse_gmt)
#' @param readable Whether output should include gene symbols
#' @param return_all Whether to remove list items with no enrichment
#' @param ... other arguments passed to the enrichment function
#' @importFrom DOSE setReadable
#' @export
multienricher_ora <- function(all_funsys=NULL, genes_list, universe=NULL, 
  organism_info, org_db = NULL, task_size=1, workers=1, pvalueCutoff = 0.05, 
  qvalueCutoff = 0.2, pAdjustMethod = "BH", kegg_file=NULL, 
  custom_sets=NULL, readable=TRUE, return_all=FALSE, ...){
  unlisted_input_flag <- FALSE
  if(! is.list(genes_list)) {
    unlisted_input_flag <- TRUE
    genes_list <- list(genes_list)
  }
  if( row.names(organism_info) != "Human" && "DO" %in% all_funsys) {
    warning("Cannot run DO with non-human organim, will be skipped")
    all_funsys <- all_funsys[! all_funsys %in% "DO"]
  }
  common_params <- list(universe = universe, pvalueCutoff = pvalueCutoff, 
    qvalueCutoff = qvalueCutoff, 
    pAdjustMethod = pAdjustMethod, ...)

  if(! is.null(custom_sets)) {
    if(is.null(names(custom_sets))) 
                stop("Custom sets enrichment object must be a named list")
    all_funsys <- c(all_funsys, names(custom_sets))
  }
  if (is.null(org_db)) {
    org_db <- get_org_db(organism_info)
  }
  enrichments_ORA <- vector("list", length(all_funsys))
  names(enrichments_ORA) <- all_funsys
  for(funsys in all_funsys) {
    specific_params <- NULL
    if (funsys %in% c("CC","BP","MF","ALL")){
      enrf <- prepare_enrichment_GO(enrichment_type="ora", subont = funsys, 
        org_db = org_db)
      specific_params <- list(OrgDb = org_db, ont = funsys)

    } else  if (funsys == "Reactome"){
     
      enrf <- prepare_enrichment_Reactome(enrichment_type="ora", 
        reactome_id = organism_info$Reactome_ID[1])
      specific_params <- list(organism = organism_info$Reactome_ID[1])

    } else if (funsys == "KEGG"){
      specific_params <- list(organism = organism_info$KeggCode[1])
      if (!is.null(kegg_file)){ 
        enrf <- enrichKEGG_user_data 
        ENRICH_DATA <- readRDS(kegg_file)
        specific_params<- c(specific_params, list(user_data = ENRICH_DATA))
      } else {
          stop("kegg_file not found or not provided. 
          It can be downloaded using download_latest_kegg_db()")
      }      
    } else if (funsys == "DO"){
      # No need for prepare enrichment function - hack not necessary
      enrf <- DOSE::enrichDO
    } else if (funsys == "DGN"){
      enrf <- DOSE::enrichDGN
    } else if (funsys == "DGNv"){
      enrf <- DOSE::enrichDGNv
    } else if (funsys == "NCG"){
      enrf <- DOSE::enrichNCG 
    } else if (funsys %in% names(custom_sets)) {
      enrf <- clusterProfiler::enricher
      specific_params <- list(TERM2GENE = custom_sets[[funsys]])
    } else {
      stop("funsys", funsys, "not recognized")
    }
    temp_workers <- workers
    if(funsys == "DO") workers <- 1 #DOSE not working in parallel, but its quick

    enriched_cats <- parallel_list(genes_list, function(l_genes){
        params_genes <- c(specific_params, common_params, list(gene = l_genes))
        enr <- do.call("enrf", params_genes)
     } , workers= workers, task_size = task_size
    )
    workers <- temp_workers
    if(return_all == FALSE) enriched_cats[sapply(enriched_cats,is.null)] <- NULL

    if(readable == TRUE) {
      enriched_cats <- lapply(enriched_cats, function(x) { 
        if(! is.null(x)) return(
          DOSE::setReadable(x, OrgDb = org_db, 
          keyType="ENTREZID")
          )
        else return(data.frame())
      })
    }
    # Remove species name from reactome enrichment
    if(funsys == "Reactome") { 
      enriched_cats <- lapply(enriched_cats, function(enr) {
        if(nrow(enr) > 0)
        enr@result[,"Description"] <- gsub("^.+\\r: ", "", 
                                                    enr@result[,"Description"])
        return(enr)
      })
    }
    if(unlisted_input_flag) enriched_cats <- enriched_cats[[1]]    
    enrichments_ORA[[funsys]] <- enriched_cats
  }
  return(enrichments_ORA)
}

add_term_sim_ora <- function(deg_enr_ora) {
  enr_with_termsim <- list()
  for(funsys in names(deg_enr_ora)) {
    enrichment <- deg_enr_ora[[funsys]]
    if(is.list(enrichment) && length(enrichment) > 1) {
      enr_with_termsim[[funsys]] <- sapply(enrichment, function(enr) {
        if(nrow(enr) > 0) {
          return(trycatch_pairwise_termsim(enr))
        } else {
          return(enr)
        }
      })
    } else {
      if(nrow(enrichment) > 0) {
        enr_with_termsim[[funsys]] <- trycatch_pairwise_termsim(enrichment)
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
#' @return enrichment object after add termsim info
#' @importFrom enrichplot pairwise_termsim
trycatch_pairwise_termsim <- function(enr){ 
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
    IDs_with_dupl_desc <- unique_desc_id$ID[
                                         duplicated(unique_desc_id$Description)]
    enr@compareClusterResult <- ccr[! ccr$ID %in% IDs_with_dupl_desc, ]

    enrichplot::pairwise_termsim(enr)
  })
  return(enr)
}


merge_clusters <- function(results_list) {
  merged_clusters <- lapply(results_list, clusterProfiler::merge_result)
}

filter_cluster_enrichment <- function(compareCluster, filter_list){
    compareCluster@compareClusterResult <- compareCluster@compareClusterResult[
              compareCluster@compareClusterResult$Cluster %in% filter_list,]
    return(compareCluster)
}

hamming_binary <- function(X, Y = NULL) {
    if (is.null(Y)) {
        D <- t(1 - X) %*% X
        D + t(D)
    } else {
        t(1 - X) %*% Y + t(X) %*% (1 - Y)
    }
} # from https://johanndejong.wordpress.com/
  #2015/10/02/faster-hamming-distance-in-r-2/

summarize_merged_ora <- function(ORA_merged, sim_thr=0.7, 
  common_name = "ancestor", pthreshold=0.1) {
  if(! any(c("BP","MF","CC") %in% names(ORA_merged))) {
    stop("Summarize only works with GO categories. Skipping")
    return(NULL)
  }
  ORA_merged <- ORA_merged[names(ORA_merged) %in% c("BP","MF","CC")]

  sum_enrichments <- summarize_categories(ORA_merged, sim_thr = sim_thr, 
    common_name = common_name)

  summ_ora_to_plot <- list()
  for (funsys in names(ORA_merged)) {
    sum_table <- sum_enrichments[[funsys]]
    sum_table <- (sum_table > pthreshold) + 0
    sum_table <- sum_table[rownames(sum_table) != "to_remove",]
    summ_enr_table <- sum_table
    rownames(summ_enr_table) <- get_GOid_term(rownames(summ_enr_table))

    summ_enr_clean <-  clean_parentals_in_matrix(sum_table, funsys)
    rownames(summ_enr_clean) <- get_GOid_term(rownames(summ_enr_clean))
    summ_enr_clean <- summ_enr_clean[rowSums(summ_enr_clean < pthreshold) != 0,]
   
    full_enr_table <- cluster_enr_to_matrix(
                                      ORA_merged[[funsys]]@compareClusterResult)
    full_enr_table <- (full_enr_table > pthreshold) + 0
    summ_ora_to_plot[[funsys]] <- list(summ_enr_table=summ_enr_table, 
      summ_enr_clean=summ_enr_clean, full_enr_table=full_enr_table)
  }

  return(summ_ora_to_plot)
}

summarize_categories <- function(all_enrichments, sim_thr = 0.7, 
  common_name = "significant"){
   enrichment_summary <- list()
       # save(list = ls(all.names = TRUE), file = "summarize_environment.RData")
  clusterized_terms <- clusterize_terms(all_enrichments, threshold = sim_thr, 
    common_name = common_name)

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
#    print("time calc_all_paths")
#     print(system.time(
    paths <- calc_all_paths(ids = unique(names(GO_parents)), 
      GO_parents = GO_parents)
#    ))
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


filter_top_categories <- function(enrichments_ORA_merged, top_c = 50){
  # save(enrichments_ORA_merged, file="enrichments_ORA_merged.RData")
  if(! is.null(top_c)) {
    for (funsys in names(enrichments_ORA_merged)){
      filtered_enrichments <- 
        enrichments_ORA_merged[[funsys]]@compareClusterResult
      if (nrow(filtered_enrichments) == 0) next 
      filtered_enrichments <- filtered_enrichments[order(
        filtered_enrichments$p.adjust, decreasing = FALSE), ]
      filtered_enrichments <- Reduce(rbind,by(
        filtered_enrichments,filtered_enrichments["Cluster"], head, n = top_c))
      filtered_terms <- unique(filtered_enrichments$Description)
      enrichments_ORA_merged[[funsys]]@compareClusterResult <- 
        filtered_enrichments
      enrichments_ORA_merged[[funsys]]@termsim <- 
        enrichments_ORA_merged[[funsys]]@termsim[filtered_terms,filtered_terms]
    }
  }
  return(enrichments_ORA_merged)
}


process_cp_list <- function(enrichments_ORA, simplify_results, 
  clean_parentals){
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
      enr_obj <- trycatch_pairwise_termsim(enr_obj)
    }                              
    enrichments_ORA_tr[[funsys]] <- enr_obj 
  }
  return(enrichments_ORA_tr)
}

#' @importFrom GO.db GOOBSOLETE
clean_GO_obsolete <- function(enr_obj) {
  obsolete_terms <- names(as.list(GO.db::GOOBSOLETE))

  if (is(enr_obj, "compareClusterResult")){
    enriched_table <- enr_obj@compareClusterResult
    enriched_table <- enriched_table[!enriched_table$ID %in% obsolete_terms,]
    enr_obj@compareClusterResult <- enriched_table
  } else if (is(enr_obj, "enrichResult")) {
    enriched_table <- enr_obj@result
    enriched_table <- enriched_table[!enriched_table$ID %in% obsolete_terms,]
    enr_obj@result <- enriched_table
  } 
  return(enr_obj)
}

#' @importFrom GO.db GOBPANCESTOR GOMFANCESTOR GOCCANCESTOR
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
  enr_obj <- clean_GO_obsolete(enr_obj)
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

  hamming_DT <- data.table::setDT(hamming_0, key="Var1")
  to_remove_lapply <- lapply(unique(hamming_DT$Var1), function(id) {
    ancs <- GO_ancestors[[id]]
    id_pairs <- hamming_DT[hamming_DT$Var1 == id,]$Var2
    ancs[ancs %in% id_pairs]
  })
  to_remove_lapply <- unique(unlist(to_remove_lapply))
  enr_obj@compareClusterResult <- 
    enrich_obj[!enrich_obj$ID %in% to_remove_lapply,]
  return(enr_obj)
}


#' @importFrom GO.db GOBPANCESTOR GOMFANCESTOR GOCCANCESTOR
clean_topGO_parentals <- function(top_GO_enr, subont, p_value_cutoff = 0.05){
  if (nrow(top_GO_enr) == 0) return(top_GO_enr)

  if (subont=="BP"){
    GO_ancestors <- GO.db::GOBPANCESTOR
  } else if (subont=="MF"){
    GO_ancestors <-  GO.db::GOMFANCESTOR
  } else if (subont=="CC"){
    GO_ancestors <- GO.db::GOCCANCESTOR
  }
  GO_ancestors <- as.list(GO_ancestors)
  GO_ancestors <- GO_ancestors[!is.na(GO_ancestors)]
  obsolete_terms <- names(as.list(GO.db::GOOBSOLETE))
  top_GO_enr <- top_GO_enr[!top_GO_enr$GO.ID %in% obsolete_terms,]
  go_to_keep <- rep(TRUE, nrow(top_GO_enr))
  names(go_to_keep) <- top_GO_enr$GO.ID

  for (go_term in top_GO_enr[top_GO_enr$p.value < p_value_cutoff,"GO.ID"]){ 
    term_stats <- top_GO_enr[top_GO_enr$GO.ID == go_term,]
    go_to_keep[names(go_to_keep) %in% GO_ancestors[[go_term]]] <- FALSE  
  }
  top_GO_enr <- top_GO_enr[go_to_keep,]
  return(top_GO_enr)
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


obtain_gene_attributes <- function(gene_attribute_file, org_db, keytype) {
  gene_attributes <- read.table(gene_attribute_file, header=TRUE)
  gene_attribute_name <- colnames(gene_attributes)[3]
  ids <- translate_ids_orgdb(ids=gene_attributes[,2], input_id = keytype,
           org_db = org_db)
  ids <- unique(ids)
  gene_attributes <- merge(gene_attributes, ids, by.x = "geneid", 
    by.y = keytype)
  gene_attributes <- gene_attributes[!is.na(gene_attributes$ENTREZID),]
  gene_attributes <- gene_attributes[,2:4]
  gene_attr_list <- split(gene_attributes, gene_attributes$cluster)
  gene_attr_list <-lapply(gene_attr_list, function(gene_attr){
    gene_map <- gene_attr[, gene_attribute_name]
    names(gene_map) <- gene_attr$ENTREZID
    return(gene_map)
  })
  return(list(gene_attribute_name = gene_attribute_name, 
    gene_attributes = gene_attr_list))
}

#' @importFrom clusterProfiler simplify
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
      enr <- trycatch_pairwise_termsim(enr)
      enrichments_for_reports[[cluster]][[funsys]] <- enr 
    }
  }
  return(enrichments_for_reports)
}



#' @importFrom ggridges geom_density_ridges
#' @importFrom ggplot2 ggplot aes theme_classic
#' @export
enrich_density <- function(enrich_result, 
                           attributes, 
                           showCategory = 30) {
  attribute <- Description <- p.adjust <- NULL

  terms_attr <- get_attr_by_terms(enrich_result, attributes)
  enr_fort <- utils::getFromNamespace("fortify.enrichResult", "enrichplot")

  enrich_result <- enr_fort(enrich_result,showCategory = showCategory,by = "p.adjust")

  enriched_dist <- merge(terms_attr, enrich_result, by ="ID", all.y = TRUE)
  plot <- ggplot2::ggplot(enriched_dist, ggplot2::aes(x = attribute, 
                                              y = Description, 
                                              fill = p.adjust)) + 
  ggridges::geom_density_ridges(jittered_points = TRUE, 
                                position = "raincloud", 
                                alpha = 0.7, 
                                scale = 0.9, 
                                quantile_lines = TRUE) + 
  ggplot2::theme_classic()

  return(plot)
}

get_attr_by_terms <-function(enrich_result, attributes){
  term_list <- split(enrich_result@result, enrich_result@result$ID)
  gene2Symbol <- enrich_result@gene2Symbol
  attr_list <- lapply(term_list, get_term_attr, 
                      attributes = attributes, 
                      gene2Symbol = gene2Symbol)
 
  term_attr <- data.table::rbindlist(attr_list, 
                                     use.names = TRUE, 
                                     idcol = "ID")
  return(as.data.frame(term_attr))
}

get_term_attr <- function(term_data, attributes, gene2Symbol) {
  g_symbols <- unlist(stringr::str_split(term_data$geneID,"/"))
  g_id <- names(gene2Symbol)[match(g_symbols, gene2Symbol)]
  g_attr <- attributes[match(g_id, names(attributes))]
  names(g_attr) <- NULL
  return(data.frame(attribute = g_attr))
}