###############################################################################
########################## FUNCTIONAL ANALYSIS LIBRARY ########################
###############################################################################

#' Obtains biological info from BiomaRt API
#' @param orthologues orthologue organism ID
#' @param id_type gene code type
#' @param mart BioMart database name you want to connect to
#' @param dataset Dataset you want to use
#' @param host host to connect to
#' @param attr Attributes you want to retrieve
#' @keywords method
#' @importFrom biomaRt useMart getBM
#' @return downloaded biomaRt response
obtain_info_from_biomaRt <- function(orthologues, 
                                     id_type, 
                                     mart, 
                                     dataset, 
                                     host, 
                                     attr){
    # require(biomaRt)

    ensembl <- biomaRt::useMart(mart,dataset=dataset, host=host)
    filt <- id_type
    val <- orthologues
    # attr <- attr
    
    # Check if already exists
    if(file.exists("query_results_temp")){
      container <- readRDS("query_results_temp")
      # Check
      if(nrow(container) != length(attr)){
        warning("Current query results does not match. Will be overwritten") 
        container <- NULL
      } else if(!all(container[,1] %in% val)){
        warning("Current query results does not match. Will be overwritten")
        container <- NULL
      } else if(all(val %in% container[,1])){
        message(paste0("All IDs are already stored at temporal query results",
                       " file. Query will not be performed"))
        return(container)
      } else if(any(val %in% container[,1])){
          # Filter already calculated 
          val <- val[-which(val %in% container[,1])]
      }
    } else{
        container <- NULL           
    }

    # Frag values
    frag_size <- 3000
    indx <- seq(1,length(val),frag_size)

    invisible(lapply(seq_along(indx),function(i){
    # for(i in seq_along(indx)){
        # Check
        if(i == length(indx)){ # Last set
            end <- length(val)
        } else{
            end <- indx[i+1] - 1
        }
        interval <- seq(indx[i],end)
        message(paste("Fragment: ",i,"/",length(indx),"  (",
                length(interval),")",sep=""))

        # Run query
        query <- biomaRt::getBM(attributes = attr, 
                                filters = filt, 
                                values = val[interval], 
                                mart = ensembl)
        # Store
        if(is.null(container)){
            container <- query
        } else if(nrow(query) > 0){
            container <- rbind(container,query)
        }

        # Check
        if(!all(val[interval] %in% container[,1])){
            miss <- val[interval]
            miss <- which(!miss %in% container[,1])
            warning(paste("There are (",length(miss),
                          ") without results from the API",sep=""))
        }

        # Save
        saveRDS(container, file = file.path("query_results_temp"))

        # Wait
        if(i != length(indx)){
            Sys.sleep(5)            
        }
    # }
    }))

    # Return
    return(container)
}


#' Translates a given gene ID
#' @param input_ids gene IDs to be translated
#' @param organism_db BiocGenerics organism database
#' @param org_var_name BiocGenerics organism translation variable name
#' @keywords translate
#' @importFrom BiocGenerics get
#' @importFrom AnnotationDbi mappedkeys
#' @return translated gene IDs
id_translation_orgdb <- function(input_ids, organism_db, org_var_name){
    # Load necessary package
    require(organism_db, character.only = TRUE)
    # Obtain target variable
    org_var <- BiocGenerics::get(org_var_name)
    # Translate ENSEMBL to Entrex IDs
    translation_table <- as.list(org_var[AnnotationDbi::mappedkeys(org_var)])
    # Convert to dataframe
    translation_table_df <- as.data.frame(do.call(rbind,
      lapply(intersect(input_ids, names(translation_table)),function(matches){
        # Obtain genes
        translated_ids <- translation_table[[matches]]            
 
        if(length(translated_ids) == 0){
            return(data.frame())
        }
        # Return info
        return(data.frame(input = rep(matches,length(translated_ids)), 
                          output = translated_ids, 
                          stringsAsFactors = FALSE))
    })))
    return(translation_table_df)
}


#' Scale a matrix using minimum-maximum method
#' @param data_matrix to be scaled
#' @param norm_by_col boolean flag: if true scaling will be performed 
#' by columns intead of by rows. Default: FALSE
#' @importFrom matrixStats rowRanges rowDiffs
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

#' Catched errors fo pairwise_termsim for special cases
#' @param enr enrichment object to be studied
#' @param num_cats number of categories to be shown
#' @return enrichment object after add termsim info
#' @importFrom enrichplot pairwise_termsim
catched_pairwise_termsim <- function(enr, num_cats = 200){
  endedWithoutERRs <- FALSE
  preparedForFortify <- FALSE
  initial_num_cats <- num_cats
  res <- enr
  # Try
  while(!endedWithoutERRs){
    tryCatch(
      # MAIN
      {
        enr <- enrichplot::pairwise_termsim(res, showCategory = num_cats)
        endedWithoutERRs <- TRUE
      },
      # CATCH
      error = function(cond){
        if(!preparedForFortify){
          res <<- prepare_for_fortify(res)
          preparedForFortify <<- TRUE
        }else{
          num_cats <<- num_cats - 20
        }
        if(num_cats <= 0){
          stop(cond)
        }
      }
    )
  }
  if(num_cats < initial_num_cats) 
    warning(paste0("Finally number of categories used has been (",
                   num_cats,") for pairwise_termsim"))
  return(enr)
}




#' Performs Over Representation Analysis (ORA) 
#' enrichment of specified ontology
#' @param genes significant genes to be used
#' @param organism target organism
#' @param keyType genes code type
#' @param pvalueCutoff p-value threshold
#' @param pAdjustMethod p-valued adjust method to be applied
#' @param ont ontology to be used. Allowed (GO_MF, GO_CC, GO_BP, KEGG, REACT)
#' @param useInternal used only for KEGG enrichment, activate internal data 
#' usage mode
#' @param qvalueCutoff q-value threshold
#' @param ENRICH_DATA optional enrichment universe already loaded
#' @param semsim flag to indicate if semantic similitud must be calculated. 
#' Necessary for emaplots
#' @keywords enrich
#' @return enrichment table obtained
#' @importFrom clusterProfiler enrichGO enrichKEGG
#' @importFrom ReactomePA enrichPathway
enrichment_ORA <- function(genes,
                           organism,
                           keyType="ENTREZID",
                           pvalueCutoff,
                           pAdjustMethod = "BH",
                           ont,
                           useInternal = FALSE, 
                           qvalueCutoff = 0.2,
                           ENRICH_DATA = NULL, 
                           semsim = TRUE){
  # @import clusterProfiler KEGG.db ReactomePA
  # Check
  if(is.numeric(ont)){
    stop("Ontology specified is a number, not a string")
  }
  # Parse onto
  if(grepl("GO",ont)){
    aux = unlist(strsplit(ont,"_"))
    ont = aux[1]
    go_subonto = aux[2]
  }
  # Take enrichment function
  if(ont == "GO"){
    # require(clusterProfiler)
    enrf <- clusterProfiler::enrichGO
    patter_to_remove <- "GO_DATA *<-"
  } else if(ont == "KEGG"){
    # require(clusterProfiler)
    enrf <- clusterProfiler::enrichKEGG
    patter_to_remove <- "KEGG_DATA *<-"
  } else if(ont == "REACT"){
    # require(ReactomePA)
    enrf <- ReactomePA::enrichPathway
    patter_to_remove <- "Reactome_DATA *<-"
  }

  # Substitute if proceed
  if(!is.null(ENRICH_DATA)){
    # Find not necessary task into code
    ltorem <- grep(patter_to_remove,body(enrf))
    # Check
    if(length(ltorem) == 0){ # Warning, task not found
      warning(paste0("ern_fun: Can not find annot task to be removed.",
                     " Regular version will be used."))
    } else{ # Remove task from code
      if(ont == "GO"){
  body(enrf)[[ltorem]] <- substitute(GO_DATA <- parent.frame()$ENRICH_DATA)
      } else if(ont == "KEGG"){
  body(enrf)[[ltorem]] <- substitute(KEGG_DATA <- parent.frame()$ENRICH_DATA)
      } else if(ont == "REACT"){
body(enrf)[[ltorem]] <- substitute(Reactome_DATA <- parent.frame()$ENRICH_DATA)
      }
    }
  }

  # Check ontology 
  if(ont == "GO"){
    enrichment <- enrf(gene          = genes,
                          OrgDb         = organism,
                          keyType       = keyType,
                          ont           = go_subonto,
                          pvalueCutoff  = pvalueCutoff,
                          pAdjustMethod = pAdjustMethod,
                          qvalueCutoff  = qvalueCutoff) 
  } else if(ont == "KEGG"){
    enrichment <- enrf(gene          = genes,
                          organism      = organism,
                          keyType       = keyType,
                          pvalueCutoff  = pvalueCutoff,
                          pAdjustMethod = pAdjustMethod,
                          use_internal_data = useInternal,
                          qvalueCutoff  = qvalueCutoff)
  } else if(ont == "REACT"){
    enrichment <- enrf(gene          = genes,
                          organism      = organism,
                          pAdjustMethod = pAdjustMethod,
                          pvalueCutoff  = pvalueCutoff,
                          qvalueCutoff  = qvalueCutoff)
  } else{
    stop("Error, ontology specified is not supported to be enriched")
  }

  if(!is.null(enrichment)){
    # enrichment <- prepare_for_fortify(enrichment) 
    if(nrow(enrichment) > 0 && semsim) 
        enrichment <- catched_pairwise_termsim(enrichment)
  }

  # Return enrichment
  return(enrichment)
}


#' Performs several CUSTOM enrichments using ORA method. 
#' You can give universes or load GMT files
#' @param custom_sets custom universes to be used.
#' @param custom_files custom files to load universes to be used
#' @param p_val_threshold p-value threshold
#' @param genes significant genes to be used
#' @param write_path optional output path where enrichment tables will be stored
#' @return list with enrichment tables obtained
#' @keywords enrich
#' @export
#' @importFrom clusterProfiler enricher
#' @importFrom utils write.table
#' @examples
#' enrich_all_customs()
enrich_all_customs <- function(custom_sets = NULL, 
                               custom_files = NULL, 
                               p_val_threshold = 0.05, 
                               genes, 
                               write_path = NULL){
  if(!is.null(custom_sets)){
    custom_set <- custom_sets
    load_files <- FALSE
  }else{
    custom_set <- custom_files
    names(custom_set) <- custom_files
    load_files <- TRUE
  }
  custom_enrichments <- lapply(seq_along(custom_set), function(i) {
      # Load info
      if(load_files){
        custom_gmt <- load_and_parse_gmt(custom_set[i])
        custom_file <- custom_set[i]
      }else{
        custom_gmt <- custom_set[[i]]
        custom_file <- names(custom_set[i])
      }
      enr <- clusterProfiler::enricher(genes, 
                                       pvalueCutoff = p_val_threshold, 
                                       TERM2GENE = custom_gmt)
      if(nrow(enr) > 0) enr <- catched_pairwise_termsim(enr)
      # # Store results
      if(!is.null(write_path)) 
          utils::write.table(enr, 
                 file=file.path(write_path,
                                paste0(basename(custom_file),"_ora_results")), 
                 quote=FALSE, 
                 col.names=TRUE, 
                 row.names = FALSE, 
                 sep="\t")      
      # Return
      return(enr)
    })
  names(custom_enrichments) <- names(custom_set)
  return(custom_enrichments)
}


#' Performs custom enrichments over cluster sets
#' @param custom_set custom unvierses
#' @param genes_in_modules list of genes (clusters)
#' @param p_val_threshold p-value threshold
#' @param cores optional parallel cores to be used. Default: 1
#' @param task_size number of elements per packages used
#' @return enrichment tables obtained
#' @keywords enrich
#' @export
#' @importFrom clusterProfiler enricher
#' @examples
#' # Will return NULL
#' enrich_clusters_with_gmt()
enrich_clusters_with_gmt <- function(custom_set = NULL, 
                                     genes_in_modules = NULL, 
                                     p_val_threshold, 
                                     cores = 1, 
                                     task_size = 1){
      if(is.null(genes_in_modules)) {
        warning("no value for genes_in_modules argument given")
        return(NULL)
      } 
      if(is.null(custom_set)){
        warning("No custom set defined")
        return(NULL)
    }
      modules_enrichment <- parallel_list(genes_in_modules, function(genesset) {
        enr <- clusterProfiler::enricher(genesset, 
                                         pvalueCutoff = p_val_threshold, 
                                         TERM2GENE = custom_set)
        if(nrow(enr) > 0) enr <- catched_pairwise_termsim(enr)
        return(enr)
      }, workers = cores, task_size = task_size)
      names(modules_enrichment) <- names(genes_in_modules)
      # Return
      return(modules_enrichment)
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


#' Performs ORA enrichment over several gene sets (clusters)
#' @param genes list of gene sets (clusters9)
#' @param organism target organism
#' @param keyType gene code type
#' @param pvalueCutoff p-value threshold
#' @param pAdjustMethod p_value adjust method to be applied
#' @param ont ontology to be used. Allowed (GO_MF, GO_CC, GO_BP, KEGG, REACT)
#' @param useInternal used only for KEGG enrichment, activate internal 
#' data usage mode
#' @param qvalueCutoff q-value threshold
#' @param ENRICH_DATA optional enrichment universe already loaded 
#' @param cores optional number of parallel cores to be used. See mcapply
#' @param task_size number of elements per packages used
#' @keywords enrich
#' @return enrichment performed
enrichment_clusters_ORA <- function(genes,
                                    organism,
                                    keyType="ENTREZID",
                                    pvalueCutoff,
                                    pAdjustMethod = "BH",
                                    ont,
                                    useInternal = FALSE, 
                                    qvalueCutoff, 
                                    ENRICH_DATA = NULL, 
                                    cores = 1, 
                                    task_size = 1){
  src_ont = ont
  if(grepl("GO",ont)){
    aux = unlist(strsplit(ont, "_"))
    ont = aux[1]
    go_subonto = aux[2]
  }

  # Prepare set
  if(is.null(ENRICH_DATA) & length(genes) > 1){
    if(ont == "GO"){
      # require(clusterProfiler)
      get_GO_data <- get_unexported_function("clusterProfiler", "get_GO_data")
      ENRICH_DATA <- get_GO_data(organism, go_subonto, keyType)
    } else if(ont == "KEGG"){
      # require(clusterProfiler)
      get_data_from_KEGG_db <- get_unexported_function("clusterProfiler", 
                                                       "get_data_from_KEGG_db")
      organismMapper <- get_unexported_function("clusterProfiler", 
                                                "organismMapper")
      ENRICH_DATA <- get_data_from_KEGG_db(organismMapper(organism))
    } else if(ont == "REACT"){
      # require(ReactomePA)
      get_Reactome_DATA <- get_unexported_function("ReactomePA", 
                                                   "get_Reactome_DATA")
      ENRICH_DATA <- get_Reactome_DATA(organism)
    }
  }

  # require(parallel)
  # Enrich per gene_set
  enrichment <- parallel_list(genes, function(setg){
    # Obtain enrichment 
    enr <- enrichment_ORA(genes         = setg,
                          organism      = organism,
                          keyType       = keyType,
                          pvalueCutoff  = pvalueCutoff,
                          pAdjustMethod = pAdjustMethod,
                          ont           = src_ont,
                          useInternal   = useInternal, 
                          qvalueCutoff  = qvalueCutoff,
                          ENRICH_DATA   = ENRICH_DATA,
                          semsim = TRUE)
    # Return
    return(enr)
  }, workers= cores, task_size = task_size)

  names(enrichment) <- names(genes)

  # Return enrichment
  return(enrichment)
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

get_translation_tables_orgdb <- function(input_gene_id, input_ids, current_organism_info) {
  org_db <- get_org_db(current_organism_info)
  if(input_gene_id == "ENTREZID") {
    input_to_entrezgene <- data.frame(input=input_ids, 
                                      ENTREZID=input_ids)
  } else {
    # Check input gene ID valid
    check_id_valid_orgdb(gene_id=input_gene_id, id_type="input", organism_info=current_organism_info, outcome_action="stop")

        input_to_entrezgene <- translate_ids_orgdb(ids=input_ids, 
        input_id=input_gene_id, org_db=org_db) 
  }
  symbol_output_available <- check_id_valid_orgdb(gene_id=input_gene_id, id_type="output", 
                                                  organism_info=current_organism_info, outcome_action="warning")
            
  if(symbol_output_available == TRUE) {
    input_to_symbol <- translate_ids_orgdb(ids=input_ids, 
    input_id=input_gene_id, output_id="SYMBOL", org_db=org_db)
  } else {
    input_to_symbol <- NULL
  }
  return(list(input_to_entrezgene = input_to_entrezgene, input_to_symbol = input_to_symbol))
}

translate_ids_orgdb <- function(ids, input_id, output_id="ENTREZID", org_db=org_db, just_output_ids=FALSE){
  possible_ids <- AnnotationDbi::columns(org_db)
  if(! input_id %in% possible_ids) 
    stop(paste(c("gene keytype must be one of the following:", possible_ids), collapse=" "))
    ids <- tryCatch(
      ids <- AnnotationDbi::select(org_db, keys=ids, column=output_id, keytype=input_id),
      error=function(cond){
            ids <- NULL
        }
    )
    ids <- ids[!is.na(ids[,2]),]
    if(just_output_ids == TRUE) {
      return(unique(ids[,2]))
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
  organism_info, gene_id="entrez", algorithm = "classic", statistic = "fisher", 
  nodeSize = 5, task_size=1, workers=1, ...){

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
     

    # Create environment variables used for initialization of the topGOdata obj
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
      resultFis <- topGO::runTest(l_GOdata, algorithm = algorithm, 
      statistic = statistic)
    }, workers = workers, task_size = task_size )
    
    if(unlisted_input_flag) enriched_cats <- enriched_cats[[1]]
    enrichments_topGO[[funsys]] <- enriched_cats
  }
  return(enrichments_topGO)
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

  common_params <- list(universe = universe, pvalueCutoff = pvalueCutoff, 
    qvalueCutoff = qvalueCutoff, 
    pAdjustMethod = pAdjustMethod, ...)

  if(! is.null(custom_sets)) {
    if(is.null(names(custom_sets))) stop("Custom sets enrichment object must be a named list")
    all_funsys <- c(all_funsys, names(custom_sets))
  }
  if (is.null(org_db)) {
    org_db <- get_org_db(organism_info)
  }
  enrichments_ORA <- vector("list", length(all_funsys))
  names(enrichments_ORA) <- all_funsys
  for(funsys in all_funsys) {
    if (funsys %in% c("CC","BP","MF")){
      enrf <- prepare_enrichment_GO(enrichment_type="ora", subont = funsys, 
        org_db = org_db)
      specific_params <- list(OrgDb = org_db, ont = funsys)

    } else  if (funsys == "Reactome"){
      enrf <- prepare_enrichment_Reactome(enrichment_type="ora", 
        reactome_id = organism_info$Reactome_ID[1])
      specific_params <- list(organism = organism_info$Reactome_ID[1])

    } else if (funsys == "KEGG"){
      enrf <- prepare_enrichment_KEGG(enrichment_type="ora", 
        kegg_file = kegg_file)
      specific_params <- list(organism = organism_info$KeggCode[1])
                       
    } else if (funsys %in% names(custom_sets)) {
      enrf <- clusterProfiler::enricher
      specific_params <- list(TERM2GENE = custom_sets[[funsys]])

    } else {
      stop("funsys", funsys, "not recognized")
    }
    enriched_cats <- parallel_list(genes_list, function(l_genes){
        params_genes <- c(specific_params, common_params, list(gene = l_genes))
        enr <- do.call("enrf", params_genes)
      }, 
      workers= workers, task_size = task_size
    )
    # save(enriched_cats, file=paste0("enriched_cats_", funsys, ".RData"))
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
        enr@result[,"Description"] <- gsub("^.+\\r: ", "", enr@result[,"Description"])
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


merge_clusters <- function(results_list) {
  merged_clusters <- lapply(results_list, clusterProfiler::merge_result)
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

filter_top_categories <- function(enrichments_ORA_merged, top_c = 50){
    # save(enrichments_ORA_merged, file="enrichments_ORA_merged.RData")
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
      enr_obj <- catched_pairwise_termsim(enr_obj, 200)
    }                              
    enrichments_ORA_tr[[funsys]] <- enr_obj 
  }
  return(enrichments_ORA_tr)
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

