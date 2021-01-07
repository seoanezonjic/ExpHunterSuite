#################################################################################
############################ FUNCTIONAL ANALYSIS LIBRARY ########################
#################################################################################

# #'
# #' @param opt list with options where folders will be taken
# #' @keywords
# #' @return
# defining_functional_hunter_subfolders <- function(opt){
#   subfolders <- c()
#   if (grepl("G", opt$functional_analysis)){
#     subfolders <- c(subfolders, 'Results_topGO')
#   }
#   if (grepl("K", opt$functional_analysis)){
#     subfolders <- c(subfolders, 'Results_KEGG')
#   }
# }

#' Function used to check if a specific result has been calculated
#' @param var_name variable name to be checked
#' @param level_depth variable depth into a list of objects
#' @keywords check
#' @return true if variable exists or false in other cases
exists_enrich_df <- function(var_name,level_depth = 3){
  # Check if it is instantiated
  if(!var_name %in% ls(envir = parent.frame(n = level_depth))){
    return(FALSE)
  }
  # Check if is a DF with info
  df <- get(var_name, envir = parent.frame(n = level_depth))
  return(check_results(df))
}


#' Generates the full path to a specific graph file
#' @param path to file
#' @param file_name file name
#' @param ext file extension
#' @keywords file
#' @return full file path and name
graph_file_name <- function(path, file_name, ext = "pdf"){
  return(paste(path,.Platform$file.sep,file_name,".",ext,sep=""))
}


#' Checks  if a graph file exists
#' @param path to file
#' @param file_name file name
#' @param ext file extension
#' @keywords file
#' @return true if file exists or false in other cases
exists_graph_file <- function(path, file_name, ext = "pdf"){
  f <- graph_file_name(path,file_name,ext)
  return(file.exists(f))
}


#' Check if a given object is a correct enrichment result object
#' @param df variable to be checked
#' @keywords check
#' @return true if it is an allowed object or false in other cases
check_results <- function(df){
  if(is.data.frame(df)){
    if(nrow(df) <= 0){
      return(FALSE)
    }
  }else if(typeof(df)=="S4"){
    if("pvalueCutoff" %in% slotNames(df)){
      if(nrow(df@result[which(df@result$p.adjust <= df@pvalueCutoff),]) < 1){
        return(FALSE)
      }
    }else if(length(get_categories(df)) < 1){
      return(FALSE)
    }
    return(TRUE)
  }else if(!is.list(df) & !typeof(df)=="S4"){
    return(FALSE)
  }else if(is.list(df)){
    if(!all(unlist(lapply(df,function(set){
      if(typeof(set)=="S4"){
        if("pvalueCutoff" %in% slotNames(set)){
          if(nrow(set@result[which(set@result$p.adjust <= set@pvalueCutoff),]) >= 1){
            return(TRUE)
          }
        }else if("params" %in% slotNames(set)){
          if(nrow(set@result[which(set@result$p.adjust <= set@params$pvalueCutoff),]) >= 1){
            return(TRUE)
          }
        }else{
          return(FALSE)
        }
        return(FALSE)
      } else if (is.null(set)) {
        return(FALSE)
      }

      return(TRUE)
    })))){ # IF
      return(FALSE)
    }
  }
  # Everything OK
  return(TRUE)
}


#' Translates a given gene ID using a dictionary
#' @param ids_to_translate set of IDs to be translated
#' @param annot_table dictionary to translate IDs
#' @keywords translate
#' @return translated IDs or NA if it's not possible to translate
translate_id <- function(ids_to_translate, annot_table){ 
#This function translates unknown transcripts IDs to known gen IDs
#ids_to_translate: unknown gene IDs
#annot_table: DF, KNOWN gene IG on first column and UNKNOWN on second
#ONE unknown ID can correspond to MANY Known IDs  
 # save(list = ls(all.names = TRUE), file = "test.RData", envir = environment())

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
obtain_info_from_biomaRt <- function(orthologues, id_type, mart, dataset, host, attr){
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
        message("All IDs are already stored at temporal query results file. Query will not be performed")
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
        message(paste("Fragment: ",i,"/",length(indx),"  (",length(interval),")",sep=""))

        # Run query
        query <- biomaRt::getBM(attributes = attr, filters = filt, values = val[interval], mart = ensembl)
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
            warning(paste("There are (",length(miss),") without results from the API",sep=""))
            # # Add empty lines
            # aux_NAs <- rep(NA,length(miss))
            # to_concat <- cbind(miss,as.data.frame(matrix(aux_NAs,length(aux_NAs), ncol(container) - 1)))
            # # Concat
            # container <- rbind(container,to_concat)
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
    translation_table_df <- as.data.frame(do.call(rbind,lapply(intersect(input_ids, names(translation_table)),function(matches){
        # Obtain genes
        translated_ids <- translation_table[[matches]]            
 
        if(length(translated_ids) == 0){
            return(data.frame())
        }
        # Return info
        return(data.frame(input = rep(matches,length(translated_ids)), output = translated_ids, stringsAsFactors = FALSE))
    })))
    return(translation_table_df)
}



#' Pipeline to perform topGO enrichment
#' @param entrez_targets genes used to enrich
#' @param entrez_universe total of allowed genes
#' @param sub_ontology GO subontology to be used
#' @param outFile output path
#' @param organism target organism
#' @param plot_graph boolean to generates a plot 
#' @keywords enrich
#' @importFrom topGO runTest GenTable showSigOfNodes
#' @importFrom BiocGenerics score
#' @importFrom utils write.table
#' @importFrom grDevices pdf dev.off
#' @return enrichment tables obtained
perform_topGO_local <- function(entrez_targets,entrez_universe,sub_ontology,outFile=NULL,organism, plot_graph = TRUE){
  #! GOFisherTest
  # require(topGO)
  
  # Prepare necessary info
  allGenes <- unique(entrez_universe) 
  all_results_table <- data.frame() 
  
  # creating contigency vector for gene identifiers
  geneList <- factor(as.integer(allGenes %in% entrez_targets)) 
  names(geneList) <- allGenes
  
  if(length(levels(geneList)) == 2){ # launch analysis
    TopGOobject <- new("topGOdata", ontology = sub_ontology, allGenes = geneList, annot = topGO::annFUN.org, mapping = organism)
    # # Possible option 2
    results_fisher  <- topGO::runTest(TopGOobject, algorithm = "classic", statistic = "fisher") 
    results_KS      <- topGO::runTest(TopGOobject, algorithm = "classic", statistic = "ks")
    results_KS_elim <- topGO::runTest(TopGOobject, algorithm = "elim", statistic = "ks")

    # Create results table
    all_results_table <- topGO::GenTable(TopGOobject, 
                                  classicFisher = results_fisher,
                                  classicKS = results_KS,
                                  elimKS = results_KS_elim,
                                  orderBy = "classicFisher",
                                  ranksOf = "classicFisher", 
                                  topNodes = length(TopGOobject@graph@nodes)) 
  }
    
  # Write info
  if(!is.null(outFile)){
    utils::write.table(all_results_table, file=outFile, quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")    
  }

  # Plot graphs
  if(plot_graph){
    grDevices::pdf(paste(outFile,"pdf",sep="."), w=11, h=8.5)
      topGO::showSigOfNodes(TopGOobject, BiocGenerics::score(results_fisher), firstSigNodes = 10, useInfo = 'all')
    grDevices::dev.off()
  }
  
  # Return info
  return(all_results_table)
}


#' Scale a matrix using minimum-maximum method
#' @param data_matrix to be scaled
#' @param transpose boolena which indicates if matrix must be transposed. Default: FALSE
#' @keywords method
#' @return scaled matrix
scale_data_matrix <- function(data_matrix, transpose = FALSE){
  if ( transpose ) {
    data_matrix <- t(data_matrix)
  } 
  scaled_counts <- apply(data_matrix, 2, function(column) {
    minimum <- min(column, na.rm=TRUE)
    maximum <- max(column, na.rm=TRUE)
    difference <- maximum - minimum
    if (difference == 0) { 
      scaled_column <- column
    } else {
      scaled_column <- (column - minimum) / difference
    }
    return(scaled_column)
  })
  if ( transpose ) {
    scaled_counts <- t(scaled_counts)
  } 
  return(scaled_counts)
}




#' Run topGO enrichment and generates target plots storing into PDF files
#' @param attr_name target column of DEG_annot_table with attributes
#' @param interesting_genenames target genes
#' @param DEG_annot_table DEgenes Hunter expression result table
#' @param ontology GO subontology to be used
#' @param graphname output file name
#' @param filter_name target column of DEG_annot_table with filtering param
#' @param output_files output path
#' @return void
#' @importFrom topGO runTest GenTable showSigOfNodes
#' @importFrom utils write.table
#' @importFrom grDevices pdf dev.off
#' @importFrom BiocGenerics score
#' @keywords enrich
perform_topGO <- function(attr_name, interesting_genenames, DEG_annot_table, ontology, graphname,filter_name, output_files){
    geneID2GO <- split(DEG_annot_table[,attr_name], DEG_annot_table[,filter_name])
    geneID2GO <- lapply(geneID2GO, unique)
    geneNames <- names(geneID2GO)
    geneList <- factor(as.integer(geneNames %in% interesting_genenames))
    names(geneList) <- geneNames

    GOdata <- new("topGOdata", ontology =ontology, allGenes = geneList, annot=topGO::annFUN.gene2GO, gene2GO = geneID2GO)
    resultFis <- topGO::runTest(GOdata, algorithm = "classic", statistic = "fisher") 

    resultKS <- topGO::runTest(GOdata, algorithm = "classic", statistic = "ks")
    resultKS.elim <- topGO::runTest(GOdata, algorithm = "elim", statistic = "ks")

    allRes <- topGO::GenTable(GOdata, classicFisher = resultFis,
        classicKS = resultKS, elimKS = resultKS.elim,
        orderBy = "elimKS", ranksOf = "classicFisher", topNodes = length(GOdata@graph@nodes))

    utils::write.table(allRes, file=file.path(output_files, "allResGOs.txt"), quote=FALSE, col.names=NA, sep="\t")

    grDevices::pdf(file.path(output_files, graphname), w=11, h=8.5)
        topGO::showSigOfNodes(GOdata, BiocGenerics::score(resultFis), firstSigNodes = 10, useInfo = 'all')
    grDevices::dev.off()

}



#' Catched errors fo pairwise_termsim for special cases
#' @param enr enrichment object to be studied
#' @return enrichment object after add termsim info
#' @importFrom enrichplot pairwise_termsim
catched_pairwise_termsim <- function(enr){
  endedWithoutERRs <- FALSE
  initial_num_cats <- 200
  num_cats <- initial_num_cats
  res <- enr
  while(!endedWithoutERRs){
    tryCatch(
      # MAIN
      {
        res <- enrichplot::pairwise_termsim(enr, showCategory = num_cats)
        endedWithoutERRs <- TRUE
      },
      # CATCH
      error = function(cond){
        num_cats <<- num_cats - 20
        if(num_cats <= 0){
          stop(cond)
        }
      }
    )
  }
  if(num_cats < initial_num_cats) warning(paste0("Finally number of categories used has been (",num_cats,") for pairwise_termsim"))
  return(res)
}




#' Performs Over Representation Analysis (ORA) enrichment of specified ontology
#' @param genes significant genes to be used
#' @param organism target organism
#' @param keyType genes code type
#' @param pvalueCutoff p-value threshold
#' @param pAdjustMethod p-valued adjust method to be applied
#' @param ont ontology to be used. Allowed (GO_MF, GO_CC, GO_BP, KEGG, REACT)
#' @param useInternal used only for KEGG enrichment, activate internal data usage mode
#' @param qvalueCutoff q-value threshold
#' @param ENRICH_DATA optional enrichment universe already loaded
#' @param semsim flag to indicate if semantic similitud must be calculated. Necessary for emaplots
#' @keywords enrich
#' @return enrichment table obtained
#' @importFrom clusterProfiler enrichGO enrichKEGG
#' @importFrom ReactomePA enrichPathway
enrichment_ORA <- function(genes,organism,keyType="ENTREZID",pvalueCutoff,pAdjustMethod = "BH",ont,useInternal = FALSE, qvalueCutoff = 0.2,ENRICH_DATA = NULL, semsim = TRUE){
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
    enr_fun <- clusterProfiler::enrichGO
    patter_to_remove <- "GO_DATA *<-"
  } else if(ont == "KEGG"){
    # require(clusterProfiler)
    enr_fun <- clusterProfiler::enrichKEGG
    patter_to_remove <- "KEGG_DATA *<-"
  } else if(ont == "REACT"){
    # require(ReactomePA)
    enr_fun <- ReactomePA::enrichPathway
    patter_to_remove <- "Reactome_DATA *<-"
  }

  # Substitute if proceed
  if(!is.null(ENRICH_DATA)){
    # Find not necessary task into code
    line_to_remove <- grep(patter_to_remove,body(enr_fun))
    # Check
    if(length(line_to_remove) == 0){ # Warning, task not found
      warning("ern_fun: Can not find annot task to be removed. Regular version will be used.")
    } else{ # Remove task from code
      if(ont == "GO"){
        body(enr_fun)[[line_to_remove]] <- substitute(GO_DATA <- parent.frame()$ENRICH_DATA)      
      } else if(ont == "KEGG"){
        body(enr_fun)[[line_to_remove]] <- substitute(KEGG_DATA <- parent.frame()$ENRICH_DATA)      
      } else if(ont == "REACT"){
        body(enr_fun)[[line_to_remove]] <- substitute(Reactome_DATA <- parent.frame()$ENRICH_DATA)              
      }
    }
  }

  # Check ontology 
  if(ont == "GO"){
    enrichment <- enr_fun(gene          = genes,
                          OrgDb         = organism,
                          keyType       = keyType,
                          ont           = go_subonto,
                          pvalueCutoff  = pvalueCutoff,
                          pAdjustMethod = pAdjustMethod,
                          qvalueCutoff  = qvalueCutoff) 
  } else if(ont == "KEGG"){
    enrichment <- enr_fun(gene          = genes,
                          organism      = organism,
                          keyType       = keyType,
                          pvalueCutoff  = pvalueCutoff,
                          pAdjustMethod = pAdjustMethod,
                          use_internal_data = useInternal,
                          qvalueCutoff  = qvalueCutoff)
  } else if(ont == "REACT"){
    enrichment <- enr_fun(gene          = genes,
                          organism      = organism,
                          pAdjustMethod = pAdjustMethod,
                          pvalueCutoff  = pvalueCutoff,
                          qvalueCutoff  = qvalueCutoff)
  } else{
    stop("Error, ontology specified is not supported to be enriched")
  }

  if(!is.null(enrichment)){
    if(nrow(enrichment) > 0 && semsim) enrichment <- catched_pairwise_termsim(enrichment)
  }

  # Return enrichment
  return(enrichment)
}


#' Performs several CUSTOM enrichments using ORA method. You can give universes or load GMT files
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
enrich_all_customs <- function(custom_sets = NULL, custom_files = NULL, p_val_threshold = 0.05, genes, write_path = NULL){
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
      enr <- clusterProfiler::enricher(genes, pvalueCutoff = p_val_threshold, TERM2GENE = custom_gmt)
      if(nrow(enr) > 0) enr <- catched_pairwise_termsim(enr)
      # # Store results
      if(!is.null(write_path)) utils::write.table(enr, file=file.path(write_path,paste0(basename(custom_file),"_ora_results")), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")      
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
enrich_clusters_with_gmt <- function(custom_set, genes_in_modules = NULL, p_val_threshold, cores = 1, task_size = 1){
      if(is.null(genes_in_modules)) {
        warning("no value for genes_in_modules argument given")
        return(NULL)
      }
      modules_enrichment <- parallel_list(genes_in_modules, function(genesset) {
        enr <- clusterProfiler::enricher(genesset, pvalueCutoff = p_val_threshold, TERM2GENE = custom_set)
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
load_and_parse_gmt <- function(gmt_file) {
    # Load file
    gmt <- readLines(con = gmt_file)
    gmt_list <- strsplit(gmt, "\t")
    parsed_gmt <- do.call(rbind, lapply(gmt_list, function(category) {
          category_name <- category[1]
          genes <-category[3:length(category)]
          parsedTerms <- data.frame(Term = category_name, Gene= genes, stringsAsFactors = FALSE)
          return(parsedTerms)

    }))
    return(parsed_gmt)
}


#' Performs Gene Set Enrichment Analysis (GSEA) using a specified ontology
#' @param geneList significant genes to be used
#' @param organism target organism
#' @param keyType gene code type
#' @param pvalueCutoff p-value threshold
#' @param pAdjustMethod p-value adjust method to be applied
#' @param ont ontology to be used. Allowed (GO_MF, GO_CC, GO_BP, KEGG, REACT)
#' @param useInternal used only for KEGG enrichment, activate internal data usage mode
#' @return enrichment performed
#' @keywords enrich
#' @importFrom clusterProfiler gseGO gseKEGG
#' @importFrom ReactomePA gsePathway
enrich_GSEA <- function(geneList,organism,keyType="ENTREZID",pvalueCutoff,pAdjustMethod = "BH",ont,useInternal = FALSE){
  # @import clusterProfiler KEGG.db ReactomePA
  # require(clusterProfiler)
  # if(useInternal)
  #   require(KEGG.db)
  # require(ReactomePA)
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
  # Check ontology 
  if(ont == "GO"){
    require(organism, character.only = TRUE)
    enrichment <- clusterProfiler::gseGO(geneList      = geneList,
                                        OrgDb         = get(organism),
                                        keyType       = keyType,
                                        ont           = go_subonto,
                                        pvalueCutoff  = pvalueCutoff,
                                        pAdjustMethod = pAdjustMethod)
  } else if(ont == "KEGG"){
    enrichment <- clusterProfiler::gseKEGG(geneList     = geneList,
                                          organism     = organism,
                                          use_internal_data = useInternal,
                                          # nPerm        = 1000,
                                          # minGSSize    = 120,
                                          pvalueCutoff = pvalueCutoff,
                                          verbose      = FALSE)
  } else if(ont == "REACT"){
    enrichment<- ReactomePA::gsePathway(geneList, 
                                        organism = organism,
                                        # exponent = 1, 
                                        # nPerm = 1000,
                                        # minGSSize = 10, 
                                        # maxGSSize = 500, 
                                        pvalueCutoff = pvalueCutoff,
                                        pAdjustMethod = pAdjustMethod)
  } else{
    stop("Error, ontology specified is not supported to be enriched")
  }

  # Return enrichment
  if (nrow(enrichment) == 0){
    return(NULL)
  } else {
    return(enrichment)
  }
}


#' Performs GSEA enrichment over several gene sets (clusters)
#' @param all_clusters list of gene sets (clusters)
#' @param organism target organism
#' @param keyType gene code type
#' @param pvalueCutoff p-value threshold
#' @param pAdjustMethod p-value adjust method to be applied
#' @param ont ontology to be used
#' @param useInternal optional KEGG param used to indicate if already downloaded KEGG database must be used or online API must be called. DEfault: FALSE
#' @keywords enrich
#' @return enrichment tables obtained
perform_GSEA_clusters <- function(all_clusters, organism, keyType, pvalueCutoff, pAdjustMethod = "BH", ont, useInternal = FALSE){
  enriched_clusters <- lapply(all_clusters, function(cl_genes) {
        # Enrich
        cl_GSEA <- enrich_GSEA(geneList = cl_genes,
                      organism = organism,
                      keyType = keyType,
                      pvalueCutoff = pvalueCutoff,
                      pAdjustMethod = pAdjustMethod,
                      ont = ont, 
                      useInternal = useInternal)
        # Return
        return(cl_GSEA)
  })
}



#' Performs ORA enrichment over several gene sets (clusters)
#' @param genes list of gene sets (clusters9)
#' @param organism target organism
#' @param keyType gene code type
#' @param pvalueCutoff p-value threshold
#' @param pAdjustMethod p_value adjust method to be applied
#' @param ont ontology to be used. Allowed (GO_MF, GO_CC, GO_BP, KEGG, REACT)
#' @param useInternal used only for KEGG enrichment, activate internal data usage mode
#' @param qvalueCutoff q-value threshold
#' @param ENRICH_DATA optional enrichment universe already loaded 
#' @param cores optional number of parallel cores to be used. See mcapply
#' @param task_size number of elements per packages used
#' @keywords enrich
#' @return enrichment performed
enrichment_clusters_ORA <- function(genes,organism,keyType="ENTREZID",pvalueCutoff,pAdjustMethod = "BH",ont,useInternal = FALSE, qvalueCutoff, ENRICH_DATA = NULL, cores = 1, task_size = 1){
  # @import clusterProfiler KEGG.db ReactomePA parallel
  # Parse onto
  # save(list = ls(all.names = TRUE), file = "test.RData", envir = environment())
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
      get_GO_data <- privateFun("clusterProfiler","get_GO_data")
      ENRICH_DATA <- get_GO_data(organism, go_subonto, keyType)
    } else if(ont == "KEGG"){
      # require(clusterProfiler)
      get_data_from_KEGG_db <- privateFun("clusterProfiler", "get_data_from_KEGG_db")
      organismMapper <- privateFun("clusterProfiler", "organismMapper")
      ENRICH_DATA <- get_data_from_KEGG_db(organismMapper(organism))
    } else if(ont == "REACT"){
      # require(ReactomePA)
      get_Reactome_DATA <- privateFun("ReactomePA", "get_Reactome_DATA")
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



#' Returns correct organism ID to be used
#' @param organism_info organism table entry
#' @param ont ontology to be used
#' @return organism ID to be used
get_organismID_byOnto <- function(organism_info, ont){
  if(grepl("GO",ont)){
    return(organism_info$Bioconductor_DB[1])
  }else if(grepl("KEGG",ont)){
    return(organism_info$KeggCode[1])
  }else if(grepl("REACT",ont)){
    return(organism_info$Reactome_ID[1])
  }else{
    return(NULL)
  }
}


#' Table with information abaut all organism available
#' @param file to be loaded. Default: internal organism table
#' @return organism table
#' @keywords method
#' @export
#' @importFrom utils read.table
#' @examples
#' ot <- get_organism_table()
get_organism_table <- function(file = file.path(find.package('ExpHunterSuite'), "external_data", "organism_table.txt")){
  return(utils::read.table(file, header = TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE, fill = NA))
}


#' Perform enrichment of ORA and/or GSEA using a set of genes given
#' @param genes significant genes to be used (vector or list of vectors - clusters)
#' @param organism_info organism table info. Must include entries depending on ontologies selected: Bioconductor_DB (GO), KeggCode (KEGG), Reactome_ID (Reactome)
#' @param keytype gene code type. If is a vector, will be used in the same order than ontology items are given
#' @param ontology string with ontologIDs to be used. Allowed: Gene Ontology (BP: b ; MF: m ; CC: c), KEGG (k), Reactome (r) 
#' @param enrichmentType enrichment technique to be used. Allowed: ORA (o) and GSEA (g)
#' @param pvalueCutoff p-value threshold
#' @param pAdjustMethod p-value adjust method to be applied
#' @param qvalueCutoff q-value threshold. ONLY USED FOR ORA.
#' @param useInternal optional KEGG internal usage flag
#' @param cores for parallel. Only used in case of genes are cluster genes
#' @param task_size number of elements per packages used
#' @param verbose activate verbose mode
#' @return list with enrichments performed
#' @keywords enrich
#' @importFrom pbapply pblapply
#' @importFrom clusterProfiler merge_result
#' @export
#' @examples
#' data(degh_output)
#' genes <- head(rownames(degh_output$raw_filter),1000) # You need to translate to entrez if you want use all ontology types
#' ontologies <- "" # Select your wanted ontologies
#' enrch <- multienricher(genes = genes, ontology = "")
multienricher <- function(genes,
                          organism_info,
                          keytype = "ENTREZID",
                          ontology = "bkr",
                          enrichmentType = "o",
                          pvalueCutoff = 0.01,
                          pAdjustMethod = "BH",
                          useInternal = TRUE,
                          qvalueCutoff = 0.02,
                          cores = 1,
                          task_size = 1,
                          verbose = FALSE){
  # Initialize
  func_results <- list()
  flags <- list()
  onts <- c()

  # Create flags
  flags$ORA <- grepl("o", enrichmentType)
  flags$GSEA <- grepl("g", enrichmentType)
  flags$clusters <- is.list(genes)
  allowed_onts <- c("b" = "GO_BP", "m" = "GO_MF", "c" = "GO_CC", "k" = "KEGG", "r" = "REACT")
  onts <- unlist(lapply(seq_along(allowed_onts),function(i){if(grepl(names(allowed_onts)[i], ontology)) return(allowed_onts[[i]])}))

  # Check genes format
  if(flags$clusters){
    func_results[["WGCNA"]] <- list()
    if(is(genes[[1]], "character")){
      genes_ora <- genes
      genes_gsea <- NULL
    }else{
      genes_ora <- lapply(genes,names)
      genes_gsea <- genes
    }
  }else if(is(genes, "character")){ # Vector with genes
    genes_ora  <- genes
    genes_gsea <- NULL 
  }else{ # or named vector with logFC
    genes_ora  <- unique(names(genes))
    genes_gsea <- genes 
  }

  apply_fun <- lapply
  if(verbose){
    apply_fun <- pbapply::pblapply
  }

  # Prepare keytypes
  if(length(keytype) > 1){
    if(length(keytype) != length(onts)) stop("Given keytypes have not the same length than ontologies specified")
    ontology_splitted <- unlist(strsplit(ontology,""))
    keytypes <- lapply(onts,function(o){
      indx_ont <- which(allowed_onts == o)
      return(keytype[which(ontology_splitted == names(allowed_onts)[indx_ont])])
    })    
  }else{
    keytypes <- lapply(onts,function(o){keytype})    
  }
  names(keytypes) <- onts

  ## ORA
  if(flags$ORA){
    if(flags$clusters){
      if(verbose) message(paste0("Performing clusters ORA enrichment for (",length(onts),") ontologies"))
      func_results[["WGCNA"]][["ORA_expanded"]] <- apply_fun(onts, function(ont){
        enrichment_clusters_ORA(genes = genes_ora,
                                organism = get_organismID_byOnto(organism_info,ont),
                                keyType = keytypes[[ont]],
                                pvalueCutoff = pvalueCutoff,
                                pAdjustMethod = pAdjustMethod,
                                ont = ont,
                                useInternal = useInternal,
                                qvalueCutoff = qvalueCutoff,
                                cores = cores,
                                task_size = task_size)
      }) 
      names(func_results$WGCNA$ORA_expanded) <- onts
    }else{
      if(verbose) message(paste0("Performing ORA enrichment for (",length(onts),") ontologies"))
      func_results[["ORA"]] <- apply_fun(onts, function(ont){
        enrichment_ORA(genes = genes_ora,
                       organism = get_organismID_byOnto(organism_info,ont),
                       keyType = keytypes[[ont]],
                       pvalueCutoff = pvalueCutoff,
                       pAdjustMethod = pAdjustMethod,
                       ont = ont,
                       useInternal = useInternal, 
                       qvalueCutoff = qvalueCutoff)
      })
      names(func_results$ORA) <- onts
    }
  }

  ## GSEA
  if(flags$GSEA && !is.null(genes_gsea)){
    if(flags$clusters){
      if(verbose) message(paste0("Performing clusters GSEA enrichment for (",length(onts),") ontologies"))
      func_results[["WGCNA"]][["GSEA_expanded"]] <- apply_fun(onts, function(ont){
        perform_GSEA_clusters(all_clusters = genes_gsea, 
                              organism = get_organismID_byOnto(organism_info,ont), 
                              keyType = keytypes[[ont]], 
                              pvalueCutoff = pvalueCutoff, 
                              pAdjustMethod = pAdjustMethod, 
                              ont = ont, 
                              useInternal = useInternal)   
      })
      names(func_results$WGCNA$GSEA_expanded) <- onts
    }else{
      if(verbose) message(paste0("Performing GSEA enrichment for (",length(onts),") ontologies"))
      func_results[["GSEA"]] <- apply_fun(onts, function(ont){
        enrich_GSEA(geneList = genes_gsea,
                    organism = get_organismID_byOnto(organism_info,ont),
                    keyType = keytypes[[ont]],
                    pvalueCutoff = pvalueCutoff,
                    pAdjustMethod = pAdjustMethod,
                    ont = ont,
                    useInternal = useInternal)
      })
      names(func_results$GSEA) <- onts
    }
  }else if(flags$GSEA){
    warning("Genes are not in GSEA format, GSEA enrichment will not be performed")
  }

  # Compact
  if(!is.null(func_results$WGCNA)){
    if(!is.null(func_results$WGCNA$ORA_expanded)){
      enrichments_ORA <- lapply(func_results$WGCNA$ORA_expanded, clusterProfiler::merge_result)
      func_results$WGCNA$ORA <- lapply(enrichments_ORA, function(res){
        if(nrow(res) > 0) res <- catched_pairwise_termsim(res)
        return(res)
      })
    }
    if(!is.null(func_results$WGCNA$GSEA_expanded)){
      func_results$WGCNA$GSEA <- lapply(func_results$WGCNA$GSEA_expanded, clusterProfiler::merge_result)
    }
  }

  # End and return
  return(func_results)
}
