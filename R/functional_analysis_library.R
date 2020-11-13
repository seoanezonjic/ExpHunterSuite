#################################################################################
############################ FUNCTIONAL ANALYSIS LIBRARY ########################
#################################################################################

#'
#' @param opt
#' @keywords
#' @return
defining_functional_hunter_subfolders <- function(opt){
  subfolders <- c()
  if (grepl("G", opt$functional_analysis)){
    subfolders <- c(subfolders, 'Results_topGO')
  }
  if (grepl("K", opt$functional_analysis)){
    subfolders <- c(subfolders, 'Results_KEGG')
  }
}

#'
#' @param var_name
#' @param level_depth
#' @keywords
#' @return
exists_enrich_df <- function(var_name,level_depth = 3){
  # Check if it is instantiated
  if(!var_name %in% ls(envir = parent.frame(n = level_depth))){
    return(FALSE)
  }
  # Check if is a DF with info
  df <- get(var_name, envir = parent.frame(n = level_depth))
  return(check_results(df))
}


#'
#' @param path
#' @param file_name
#' @param ext
#' @keywords
#' @return
graph_file_name <- function(path, file_name, ext = "pdf"){
  return(paste(path,.Platform$file.sep,file_name,".",ext,sep=""))
}


#'
#' @param path
#' @param file_name
#' @param ext
#' @keywords
#' @return
exists_graph_file <- function(path, file_name, ext = "pdf"){
  f <- graph_file_name(path,file_name,ext)
  return(file.exists(f))
}


#'
#' @param df
#' @keywords
#' @return
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


#'
#' @param ids_to_translate
#' @param annot_table
#' @keywords
#' @return
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


#'
#' @param orthologues
#' @param id_type
#' @param mart
#' @param dataset
#' @param host
#' @param attr
#' @keywords
#' @return
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


#'
#' @param custom_sets
#' @param custom_files
#' @param p_val_threshold
#' @param likely_degs_entrez
#' @param write
#' @keywords
#' @return
enrich_all_customs <- function(custom_sets = NULL, custom_files = NULL, p_val_threshold, likely_degs_entrez, write_res = TRUE){
  if(is.null(custom_gmt)){
    custom_set <- custom_sets
    load_files <- FALSE
  }else{
    custom_set <- custom_files
    names(custom_set) <- custom_files
    load_files <- TRUE
  }
  custom_enrichments <- lapply(custom_set, function(custom_gmt) {
      # Load info
      if(load_files){
        c_terms <- unlist(read.table(file = custom_gmt, sep = "\n", header = FALSE, stringsAsFactors = FALSE))
      }else{
        c_terms <- custom_gmt
      }
      # Split all
      c_terms <- as.data.frame(do.call(rbind,lapply(c_terms,function(GeneTerms) {
        aux <- unlist(strsplit(GeneTerms,"\t"))
        return(data.frame(Term = rep(aux[1],length(aux)-2),
                  Gene = tail(aux,-2),
                  stringsAsFactors = FALSE))
      })))
      enr <- clusterProfiler::enricher(likely_degs_entrez, pvalueCutoff = p_val_threshold, TERM2GENE = c_terms)
      # # Store results
      # TODO: using external variables which has not been given by parameters? Uh someone has been a bad boy. Please Pepe/Jim be careful with that
      # TODO: stay until targets_functional script is refactored. then REMOVE WRITE FLAG 
      if(write_res) write.table(enr, file=file.path(paths$root, paste0(basename(names(custom_gmt)[1]),"_ora_results")), quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")      
      # Return
      return(list(File = names(custom_gmt)[1],
            Result = enr))
    })
  return(custom_enrichments)
}


#'
#' @param input_ids
#' @param organism_db
#' @param org_var_name
#' @keywords
#' @return
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



#'
#' @param entrez_targets
#' @param entrez_universe
#' @param sub_ontology
#' @keywords
#' @return
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
    write.table(all_results_table, file=outFile, quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")    
  }

  # Plot graphs
  if(plot_graph){
    pdf(paste(outFile,"pdf",sep="."), w=11, h=8.5)
      topGO::showSigOfNodes(TopGOobject, BiocGenerics::score(results_fisher), firstSigNodes = 10, useInfo = 'all')
    dev.off()
  }
  
  # Return info
  return(all_results_table)
}


#'
#' @param data_matrix
#' @param transpose
#' @keywords
#' @return
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




#'
#' @param attr_name
#' @param interesting_genenames
#' @param DEG_annot_table
#' @param ontology
#' @param graphname
#' @param filter_name
#' @param output_files
#' @keywords
#' @return
perform_topGO <- function(attr_name, interesting_genenames, DEG_annot_table, ontology, graphname,filter_name, output_files){
    geneID2GO <- split(DEG_annot_table[,attr_name], DEG_annot_table[,filter_name])
    geneID2GO <- lapply(geneID2GO, unique)
    geneNames <- names(geneID2GO)
    geneList <- factor(as.integer(geneNames %in% interesting_genenames))
    names(geneList) <- geneNames

    GOdata <- topGO::new("topGOdata", ontology =ontology, allGenes = geneList, annot=annFUN.gene2GO, gene2GO = geneID2GO)
    resultFis <- topGO::runTest(GOdata, algorithm = "classic", statistic = "fisher") 

    resultKS <- topGO::runTest(GOdata, algorithm = "classic", statistic = "ks")
    resultKS.elim <- topGO::runTest(GOdata, algorithm = "elim", statistic = "ks")

    allRes <- topGO::GenTable(GOdata, classicFisher = resultFis,
        classicKS = resultKS, elimKS = resultKS.elim,
        orderBy = "elimKS", ranksOf = "classicFisher", topNodes = length(GOdata@graph@nodes))

    write.table(allRes, file=file.path(output_files, "allResGOs.txt"), quote=FALSE, col.names=NA, sep="\t")

    pdf(file.path(paths, graphname), w=11, h=8.5)
        topGO::showSigOfNodes(GOdata, BiocGenerics::score(resultFis), firstSigNodes = 10, useInfo = 'all')
    dev.off()

}



######################### FUNCTIONAL ANALYSIS WITH KEGGREST ##################3
# obtain_pathways_gene_pval <- function(raw_filter, functional_parameters){
#     pathway_gene <- NULL
#     pathway_pval <- NULL
    
#     if(file.exists("pathway_gene_temp") & file.exists("pathway_pval_temp")){
#       pathway_gene <- readRDS("pathway_gene_temp")
#       pathway_pval <- readRDS("pathway_pval_temp")
#     } else {
#       find_interesting_pathways(raw_filter, functional_parameters)
#     }
#     return(list(pathway_gene, pathway_pval))
# }


# calculate_pvalue <- function(genes_in_pathway, genes_of_interest, total_genes) {
#     white_balls_drawn <- length(intersect(genes_of_interest, genes_in_pathway))
#     white_balls_in_urn <- length(genes_in_pathway)
#     total_balls_in_urn <- total_genes
#     black_balls_in_urn <- total_balls_in_urn - white_balls_in_urn
#     total_balls_drawn_from_urn <- length(genes_of_interest)
#     pvalue <-dhyper(
#             white_balls_drawn,
#             white_balls_in_urn,
#             black_balls_in_urn,
#             total_balls_drawn_from_urn
#     )
#     return(pvalue)
# }


# getting_number_geneIDs <- function(mart, dataset, host, biomaRt_filter){
#     require(biomaRt)

#     ensembl <- useMart(biomart=mart, dataset=dataset, host=host)
#     ensemblorganism <- useDataset(dataset, mart=ensembl)
#     genes <- getBM(attributes = c(biomaRt_filter), mart=ensemblorganism)
#     genenames <- genes[[1]]
#     total_genes <- length(genenames)
#     return(total_genes)
# }

# generate_FA_report <- function(){
#   template_path_functional_report <- file.path(main_path_script, 'templates', 'generate_functional_report.tex')
#   current_template_path_functional <- file.path(paths$root, 'generate_functional_report.tex')
#   latex_file_path_FA <- file.path(paths$root, 'Functional_analysis_report.tex')

#   file.copy(template_path_functional_report, paths$root, overwrite = TRUE)
#   opts_chunk$set(echo=FALSE)
#   knit(current_template_path_functional , output = latex_file_path_FA)

#   cmd <- paste('pdflatex', latex_file_path_FA, sep=' ')
#   system(cmd)
# }

######################### FUNCTIONAL ANALYSIS WITH KEGGREST ##################3




#'
#' @param genes 
#' @param organism 
#' @param keyType 
#' @param pvalueCutoff 
#' @param pAdjustMethod 
#' @param ont ontology to be used. Allowed [GO_MF,GO_CC,GO_BP,KEGG,REACT]
#' @param useInternal used only for KEGG enrichment, activate internal data usage mode
#' @param qvalueCutoff 
#' @param ENRICH_DATA 
#' @keywords
#' @return enrichment performed
enrichment_ORA <- function(genes,organism,keyType="ENTREZID",pvalueCutoff,pAdjustMethod = "BH",ont,useInternal = FALSE, qvalueCutoff = 0.2,ENRICH_DATA = NULL){
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

  # Return enrichment
  return(enrichment)
}


#'
#' @param gmt
#' @param genes_in_modules
#' @param pthreshold
#' @param cores
#' @keywords
#' @return
enrich_clusters_with_gmt <- function(gmt, genes_in_modules, pthreshold,cores){
      # Enrich
      modules_enrichment <- parallel::mclapply(genes_in_modules, function(genesset) {
        clusterProfiler::enricher(genesset, pvalueCutoff = pthreshold, TERM2GENE = gmt)
      },mc.cores = cores)
      names(modules_enrichment) <- names(genes_in_modules)
      # Return
      return(modules_enrichment)
}

#'
#' @param gmt_file
#' @keywords
#' @return
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


#'
#' @param genes
#' @param organism  
#' @param keyType  
#' @param pvalueCutoff  
#' @param pAdjustMethod  
#' @param ont ontology to be used. Allowed [GO_MF,GO_CC,GO_BP,KEGG,REACT]
#' @param useInternal used only for KEGG enrichment, activate internal data usage mode
#' @return enrichment performed
enrich_GSEA <- function(geneList,organism,keyType="ENTREZID",pvalueCutoff,pAdjustMethod = "BH",ont,useInternal = FALSE){
  #' @import clusterProfiler KEGG.db ReactomePA
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
    enrichment <- clusterProfiler::gseGO(geneList      = geneList,
                                        OrgDb         = organism,
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


#'
#' @param all_clusters
#' @param organism
#' @param keyType
#' @param pvalueCutoff
#' @param pAdjustMethod
#' @param ont
#' @param useInternal
#' @keywords
#' @return
perform_GSEA_clusters <- function(all_clusters, organism, keyType, pvalueCutoff, pAdjustMethod = "BH", ont, useInternal){
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





#'
#' @param genes 
#' @param organism 
#' @param keyType 
#' @param pvalueCutoff 
#' @param pAdjustMethod 
#' @param ont ontology to be used. Allowed [GO_MF,GO_CC,GO_BP,KEGG,REACT]
#' @param useInternal used only for KEGG enrichment, activate internal data usage mode
#' @param qvalueCutoff 
#' @param ENRICH_DATA 
#' @param mc.cores 
#' @param mc.preschedule 
#' @keywords
#' @return enrichment performed
enrichment_clusters_ORA <- function(genes,organism,keyType="ENTREZID",pvalueCutoff,pAdjustMethod = "BH",ont,useInternal = FALSE, qvalueCutoff, ENRICH_DATA = NULL, mc.cores = 1, mc.preschedule
 = TRUE){
  #' @import clusterProfiler KEGG.db ReactomePA parallel
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
      ENRICH_DATA <- clusterProfiler:::get_GO_data(organism, go_subonto, keyType)
    } else if(ont == "KEGG"){
      # require(clusterProfiler)
      ENRICH_DATA <- clusterProfiler:::get_data_from_KEGG_db(clusterProfiler:::organismMapper(organism))
    } else if(ont == "REACT"){
      # require(ReactomePA)
      ENRICH_DATA <- ReactomePA:::get_Reactome_DATA(organism)
    }
  }

  # require(parallel)
  # Enrich per gene_set
  enrichment <- parallel::mclapply(genes,function(setg){
    # Obtain enrichment 
    enr <- enrichment_ORA(genes         = setg,
                          organism      = organism,
                          keyType       = keyType,
                          pvalueCutoff  = pvalueCutoff,
                          pAdjustMethod = pAdjustMethod,
                          ont           = src_ont,
                          useInternal   = useInternal, 
                          qvalueCutoff  = qvalueCutoff,
                          ENRICH_DATA   = ENRICH_DATA)
    # Return
    return(enr)
  }, mc.cores = mc.cores, mc.preschedule = mc.preschedule)

  names(enrichment) <- names(genes)

  # Return enrichment
  return(enrichment)
}
