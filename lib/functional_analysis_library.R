#################################################################################
############################ FUNCTIONAL ANALYSIS LIBRARY ########################
#################################################################################
defining_functional_hunter_subfolders <- function(opt){
  subfolders <- c()
  if (grepl("G", opt$functional_analysis)){
    subfolders <- c(subfolders, 'Results_topGO')
  }
  if (grepl("K", opt$functional_analysis)){
    subfolders <- c(subfolders, 'Results_KEGG')
  }
}


obtain_info_from_biomaRt <- function(orthologues, id_type, mart, dataset, host, attr){
    require(biomaRt)

    ensembl <- useMart(mart,dataset=dataset, host=host)
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
      }else if(!all(container[,1] %in% val)){
        warning("Current query results does not match. Will be overwritten")
        container <- NULL
      }else if(all(val %in% container[,1])){
        message("All IDs are already stored at temporal query results file. Query will not be performed")
        return(container)
      }else if(any(val %in% container[,1])){
          # Filter already calculated 
          val <- val[-which(val %in% container[,1])]
      }
    }else{
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
        }else{
            end <- indx[i+1] - 1
        }
        interval <- seq(indx[i],end)
        message(paste("Fragment: ",i,"/",length(indx),"  (",length(interval),")",sep=""))

        # Run query
        query <- getBM(attributes = attr, filters = filt, values = val[interval], mart = ensembl)
        # Store
        if(is.null(container)){
            container <- query
        }else if(nrow(query) > 0){
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




ensembl_to_entrez <- function(ensembl_ids,organism_db, organism_var){
    # Load necessary package
    require(organism_db, character.only = TRUE)
    # Obtain target variable
    aux_var <- get(organism_var)
    # Translate ENSEMBL to Entrex IDs
    ensembl2entrez <- as.list(aux_var[mappedkeys(aux_var)])
    # Convert to dataframe
    ensembl2entrez_df <- as.data.frame(do.call(rbind,lapply(intersect(ensembl_ids,names(ensembl2entrez)),function(ensembl){
        # Obtain genes
        genes <- ensembl2entrez[[ensembl]]            
 
        if(length(genes) == 0){
            return(data.frame())
        }
        # Return info
        return(data.frame(ENSEMBL = rep(ensembl,length(genes)),ENTREZ = genes, stringsAsFactors = FALSE))
    })))
    return(ensembl2entrez_df)
}




#' @param entrez_targets
#' @param entrez_universe
#' @param sub_ontology
perform_GSEA_analysis_local <- function(entrez_targets,entrez_universe,sub_ontology,outFile=NULL,organism, plot_graph = TRUE){
  
  #! GOFisherTest
  
  require(topGO)
  
  # Prepare necessary info
  allGenes <- unique(entrez_universe) 
  all_results_table <- data.frame() 
  
  # creating contigency vector for gene identifiers
  geneList <- factor(as.integer(allGenes %in% entrez_targets)) 
  names(geneList) <- allGenes
  
  if(length(levels(geneList)) == 2){ # launch analysis
    TopGOobject <- new("topGOdata", ontology = sub_ontology, allGenes = geneList, annot = annFUN.org, mapping = organism)
    # # Possible option 2
    results_fisher  <- runTest(TopGOobject, algorithm = "classic", statistic = "fisher") 
    results_KS      <- runTest(TopGOobject, algorithm = "classic", statistic = "ks")
    results_KS_elim <- runTest(TopGOobject, algorithm = "elim", statistic = "ks")

    # Create results table
    all_results_table <- GenTable(TopGOobject, 
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
      showSigOfNodes(TopGOobject, score(results_fisher), firstSigNodes = 10, useInfo = 'all')
    dev.off()
  }
  
  # Return info
  return(all_results_table)
}






perform_GSEA_analysis <- function(attr_name, interesting_genenames, DEG_annot_table, ontology, graphname,filter_name){
    geneID2GO <- split(DEG_annot_table[,attr_name], DEG_annot_table[,filter_name])
    geneID2GO <- lapply(geneID2GO, unique)
    geneNames <- names(geneID2GO)
    geneList <- factor(as.integer(geneNames %in% interesting_genenames))
    names(geneList) <- geneNames

    GOdata <- new("topGOdata", ontology =ontology, allGenes = geneList, annot=annFUN.gene2GO, gene2GO = geneID2GO)
    resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher") 

    resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
    resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

    allRes <- GenTable(GOdata, classicFisher = resultFis,
        classicKS = resultKS, elimKS = resultKS.elim,
        orderBy = "elimKS", ranksOf = "classicFisher", topNodes = length(GOdata@graph@nodes))

    write.table(allRes, file=file.path(paths$root, "allResGOs.txt"), quote=F, col.names=NA, sep="\t")

    pdf(file.path(paths, graphname), w=11, h=8.5)
        showSigOfNodes(GOdata, score(resultFis), firstSigNodes = 10, useInfo = 'all')
    dev.off()

}



######################### FUNCTIONAL ANALYSIS WITH KEGGREST ##################3

obtain_pathways_gene_pval <- function(raw_filter, functional_parameters){
    pathway_gene <- NULL
    pathway_pval <- NULL
    
    if(file.exists("pathway_gene_temp") & file.exists("pathway_pval_temp")){
      pathway_gene <- readRDS("pathway_gene_temp")
      pathway_pval <- readRDS("pathway_pval_temp")
    }else {
      find_interesting_pathways(raw_filter, functional_parameters)
    }
    return(list(pathway_gene, pathway_pval))
}


calculate_pvalue <- function(genes_in_pathway, genes_of_interest, total_genes) {
    white_balls_drawn <- length(intersect(genes_of_interest, genes_in_pathway))
    white_balls_in_urn <- length(genes_in_pathway)
    total_balls_in_urn <- total_genes
    black_balls_in_urn <- total_balls_in_urn - white_balls_in_urn
    total_balls_drawn_from_urn <- length(genes_of_interest)
    pvalue <-dhyper(
            white_balls_drawn,
            white_balls_in_urn,
            black_balls_in_urn,
            total_balls_drawn_from_urn
    )
    return(pvalue)
}


getting_number_geneIDs <- function(mart, dataset, host, biomaRt_filter){
    require(biomaRt)

    ensembl <- useMart(biomart=mart, dataset=dataset, host=host)
    ensemblorganism <- useDataset(dataset, mart=ensembl)
    genes <- getBM(attributes = c(biomaRt_filter), mart=ensemblorganism)
    genenames <- genes[[1]]
    total_genes <- length(genenames)
    return(total_genes)
}


generate_FA_report <- function(){
  template_path_functional_report <- file.path(main_path_script, 'templates', 'generate_functional_report.tex')
  current_template_path_functional <- file.path(paths$root, 'generate_functional_report.tex')
  latex_file_path_FA <- file.path(paths$root, 'Functional_analysis_report.tex')

  file.copy(template_path_functional_report, paths$root, overwrite = TRUE)
  opts_chunk$set(echo=FALSE)
  knit(current_template_path_functional , output = latex_file_path_FA)

  cmd <- paste('pdflatex', latex_file_path_FA, sep=' ')
  system(cmd)
}



#'
#' @param genes :: 
#' @param organism :: 
#' @param keyType :: 
#' @param pvalueCutoff :: 
#' @param pAdjustMethod :: 
#' @param ont :: ontology to be used. Allowed [GO_MF,GO_CC,GO_BP,KEGG,REACT]
#' @param useInternal :: used only for KEGG enrichment, activate internal data usage mode
#' @return enrichment performed
#' @import clusterProfiler, KEGG.db, ReactomePA
enrichment_ORA <- function(genes,organism,keyType="ENTREZID",pvalueCutoff,pAdjustMethod = "BH",ont,useInternal = FALSE){
  require(clusterProfiler)
  if(useInternal)
    require(KEGG.db)
  require(ReactomePA)
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
    enrichment <- enrichGO(gene          = genes,
                           OrgDb         = organism,
                           keyType       = keyType,
                           ont           = go_subonto,
                           pvalueCutoff  = pvalueCutoff,
                           pAdjustMethod = pAdjustMethod) 
  }else if(ont == "KEGG"){
    enrichment <- enrichKEGG(gene          = genes,
                             organism      = organism,
                             keyType       = keyType,
                             pvalueCutoff  = pvalueCutoff,
                             pAdjustMethod = pAdjustMethod,
                             use_internal_data = useInternal)
  }else if(ont == "REACT"){
    enrichment <- enrichPathway(genes,
                                organism = organism,
                                pAdjustMethod = pAdjustMethod,
                                pvalueCutoff = pvalueCutoff)
  }else{
    stop("Error, ontology specified is not supported to be enriched")
  }

  # Return enrichment
  return(enrichment)
}


#'
#' @param genes :: 
#' @param organism :: 
#' @param keyType :: 
#' @param pvalueCutoff :: 
#' @param pAdjustMethod :: 
#' @param ont :: ontology to be used. Allowed [GO_MF,GO_CC,GO_BP,KEGG,REACT]
#' @param useInternal :: used only for KEGG enrichment, activate internal data usage mode
#' @return enrichment performed
#' @import clusterProfiler, KEGG.db, ReactomePA
enrichment_GSEA <- function(geneList,organism,keyType="ENTREZID",pvalueCutoff,pAdjustMethod = "BH",ont,useInternal = FALSE){
require(clusterProfiler)
  if(useInternal)
    require(KEGG.db)
  require(ReactomePA)
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
    enrichment <- gseGO(geneList      = geneList,
                        OrgDb         = organism,
                        keyType       = keyType,
                        ont           = go_subonto,
                        pvalueCutoff  = pvalueCutoff,
                        pAdjustMethod = pAdjustMethod)
  }else if(ont == "KEGG"){
    enrichment <- gseKEGG(geneList     = geneList,
                          organism     = organism,
                          use_internal_data = useInternal,
                          # nPerm        = 1000,
                          # minGSSize    = 120,
                          pvalueCutoff = pvalueCutoff,
                          verbose      = FALSE)
  }else if(ont == "REACT"){
    enrichment<- gsePathway(geneList, 
                            organism = organism,
                            # exponent = 1, 
                            # nPerm = 1000,
                            # minGSSize = 10, 
                            # maxGSSize = 500, 
                            pvalueCutoff = pvalueCutoff,
                            pAdjustMethod = pAdjustMethod)
  }else{
    stop("Error, ontology specified is not supported to be enriched")
  }

  # Return enrichment
  return(enrichment)
}