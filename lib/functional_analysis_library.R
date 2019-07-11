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




ensembl_to_entrez <- function(ensembl_ids,organism){
    # require(org.Hs.eg.db)
    require(organism, character.only = TRUE)
    # Obtain target variable
    aux_var_name <- paste(paste(head(unlist(strsplit(organism,"\\.")),-1),collapse="."),"ENSEMBL2EG",sep="")
    aux_var <- get(aux_var_name)
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




# human_ENSEMBL_to_GO <- function(ensembl_ids,sub_ontology){
#     require(org.Hs.eg.db)
#     require(topGO)

#     # Translate ENSEMBL to Entrex IDs
#     ensembl2entrez <- as.list(org.Hs.egENSEMBL2EG[mappedkeys(org.Hs.egENSEMBL2EG)])

#     # Obtain related GO terms
#     go2entrez <- annFUN.org(sub_ontology, mapping = "org.Hs.eg.db")
#     entrez2go <- as.data.frame(do.call(rbind,lapply(seq_along(go2entrez),function(i){
#         go_term <- names(go2entrez)[i]
#         genes <- go2entrez[[i]]
#         return(data.frame(Gene = genes, GO = rep(go_term,length(genes)), stringsAsFactors = FALSE))    
#     })))

#     # Revert GO-Entrex to GO-ENSEMBL
#     ensembl2go <- as.data.frame(do.call(rbind,lapply(ensembl_ids,function(ensembl){
#         # Obtain related entrez genes
#         genes <- ensembl2entrez[ensembl]
#         if(length(genes) == 0){
#             return(data.frame())
#         }
#         # Obtain related go
#         info <- as.data.frame(do.call(rbind,lapply(genes,function(gene){
#             go_terms <- entrez2go$GO[which(entrez2go$Gene == gene)]
#             # Return info
#             return(data.frame(Gene = rep(ensembl,length(go_terms)), 
#                               GO = go_terms, 
#                               Entrez = rep(gene,length(go_terms)), stringsAsFactors = FALSE))
#         })))
#         return(info)  
#     })))

#     # Add wanted names
#     colnames(ensembl2go) <- c("ensembl_gene_id", "go_id", "entrezgene")

#     # Return info
#     return(ensembl2go)
# }




# getting_information_with_BiomaRt <- function(orthologues, id_type, mart, dataset, host, attr){
#     require(biomaRt)
#     # Load
#     ensembl <- useMart(mart,dataset=dataset, host=host)
#     query <- getBM(attributes = attr, filters = id_type, values = orthologues, mart = ensembl)
# }





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
    # Prepare metrics
    # classic <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher_Test")
    # KS <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
    # KS_elim <- new("elimScore", testStatistic = GOKSTest, name = "KS tests")
    # # Perform metrics
    # results_fisher  <- getSigGroups(TopGOobject, classic)
    # results_KS      <- getSigGroups(TopGOobject,KS)
    # results_KS_elim <- getSigGroups(TopGOobject,KS_elim)
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


# #### OBSOLETE
# find_interesting_pathways <- function(biomaRt_organism_info, genes_of_interest, filter_name){
#     # Load necessary packages
#     require(KEGGREST)

#     # Prepare containers
#     pathway_pval <- list()
#     complete_pathway_df <- NULL
#     EC_number <- character(0)
#     path_name <- character(0)
#     entrez_ID <- character(0)

#     pathway_list <- keggList("pathway", as.character(biomaRt_organism_info[,"KeggCode"]))
#     pathway_codes <- sub("path:", "", names(pathway_list))

#     total_genes <- getting_number_geneIDs(mart = as.character(biomaRt_organism_info[,"Mart"]),
#                                           dataset = as.character(biomaRt_organism_info[,"Dataset"]), 
#                                           host = organism_host, 
#                                           biomaRt_filter = filter_name)
#     total_pathways <- length(pathway_codes)

#     for(i in c(1:total_pathways)){
#         left <- total_pathways - i
#         range <- 9
#         if(left < range){
#             range <- left
#         }
#         query_genespathway <- keggGet(pathway_codes[c(i:i+range)])

#         for(i in 1:10){ #Processing package retrieved from kegg
#             genes_in_pathway <- query_genespathway[[1]]$GENE[c(TRUE, FALSE)]
#             pVal_for_pathway <- calculate_pvalue(genes_in_pathway, genes_of_interest, total_genes)
#             if(pVal_for_pathway < 0.1){                
#                 path_name <- query_genespathway[[1]]$ENTRY
#                 EC_list <- keggGet(path_name)[[1]]$GENE
#                 if (!is.null(EC_list)){
#                     number_genes <- c(1:length(genes_in_pathway))
#                     even_number_v <- number_genes[number_genes%%2==0]
#                     odd_number_v <- number_genes[number_genes%%2!=0]
#                     real_genes_number <- length(number_genes)/2
#                     if (length(even_number_v != 0)){
#                         for (i in c(1:real_genes_number)){
#                             line_EC <- as.vector(unlist(EC_list[[even_number_v[i]]]))
#                             EC_number <- str_match(line_EC, "(EC\\:[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+)")
#                             EC_number <- EC_number[1]
#                             if (length(EC_number) != 0){
#                                 entrez_ID <- unlist(EC_list[[odd_number_v[i]]])
#                                 complete_pathway_df <- rbind(complete_pathway_df, data.frame(pathway = path_name, 
#                                 entrez = entrez_ID, EC = EC_number, row.names = NULL))
#                             } 
#                         }
#                     }
#                 }
#             }
#         }
#     } 
#     return(complete_pathway_df)
# }
# 
# 
# 
# visualizing_KEGG_pathways <- function(name, pathway_id_EC){
#     kegg_info <- pathway_id_EC 

#     report = paste('<html>',
#             '<head>',
#             '<title>Pathway table</title>',
#             '<style>',
#             'azul {color:rgb(255,0,0);}',
#             '</style>',
#             '</head>',
#             '<body bgcolor="#FFFFFF">',
#             '<center>',
#             '<table border="2" cellspacing="0" cellpadding="2">',
#             '<tr><th>PATHWAY</th></tr>',
#             sep="\n")
#     #table row generating

#     curret_pathway <- ''
#     current_entrez <- c()
#     current_ec <- c()
#     for (i in 1:nrow(kegg_info)) {  
#             record = kegg_info[i,]
#             if(record$pathway != curret_pathway ||  i == nrow(kegg_info)){
#                     # Add row to report table
#                     if(i > 1){
#                             pathway_name <- '' 
#                             row_table <- paste('<tr><td><a href="http://www.kegg.jp/kegg-bin/show_pathway?',
#                                                     curret_pathway,
#                                                     '/default%3d',
#                                                     'red',
#                                                     '/',
#                                                     paste(unique(current_ec), collapse='%09/,'),
#                                                     '">',
#                                                     curret_pathway,
#                                                     '</a></td></tr>',
#                                                     sep='')

#                             report <- paste(report, row_table, sep="\n")
#                     }
#                     # Reset pathway to new path
#                     current_entrez <- c(record$entrez)
#                     current_ec  <- c(record$EC)
#             }else{
#                     current_entrez <- c(current_entrez, as.character(record$entrez)) #Fields are factor so we use as.character
#                     current_ec <- c(current_ec, as.character(record$EC))
#             }
#             curret_pathway <- record$pathway
#     }


#     #footer html
#     report = paste(report, '</table>',
#             '</center>',
#             '</body>',
#             '</html>',
#             sep="\n")

#     write(report, file = file.path(paths, name), sep='')
# }




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

