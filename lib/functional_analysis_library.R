#################################################################################
############################ FUNCTIONAL ANALYSIS LIBRARY ########################
#################################################################################

obtain_info_from_biomaRt <- function(orthologues, id_type, mart, dataset, host, attr){
    query <- NULL
    if(file.exists("query_results_temp")){
      query <- readRDS("query_results_temp")
    }else {
      query <- getting_information_with_BiomaRt(orthologues, id_type, mart, dataset, host, attr)
    }
    return(query)
}


getting_information_with_BiomaRt <- function(orthologues, id_type, mart, dataset, host, attr){
    ensembl <- useMart(mart,dataset=dataset, host=host)
    filt <- c(opt$biomaRt_filter)
    val <- orthologues
    attr <- attr
    query <- getBM(attributes = attr, filters = filt, values = val, mart = ensembl)
}


perform_GSEA_analysis <- function(interesting_genenames, DEG_annot_table, ontology, graphname){
    geneID2GO <- split(DEG_annot_table$go_accession, DEG_annot_table[,opt$biomaRt_filter])

    geneID2GO <- lapply(geneID2GO, unique)
    geneNames <- names(geneID2GO)
    geneList <- factor(as.integer(geneNames %in% interesting_genenames))
    names(geneList) <- geneNames

    GOdata <- new("topGOdata", ontology =ontology, allGenes = geneList, annot=annFUN.gene2GO, gene2GO = geneID2GO)
    resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher") 

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
    ensembl <- useMart(biomart=mart, dataset=dataset, host=host)
    ensemblorganism <- useDataset(dataset, mart=ensembl)
    genes <- getBM(attributes = c(biomaRt_filter), mart=ensemblorganism)
    genenames <- genes[[1]]
    total_genes <- length(genenames)
    return(total_genes)
}


find_interesting_pathways <- function(biomaRt_organism_info, genes_of_interest){
    print(length(genes_of_interest))
    pathway_pval <- list()
    complete_pathway_df <- NULL
    EC_number<- character(0)
    path_name <- character(0)
    entrez_ID <- character(0)

    pathway_list <- keggList("pathway", as.character(biomaRt_organism_info[,"KeggCode"]))
    pathway_codes <- sub("path:", "", names(pathway_list))

    total_genes <- getting_number_geneIDs(mart = as.character(biomaRt_organism_info[,"Mart"]),as.character(biomaRt_organism_info[,"Dataset"]), organism_host, opt$biomaRt_filter)
    total_pathways <- length(pathway_codes)

    for(i in seq(1, total_pathways, by=10 )){ # We get paths from kegg by packages of 10
        left <- total_pathways - i
        range <- 9
        if(left < range){
            range <- left
        }
        query_genespathway <- keggGet(pathway_codes[c(i:i+range)])

        for(i in 1:10){ #Processing package retrieved from kegg
            genes_in_pathway <- query_genespathway[[1]]$GENE[c(TRUE, FALSE)]
            pVal_for_pathway <- calculate_pvalue(genes_in_pathway, genes_of_interest, total_genes)
            if(pVal_for_pathway < 0.05){                
                path_name <- query_genespathway[[1]]$ENTRY
                EC_list <- keggGet(path_name)[[1]]$GENE
                number_genes <- c(1:length(genes_in_pathway))
                even_number_v <- number_genes[number_genes%%2==0]
                odd_number_v <- number_genes[number_genes%%2!=0]
                real_genes_number <- length(number_genes)/2
                for (i in c(1:real_genes_number)){
                    line_EC <- as.vector(unlist(EC_list[[even_number_v[i]]]))
                    EC_number <- str_match(line_EC, "(EC\\:[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+)")
                    EC_number <- EC_number[1]
                    if (length(EC_number) != 0){
                       entrez_ID <- unlist(EC_list[[odd_number_v[i]]])
                       complete_pathway_df <- rbind(complete_pathway_df, data.frame(pathway = path_name, 
                       entrez = entrez_ID, EC = EC_number, row.names = NULL))
                    } 
                }
            }
        }
    } 
    return(complete_pathway_df)
}



visualizing_KEGG_pathways <- function(name, pathway_id_EC){
kegg_info <- pathway_id_EC 

report = paste('<html>',
        '<head>',
        '<title>Pathway table</title>',
        '<style>',
        'azul {color:rgb(255,0,0);}',
        '</style>',
        '</head>',
        '<body bgcolor="#FFFFFF">',
        '<center>',
        '<table border="2" cellspacing="0" cellpadding="2">',
        '<tr><th>PATHWAY</th></tr>',
        sep="\n")
#table row generating

curret_pathway <- ''
current_entrez <- c()
current_ec <- c()
for (i in 1:nrow(kegg_info)) {  
        record = kegg_info[i,]
        if(record$pathway != curret_pathway ||  i == nrow(kegg_info)){
                # Add row to report table
                if(i > 1){
                        pathway_name <- '' # Habría que definirlo, por ahora no lo uso
                        row_table <- paste('<tr><td><a href="http://www.kegg.jp/kegg-bin/show_pathway?',
                                                curret_pathway,
                                                '/default%3d',
                                                'red',
                                                '/',
                                                paste(unique(current_ec), collapse='%09/,'),
                                                '">',
                                                curret_pathway, # esto cambiar por pathway_name cuando se pueda
                                                '</a></td></tr>',
                                                sep='')

                        report <- paste(report, row_table, sep="\n")
                }
                # Reset pathway to new path
                current_entrez <- c(record$entrez)
                current_ec  <- c(record$EC)
        }else{
                current_entrez <- c(current_entrez, as.character(record$entrez)) #Fields are factor so we use as.character
                current_ec <- c(current_ec, as.character(record$EC))
        }
        curret_pathway <- record$pathway
}


#footer html
report = paste(report, '</table>',
        '</center>',
        '</body>',
        '</html>',
        sep="\n")

write(report, file = file.path(paths, name), sep='')
}

