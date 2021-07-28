#!/usr/bin/env Rscript

options(warn=1)
if( Sys.getenv('DEGHUNTER_MODE') == 'DEVELOPMENT' ){
  # Obtain this script directory
  full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile), 
                 error=function(e) # works when using R CMD
                normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', 
                  commandArgs())], '='))[2]))
  main_path_script <- dirname(full.fpath)
  root_path <- file.path(main_path_script, '..', '..')
  # Load custom libraries
  custom_libraries <- c('general_functions.R', 
    'functional_analysis_library.R','plotting_functions.R')
  for (lib in custom_libraries){
    source(file.path(root_path, 'R', lib))
  }
  template_folder <- file.path(root_path, 'inst', 'templates')
  organisms_table_file <- file.path(root_path, "inst", "external_data", 
        "organism_table.txt")
}else{
  require('ExpHunterSuite')
  root_path <- find.package('ExpHunterSuite')
  template_folder <- file.path(root_path, 'templates')
  organisms_table_file <- file.path(root_path, "external_data", 
        "organism_table.txt")
}

col <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
         "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
         "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
         "#8A7C64", "#599861")


convert_ids_to_entrez <- function(ids, gene_keytype, org_db){
  possible_ids <- AnnotationDbi::columns(org_db)
  if(! gene_keytype %in% possible_ids) 
    stop(paste(c("gene keytype must be one of the following:", possible_ids), collapse=" "))
  ids <- tryCatch(
    ids <- AnnotationDbi::mapIds(org_db, keys=ids, column="ENTREZID", keytype=gene_keytype),
    error=function(cond){
      ids <- NULL
    }
  )
  return(ids[!is.na(ids)])
}


get_organism_id <- function(organism_info, funsys){
  if(funsys %in% c("BP","CC","MF")){
    org_id <- organism_info$Bioconductor_DB[1]
  }else if(funsys == "KEGG"){
    org_id <- organism_info$KeggCode[1]
  }else if(funsys == "REACT"){
    org_id <- organism_info$Reactome_ID[1]
  }
  return(org_id)
}



get_org_db <- function(current_organism_info){
  org_db <- current_organism_info$Bioconductor_DB[1]
  org_db <- eval(parse(text = paste0(org_db,"::",org_db)))
  return(org_db)
}


multienricher_2 <- function(funsys, cluster_genes_list, organism_info,org_db = NULL,task_size, workers, pvalcutoff = 0.05, qvalcutoff = 0.02, readable = TRUE){
  enrichments_ORA <- list()
  if (is.null(org_db)){
    org_db <- get_org_db(current_organism_info)
  }

  for(funsys in all_funsys) {
    if (funsys %in% c("CC","BP","MF")){
      enrf <- clusterProfiler::enrichGO
      get_enr_data <- get("get_GO_data", envir = asNamespace("clusterProfiler"),inherits = FALSE)
      pattern_to_remove <- "GO_DATA *<-"
      ENRICH_DATA <- get_enr_data(org_db, funsys, "ENTREZID")
      ltorem <- grep(pattern_to_remove, body(enrf))
      body(enrf)[[ltorem]] <- substitute(GO_DATA <- ENRICH_DATA)

    } else  if (funsys == "REACT"){
    
    }else if (funsys == "KEGG"){

    }
    params <- list(OrgDb = org_db,
                   pAdjustMethod = "BH",minGSSize = 5, ont = funsys,
                   pvalueCutoff  = pvalcutoff, qvalueCutoff = qvalcutoff, readable = readable) #Now Param is configured Only for GO
    enriched_cats <- parallel_list(cluster_genes_list, function(cl_genes){

      print(str(cl_genes))
       params_cl <- c(params, list(gene = cl_genes) )
       enr <- do.call("enrf", params_cl)
      }, 
      workers= workers, task_size = task_size
    )
    enrichments_ORA[[funsys]] <- enriched_cats
  }
  return(enrichments_ORA)
}

parse_cluster_results <- function(enrichments_ORA, simplify_results = TRUE, clean_parentals = FALSE){
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


clean_all_parentals <- function(enr_obj, subont){
  GO_ancestors <-  utils::getAnywhere(paste0("GO",subont,"ANCESTOR"))
  GO_ancestors <- as.list(GO_ancestors$objs[[1]])
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
  for(pair in seq(nrow(hamming_0))){
    possible_parentals <- unlist(hamming_0[pair, c(1,2)])
    parental <- are_parentals(possible_parentals[1], possible_parentals[2], GO_ancestors)
    if (!is.na(parental)){
      terms_to_discard  <- c(terms_to_discard,parental)
    }
  }
  enr_obj@compareClusterResult <- enrich_obj[
             enrich_obj$ID %in% terms_to_discard,]
  return(enr_obj)
}

are_parentals <- function(term_A, term_B, GO_ancestors){
  parental <- NA
   if (term_A %in% GO_ancestors[[term_B]]){
     parental <- term_A
  } else if (term_B %in% GO_ancestors[[term_A]]) {
    parental <- term_B
  }
  return(parental)
}

hamming_binary <- function(X, Y = NULL) {
    if (is.null(Y)) {
        D <- t(1 - X) %*% X
        D + t(D)
    } else {
        t(1 - X) %*% Y + t(X) %*% (1 - Y)
    }
} #taken from https://johanndejong.wordpress.com/2015/10/02/faster-hamming-distance-in-r-2/

write_fun_enrichments <- function(enrichments_ORA, output_path, all_funsys){
  for(funsys in all_funsys) {
    enriched_cats <- enrichments_ORA[[funsys]]
    enriched_cats_dfs <- lapply(enriched_cats, data.frame)
    enriched_cats_bound <- data.table::rbindlist(enriched_cats_dfs, use.names= TRUE, idcol= "Cluster_ID" )
    if (nrow(enriched_cats_bound) == 0) next 
    utils::write.table(enriched_cats_bound, 
                       file=file.path(output_path, paste0("enrichment_",funsys,".csv")),
                       quote=FALSE, col.names=TRUE, row.names = FALSE, sep="\t")
  }
}

parse_results_for_report <- function(enrichments, simplify_results = FALSE){
  enrichments_for_reports <- list()
  for (funsys in names(enrichments)){
    for (cluster in names(enrichments[[funsys]])){
      if (is.null(enrichments_for_reports[[cluster]]))
      enrichments_for_reports[[cluster]] <- list()
      enr <- enrichments[[funsys]][[cluster]]
      if (simplify_results)
        enr <- clusterProfiler::simplify(enr) 
      if (length(enr$Description) < 3 ) next
      enr <- catched_pairwise_termsim(enr, 200)
      enrichments_for_reports[[cluster]][[funsys]] <- enr 
    }
  }
  return(enrichments_for_reports)
}


parse_mappings <- function(gene_mapping, org_db, keytype){
  cl_genes_mapping <- read.table(gene_mapping, header=TRUE)
  gane_mapping_name <- colnames(cl_genes_mapping)[3]
  ids <- AnnotationDbi::select(org_db , keys=cl_genes_mapping[,2], column="ENTREZID", keytype=keytype)
  ids <- unique(ids)
  cl_genes_mapping <- merge(cl_genes_mapping, ids, by.x = "geneid", by.y = keytype)
  cl_genes_mapping <- cl_genes_mapping[!is.na(cl_genes_mapping$ENTREZID),]
  gene_mapping <- cl_genes_mapping[,2:4]
  gene_mapping <- split(gene_mapping, gene_mapping$cluster)
  gene_mapping <-lapply(gene_mapping, function(cluster_mapping){
    gene_map <- cluster_mapping[,gane_mapping_name]
    names(gene_map) <- cluster_mapping$ENTREZID
    return(gene_map)
  })
  return(list(gane_mapping_name = gane_mapping_name, gene_mapping = gene_mapping))
}

write_func_cluster_report <- function(enrichments_for_reports, output_path, gene_mapping){
  temp_path <- file.path(output_path, "temp")
  dir.create(temp_path) #perform parallel
  for (cluster in names(enrichments_for_reports)){
      # cluster_enrichments <- enrichments_for_reports[[cluster]]
      geneList <- NULL
      if (!is.null(gene_mapping)){ 
        geneList <- gene_mapping[[cluster]]
      }
      outfile <- file.path(output_path, paste0(cluster, "_func_report.html"))
      temp_path_cl <- file.path(temp_path, paste0(cluster,"_temp"))
      rmarkdown::render(file.path(template_folder, 
                  'clusters_to_enrichment.Rmd'),output_file = outfile, 
              intermediates_dir = temp_path_cl)
  }
  unlink(temp_path, recursive = TRUE)
}


option_list <- list(
  optparse::make_option(c("-i", "--input_file"), type="character", default=NULL,
                        help="2 columns - cluster and comma separated gene ids"),
  optparse::make_option(c("-w", "--workers"), type="integer", default=1,
                        help="number of processes for parallel execution. Default=%default"),
  optparse::make_option(c("-p", "--pvalcutoff"), type="double", default=0.05,
                        help="Cutoff for P value and adjusted P value for enrichments. Default=%default"),
  optparse::make_option(c("-q", "--qvalcutoff"), type="double", default=0.02,
                        help="Cutoff for Q value for enrichments. Default=%default"),
  optparse::make_option(c("-t", "--task_size"), type="integer", default=1,
                        help="number of clusters per task. Default=%default"),
  optparse::make_option(c("-F", "--force"), type="logical", default=FALSE, 
                        action = "store_true", help="Ignore temporal files"),
  optparse::make_option(c("-f", "--funsys"), type="character", default="BP,MF,CC", 
                        help="Funsys to execute: MF => GO Molecular Function, BP => GO Biological Process, CC => GO Celular Coponent. Default=%default"),
  optparse::make_option(c("-c", "--clean_parentals"), type="logical", default=FALSE, 
                        action = "store_true", help="Clean parentals GO terms that appears on the same clusters than child."),
    optparse::make_option(c("-O", "--model_organism"), type="character", default="Human", 
                        help="Model organism. Human or Mouse"),
  optparse::make_option(c("-d", "--description_file"), type="character", default=NULL,
                        help="Markdown file describing of the enriched clusters."),
  optparse::make_option(c("-k", "--gene_keytype"), type="character", default="ENTREZID",
                        help="What identifier is being used for the genes in the clusters?. Default=%default"),
  optparse::make_option(c("-g", "--gene_mapping"), type="character", default=NULL,
                        help="3 columns tabular file- Cluster - InputGeneID - NumericGeneMapping. Header must be indicated as cluster - geneid - [numeric_mapping]"),
  optparse::make_option(c("-o", "--output_file"), type="character", default="results",
                        help="Define the output path.")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

##################################### INITIALIZE ##
library(clusterProfiler)
output_path <- paste0(opt$output_file, "_functional_enrichment")
dir.create(output_path)
output_path <- normalizePath(output_path)
temp_file <- file.path(output_path, "enr_tmp.RData")
all_funsys <- c("MF", "CC", "BP") 
n_category <- 30
organisms_table <- get_organism_table(organisms_table_file)
current_organism_info <- organisms_table[rownames(organisms_table) %in% opt$model_organism,]
gene_mapping <- NULL
org_db <- get_org_db(current_organism_info)
#################################### MAIN ##


if (!is.null(opt$gene_mapping)){
  parsed_mappings <- parse_mappings(opt$gene_mapping, org_db, opt$gene_keytype)
  gene_mapping <- parsed_mappings[["gene_mapping"]]
  gane_mapping_name <-  parsed_mappings[["gane_mapping_name"]]
}

if (!file.exists(temp_file) || opt$force) {

  cluster_genes <- read.table(opt$input_file, header=FALSE)
  cluster_genes_list <- strsplit(cluster_genes[,2], ",")
  print(str(cluster_genes_list))
  if(opt$gene_keytype != "ENTREZID") {
    cluster_genes_list <- lapply(cluster_genes_list, function(x){
                                      convert_ids_to_entrez(ids=x, 
                                                            gene_keytype=opt$gene_keytype,
                                                            org_db = org_db)}) 
  }
    print(str(cluster_genes_list))

  names(cluster_genes_list) <- cluster_genes[,1]
  print(str(cluster_genes_list))

  enrichments_ORA <- multienricher_2(funsys =  all_funsys, 
                                  cluster_genes_list =  cluster_genes_list, 
                                  task_size = opt$task_size, 
                                  org_db= org_db,
                                  workers = opt$workers, 
                                  pvalcutoff =  opt$pvalcutoff, 
                                  qvalcutoff = opt$qvalcutoff)
  save(enrichments_ORA, file = temp_file)

} else {
  load(temp_file)
}

enrichments_for_reports <- parse_results_for_report(enrichments_ORA)
write_fun_enrichments(enrichments_ORA, output_path, all_funsys)
enrichments_ORA_merged <- parse_cluster_results(enrichments_ORA, simplify_results = TRUE, clean_parentals = opt$clean_parentals)

for (funsys in names(enrichments_ORA_merged)){
  if (length(unique(enrichments_ORA_merged[[funsys]]@compareClusterResult$Description)) < 2 ) next

  pp <- enrichplot::emapplot(enrichments_ORA_merged[[funsys]], showCategory= n_category, pie="Count", layout = "kk")# + ggplot2::scale_fill_manual(values = col)

  ggplot2::ggsave(filename = file.path(output_path,paste0("emaplot_",funsys,"_",opt$output_file,".png")), pp, width = 30, height = 30, dpi = 300, units = "cm", device='png')

  pp <- enrichplot::dotplot(enrichments_ORA_merged[[funsys]], showCategory= n_category)
  ggplot2::ggsave(filename = file.path(output_path,paste0("dotplot_",funsys,"_",opt$output_file,".png")), pp, width = 60, height = 40, dpi = 300, units = "cm", device='png')

}
write_func_cluster_report(enrichments_for_reports,output_path,gene_mapping)

