#!/usr/bin/env Rscript

# options(warn=1)
# if( Sys.getenv('DEGHUNTER_MODE') == 'DEVELOPMENT' ){
#   # Obtain this script directory
#   full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile), 
#                  error=function(e) # works when using R CMD
#                 normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', 
#                   commandArgs())], '='))[2]))
#   main_path_script <- dirname(full.fpath)
#   root_path <- file.path(main_path_script, '..', '..')
#   # Load custom libraries
#   custom_libraries <- c('io_handling.R', 
#     'functional_analysis_library.R')
#   for (lib in custom_libraries){
#     source(file.path(root_path, 'R', lib))
#   }
#   template_folder <- file.path(root_path, 'inst', 'templates')
# }else{
#   require('ExpHunterSuite')
#   root_path <- find.package('ExpHunterSuite')
#   template_folder <- file.path(root_path, 'templates')
# }


col <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
         "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
         "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
         "#8A7C64", "#599861")


convert_ids_to_entrez <- function(ids, gene_keytype){
  possible_ids <- columns(org.Hs.eg.db)
  if(! gene_keytype %in% possible_ids) 
    stop(paste(c("gene keytype must be one of the following:", possible_ids), collapse=" "))
  ids <- tryCatch(
    ids <- mapIds(org.Hs.eg.db, keys=ids, column="ENTREZID", keytype=gene_keytype),
    error=function(cond){
      ids <- NULL
    }
  )
  return(ids[!is.na(ids)])
}

multienricher <- function(funsys, cluster_genes_list, task_size, workers, pvalcutoff = 0.05, qvalcutoff = 0.02, readable = TRUE){
  enrichments_ORA <- list()
  for(funsys in all_funsys) {
    ENRICH_DATA <- clusterProfiler:::get_GO_data(org.Hs.eg.db::org.Hs.eg.db, funsys, "ENTREZID")
    enrf <- clusterProfiler::enrichGO
    patter_to_remove <- "GO_DATA *<-"
    ltorem <- grep(patter_to_remove, body(enrf))
    body(enrf)[[ltorem]] <- substitute(GO_DATA <- ENRICH_DATA)
    enriched_cats <- ExpHunterSuite:::parallel_list(cluster_genes_list, function(cl_genes){
       enr <- enrf(gene = cl_genes,
                   OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                   pAdjustMethod = "BH", ont = funsys,
                   pvalueCutoff  = pvalcutoff, qvalueCutoff = qvalcutoff, readable = readable)
      }, 
      workers= workers, task_size = task_size
    )
    enrichments_ORA[[funsys]] <- enriched_cats
  }
  return(enrichments_ORA)
}

parse_cluster_results <- function(enrichments_ORA, simplify_results = TRUE){
  enrichments_ORA_tr <- list()
  for (funsys in names(enrichments_ORA)){
    enr_obj <- clusterProfiler::merge_result(enrichments_ORA[[funsys]])
    if(nrow(enr_obj@compareClusterResult) > 0){
      if (funsys %in% c("MF", "CC", "BP") && simplify_results){
        enr_obj@fun <- "enrichGO"
        # enr_obj <- clusterProfiler::simplify(enr_obj) 
      } 
      enr_obj <- ExpHunterSuite:::catched_pairwise_termsim(enr_obj, 200)
    }                              
    enrichments_ORA_tr[[funsys]] <- enr_obj 
  }
  return(enrichments_ORA_tr)
}

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

option_list <- list(
  optparse::make_option(c("-i", "--input_file"), type="character", default=NULL,
                        help="2 columns - cluster and comma separated gene ids"),
  optparse::make_option(c("-w", "--workers"), type="integer", default=1,
                        help="number of processes for parallel execution"),
  optparse::make_option(c("-p", "--pvalcutoff"), type="double", default=0.05,
                        help="Cutoff for P value and adjusted P value for enrichments"),
  optparse::make_option(c("-q", "--qvalcutoff"), type="double", default=0.02,
                        help="Cutoff for Q value for enrichments"),
  optparse::make_option(c("-t", "--task_size"), type="integer", default=1,
                        help="number of clusters per task"),
  optparse::make_option(c("-F", "--force"), type="logical", default=FALSE, 
                        action = "store_true", help="Ignore temporal files"),
  optparse::make_option(c("-k", "--gene_keytype"), type="character", default="ENTREZID",
                        help="What identifier is being used for the genes in the clusters?"),
  optparse::make_option(c("-g", "--gene_mapping"), type="character", default=NULL,
                        help="3 columns tabular file- Cluster - InputGeneID - NumericGeneMapping. Header must be indicated as cluster - geneid - [numeric_mapping]"),
  optparse::make_option(c("-o", "--output_file"), type="character", default="results",
                        help="Define the output path.")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

##################################### INITIALIZE ##
library(optparse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
output_path <- paste0(opt$output_file, "_functional_enrichment")
dir.create(output_path)
output_path <- normalizePath(output_path)
temp_file <- file.path(output_path, "enr_tmp.RData")
all_funsys <- c("MF", "CC", "BP") 
n_category <- 30
#################################### MAIN ##

if (!is.null(opt$gene_mapping)){
  cl_genes_mapping <- read.table(opt$gene_mapping, header=TRUE)
  gane_mapping_name <- colnames(cl_genes_mapping)[3]
  ids <- select(org.Hs.eg.db, keys=cl_genes_mapping[,2], column="ENTREZID", keytype=opt$gene_keytype)
  ids <- unique(ids)
  cl_genes_mapping <- merge(cl_genes_mapping, ids, by.x = "geneid", by.y = opt$gene_keytype)
  cl_genes_mapping <- cl_genes_mapping[!is.na(cl_genes_mapping$ENTREZID),]
  gene_mapping <- cl_genes_mapping[,2:4]
  gene_mapping <- split(gene_mapping, gene_mapping$cluster)
  gene_mapping <-lapply(gene_mapping, function(cluster_mapping){
    gene_map <- cluster_mapping[,gane_mapping_name]
    names(gene_map) <- cluster_mapping$ENTREZID
    return(gene_map)
  })

}



if (!file.exists(temp_file) || opt$force) {

  cluster_genes <- read.table(opt$input_file, header=FALSE)
  cluster_genes_list <- strsplit(cluster_genes[,2], ",")
  if(opt$gene_keytype != "ENTREZID") {
    cluster_genes_list <- sapply(cluster_genes_list, function(x){
                                      convert_ids_to_entrez(ids=x, 
                                                            gene_keytype=opt$gene_keytype)}) 
  }
  names(cluster_genes_list) <- cluster_genes[,1]

  enrichments_ORA <- multienricher(funsys =  all_funsys, 
                                  cluster_genes_list =  cluster_genes_list, 
                                  task_size = opt$task_size, 
                                  workers = opt$workers, 
                                  pvalcutoff =  opt$pvalcutoff, 
                                  qvalcutoff = opt$qvalcutoff)
  save(enrichments_ORA, file = temp_file)

} else {
  load(temp_file)
}

for (funsys in names(enrichments_ORA)){
print(funsys)
  for (cluster in names(enrichments_ORA[[funsys]])){
    cl_path <- file.path(output_path, paste0(cluster,"_cl_enr"))
    if (!dir.exists(cl_path)){
      dir.create(cl_path, recursive = TRUE)
    } 
    enr <-enrichments_ORA[[funsys]][[cluster]]
    print(cluster)
    # save(enr, file = "test.RData")
    if (length(enr$Description) < 3 ) next
    enr <- ExpHunterSuite:::catched_pairwise_termsim(enr, 200)
      pp <- enrichplot::emapplot(enr, showCategory= n_category, layout = "kk")    
    ggplot2::ggsave(filename = file.path(cl_path,paste0(cluster,"_emaplot_",funsys,"_",opt$output_file,".png")), pp, width = 30, height = 30, dpi = 300, units = "cm", device='png')
    if (!is.null(opt$gene_mapping)){
  print(str(gene_mapping[[cluster]]))

      pp <- enrichplot::cnetplot(enr, showCategory= n_category, foldChange = gene_mapping[[cluster]]) + ggplot2::scale_colour_gradient(name = gane_mapping_name, low = "#AFCAFF",breaks=unique(gene_mapping[[cluster]]),high = "#00359C") 
    } else {
      pp <- enrichplot::cnetplot(enr, showCategory= n_category)
    }
    ggplot2::ggsave(filename = file.path(cl_path,paste0(cluster,"_cnetplot_",funsys,"_",opt$output_file,".png")), pp, width = 30, height = 30, dpi = 300, units = "cm", device='png')

  }
}
    q()   
write_fun_enrichments(enrichments_ORA, output_path, all_funsys)
enrichments_ORA <- parse_cluster_results(enrichments_ORA, simplify_results = TRUE)

for (funsys in names(enrichments_ORA)){
  if (length(unique(enrichments_ORA[[funsys]]@compareClusterResult$Description)) < 2 ) next

  pp <- enrichplot::emapplot(enrichments_ORA[[funsys]], showCategory= n_category, pie="Count", layout = "kk")# + ggplot2::scale_fill_manual(values = col)

  ggplot2::ggsave(filename = file.path(output_path,paste0("emaplot_",funsys,"_",opt$output_file,".png")), pp, width = 30, height = 30, dpi = 300, units = "cm", device='png')

  pp <- enrichplot::dotplot(enrichments_ORA[[funsys]], showCategory= n_category)
  ggplot2::ggsave(filename = file.path(output_path,paste0("dotplot_",funsys,"_",opt$output_file,".png")), pp, width = 60, height = 40, dpi = 300, units = "cm", device='png')

}
   

# organisms_table <- ExpHunterSuite::get_organism_table()
# current_organism_info <- subset(organisms_table, rownames(organisms_table) == "Human")
# system.time(
# enriched_cats4 <- ExpHunterSuite::multienricher(cluster_genes_list, organism_info = current_organism_info
# , ontology="m", pvalueCutoff = 0.05, pAdjustMethod = "BH")
# )
