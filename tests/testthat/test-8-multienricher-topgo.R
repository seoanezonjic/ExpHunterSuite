# test_that("Multienricher gsea gives same results as running packets separately", {

#   funsys <- "BP"
#   org_db <- "org.Mm.eg.db"
#   gene_id <- "entrez"

#   entrez_to_enrich <- c("13853", "16819", "12266", "497097", "19062")
#   entrez_to_enrich_kegg <- c("16176", "18049", "18211", "193034", "277328")

#   organisms_table <- get_organism_table()
#   current_organism_info <- subset(organisms_table, 
#                           rownames(organisms_table) == "Mouse")

# # start_time <- proc.time()
#   new_res_tg <- ExpHunterSuite::multienricher_topGO(all_funsys=funsys, 
#     genes_list=list(entrez_to_enrich_kegg, entrez_to_enrich),
#     organism_info=current_organism_info)
# # end_time <- proc.time()
# # pack_time_taken <- end_time - start_time
# # cat(pack_time_taken, file="~/new_topgo_time_taken_5.txt")

#   library(topGO)
#   funsys <- "BP"
#   org_db <- "org.Mm.eg.db"
#   gene_id <- "entrez"

#   go_to_genes <- topGO::annFUN.org(funsys, mapping = org_db, ID = gene_id)

#   universe <- unique(unlist(go_to_genes))
#   geneList <- factor(as.integer(universe %in% entrez_to_enrich_kegg))
#   names(geneList) <- universe



# # start_time <- proc.time()
# # for(i in 1:5) {
#   GOdata <- new("topGOdata",
#                    ontology = funsys,
#                    allGenes = geneList,
#                    nodeSize = 5,
#                    mapping = org_db,
#                    annotationFun = annFUN.org,
#                    ID = gene_id)

# #   resultFis <- topGO::runTest(GOdata, algorithm = "classic", statistic = "fisher")
# # }
# end_time <- proc.time()
# pack_time_taken <- end_time - start_time
# cat(pack_time_taken, file="~/pack_topgo_time_taken_5.txt")


#   save(list = ls(all.names = TRUE), file = "~/environment_test8.RData")

#   testthat::expect_identical(object=resultFis, expected=new_res_tg[["BP"]][[2]])


# })




test_that("Multienricher gsea gives same results as running packets separately with universe", {

  funsys <- "BP"
  org_db <- "org.Mm.eg.db"
  gene_id <- "entrez"

 # entrez_to_enrich <- c("13853", "16819", "12266", "497097", "19062")
 # entrez_to_enrich_kegg <- c("16176", "18049", "18211", "193034", "277328")

  genes <- c("11465", "11472", "12154", "12335", "12372", 
    "12373")

  universe <- c("15118", "16873", "12483", "15171", "67490", "18787", "11911", 
    "258390", "109019", "381484", "14377", "56706", "231474", "54427", 
    "192176", "70427", "20293", "13559", "66508", "14804", "19224", 
    "381290", "12909", "54633", "619846", "22234", "22160", "12464", 
    "101497", "16193", "258952", "503491", "208043", "12034", "629378", 
    "20355", "19106", "73679", "19824", "18596", "230857", "84505", 
    "14173", "11813", "19229", "13555", "68268", "229521", "12558", 
    "16869")

  universe <- unique(c(genes, universe))

  organisms_table <- get_organism_table()
  current_organism_info <- subset(organisms_table, 
                          rownames(organisms_table) == "Mouse")

# start_time <- proc.time()
  new_res_tg <- ExpHunterSuite::multienricher_topGO(all_funsys=funsys, universe=universe,
    genes_list=list(genes, genes),
    organism_info=current_organism_info)
# end_time <- proc.time()
# pack_time_taken <- end_time - start_time
# cat(pack_time_taken, file="~/new_topgo_time_taken_5.txt")

  library(topGO)
  funsys <- "BP"
  org_db <- "org.Mm.eg.db"
  gene_id <- "entrez"

  #go_to_genes <- topGO::annFUN.org(funsys, mapping = org_db, ID = gene_id)

  geneList <- factor(as.integer(universe %in% genes))
  names(geneList) <- universe



# start_time <- proc.time()
# for(i in 1:5) {
  GOdata <- new("topGOdata",
                   ontology = funsys,
                   allGenes = geneList,
                   nodeSize = 5,
                   mapping = org_db,
                   annotationFun = annFUN.org,
                   ID = gene_id)

  resultFis <- topGO::runTest(GOdata, algorithm = "classic", statistic = "fisher")
# }
# end_time <- proc.time()
# pack_time_taken <- end_time - start_time
# cat(pack_time_taken, file="~/pack_topgo_time_taken_5.txt")


  # save(list = ls(all.names = TRUE), file = "~/environment_test8.RData")

  testthat::expect_identical(object=resultFis, expected=new_res_tg[["BP"]][[2]])


})