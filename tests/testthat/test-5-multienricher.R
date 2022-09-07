  test_that("Multienricher gives same results as running packets separately", {

  entrez_to_enrich_GO_react <- c("13853", "16819", "12266", "497097", "19062")
  entrez_to_enrich_kegg <- c("16176", "18049", "18211", "193034", "277328")

  packet_res_GO <- clusterProfiler::enrichGO(gene=entrez_to_enrich_GO_react, 
    OrgDb= "org.Mm.eg.db", 
    pAdjustMethod = "BH",
    ont = "BP",
    pvalueCutoff  = 0.05, qvalueCutoff = 0.2, 
    readable = FALSE)

  packet_res_reactome <- ReactomePA::enrichPathway(gene=entrez_to_enrich_GO_react, 
    organism      = "mouse",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.2,
    readable = FALSE)

  packet_res_GO <- DOSE::setReadable(packet_res_GO, OrgDb = "org.Mm.eg.db", 
        keyType="ENTREZID")

  packet_res_reactome <- DOSE::setReadable(packet_res_reactome, OrgDb = "org.Mm.eg.db", 
        keyType="ENTREZID")

  organisms_table <- get_organism_table()
  current_organism_info <- subset(organisms_table, 
                          rownames(organisms_table) == "Mouse")
 
#  new_res_all <- ExpHunterSuite::multienricher_ora(all_funsys=c("BP", "Reactome", "KEGG"),
  new_res_all <- ExpHunterSuite::multienricher_ora(all_funsys=c("BP", "Reactome"),
    genes_list=list(entrez_to_enrich_GO_react, entrez_to_enrich_kegg), 
    organism_info = current_organism_info,
#    kegg_file=kegg_file_mouse,
    pvalueCutoff = 0.05, qvalueCutoff = 0.2
  )


  testthat::expect_identical(object=packet_res_GO, expected=new_res_all[["BP"]][[1]])
})

test_that("Multienricher gives same results as running packets separately using universe", {

  organisms_table <- get_organism_table()
  current_organism_info <- subset(organisms_table, 
                           rownames(organisms_table) == "Mouse")

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
 
#  new_res_all <- ExpHunterSuite::multienricher_ora(all_funsys=c("BP", "Reactome", "KEGG"),
  new_res_all <- ExpHunterSuite::multienricher_ora(all_funsys=c("BP", "Reactome"),
    universe=c(genes,universe),
    genes_list=list(genes), 
    organism_info = current_organism_info,
#    kegg_file=kegg_file_mouse, 
    minGSSize=5,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE
  )

  packet_res_GO <- clusterProfiler::enrichGO(gene=genes,
    universe=c(genes,universe),
    OrgDb= "org.Mm.eg.db", minGSSize=5,
    pAdjustMethod = "BH",
    ont = "BP",
    pvalueCutoff  = 0.05, qvalueCutoff = 0.2, 
    readable = TRUE)

  packet_res_GO <- DOSE::setReadable(packet_res_GO, OrgDb = "org.Mm.eg.db", 
        keyType="ENTREZID")


  packet_res_reactome <- ReactomePA::enrichPathway(gene=genes, 
    universe=c(genes,universe),
    organism      = "mouse", minGSSize=5,
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.2,
    readable = TRUE)
  packet_res_reactome <- DOSE::setReadable(packet_res_reactome, OrgDb = "org.Mm.eg.db", 
        keyType="ENTREZID")

  testthat::expect_identical(object=packet_res_GO, expected=new_res_all[["BP"]][[1]])

})



