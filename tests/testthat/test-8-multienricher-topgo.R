
test_that("multiple gene list Multienricher topGO", {

  funsys <- c("BP", "MF")
  org_db <- "org.Mm.eg.db"
  gene_id <- "entrez"

  genes <- list(a = c("11465", "11472", "12154", "12335", "12372", "12373"),
                b = c("13853", "16819", "12266", "497097", "19062"))

  universe <- c("15118", "16873", "12483", "15171", "67490", "18787", "11911", 
    "258390", "109019", "381484", "14377", "56706", "231474", "54427", 
    "192176", "70427", "20293", "13559", "66508", "14804", "19224", 
    "381290", "12909", "54633", "619846", "22234", "22160", "12464", 
    "101497", "16193", "258952", "503491", "208043", "12034", "629378", 
    "20355", "19106", "73679", "19824", "18596", "230857", "84505", 
    "14173", "11813", "19229", "13555", "68268", "229521", "12558", 
    "16869")

  universe <- unique(c(unlist(genes), universe))

  organisms_table <- get_organism_table()
  current_organism_info <- subset(organisms_table, 
                          rownames(organisms_table) == "Mouse")

  new_res_tg <- multienricher_topGO(all_funsys=funsys, universe=universe,
    genes_list=genes,
    organism_info=current_organism_info)


  library(topGO)
  funsys <- "BP"
  org_db <- "org.Mm.eg.db"
  gene_id <- "entrez"


  ### PREPARING EXPECTED 1
  geneList <- factor(as.integer(universe %in% genes[[1]]))
  names(geneList) <- universe

  GOdata <- new("topGOdata",
                   ontology = funsys,
                   allGenes = geneList,
                   nodeSize = 5,
                   mapping = org_db,
                   annotationFun = annFUN.org,
                   geneSel =  get_default_geneSel("fisher"),
                   ID = gene_id)
  resultFis_1 <- topGO::runTest(GOdata, algorithm = "classic", statistic = "fisher")

  restable_1 <- topGO::GenTable(GOdata, 
                                 p.value = resultFis_1,
                                 topNodes = length(resultFis_1@score), 
                                 format.FUN = function(x, dig, eps){x})
  restable_1$fdr <- p.adjust(restable_1$p.value, method = "BH", n = length(restable_1$p.value))


  ### PREPARING EXPECTED 2
  geneList2 <- factor(as.integer(universe %in% genes[[2]]))
  names(geneList2) <- universe
  GOdata2 <- topGO::updateGenes(object = GOdata, geneList = geneList2)
  resultFis_2 <- topGO::runTest(GOdata2, algorithm = "classic", statistic = "fisher")
  restable_2 <- topGO::GenTable(GOdata2, 
                                 p.value = resultFis_2,
                                 topNodes = length(resultFis_2@score), 
                                 format.FUN = function(x, dig, eps){x})
  restable_2$fdr <- p.adjust(restable_2$p.value, method = "BH", n = length(restable_2$p.value))



  testthat::expect_identical(object=restable_1, expected=new_res_tg[["BP"]][[1]])
  testthat::expect_identical(object=restable_2, expected=new_res_tg[["BP"]][[2]])

})



test_that("unique gene list Multienricher topGO", {

  funsys <- c("BP", "MF")
  org_db <- "org.Mm.eg.db"
  gene_id <- "entrez"

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

  new_res_tg <- multienricher_topGO(all_funsys=funsys, universe=universe,
    genes_list=genes,
    organism_info=current_organism_info)


  library(topGO)
  funsys <- "BP"
  org_db <- "org.Mm.eg.db"
  gene_id <- "entrez"


  geneList <- factor(as.integer(universe %in% genes))
  names(geneList) <- universe

  GOdata <- new("topGOdata",
                   ontology = funsys,
                   allGenes = geneList,
                   nodeSize = 5,
                   mapping = org_db,
                   annotationFun = annFUN.org,
                   ID = gene_id)

  resultFis <- topGO::runTest(GOdata, algorithm = "classic", statistic = "fisher")
  restable <- topGO::GenTable(GOdata, 
                               p.value = resultFis,
                               topNodes = length(resultFis@score), 
                               format.FUN = function(x, dig, eps){x})
   restable$fdr <- p.adjust(restable$p.value, method = "BH", n = length(restable$p.value))

  testthat::expect_identical(object=restable, expected=new_res_tg[["BP"]])

})


test_that("multi_topGOTest BP", {


  funsys <-"BP"
  org_db <- "org.Mm.eg.db"
  gene_id <- "entrez"

  genes <- list(c("11465", "11472", "12154", "12335", "12372", 
    "12373"))

  universe <- c("15118", "16873", "12483", "15171", "67490", "18787", "11911", 
    "258390", "109019", "381484", "14377", "56706", "231474", "54427", 
    "192176", "70427", "20293", "13559", "66508", "14804", "19224", 
    "381290", "12909", "54633", "619846", "22234", "22160", "12464", 
    "101497", "16193", "258952", "503491", "208043", "12034", "629378", 
    "20355", "19106", "73679", "19824", "18596", "230857", "84505", 
    "14173", "11813", "19229", "13555", "68268", "229521", "12558", 
    "16869")

  universe <- unique(c(unlist(genes), universe))

  organisms_table <- get_organism_table()
  current_organism_info <- subset(organisms_table, 
                          rownames(organisms_table) == "Mouse")

  geneList <- factor(as.integer(universe %in% unlist(genes)))
  names(geneList) <- universe

  GOdata <- new("topGOdata",
                   ontology = funsys,
                   allGenes = geneList,
                   nodeSize = 5,
                   mapping = org_db,
                   annotationFun = annFUN.org,
                   ID = gene_id)
  new_res_tg <- multi_topGOTest(funsys=funsys, universe=universe,
                                genes_list=genes, org_db = org_db )

  resultFis <- topGO::runTest(GOdata, algorithm = "classic", statistic = "fisher")
  restable <- topGO::GenTable(GOdata, 
                             p.value = resultFis,
                             topNodes = length(resultFis@score), 
                             format.FUN = function(x, dig, eps){x})
  restable$fdr <- p.adjust(restable$p.value, method = "BH", n = length(restable$p.value))

  testthat::expect_identical(object=restable, expected=new_res_tg)

})


test_that("multi_topGOTest BP ks", {


  funsys <-"BP"
  org_db <- "org.Mm.eg.db"
  gene_id <- "entrez"

  genes <- list(c(1,2,3,4,5,6))
  names(genes[[1]]) <-  c("11465", "11472", "12154", "12335", "12372", 
    "12373")



  organisms_table <- get_organism_table()
  current_organism_info <- subset(organisms_table, 
                          rownames(organisms_table) == "Mouse")
  GOdata <- new("topGOdata",
                   ontology = funsys,
                   allGenes = genes[[1]],
                   nodeSize = 5,
                   mapping = org_db,
                   geneSel = function(x){x},
                   annotationFun = annFUN.org,
                   ID = gene_id)
  new_res_tg <- multi_topGOTest(funsys=funsys, universe=NULL,
                                genes_list=genes, org_db = org_db, statistic = "ks")

  resultKS <- topGO::runTest(GOdata, algorithm = "classic", statistic = "ks")
  restable <- topGO::GenTable(GOdata, 
                             p.value = resultKS,
                             topNodes = length(resultKS@score), 
                             format.FUN = function(x, dig, eps){x})
  restable$fdr <- p.adjust(restable$p.value, method = "BH", n = length(restable$p.value))
  testthat::expect_identical(object=restable, expected=new_res_tg)

})