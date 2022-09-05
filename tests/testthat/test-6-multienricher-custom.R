test_that("Multienricher custom category gives same results as old and packet implementations", {

	genes <- c("GeneB", "GeneC", "GeneE","GeneF")

  gmt_file_1 <- system.file("extData", "toy_categories_1.gmt", package = "ExpHunterSuite")
  gmt_file_2 <- system.file("extData", "toy_categories_2.gmt", package = "ExpHunterSuite")

  custom_term2gene_1 <- ExpHunterSuite::load_and_parse_gmt(gmt_file_1)
  custom_term2gene_2 <- ExpHunterSuite::load_and_parse_gmt(gmt_file_2)

  organisms_table <- get_organism_table()
  current_organism_info <- subset(organisms_table, 
                        rownames(organisms_table) == "Mouse")

  packet_res_custom_1 <- clusterProfiler::enricher(genes, 
                            pvalueCutoff = 0.2,
                            pAdjustMethod = "BH",
                   		      qvalueCutoff  = 0.2,
                            minGSSize = 3,
                            TERM2GENE = custom_term2gene_1)

  # deliberately chosen to return no enrichment
  packet_res_custom_2 <- clusterProfiler::enricher(genes, 
                            pvalueCutoff = 0.2,
                            pAdjustMethod = "BH",
                            qvalueCutoff  = 0.2,
                            minGSSize = 3,
                            TERM2GENE = custom_term2gene_2)

 res_all <- ExpHunterSuite::multienricher_ora( 
    genes_list=genes, minGSSize = 3, organism_info = current_organism_info,
  	custom_sets = list(custom_term2gene_1 = custom_term2gene_1, custom_term2gene_2 = custom_term2gene_2),
    qvalueCutoff  = 0.2, pvalueCutoff = 0.2, readable=FALSE, return_all=TRUE)

  testthat::expect_equal(packet_res_custom_1, res_all[["custom_term2gene_1"]])
  testthat::expect_equal(packet_res_custom_2, res_all[["custom_term2gene_2"]])

})