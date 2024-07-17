test_that("Run c2e main", {
  organisms_table <- get_organism_table()
  current_organism_info <- organisms_table[rownames(organisms_table) %in% "Mouse",]
  org_db <- get_org_db(current_organism_info)

  input_file <- system.file("extData", "cluster_genes.txt", package = "ExpHunterSuite")
  enr_lists <- main_clusters_to_enrichment(input_file, org_db=org_db, 
    current_organism_info=current_organism_info, gene_keytype="ENSEMBL")


   testthat::expect_equal(
    as.data.frame(enr_lists$enrichments_ORA$MF$two)[1, "ID"], 
    "GO:0004866")
})