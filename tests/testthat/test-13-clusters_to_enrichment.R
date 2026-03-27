test_that("Run c2e main", {
  organisms_table <- get_organism_table()
  current_organism_info <- organisms_table[rownames(organisms_table) %in% "Mouse",]
  org_db <- get_org_db(current_organism_info)

  input_file <- system.file("extData", "cluster_genes.txt", package = "ExpHunterSuite")
  enr_lists <- main_clusters_to_enrichment(input_file, org_db=org_db, 
    current_organism_info=current_organism_info, gene_keytype="ENSEMBL")
  top_cat <- as.data.frame(enr_lists$enrichments_ORA$MF$two)[1, "ID"]
  
  testthat::expect_equal(top_cat, "GO:0004866")

  
  #compC <- clusterProfiler::merge_result(enr_lists$enrichments_ORA$MF)
  compC <- process_cp_list(enr_lists$enrichments_ORA,
    simplify_results = FALSE, clean_parentals = FALSE)
  showCats <- calc_showCat_compareCluster(compC$MF@compareClusterResult,
   10, 20)
  testthat::expect_equal(showCats, 10)
})

