test_that("main functional enrichment function works Reactome", {

  precomp_degh_res_file <- system.file("extData", "testdata", "precomp_expression_results.RData", 
                           package="ExpHunterSuite")
  load(precomp_degh_res_file)
  organisms_table_file <- system.file("external_data", "organism_table.txt", 
                          package="ExpHunterSuite")
#  organisms_table <- get_organism_table(organisms_table_file)

  fh_out <- main_functional_hunter( #Perform enrichment analysis
         precomp_degh_out,
         'Mouse', # Use specified organism database 
#         organisms_table = organisms_table,
         enrich_dbs = c("MF", "BP","Reactome"), # Enrichment analysis for GO, KEGG and Reactome
         enrich_methods = "ORA",
  )

  testthat::expect_equal(nrow(as.data.frame(fh_out$ORA$MF)), 22)
  testthat::expect_equal(as.data.frame(fh_out$ORA$BP)[1:2,"ID"], c("GO:0043277", "GO:0006911"))
})
