test_that("main functional enrichment function works Reactome", {

  precomp_degh_res_file <- system.file("extData", "testdata", "precomp_expression_results.RData", 
                           package="ExpHunterSuite")
  load(precomp_degh_res_file)
  organisms_table_file <- system.file("external_data", "organism_table.txt", 
                          package="ExpHunterSuite")
  fh_out <- main_functional_hunter( #Perform enrichment analysis
         precomp_degh_out,
         'Mouse', # Use specified organism database 
         enrich_dbs = c("Reactome"), # Enrichment for Reactome
         enrich_methods = "ORA",
  )
  
  testthat::expect_no_error(testthat::expect_no_error(nrow(as.data.frame(fh_out$ORA$MF))) > 0)
})
