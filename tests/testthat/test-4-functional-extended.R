test_that("main functional enrichment function works GO and Reactome - mouse", {

  precomp_degh_res_file <- system.file("extData", "testdata", "precomp_expression_results.RData", 
                           package="ExpHunterSuite")
  load(precomp_degh_res_file)

  organisms_table_file <- system.file("external_data", "organism_table.txt", 
                          package="ExpHunterSuite")
  organisms_table <- get_organism_table(organisms_table_file)

  fh_ext_out <- main_functional_hunter( #Perform enrichment analysis
         precomp_degh_out,
         'Mouse', # Use specified organism database 
         organisms_table = organisms_table,
         enrich_dbs = c("MF", "BP","Reactome"), # Enrichment analysis for GO, KEGG and Reactome
         enrich_methods = "ORA",
  )
  expect_equal(nrow(fh_ext_out$ORA$Reactome), 12)
})
