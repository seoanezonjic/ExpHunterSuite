test_that("main functional enrichment function works Reactome", {

  precomp_degh_res_file <- system.file("extData", "testdata", "precomp_expression_results.RData", 
                           package="ExpHunterSuite")
  load(precomp_degh_res_file)
  organisms_table_file <- system.file("external_data", "organism_table.txt", 
                          package="ExpHunterSuite")
#  organisms_table <- get_organism_table(organisms_table_file)

  fh_out <- functional_hunter( #Perform enrichment analysis
         precomp_degh_out,
         'Mouse', # Use specified organism database 
#         organisms_table = organisms_table,
         func_annot_db = "R", # Enrichment analysis for GO, KEGG and Reactome
         GO_subont = "BMC",
         analysis_type= "o" # Use overepresentation analysis only (Not GSEA)
  )

  # IMPORTANT: SECTION TO CREATE THE GROUND TRUTH OUTPUT TO PERFORM THE TEST
  # Only uncomment this section when you need to regenerate it, i.e. new version of Bioc
  # precomp_fh_out <- fh_out
  # save(precomp_fh_out, file="../../../ExpHunterSuite/inst/extData/testdata/precomp_fh_out.RData")
  
  precomp_fh_file <- system.file("extData", "testdata", "precomp_fh_out.RData", 
                          package="ExpHunterSuite")
  load(precomp_fh_file)
  expect_equivalent(fh_out, precomp_fh_out)
    
  # save(list = ls(all.names = TRUE), file = "~/environment_test2.RData")
})