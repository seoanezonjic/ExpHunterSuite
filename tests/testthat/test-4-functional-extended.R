test_that("main functional enrichment function works GO and Reactome - mouse", {

  precomp_degh_res_file <- system.file("extData", "testdata", "precomp_expression_results.RData", 
                           package="ExpHunterSuite")
  load(precomp_degh_res_file)

  organisms_table_file <- system.file("external_data", "organism_table.txt", 
                          package="ExpHunterSuite")
  organisms_table <- get_organism_table(organisms_table_file)

  fh_ext_out <- functional_hunter( #Perform enrichment analysis
         precomp_degh_out,
         'Mouse', # Use specified organism database 
         organisms_table = organisms_table,
         func_annot_db = "gR", # Enrichment analysis for GO, KEGG and Reactome
         GO_subont = "BMC",
         analysis_type= "o" # Use overepresentation analysis only (Not GSEA)
  )

  # IMPORTANT: SECTION TO CREATE THE GROUND TRUTH OUTPUT TO PERFORM THE TEST
  # Only uncomment this section when you need to regenerate it, i.e. new version of Bioc
  #precomp_fh_ext_out <- fh_ext_out
  #save(precomp_fh_ext_out, file="../../../ExpHunterSuite/inst/extData/testdata/precomp_fh_ext_out.RData")
  # For the future - when we can't compare the objects directly - check the tables of enriched functions  
  #go_cc_res <- as.data.frame(fh_ext_out$ORA[["GO_CC"]])
  #go_react_res <- as.data.frame(fh_ext_out$ORA[["REACT"]])
  #save(go_cc_res, go_react_res, file="../../../ExpHunterSuite/inst/extData/testdata/fh_ext_enrich_tables.RData")

  precomp_fh_ext_file <- system.file("extData", "testdata", "precomp_fh_ext_out.RData", 
                          package="ExpHunterSuite")
  load(precomp_fh_ext_file)

  # save(list = ls(all.names = TRUE), file = "~/environment_test4.RData")
  expect_equivalent(fh_ext_out, precomp_fh_ext_out)

})

# test_that("main functional enrichment function works when no DEGs", {

#   precomp_degh_res_file <- system.file("extData", "precomp_expression_results_no_changes.RData", 
#                            package="ExpHunterSuite")
#   load(precomp_degh_res_file)

#   organisms_table_file <- system.file("external_data", "organism_table.txt", 
#                           package="ExpHunterSuite")
#   organisms_table <- get_organism_table(organisms_table_file)

#   fh_ext_out <- functional_hunter( #Perform enrichment analysis
#          precomp_degh_out,
#          'Mouse', # Use specified organism database 
#          organisms_table = organisms_table,
#          func_annot_db = "gR", # Enrichment analysis for GO, KEGG and Reactome
#          GO_subont = "M",
#          analysis_type= "o" # Use overepresentation analysis only (Not GSEA)
#   )
#   # IMPORTANT: SECTION TO CREATE THE GROUND TRUTH OUTPUT TO PERFORM THE TEST
#   # Only uncomment this section when you need to regenerate it, i.e. new version of Bioc
#   #precomp_fh_ext_out <- fh_ext_out
#   #save(precomp_fh_ext_out, file="../../../ExpHunterSuite/inst/extData/precomp_fh_no_changes_out.RData")
#   precomp_fh_ext_file <- system.file("extData", "precomp_fh_no_changes_out.RData", 
#                           package="ExpHunterSuite")
#   load(precomp_fh_ext_file)
#   save(list = ls(all.names = TRUE), file = "~/environment_fh_no_degs.RData")


#   expect_equivalent(fh_ext_out, precomp_fh_ext_out)

# })