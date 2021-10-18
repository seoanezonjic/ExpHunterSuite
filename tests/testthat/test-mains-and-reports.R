test_that("main expression function works", {
  toc_file <- system.file("extData", "table_of_counts.txt", 
  	                      package="ExpHunterSuite")
  print(toc_file)
  toc <- read.table(toc_file, 
                    header=TRUE, row.names=1, sep="\t")
  
  target_file <- system.file("extData", "target.txt", 
  	                         package="ExpHunterSuite")

  #target <- read.table()
  target <- target_generation(from_file=target_file)
  print(target[1:2,])
  print(as.character(target$sample[target$treat == "Treat"]))

  degh_out <- main_degenes_Hunter(raw=toc, 
                                  target=target,
                                  modules="DW")
  # IMPORTANT: SECTION TO CREATE THE GROUND TRUTH OUTPUT TO PERFORM THE TEST
  # Only uncomment this section when you need to regenerate it, i.e. new version of Bioc
  #precomp_degh_out <- degh_out
  #save(precomp_degh_out, file="../../../ExpHunterSuite/inst/extData/precomp_expression_results.RData")
  #save(list = ls(all.names = TRUE), file = "~/environment1.RData")

  precomp_degh_res_file <- system.file("extData", "precomp_expression_results.RData", 
  	                      package="ExpHunterSuite")
  load(precomp_degh_res_file)

  expect_equivalent(degh_out, precomp_degh_out)
})

test_that("writing expression report works", {
  precomp_degh_res_file <- system.file("extData", "precomp_expression_results.RData", 
                          package="ExpHunterSuite")
  load(precomp_degh_res_file)
  write_expression_report(exp_results=precomp_degh_out)


})

test_that("main functional enrichment function works", {

  precomp_degh_res_file <- system.file("extData", "precomp_expression_results.RData", 
                           package="ExpHunterSuite")
  load(precomp_degh_res_file)

  fh_out <- functional_hunter( #Perform enrichment analysis
         precomp_degh_out,
         'Mouse', # Use specified organism database 
         func_annot_db = "gR", # Enrichment analysis for GO, KEGG and Reactome
         GO_subont = "BMC",
         analysis_type= "o" # Use overepresentation analysis only (Not GSEA)
  )

  # IMPORTANT: SECTION TO CREATE THE GROUND TRUTH OUTPUT TO PERFORM THE TEST
  # Only uncomment this section when you need to regenerate it, i.e. new version of Bioc
  precomp_fh_out <- fh_out
  save(precomp_fh_out, file="../../../ExpHunterSuite/inst/extData/precomp_fh_out.RData")

  expect_equivalent(fh_out, precomp_fh_out)

})


