test_that("main expression function works with 0 DEGs", {
  toc_file <- system.file("extData", "table_of_counts_no_changes.txt", 
  	                      package="ExpHunterSuite")
  toc <- read.table(toc_file, 
                    header=TRUE, row.names=1, sep="\t")

  target_file <- system.file("extData", "target.txt", 
  	                         package="ExpHunterSuite")

  #target <- read.table()
  target <- target_generation(from_file=target_file)

  degh_out <- main_degenes_Hunter(raw=toc, 
                                  target=target,
                                  modules="D",
                                  minpack_common=1)
  # IMPORTANT: SECTION TO CREATE THE GROUND TRUTH OUTPUT TO PERFORM THE TEST
  # Only uncomment this section when you need to regenerate it, i.e. new version of Bioc
  # precomp_degh_out <- degh_out
  # save(precomp_degh_out, file="../../../ExpHunterSuite/inst/extData/testdata/precomp_expression_results_no_changes.RData")

  precomp_degh_res_file <- system.file("extData", "testdata", "precomp_expression_results_no_changes.RData", 
  	                      package="ExpHunterSuite")
  load(precomp_degh_res_file)

  testthat::expect_equivalent(degh_out, precomp_degh_out)
  save(list = ls(all.names = TRUE), file = "~/environment_test3.RData")

})