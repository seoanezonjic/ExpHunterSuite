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
  precomp_degh_out <- degh_out
  save(precomp_degh_out, file="../../../ExpHunterSuite/inst/extData/precomp_expression_results_no_changes.RData")

  # Code to obtain a small table of counts
  #degs <- row.names(degh_out$DE_all_genes)[degh_out$DE_all_genes$genes_tag == "PREVALENT_DEG"][1:6]
  #toc_no_changes <-toc[1:15, ]
  #toc_no_changes_file <- "../../../ExpHunterSuite/inst/extData/table_of_counts_no_changes.txt"
  #write.table(toc_no_changes, file=toc_no_changes_file, quote=FALSE, sep="\t")

  precomp_degh_res_file <- system.file("extData", "precomp_expression_results_no_changes.RData", 
  	                      package="ExpHunterSuite")
  load(precomp_degh_res_file)
  save(list = ls(all.names = TRUE), file = "~/environment_exp_0_changes.RData")

  testthat::expect_equivalent(degh_out, precomp_degh_out)

  save(list = ls(all.names = TRUE), file = "~/environment_test3.RData")

})