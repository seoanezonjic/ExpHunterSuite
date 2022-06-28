test_that("main expression function works with 0 DEGs", {
  toc_file <- system.file("extData", "table_of_counts_no_changes.txt", 
  	                      package="ExpHunterSuite")
  toc <- read.table(toc_file, 
                    header=TRUE, row.names=1, sep="\t")

  target_file <- system.file("extData", "target.txt", 
  	                         package="ExpHunterSuite")

  target <- target_generation(from_file=target_file)

  degh_out <- main_degenes_Hunter(raw=toc, 
                                  target=target,
                                  modules="D",
                                  minpack_common=1)

 testthat::expect_equivalent(row.names(degh_out$DE_all_genes)[1:3], 
  c("ENSMUSG00000051951", "ENSMUSG00000102348", "ENSMUSG00000103201"))

})
