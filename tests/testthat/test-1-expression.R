test_that("main expression function works", {
  toc_file <- system.file("extData", "table_of_counts_min.txt", 
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
  #precomp_degh_out <- degh_out
  #save(precomp_degh_out, file="../../../ExpHunterSuite/inst/extData/testdata/precomp_expression_results.RData")

  # Code to obtain a small table of counts
  #degs <- row.names(degh_out$DE_all_genes)[degh_out$DE_all_genes$genes_tag == "PREVALENT_DEG"][1:6]
  #toc_min <-toc[unique(c(row.names(toc)[1:10],degs)),]
  #toc_min_file <- "../../../ExpHunterSuite/inst/extData/table_of_counts_min.txt"
  #write.table(toc_min, file=toc_min_file, quote=FALSE, sep="\t")

  precomp_degh_res_file <- system.file("extData", "testdata", "precomp_expression_results.RData", 
  	                      package="ExpHunterSuite")
  print(precomp_degh_res_file)
  load(precomp_degh_res_file)
  # save(list = ls(all.names = TRUE), file = "~/environment1.RData")

  testthat::expect_equivalent(degh_out, precomp_degh_out)
})

# test_that("Correlation function works", {
#   toc_file <- system.file("extData", "table_of_counts.txt", 
#                           package="ExpHunterSuite")
#   toc <- read.table(toc_file, header=TRUE, row.names=1, sep="\t")
#   target_file <- system.file("extData", "target.txt", 
#                              package="ExpHunterSuite")

#   target <- target_generation(from_file=target_file)
  
#   coexp_out <- main_degenes_Hunter(raw=toc, 
#                                   target=target,
#                                   modules="DW")
#   # IMPORTANT: SECTION TO CREATE THE GROUND TRUTH OUTPUT TO PERFORM THE TEST
#   # You must uncomment this section the first time you run the test
#   #precomp_coexp_out <- coexp_out
#   #save(precomp_coexp_out, file="../../../ExpHunterSuite/inst/extData/precomp_coexpression_results.RData")

#   precomp_coexp_res_file <- system.file("extData", "precomp_coexpression_results.RData", 
#                          package="ExpHunterSuite")
#   load(precomp_coexp_res_file)

#   expect_equivalent(coexp_out, precomp_coexp_out)
# })

test_that("writing expression report works", {
  precomp_degh_res_file <- system.file("extData", "testdata", "precomp_expression_results.RData", 
                          package="ExpHunterSuite")
  load(precomp_degh_res_file)
  # write_expression_report(exp_results=precomp_degh_out)


})



