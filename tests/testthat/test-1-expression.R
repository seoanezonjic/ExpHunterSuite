test_that("main expression function works", {
  temp_folder <- file.path(find.package('ExpHunterSuite'), "tests", "tmp_folder","1-expression")
  toc_file <- system.file("extData", "table_of_counts_min.txt", 
  	                      package="ExpHunterSuite")
  toc <- read.table(toc_file, 
                    header=TRUE, row.names=1, sep="\t")
  
  target_file <- system.file("extData", "target.txt", 
  	                         package="ExpHunterSuite")

  target <- target_generation(from_file=target_file)

  degh_out <- main_degenes_Hunter(raw=toc, 
                                  target=target,
                                  modules="D",
                                  minpack_common=1,
                                  output_files = temp_folder)
  testthat::expect_equal(row.names(degh_out$DE_all_genes)[1:3], c("ENSMUSG00000055493", "ENSMUSG00000026822", "ENSMUSG00000024164"))
  testthat::expect_equal(degh_out$DE_all_genes[1:3, "FDR_DESeq2"],  c(1.105994e-186, 3.362235e-44, 6.545641e-29))

  write_expression_report(exp_results=degh_out, output_files = temp_folder, template_folder = file.path(find.package('ExpHunterSuite'), "inst", "templates"))
  testthat::expect(file.exists(file.path(temp_folder, "DEG_report.html")), failure_message="DEG_report.html has not been created")
  unlink(temp_folder, recursive = TRUE)
})