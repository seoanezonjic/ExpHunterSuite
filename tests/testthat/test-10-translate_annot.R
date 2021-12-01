test_that("new main functional enrichment with translation", {

  precomp_degh_res_file <- system.file("extData", "testdata", "precomp_expression_results.RData", 
                           package="ExpHunterSuite")
  load(precomp_degh_res_file)
  
  fh_out_new <- main_functional_hunter( #Perform enrichment analysis
    hunter_results = precomp_degh_out,
    model_organism = 'Mouse',
    annot_table = NULL,
    organisms_table = get_organism_table(),
    input_gene_id = "E",
    func_annot_db = "gR",
    GO_subont = "BMC",
    analysis_type = "o", #g
    #enrich_annot_sources = c("BP", "MF"),
    #enrich_methods = c("ora", "gsea", "topGO")
    custom = NULL,
    remote = "",
    save_query = FALSE,
    pthreshold = 0.1,
    qthreshold = 0.2,
    cores = 1,
    task_size = 1,
    output_files = "results",
    fc_colname = "mean_logFCs"
  )

  fh_out <- functional_hunter( #Perform enrichment analysis
    hunter_results = precomp_degh_out,
    model_organism = 'Mouse',
    annot_table = NULL,
    organisms_table = get_organism_table(),
    input_gene_id = "E",
    func_annot_db = "gR",
    GO_subont = "BMC",
    custom = NULL,
    analysis_type = "o", #g
    remote = "",
    save_query = FALSE,
    pthreshold = 0.1,
    qthreshold = 0.2,
    cores = 1,
    task_size = 1,
    output_files = "results",
    fc_colname = "mean_logFCs"
  )



  save(list = ls(all.names = TRUE), file = "~/environment_test9.RData")
  expect_equivalent(fh_out$ORA$GO_BP, fh_out_new$ORA$GO_BP)

  # $ORA$GO_BP
  # $ORA$GO_MF

})