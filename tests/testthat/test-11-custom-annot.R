# test_that("custom functional enrichment function works ora", {

#   precomp_degh_res_file <- system.file("extData", "testdata", "precomp_expression_results.RData", 
#                            package="ExpHunterSuite")
#   load(precomp_degh_res_file)

#   gmt_file_1 <- system.file("extData", "toy_categories_test.gmt", package = "ExpHunterSuite")

#   custom_term2gene_1 <- ExpHunterSuite::load_and_parse_gmt(gmt_file_1)
#   custom_term2gene_2 <- ExpHunterSuite::load_and_parse_gmt(gmt_file_1)


#   custom_gmts <- list(gmt_file_1 = custom_term2gene_1, gmt_file_2=custom_term2gene_2)
#   print(precomp_degh_out$DE_all_genes)  
   

#   fh_out_new <- main_functional_hunter( #Perform enrichment analysis
#     hunter_results = precomp_degh_out,
#     model_organism = 'Mouse',
#     annot_table = NULL,
#     organisms_table = get_organism_table(),
#     input_gene_id = "ENSEMBL",
#     #func_annot_db = "gR",
#     #GO_subont = "C",
#     #analysis_type = "o", #g
#     enrich_dbs = c("MF", "Reactome"),
#     #enrich_methods = c("ora", "gsea", "topGO")
#     enrich_methods = "ORA",
#     custom = custom_gmts,
#     save_query = FALSE,
#     pthreshold = 0.9,
#     qthreshold = 0.9,
#     cores = 1,
#     task_size = 1,
#     output_files = "results",
#     fc_colname = "mean_logFCs"
#   )



#   fh_out <- functional_hunter( #Perform enrichment analysis
#     hunter_results = precomp_degh_out,
#     model_organism = 'Mouse',
#     annot_table = NULL,
#     organisms_table = get_organism_table(),
#     input_gene_id = "E",
#     func_annot_db = "gR",
#     GO_subont = "M",
#     custom = custom_gmts,
#     analysis_type = "o", #g
#     remote = "",
#     save_query = FALSE,
#     pthreshold = 0.9,
#     qthreshold = 0.9,
#     cores = 1,
#     task_size = 1,
#     output_files = "results",
#     fc_colname = "mean_logFCs"
#   )

#   save(list = ls(all.names = TRUE), file = "~/environment_test11.RData")
#   testthat::expect_equal(fh_out_new$CUSTOM$gmt_file_1, fh_out$CUSTOM$gmt_file_1)

# })


test_that("Custom functional enrichment function works ora with cluster", {

  precomp_degh_res_file <- system.file("extData", "testdata", "precomp_expression_results.RData", 
                           package="ExpHunterSuite")
  load(precomp_degh_res_file)

  precomp_degh_out$DE_all_genes$genes_tag[1:8] <- "PREVALENT_DEG"

  precomp_degh_out$DE_all_genes <- cbind(precomp_degh_out$DE_all_genes, 
    Cluster_ID = c(1,1,1,1, 2,0,2,2, sample(0:2, 8, replace=TRUE)))
   

  gmt_file_1 <- system.file("extData", "toy_categories_test.gmt", package = "ExpHunterSuite")

  custom_term2gene_1 <- ExpHunterSuite::load_and_parse_gmt(gmt_file_1)
  custom_term2gene_2 <- ExpHunterSuite::load_and_parse_gmt(gmt_file_1)


  custom_gmts <- list(gmt_file_1 = custom_term2gene_1, gmt_file_2=custom_term2gene_2)
  print(precomp_degh_out$DE_all_genes)  
   

  fh_out_new <- main_functional_hunter( #Perform enrichment analysis
    hunter_results = precomp_degh_out,
    model_organism = 'Mouse',
    annot_table = NULL,
    organisms_table = get_organism_table(),
    input_gene_id = "ENSEMBL",
    #func_annot_db = "gR",
    #GO_subont = "C",
    #analysis_type = "o", #g
    enrich_dbs = NULL,
    #enrich_methods = c("ora", "gsea", "topGO")
    enrich_methods = "ORA",
    custom = custom_gmts,
    save_query = FALSE,
    pthreshold = 0.9,
    qthreshold = 0.2,
    cores = 1,
    task_size = 1,
    output_files = "results",
    fc_colname = "mean_logFCs"
  )

  # save(list = ls(all.names = TRUE), file = "~/environment_test11.RData")


  fh_out <- functional_hunter( #Perform enrichment analysis
    hunter_results = precomp_degh_out,
    model_organism = 'Mouse',
    annot_table = NULL,
    organisms_table = get_organism_table(),
    input_gene_id = "E",
    func_annot_db = "gR",
    GO_subont = "M",
    custom = custom_gmts,
    analysis_type = "o", #g
    remote = "",
    save_query = FALSE,
    pthreshold = 0.9,
    qthreshold = 0.2,
    cores = 1,
    task_size = 1,
    output_files = "results",
    fc_colname = "mean_logFCs"
  )

  # save(list = ls(all.names = TRUE), file = "~/environment_test11.RData")
  # Can't test whol thing as previous implementation had genes named NA 
  testthat::expect_equal(as.data.frame(fh_out_new$ORA$gmt_file_1), as.data.frame(fh_out$CUSTOM$gmt_file_1))

})