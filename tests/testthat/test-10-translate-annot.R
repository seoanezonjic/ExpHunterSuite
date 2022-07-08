test_that("new main functional enrichment with translation", {

  annot_file <- system.file("extData", "ensembl_symbol_annot_external.txt", 
                           package="ExpHunterSuite")
  annot_table <-  read.table(annot_file, header=FALSE, row.names=NULL, 
        sep="\t", stringsAsFactors = FALSE, quote = "")

  ext_data_file <- system.file("extData", "table_of_counts_external.txt", 
                          package="ExpHunterSuite")
  ext_data <- read.table(ext_data_file, 
                    header=TRUE, row.names=1, sep="\t")
  degh_out <- main_degenes_Hunter(external_DEA_data=ext_data,
                                  modules="F",
                                  lfc=4)


  # save(list = ls(all.names = TRUE), file = "~/environment_test10.RData")

  fh_out_new <- main_functional_hunter( #Perform enrichment analysis
    hunter_results = degh_out,
    model_organism = 'Zebrafish',
    annot_table = annot_table,
    organisms_table = get_organism_table(),
    input_gene_id = "ENSEMBL",
    #func_annot_db = "gR",
    #GO_subont = "BMC",
    #analysis_type = "o", #g
    enrich_dbs = c("BP", "MF", "Reactome"),
    enrich_methods = c("ORA"),
    custom = NULL,
    pthreshold = 0.1,
    qthreshold = 0.2,
    cores = 1,
    task_size = 1,
    output_files = "results",
    fc_colname = "mean_logFCs"
  )

  testthat::expect_equal(nrow(fh_out_new$ORA$BP), 51)
})