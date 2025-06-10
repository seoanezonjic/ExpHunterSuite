# Empty test
# test_that("Load and translate GMT file", {
#   organisms_table <- get_organism_table()
#   current_organism_info <- organisms_table[rownames(organisms_table) %in% "Mouse",]
#   org_db <- get_org_db(current_organism_info)

#   custom_file <- system.file("extData", "cytokine.gmt", package = "ExpHunterSuite")
#   custom_gmt <- ExpHunterSuite::load_and_parse_gmt(custom_file)
#   tr_gmt <- translate_gmt(custom_gmt, "SYMBOL", org_db)
# })