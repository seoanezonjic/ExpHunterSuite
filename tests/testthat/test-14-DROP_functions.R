test_that("preprocess_gtf works as intended", {
  gtf <- system.file("extData/testdata", "gencode.v45.toy.annotation.gtf",
                     package = "ExpHunterSuite")
  txdb <- suppressMessages(txdbmaker::makeTxDbFromGFF(gtf))
  GenomeInfoDb::keepStandardChromosomes(txdb)
  count_ranges <- GenomicFeatures::exonsBy(txdb, by = "gene")
  gene_name_mapping <- .map_genes(gtf)
  count_ranges@metadata$genomeInfo$`Creation time` <- ""
  # Replace with GRanges object, comparing environments is messier
  expected_output <- list(txdb = GenomicFeatures::transcripts(txdb),
                          count_ranges = count_ranges,
                          gene_name_mapping = gene_name_mapping)
  expected_output$txdb@metadata$genomeInfo$`Creation time` <- ""
  actual_output <- suppressMessages(preprocess_gtf(gtf))
  actual_output$txdb <- GenomicFeatures::transcripts(actual_output$txdb)
  actual_output$txdb@metadata$genomeInfo$`Creation time` <- ""
  actual_output$count_ranges@metadata$genomeInfo$`Creation time` <- ""
  testthat::expect_equal(actual_output, expected_output)
})

test_that("map_genes works as intended", {
  gtf <- system.file("extData/testdata", "gencode.v45.toy.annotation.gtf",
                     package = "ExpHunterSuite")
  actual_output <- .map_genes(gtf)
  gtf_name <- basename(tools::file_path_sans_ext(gtf))
  gtf_df <- as.data.frame(rtracklayer::import(gtf))
  if (!"gene_name" %in% colnames(gtf_df)) {
    gtf_df$gene_name <- gtf_df$gene_id
  }
  gtf_df <- gtf_df[gtf_df$type == 'gene',]
  if('gene_biotype' %in% colnames(gtf_df)) {
    colnames(gtf_df)[colnames(gtf_df) == "gene_biotype"] <- "gene_type"
  }
  # Subset to the following columns only
  columns <- c('seqnames', 'start', 'end', 'strand', 'gene_id',
               'gene_name', 'gene_type', 'gene_status')
  columns <- intersect(columns, colnames(gtf_df))
  gtf_df <- gtf_df[, columns]
  expected_output <- gtf_df
  testthat::expect_equal(actual_output, expected_output)
})

test_that(".split_string_by_char works as intended", {
  string <- "Lorem ipsum dolor sit_amet"
  substr_1 <- .split_string_by_char(string, "_", 1)
  substr_2 <- .split_string_by_char(string, "_", 2)
  expected_substr_1 <- "Lorem ipsum dolor sit"
  expected_substr_2 <- "amet"
  testthat::expect_equal(substr_1, expected_substr_1)
  testthat::expect_equal(substr_2, expected_substr_2)
})

test_that(".split_string_by_char edge cases", {
  string <- "._?!#$"
  testthat::expect_warning(.split_string_by_char(string, "(", 1),
                           "character not present")
  testthat::expect_warning(.split_string_by_char(string, ".", 0),
                           "out-of-bounds")
  testthat::expect_warning(.split_string_by_char(string, ".", 200),
                           "out-of-bounds")
})
