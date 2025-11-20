
markers_df <- data.frame(samples = 1:4, seurat_clusters = 1:4,
                         gene = LETTERS[1:4])
markers_df$p_val_adj <- rep(1e-5, nrow(markers_df))
markers_df$avg_log2FC <- rep(1, nrow(markers_df))

test_that("match_cell_types, simple case", {
  cell_annotation <- data.frame(markers = c("A", "B", "C"),
               type = c("Type1", "Type2", "Type3"))
  expected_df <- markers_df
  expected_df$cell_type <- c("1. Type1", "2. Type2", "3. Type3", "4. Unknown")
  expected_df <- expected_df[order(expected_df$seurat_clusters), ]
  output_df <- match_cell_types(markers_df = markers_df,
                  cell_annotation = cell_annotation)$stats_table
  expect_equal(output_df, expected_df)
})

test_that("match_cell_types, unsorted simple case", {
  test_markers_df <- markers_df
  test_markers_df$samples <- 4:1
  test_markers_df$seurat_clusters <- c(2, 3, 1, 4)
  test_markers_df$gene <- c("A", "C", "B", "D")
  cell_annotation <- data.frame(markers = c("C", "B", "A"),
                                type = c("Type3", "Type2", "Type1"))
  expected_df <- test_markers_df
  expected_df$cell_type <- c("2. Type1", "3. Type3", "1. Type2", "4. Unknown")
  expected_df <- expected_df[order(expected_df$seurat_clusters), ]
  output_df <- match_cell_types(test_markers_df, cell_annotation)$stats_table
  expect_equal(output_df, expected_df)
})

test_that("match_cell_types, complex case", {
  test_markers_df <- data.frame(samples = seq(1, 50),
                       seurat_clusters = c(rep(1, 20), rep(2, 20), rep (3, 10)))
  test_markers_df$p_val_adj <- rep(1e-5, nrow(test_markers_df))
  test_markers_df$avg_log2FC <- rep(1, nrow(test_markers_df))
  genes <- paste0("gene", seq(1, 50))
  test_markers_df$gene <- genes
  types <- c(rep("type1", 20), rep("type2", 13), rep("type3", 17))
  cell_annotation <- data.frame(markers = genes, type = types)
  expected_df <- test_markers_df
  expected_df$cell_type <-  c(rep("1. type1", 20), rep("2. type2", 20),
                rep("3. type3", 10))
  expected_df <- expected_df[order(expected_df$seurat_clusters), ]
  output_df <- match_cell_types(test_markers_df, cell_annotation)$stats_table
  expect_equal(output_df, expected_df)
})

test_that("match_cell_types, multiple clusters assigned the same cell type", {
  test_markers_df <- markers_df
  genes <- paste0("gene", seq(3))
  genes <- c(genes, paste0("gene", 1))
  test_markers_df$gene <- genes
  types <- c("type1", "type2", "type3")
  cell_annotation <- data.frame(markers = genes[1:12], type = types)
  expected_df <- test_markers_df
  expected_df$cell_type <- c("1. type1 (a)", "2. type2", "3. type3",
                             "4. type1 (b)")
  output_df <- match_cell_types(test_markers_df, cell_annotation)$stats_table
  expect_equal(output_df, expected_df)
})

test_that("match_cell_types, draw between more than one cell type", {
  test_markers_df <- markers_df
  genes <- paste0("gene", c(1, 1, 2))
  genes <- c(genes, paste0("gene", 1))
  test_markers_df$gene <- genes
  types <- c("type1", "type2", "type3")
  cell_annotation <- data.frame(markers = genes[1:12], type = types)
  expected_df <- test_markers_df
  expected_df$cell_type <- c("1. type1 / type2 (a)", "2. type1 / type2 (b)",
                             "3. type3", "4. type1 / type2 (c)")
  output_df <- match_cell_types(test_markers_df, cell_annotation)$stats_table
  expect_equal(output_df, expected_df)
})

test_that("match_cell_types, cluster with no significant markers", {
  test_markers_df <- markers_df
  test_markers_df$p_val_adj[1] <- 1
  test_markers_df$seurat_clusters <- 4:1
  genes <- LETTERS[1:3]
  types <- c("type1", "type2", "type3")
  cell_annotation <- data.frame(markers = genes, type = types)
  expected_df <- test_markers_df[order(test_markers_df$seurat_clusters), ]
  expected_df$gene[4] <- "None"
  expected_df$cell_type <- c("1. Unknown (a)", "2. type3",
                             "3. type2", "4. Unknown (b)")
  expect_warning(match_cell_types(test_markers_df, cell_annotation),
                 "Cluster 4 contains no significant markers")
  output_df <- suppressWarnings(match_cell_types(test_markers_df,
                              cell_annotation))$stats_table
  expect_equal(output_df, expected_df)
})

test_that("match_cell_types where one cluster has no markers", {
  test_markers_df <- markers_df[-3, ]
  genes <- LETTERS[1:3]
  cell_annotation <- data.frame(markers = LETTERS[1:3],
                                type = c("type1", "type2", "type3"))
  expected_df <- data.frame(samples = c(1, 2, 4), seurat_clusters = c(1, 2, 4),
                            gene = c("A", "B", "D"), p_val_adj = 1e-05,
                            avg_log2FC = 1, cell_type = c("1. type1",
                            "2. type2", "4. Unknown"))
  output_df <- match_cell_types(markers = test_markers_df,
                                cell_annotation = cell_annotation)$stats_table
  rownames(output_df) <- NULL
  expect_equal(output_df, expected_df)
})

test_that("collapse_markers base use case", {
  markers_1 <- data.frame(a_avg_log2FC = rep(1, 3), b_avg_log2FC = rep (-1, 3))
  markers_2 <- data.frame(a_avg_log2FC = rep(-2, 3), b_avg_log2FC = rep (2, 3))
  markers_list <- list(cluster_1 = markers_1, cluster_2 = markers_2)
  output_df <- collapse_markers(markers_list)
  expected_df <- data.frame(gene = as.character(rep(1:3, 2)),
                            a_avg_log2FC = c(rep(1, 3), rep(-2, 3)),
                            b_avg_log2FC = c(rep(-1, 3), rep (2, 3)),
                            avg_log2FC = rep(0, 6))
  expect_equal(output_df, expected_df)
})

data(pbmc_tiny)
test_pbmc <- pbmc_tiny
test_pbmc@meta.data$groups <- "g1"
test_pbmc@meta.data$groups[7:15] <- "g2"
test_pbmc@meta.data$seurat_clusters <- 0
test_pbmc@meta.data$seurat_clusters[1] <- 1
test_pbmc@meta.data$seurat_clusters[15] <- 1
test_pbmc@meta.data$cell_types <- "typeA"
test_pbmc@meta.data$cell_types[3:5] <- "typeB"
test_pbmc@meta.data$cell_types[15] <- "typeB"

test_that(".has_exclusive_idents works in base case", {
  expected_warning <- "Defaulting to general marker analysis"
  expect_warning(.has_exclusive_idents(seu = test_pbmc, cond = "groups",
                                  idents = "seurat_clusters"), expected_warning)
  expect_true(suppressWarnings(.has_exclusive_idents(seu = test_pbmc,
                                  idents = "seurat_clusters", cond = "groups")))
})

test_that(".has_exclusive_idents works with alternate ident", {
  expected_warning <- "Defaulting to general marker analysis"
  expect_warning(.has_exclusive_idents(seu = test_pbmc, cond = "groups",
                                  idents = "cell_types"), expected_warning)
  expect_true(suppressWarnings(.has_exclusive_idents(seu = test_pbmc,
                                  idents = "cell_types", cond = "groups")))
})

test_that(".has_exclusive_idents can handle more than one positive", {
  expected_warning <- paste0("g1-1, g2-1")
  expect_warning(.has_exclusive_idents(seu = test_pbmc, cond = "groups",
                                  idents = "seurat_clusters"), expected_warning)
  expect_true(suppressWarnings(.has_exclusive_idents(seu = test_pbmc,
                                  idents = "seurat_clusters", cond = "groups")))
})

test_that(".has_exclusive_idents can predict pairs that should appear but
           do not (behaves correctly for strictly exclusive idents)", {
  test_pbmc <- pbmc_tiny
  test_pbmc@meta.data$groups <- "g1"
  test_pbmc@meta.data$groups[7:15] <- "g2"
  test_pbmc@meta.data$seurat_clusters <- 0
  test_pbmc@meta.data$seurat_clusters[7:15] <- 1
  expected_warning <- paste0("g2-0, g1-1")
  expect_warning(.has_exclusive_idents(seu = test_pbmc, cond = "groups",
                                  idents = "seurat_clusters"), expected_warning)
  expect_true(suppressWarnings(.has_exclusive_idents(seu = test_pbmc,
                                  idents = "seurat_clusters", cond = "groups")))
})

test_that(".has_exclusive_idents can identify that no exclusive idents
           exist", {
  test_pbmc <- pbmc_tiny
  test_pbmc@meta.data$seurat_clusters <- 0
  test_pbmc@meta.data$seurat_clusters[c(8:15)] <- 1
  expect_false(.has_exclusive_idents(seu = test_pbmc,
                                  idents = "seurat_clusters", cond = "groups"))
})

test_that("get_clusters_distribution properly calculates percentages", {
  test_pbmc <- pbmc_tiny
  test_pbmc@meta.data$sample <- c(rep("A", 5), rep("B", 5), rep("C", 5))
  test_pbmc@meta.data$seurat_clusters <- c(rep(1, 3), rep(2, 10), rep(3, 2))
  expected_pct <- matrix(nrow = 3, ncol = 3)
  expected_pct[1, ] <- c(0.6, 0.4, 0)
  expected_pct[2, ] <- c(0, 1.0, 0)
  expected_pct[3, ] <- c(0, 0.6, 0.4)
  expected_pct <- data.frame(expected_pct) * 100
  colnames(expected_pct) <- 1:3
  rownames(expected_pct) <- c("A", "B", "C")
  output_pct <- suppressMessages(get_clusters_distribution(test_pbmc, 2))
  expect_equal(output_pct, expected_pct)
})

query <- c("PPBP", "IGLL5", "VDAC3")

test_that("get_query_distribution properly sums expression levels in samples", {
  test_pbmc <- pbmc_tiny
  test_pbmc@meta.data$sample <- c(rep("A", 5), rep("B", 5), rep("C", 5))
  expected_exp <- data.frame(matrix(nrow = 3, ncol = 3))
  rownames(expected_exp) <- c("A", "B", "C")
  expected_exp[, 1] <- c(4.75, 0, 3.98)
  expected_exp[, 2] <- c(0, 0, 12.4)
  expected_exp[, 3] <- c(4.38, 9.67, 7.33)
  colnames(expected_exp) <- query
  output_exp <- get_query_distribution(seu = test_pbmc, query = query,
                sigfig = 3, layer = "data")
  expect_equal(output_exp, expected_exp)
})

test_pbmc <- pbmc_tiny
test_pbmc@meta.data$sample <- c(rep("A", 5), rep("B", 5), rep("C", 5))
test_pbmc@meta.data$seurat_clusters <- rep(0:4, 3)
cell_types <- c("0. typeA", "1. typeB", "2. typeC", "3. typeD", "4. typeE")
test_pbmc@meta.data$cell_type <- rep(cell_types, 3)

test_that("get_query_pct works in simple case", {
  expected_df <- matrix(nrow = 3, ncol = 3)
  rownames(expected_df) <- c("A", "B", "C")
  expected_df[1, ] <- c(20, 0, 20)
  expected_df[2, ] <- c(0, 0, 40)
  expected_df[3, ] <- c(20, 40, 20)
  colnames(expected_df) <- query
  output_df <- suppressMessages(get_query_pct(seu = test_pbmc, query = query,
                              by = "sample", layer = "counts", min_counts = 0))
  expect_equal(output_df, expected_df)
})

test_that("get_query_pct can apply minimum counts filter", {
  expected_df <- matrix(nrow = 3, ncol = 3)
  rownames(expected_df) <- c("A", "B", "C")
  expected_df[1, ] <- c(rep(0, 3))
  expected_df[2, ] <- c(0, 0, 0)
  expected_df[3, ] <- c(0, 20, 20)
  colnames(expected_df) <- query
  output_df <- suppressMessages(get_query_pct(seu = test_pbmc, query = query,
                              by = "sample", layer = "counts", min_counts = 10))
  expect_equal(output_df, expected_df)
})

test_that("get_query_pct gives warning if any query genes are not found", {
  missing_query <- c(query, "NOEXPA", "NOEXPB")
  expected_df <- matrix(nrow = 3, ncol = 3)
  rownames(expected_df) <- c("A", "B", "C")
  expected_df[1, ] <- c(20, 0, 20)
  expected_df[2, ] <- c(0, 0, 40)
  expected_df[3, ] <- c(20, 40, 20)
  colnames(expected_df) <- query
  warnings <- capture_warnings(suppressMessages(get_query_pct(seu = test_pbmc,
      query = missing_query, by = "sample", layer = "counts", min_counts = 0)))
  expect_match(warnings, "NOEXPA, NOEXPB")
  output_df <- suppressWarnings(suppressMessages(get_query_pct(seu = test_pbmc,
      query = missing_query, by = "sample", layer = "counts", min_counts = 0)))
  expect_equal(output_df, expected_df)
})

test_that("get_query_pct works with query of length one", {
  single_query <- "PPBP"
  expected_df <- matrix(nrow = 3, ncol = 1)
  rownames(expected_df) <- c("A", "B", "C")
  expected_df[, 1] <- c(20, 0, 20)
  colnames(expected_df) <- single_query
  output_df <- suppressMessages(get_query_pct(seu = test_pbmc, by = "sample",
                      query = single_query, layer = "counts", min_counts = 0))
  expect_equal(output_df, expected_df)
})

test_that("get_query_pct works with alternate 'by' arguments", {
  single_query <- "PPBP"
  expected_df <- matrix(nrow = 5, ncol = 1)
  rownames(expected_df) <- c(paste0(0:4, ". type", LETTERS[1:5]))
  expected_df[, 1] <- c(0, 33, 33, 0, 0)
  colnames(expected_df) <- single_query
  output_df <- suppressMessages(get_query_pct(seu = test_pbmc, by = "cell_type",
                      query = single_query, layer = "counts", min_counts = 0))
  testthat::expect_equal(output_df, expected_df)
})

test_that("get_query_pct works with 'by' argument of length 2", {
  pbmc_updated <- Seurat::CreateSeuratObject(counts = pbmc_tiny$RNA$counts)
  pbmc_updated@meta.data$sample <- c(rep("A", 5), rep("B", 5), rep("C", 5))
  pbmc_updated@meta.data$seurat_clusters <- rep(0:4, 3)
  dfA <- matrix(data = 0, nrow = 5, ncol = 3)
  colnames(dfA) <- query
  rownames(dfA) <- 0:4
  dfB <- dfC <- dfA
  dfA[, 1] <- c(0, 0, 100, 0, 0)
  dfA[, 3] <- c(0, 0, 0, 100, 0)
  dfB[, 3] <- c(0, 100, 0, 0, 100)
  dfC[, 1] <- c(0, 100, 0, 0, 0)
  dfC[, 2] <- c(100, 0, 100, 0, 0)
  dfC[, 3] <- c(0, 100, 0, 0, 0)
  expected_list <- list(A = dfA, B = dfB, C = dfC)
  output_list <- suppressMessages(get_query_pct(seu = pbmc_updated,
                  query =  query, by = c("sample", "seurat_clusters"),
                  layer = "counts", min_counts = 0))
  expect_equal(output_list, expected_list)
})

test_pbmc <- pbmc_tiny
test_pbmc@meta.data$groups <- "g1"
test_pbmc@meta.data$groups[7:15] <- "g2"
test_pbmc@meta.data$seurat_clusters <- 0
test_pbmc@meta.data$seurat_clusters[3:5] <- 1
test_pbmc@meta.data$seurat_clusters[15] <- 1
test_pbmc@meta.data$cell_types <- test_pbmc@meta.data$seurat_clusters

test_that("get_sc_markers works as intended, DEG TRUE", {
  markers_test_pbmc <- test_pbmc
  markers_test_pbmc$seurat_clusters[11:15] <- 1
  output <- suppressMessages(get_sc_markers(seu = markers_test_pbmc,
    cond = "groups", DEG = TRUE, verbose = FALSE, subset_by = "seurat_clusters",
    simple_DEG_pct = TRUE))$markers
  markers_0 <- data.frame(p_val = 0.27, avg_log2FC = -7.6, pct.1 = 0, pct.2 = 0.5,
                          p_val_adj = 1, gene = "VDAC3")
  rownames(markers_0) <- "VDAC3"
  markers_1 <- data.frame(p_val = c(0.30, 0.33, 0.70, 1),
           avg_log2FC = c(7.1, -9.9, 1.9, -3.5), pct.1 = c(0.33, 0, 0.33, 0.33),
           pct.2 = c(0, 0.4, 0.2, 0.2), p_val_adj = 1,
           gene = c("GNLY", "IGLL5", "PPBP", "VDAC3"))
  rownames(markers_1) <- c("GNLY", "IGLL5", "PPBP", "VDAC3")
  markers_global <- data.frame(p_val = c(0.27, 0.28, 0.41, 0.77),
    avg_log2FC = c(-10, 6.9, -3.9, 1.7), pct.1 = c(0, rep(0.17, 3)),
    pct.2 = c(0.22, 0, 0.33, 0.11), p_val_adj = 1,
    gene = c("IGLL5", "GNLY", "VDAC3", "PPBP"))
  rownames(markers_global) <- c("IGLL5", "GNLY", "VDAC3", "PPBP")
  expected <- list(`0` = markers_0, `1` = markers_1, global = markers_global)
  testthat::expect_equal(output, expected)
})

test_that("get_sc_markers works as intended with simple_DEG_pct set to FALSE", {
  test_pbmc <- test_pbmc
  test_pbmc$seurat_clusters[11:15] <- 1
  output <- suppressMessages(suppressWarnings(get_sc_markers(seu = test_pbmc,
    cond = "groups", DEG = TRUE, verbose = FALSE, subset_by = "seurat_clusters",
    simple_DEG_pct = FALSE)))$markers
  markers_0 <- data.frame(FALSE)
  markers_1 <- data.frame(p_val = c(0.70, 1), avg_log2FC = c(1.9, -3.5),
                          pct.1 = c(0.33, 0.33), pct.2 = c(0.2, 0.2),
                          p_val_adj = 1, gene = c("PPBP", "VDAC3"))
  rownames(markers_1) <- c("PPBP", "VDAC3")
  markers_global <- data.frame(p_val = c(0.41, 0.77),
    avg_log2FC = c(-3.9, 1.7), pct.1 = c(rep(0.17, 2)),
    pct.2 = c(0.33, 0.11), p_val_adj = 1,
    gene = c("VDAC3", "PPBP"))
  rownames(markers_global) <- c("VDAC3", "PPBP")
  expected <- list(`0` = markers_0, `1` = markers_1, global = markers_global)
  testthat::expect_equal(output, expected)
})

test_that("get_sc_markers works as intended, DEG FALSE", {
  markers_test_pbmc <- pbmc_tiny
  markers_test_pbmc$seurat_clusters <- 0
  markers_test_pbmc$seurat_clusters[c(3:5, 6:8)] <- 1
  Seurat::Idents(markers_test_pbmc) <- markers_test_pbmc$seurat_clusters
  output <- suppressMessages(get_sc_markers(seu = markers_test_pbmc,
                 cond = "groups", DEG = FALSE, min.pct = 0, verbose = TRUE,
                 subset_by = "seurat_clusters", simple_DEG_pct = TRUE))$markers
  markers_0 <- data.frame(g1_p_val = c(0.56, 0.85), g1_avg_log2FC = c(5.3, 2.9),
                          g1_pct.1 = c(0.25, 0.50), g1_pct.2 = c(0, 0.33),
                          g1_p_val_adj = 1, g2_p_val = 0.3,
                          g2_avg_log2FC = -c(7.6, 7.1), g2_pct.1 = 0,
                          g2_pct.2 = 0.33, g2_p_val_adj = 1,
                          max_pval = c(0.56, 0.85), minimump_p_val = 0.51,
                          seurat_clusters = "0")
  markers_1 <- data.frame(g1_p_val = c(0.56, 0.85),g1_avg_log2FC = -c(5.3, 2.9),
                          g1_pct.1 = c(0, 0.33), g1_pct.2 = c(0.25, 0.50),
                          g1_p_val_adj = 1, g2_p_val = 0.3,
                          g2_avg_log2FC = c(7.6, 7.1), g2_pct.1 = 0.33,
                          g2_pct.2 = 0, g2_p_val_adj = 1,
                          max_pval = c(0.56, 0.85), minimump_p_val = 0.51,
                          seurat_clusters = "1")
  rownames(markers_0) <- rownames(markers_1) <- c("PPBP", "VDAC3")
  expected <- list(`0` = markers_0, `1` = markers_1)
  testthat::expect_equal(output, expected)
})

test_that("get_sc_markers skips exclusive clusters in DEG analysis", {
  expected_warning <- paste0("Cluster 2 contains fewer than three cells for ",
                "condition\\(s\\) 'g2'. Skipping DEG analysis")
  suppressMessages(expect_warning(get_sc_markers(seu = test_pbmc,
                 cond = "groups", DEG = TRUE, verbose = FALSE,
                 subset_by = "seurat_clusters", simple_DEG_pct = TRUE),
                 expected_warning))
  expect_false(suppressWarnings(suppressMessages(get_sc_markers(seu = test_pbmc,
    cond = "groups", DEG = TRUE, verbose = FALSE, subset_by = "seurat_clusters",
    simple_DEG_pct = TRUE)$markers[[2]][[1]])))
})

test_that("get_sc_markers skips exclusive clusters in DEG analysis, alternate
           idents", {
  expected_warning <- paste0("Cluster 2 contains fewer than three cells for ",
                "condition\\(s\\) 'g2'. Skipping DEG analysis")
  suppressMessages(expect_warning(get_sc_markers(seu = test_pbmc,
                 cond = "groups", DEG = TRUE, verbose = FALSE,
                 subset_by = "cell_types", simple_DEG_pct = TRUE),
                 expected_warning))
  expect_false(suppressWarnings(suppressMessages(get_sc_markers(seu = test_pbmc,
         cond = "groups", DEG = TRUE, verbose = FALSE, subset_by = "cell_types",
         simple_DEG_pct = TRUE)$markers[[2]][[1]])))
})

test_that("get_sc_markers properly applies min_pct filter", {
  local_pbmc <- test_pbmc
  clusters_to_remove <- local_pbmc@meta.data$seurat_clusters == "1"
  local_pbmc@meta.data <- local_pbmc@meta.data[!clusters_to_remove,]
  expected_warning <- paste0("Cluster 2 contains fewer than three cells for ",
                "condition(s) 'g2'")
  suppressMessages(get_sc_markers(seu = local_pbmc, cond = "groups", DEG = TRUE,
    verbose = FALSE, subset_by = "cell_types", min.pct = 0.13,
    simple_DEG_pct = TRUE))
  expect_false(suppressWarnings(suppressMessages(get_sc_markers(seu = test_pbmc,
               cond = "groups", DEG = TRUE, verbose = FALSE,
               subset_by = "cell_types")$markers[[2]][[1]])))
})

test_that("rename_clusters simply assigns names to clusters", {
  test_pbmc$seurat_clusters <- 1:15
  Seurat::Idents(test_pbmc) <- test_pbmc$seurat_clusters
  new_clusters <- paste0("Type", LETTERS[1:15])
  annotated <- suppressWarnings(rename_clusters(test_pbmc, new_clusters))
  output <- annotated@meta.data$cell_type
  expect_equal(as.character(output), new_clusters)
})

test_pbmc <- pbmc_tiny
test_pbmc@meta.data$sample <- c(rep("A", 5), rep("B", 5), rep("C", 5))
test_pbmc@meta.data$seurat_clusters <- c(rep(0, 5), rep (1, 5), rep (2, 5))

test_that("breakdown_query works in simple case", {
  expected_df <- c(0.130, 0.130, 0.270)
  names(expected_df) <- query
  output_df <- signif(breakdown_query(test_pbmc, query, layer = "data",
                                      min_counts = 0), 2)
  expect_equal(output_df, expected_df)
})

test_that("breakdown_query works with counts matrix", {
  expected_df <- c(0.130, 0.130, 0.270)
  names(expected_df) <- query
  counts_matrix <- Seurat::GetAssayData(test_pbmc, layer = "data")
  output_df <- signif(breakdown_query(counts_matrix, query, layer = "data",
                                      min_counts = 0), 2)
  expect_equal(output_df, expected_df)
})

test_that("breakdown_query gives warning if any query genes are not found", {
  missing_query <- c(query, "NOEXPA", "NOEXPB")
  expected_df <- c(0.130, 0.130, 0.270)
  names(expected_df) <- query
  expect_warning(breakdown_query(test_pbmc, missing_query), "NOEXPA, NOEXPB")
  output_df <- suppressWarnings(signif(breakdown_query(test_pbmc,
                            missing_query, layer = "data", min_counts = 0), 2))
  expect_equal(output_df, expected_df)
})

test_that("breakdown_query works with query of length one", {
  single_query <- c("PPBP")
  expected_df <- 0.13
  names(expected_df) <- single_query
  output_df <- suppressMessages(signif(breakdown_query(test_pbmc, single_query,
                                layer = "data", min_counts = 0), 2))
  expect_equal(output_df, expected_df)
})

test_pbmc <- pbmc_tiny
test_pbmc@meta.data$sample <- c(rep("A", 5), rep("B", 5), rep("C", 5))
test_pbmc@meta.data$seurat_clusters <- rep(0:4, 3)
cell_types <- c("0. typeA", "1. typeB", "2. typeC", "3. typeD", "4. typeE")
test_pbmc@meta.data$cell_type <- rep(cell_types, 3)
query <- c("MS4A1", "CD79A", "HLA-DRB5")
counts <- GetAssayData(test_pbmc)
expected_vector <- c("PPBP", "VDAC3", "IGLL5")
expected_big <- c("PPBP", "IGLL5", "VDAC3", "CD1C", "AKR1C3", "PF4", "MYL9",
                  "GNLY", "TREML1", "CA2")
test_pbmc <- subset(test_pbmc, features = rownames(counts))

test_that("get_top_genes works as expected", {
  expect_equal(suppressMessages(get_top_genes(test_pbmc, top = 1,
    min_counts = 0)), expected_vector)
  expect_equal(suppressMessages(get_top_genes(test_pbmc, top = 10,
    min_counts = 0)), expected_big)
  expect_error(suppressMessages(get_top_genes(test_pbmc, top = 0,
    min_counts = 0)), "greater than 1, was 0")
})

test_that("get_top_genes works even if N is greater than number of expressed
       genes", {
  expect_equal(suppressMessages(get_top_genes(test_pbmc, top = 10e99999)),
               expected_big)
})

test_pbmc <- pbmc_tiny
test_pbmc$cell_type <- "g1"
test_pbmc$cell_type[8:15] <- "g2"
DEG_g1 <- data.frame(avg_log2FC = c(0, 0.2, 1, -0.1),
                       p_val_adj = c(0, 0, 0, 0),
                       gene = c("PPBP", "IGLL5", "VDAC3", "GNLY"))
rownames(DEG_g1) <- DEG_g1$gene
DEG_g2 <- data.frame(avg_log2FC = c(1, 0.2, 0.2, -0.4),
                       p_val_adj = c(0, 1, 0, 1),
                       gene = c("PPBP", "IGLL5", "VDAC3", "GNLY"))
rownames(DEG_g2) <- DEG_g2$gene
DEG_list <- list(global = "This should not exist", g1 = DEG_g1, g2 = DEG_g2)

test_that("get_fc_vs_ncells works as intended", {
  output <- get_fc_vs_ncells(seu = test_pbmc, DEG_list = DEG_list,
                             min_avg_log2FC = 0.2, p_val_cutoff = 0.01,
                             layer = "data")
  expected_DEGs <- data.frame(PPBP = c(0, 1), IGLL5 = c(0.2, 0.0),
                              VDAC3 = c(1.0, 0.2))
  expected_ncells <- data.frame(PPBP = c(14, 12), IGLL5 = c(0, 25),
                                VDAC3 = c(29, 25))
  rownames(expected_DEGs) <- c("g1", "g2")
  rownames(expected_ncells) <- c("g1", "g2")
  expected <- list(DEG_df = expected_DEGs, ncell_df = expected_ncells)
  expect_equal(output, expected)
})

test_that("get_fc_vs_ncells works with scaled data (dense matrix)", {
  output <- get_fc_vs_ncells(seu = test_pbmc, DEG_list = DEG_list,
                             min_avg_log2FC = 0.2, p_val_cutoff = 0.01,
                             layer = "scale.data")
  expected_DEGs <- data.frame(PPBP = c(0, 1), IGLL5 = c(0.2, 0.0),
                              VDAC3 = c(1.0, 0.2))
  expected_ncells <- data.frame(matrix(100, 2, 3))
  colnames(expected_ncells) <- c("PPBP", "IGLL5", "VDAC3")
  rownames(expected_DEGs) <- rownames(expected_ncells) <- c("g1", "g2")
  expected <- list(DEG_df = expected_DEGs, ncell_df = expected_ncells)
  expect_equal(output, expected)
})

test_that("get_fc_vs_ncells works with DEG tables containing different genes", {
  test_DEG_list <- DEG_list
  test_DEG_list[[2]] <- test_DEG_list[[2]][1:2, ]
  test_DEG_list[[3]] <- test_DEG_list[[3]][3:4, ]
  output <- get_fc_vs_ncells(seu = test_pbmc, DEG_list = test_DEG_list,
                             min_avg_log2FC = 0.2, p_val_cutoff = 0.01,
                             layer = "data")
  expected_DEGs <- data.frame(IGLL5 = c(0.2, 0.0), VDAC3 = c(0.0, 0.2))
  expected_ncells <- data.frame(IGLL5 = c(0, 25), VDAC3 = c(29, 25))
  rownames(expected_DEGs) <- c("g1", "g2")
  rownames(expected_ncells) <- c("g1", "g2")
  expected <- list(DEG_df = expected_DEGs, ncell_df = expected_ncells)
  expect_equal(output, expected)
})

test_that("get_fc_vs_ncells works with target list", {
  output <- suppressMessages(get_fc_vs_ncells(seu = test_pbmc, layer = "data",
                min_avg_log2FC = Inf, p_val_cutoff = -Inf, DEG_list = DEG_list,
                query = c("PPBP", "VDAC3")))
  expected_DEGs <- data.frame(PPBP = 0:1, VDAC3 = c(1.0, 0.2))
  expected_ncells <- data.frame(PPBP = c(14, 12), VDAC3 = c(29, 25))
  rownames(expected_DEGs) <- c("g1", "g2")
  rownames(expected_ncells) <- c("g1", "g2")
  expected <- list(DEG_df = expected_DEGs, ncell_df = expected_ncells)
  expect_equal(output, expected)
})

test_that("get_fc_vs_ncells can handle no target genes being present in seurat
         object", {
  ## This test could use some work, the function is being ran twice. There
  ## should be a way to capture the warning while still saving output.
  expect_warning(suppressMessages(get_fc_vs_ncells(seu = test_pbmc,
                 DEG_list = DEG_list, min_avg_log2FC = Inf, p_val_cutoff = -Inf,
                 query = c("None", "Zilch"), layer = "data")),
                 "None of the target genes are expressed in seurat object.")
  output <- suppressWarnings(suppressMessages(get_fc_vs_ncells(seu = test_pbmc,
                             DEG_list = DEG_list, min_avg_log2FC = Inf,
                             p_val_cutoff = -Inf,
                             query = c("None", "Zilch"))))
  expect_null(output)
})

test_that("get_fc_vs_ncells can handle some target genes not being present in
         seurat object", {
  expect_warning(suppressMessages(get_fc_vs_ncells(seu = test_pbmc,
                 DEG_list = DEG_list, min_avg_log2FC = Inf, p_val_cutoff = -Inf,
                 query = c("None", "Zilch", "PPBP"),
                 layer = "data")), c('"None", "Zilch"'))
  output <- suppressMessages(suppressWarnings(get_fc_vs_ncells(seu = test_pbmc,
                            DEG_list = DEG_list, min_avg_log2FC = Inf,
                            p_val_cutoff = -Inf, layer = "data",
                            query = c("None", "Zilch", "PPBP"))))
  expected_DEGs <- data.frame(PPBP = 1)
  rownames(expected_DEGs) <- c("g2")
  expected_ncells <- data.frame(PPBP = 12)
  rownames(expected_ncells) <- c("g2")
  expected <- list(DEG_df = expected_DEGs, ncell_df = expected_ncells)
  expect_equal(output, expected)
})

test_that("get_fc_vs_ncells can handle FALSE and NULL values in DEG_list", {
  test_DEG_list <- DEG_list
  test_DEG_list$g3 <- data.frame(FALSE)
  output <- suppressMessages(get_fc_vs_ncells(seu = test_pbmc,
            DEG_list = test_DEG_list, min_avg_log2FC = Inf, p_val_cutoff = -Inf,
            query = "PPBP", layer = "data"))
  expected_DEGs <- data.frame(PPBP = 1)
  rownames(expected_DEGs) <- c("g2")
  expected_ncells <- data.frame(PPBP = 12)
  rownames(expected_ncells) <- c("g2")
  expected <- list(DEG_df = expected_DEGs, ncell_df = expected_ncells)
  expect_equal(output, expected)
})

test_that("get_fc_vs_ncells returns NULL if all non-global values in DEG_list
           are FALSE", {
  test_DEG_list <- DEG_list
  test_DEG_list$g1 <- test_DEG_list$g2 <- test_DEG_list$g3 <- data.frame(FALSE)
  warnings <- capture_warnings(get_fc_vs_ncells(seu = test_pbmc,
                DEG_list = test_DEG_list, min_avg_log2FC = 0, p_val_cutoff = 1,
                layer = "data"))
  expect_match(warnings, "No per-identity DEG analysis present in DEG list.")
  expect_null(suppressWarnings(get_fc_vs_ncells(seu = test_pbmc,
               DEG_list = test_DEG_list, min_avg_log2FC = 0, p_val_cutoff = 1)))
})

test_that("get_fc_vs_ncells returns NULL if input is NULL, regardless of
           global", {
  test_DEG_list <- NULL
  warnings <- capture_warnings(get_fc_vs_ncells(seu = test_pbmc,
    DEG_list = test_DEG_list, min_avg_log2FC = Inf, p_val_cutoff = -Inf,
    layer = "data"))
  expect_match(warnings, "Empty DEG_list provided")
  expect_null(suppressWarnings(get_fc_vs_ncells(seu = test_pbmc,
    DEG_list = test_DEG_list, min_avg_log2FC = Inf, p_val_cutoff = -Inf,
    layer = "data")))
})

test_that("test check_sc_input, integrate and sketch to TRUE, SingleR_ref,
           reduce to FALSE", {
  output <- suppressMessages(check_sc_input(integrate = TRUE, sketch = TRUE,
   fine.tune = TRUE, aggr.ref = FALSE, SingleR_ref = "Non_NULL", reduce = FALSE))
  expected <- list(integrate = FALSE, sketch = FALSE, aggr.ref = FALSE,
                   fine.tune = TRUE)
  expect_equal(output, expected)
})

test_that("test check_sc_input, integrate and sketch to FALSE, SingleR_ref,
           reduce to TRUE", {
  output <- suppressMessages(check_sc_input(integrate = TRUE, sketch = TRUE,
   fine.tune = TRUE, aggr.ref = FALSE, SingleR_ref = "Non_NULL", reduce = TRUE))
  expected <- list(integrate = FALSE, sketch = FALSE, aggr.ref = TRUE,
                   fine.tune = FALSE)
  expect_equal(output, expected)
})

test_that("test check_sc_input, integrate and sketch to TRUE, no SingleR_ref,
           reduce to FALSE", {
  output <- check_sc_input(integrate = TRUE, sketch = TRUE, fine.tune = TRUE,
                           aggr.ref = FALSE, SingleR_ref = NULL, reduce = FALSE)
  expected <- list(integrate = TRUE, sketch = TRUE, aggr.ref = FALSE,
                   fine.tune = TRUE)
  expect_equal(output, expected)
})

test_that("test process_sc_params, annotation mode", {
  params <- list(name = "test", doublet_file = "", filter = TRUE, mincells = 1,
                 minfeats = 1, minfeats = 1, minqcfeats = 1, percentmt = 5,
                 normalmethod = "LogNormalize", scalefactor = 1e5, hvgs = 2e3,
                 ndims = 10, resolution = 0.33, p_adj_cutoff = 1,
                 verbose = FALSE, reduce = FALSE, samples_to_integrate = "",
                 integrate = TRUE, int_method = "RPCA", filter_dataset = "",
                 sketch = TRUE, sketch_pct = 0.25, force_ncells = "",
                 sketch_method = "LeverageScore", extra_columns = "one;two",
                 k_weight = 100, genome = "hg38", exp_design = "",
                 subset_by = "", cpu = "2", input = NULL,
                 suffix = "suffix", imported_counts = "",
                 cluster_annotation = "", cell_annotation = "",
                 SingleR_ref = "SingleR_ref", ref_version = "1",
                 ref_label = "cell_type", ref_de_method = "wilcox", ref_n = 15,
                 ref_filter = "", target_genes = "gene1;gene2", meta_file = "")
  output <- suppressMessages(process_sc_params(params = params,
                                               mode = "annotation"))
  output$opt <- output$opt[order(names(output$opt))]
  expected <- list(opt = params, doublet_list = NULL)
  expected$opt$extra_columns <- c("one", "two")
  expected$opt$top_N <- NA
  expected$out_suffix <- "annotation_report.html"
  expected$opt$target_genes <- ""
  expected$opt$filter_dataset <- NULL
  expected$opt$ref_filter <- NULL
  expected$opt <- expected$opt[order(names(expected$opt))]
  expected$opt$meta_file <- ""
  expected$opt$recalc_query <- FALSE
  output$opt <- output$opt[order(names(output$opt))]
  expected$opt <- expected$opt[order(names(expected$opt))]
  expect_equal(output, expected)
})

test_that("test process_sc_params removes duplicate target genes", {
  params <- list(name = "test", p_adj_cutoff = 1, verbose = FALSE,
                 DE_method = "", subset_by = "genotype;time", cpu = "2",
                 min_avg_log2FC = 1, target_genes = "TREM2;TREM2;DAP12",
                 mincells = 1, top_N = Inf, minpack_common = 0)
  expect_warning(suppressMessages(process_sc_params(params = params,
    mode = "DEG")), "Removing duplicates of gene TREM2 from query.")
  output <- suppressMessages(suppressWarnings(process_sc_params(params = params,
              mode = "DEG")))
  expect_equal(output$opt$target_genes, c("TREM2", "DAP12"))
  })

test_that("test process_sc_params, DEG mode", {
  params <- list(name = "test", p_adj_cutoff = 1, verbose = FALSE,
                 DE_method = "", minpack_common = 0, cpu = "2", mincells = 1,
                 subset_by = "genotype;time", min_avg_log2FC = 1,
                 target_genes = "", top_N = Inf)
  output <- suppressMessages(process_sc_params(params = params, mode = "DEG"))
  expected <- list(opt = params, doublet_list = NULL)
  new_opt <- list(cell_annotation = "", cluster_annotation = "",
             doublet_file = "", exp_design = "", extra_columns = "",
             integrate = FALSE, samples_to_integrate = "")
  expected$opt <- c(expected$opt, new_opt)
  expected$opt$subset_by <- ""
  expected$out_suffix <- "sample_annotation_report.html"
  expected$opt$imported_counts <- ""
  expected$opt$meta_file <- ""
  expected$opt$top_N <- NA
  expected$opt$recalc_query <- FALSE
  expected$opt$DE_method <- "wilcox"
  expected$opt$minpack_common <- 1
  output$opt <- output$opt[order(names(output$opt))]
  expected$opt <- expected$opt[order(names(expected$opt))]
  expect_equal(output, expected)
})

test_that("test process_sc_params, query mode", {
  params <- list(name = "test", p_adj_cutoff = 1, verbose = FALSE,
               subset_by = "genotype;time", cpu = "2", min_avg_log2FC = 1,
               target_genes = "gene1;gene2", mincells = 1, extra_columns = "")
  output <- suppressMessages(process_sc_params(params = params, mode = "query"))
  expected <- list(opt = params, doublet_list = NULL)
  new_opt <- list(cell_annotation = "", cluster_annotation = "",
             doublet_file = "", exp_design = "",
             integrate = FALSE, samples_to_integrate = "")
  expected$opt <- c(expected$opt, new_opt)
  expected$opt$subset_by <- ""
  expected$opt$integrate <- FALSE
  expected$opt$extra_columns <- ""
  expected$out_suffix <- "sample_annotation_report.html"
  expected$opt$target_genes <- c("gene1", "gene2")
  expected$opt$ref_de_method <- NULL
  expected$opt$ref_n <- NULL
  expected$opt$filter_dataset <- NULL
  expected$opt$ref_filter <- NULL
  expected$opt$imported_counts <- ""
  expected$opt$meta_file <- ""
  expected$opt$top_N <- NA
  expected$opt$recalc_query <- FALSE
  output$opt <- output$opt[order(names(output$opt))]
  expected$opt <- expected$opt[order(names(expected$opt))]
  expect_equal(output, expected)
})

test_that("test .add_target_info, simple case", {
  DEG_target <- data.frame(sample = c("A", "B"), treat = c("Ctrl", "Treat"))
  test_pbmc <- pbmc_tiny
  test_pbmc$sample <- c(rep("A", 7), rep("B", 8))
  expected <- c(rep("Ctrl", 7), rep("Treat", 8))
  names(expected) <- colnames(test_pbmc)
  output <- .add_target_info(seu = test_pbmc, DEG_target = DEG_target)$deg_group
  expect_equal(output, expected)
})

test_that("test .add_target_info, complex case", {
  DEG_target <- data.frame(sample = c("A", "B"), treat = c("Treat", "Ctrl"))
  test_pbmc <- pbmc_tiny
  test_pbmc$sample <- c("A", "B", "B", "A", "A", "B", "A", "B", "A", "B", "B",
                        "B", "A", "B", "A")
  expected <- c("Treat", "Ctrl", "Ctrl", "Treat", "Treat", "Ctrl", "Treat",
                "Ctrl", "Treat", "Ctrl", "Ctrl", "Ctrl", "Treat", "Ctrl",
                "Treat")
  names(expected) <- colnames(test_pbmc)
  output <- .add_target_info(seu = test_pbmc, DEG_target = DEG_target)$deg_group
  expect_equal(output, expected)
})

test_that("test .get_matrices, simple case", {
  test_pbmc <- pbmc_tiny
  test_pbmc$cell_type <- "g1"
  test_pbmc$cell_type[8:15] <- "g2"
  test_DEGs <- DEG_list
  test_DEGs$global <- NULL
  m1 <- data.frame(matrix(0, nrow = 3, ncol = 7))
  m2 <- data.frame(matrix(0, nrow = 3, ncol = 8))
  rownames(m1) <- c("PPBP", "IGLL5", "VDAC3")
  colnames(m1) <- c("ATGCCAGAACGACT", "CATGGCCTGTGCAT", "GAACCTGATGAACC",
                    "TGACTGGATTCTCA", "AGTCAGACTGCACA", "TCTGATACACGTGT",
                    "TGGTATCTAAACAG") 
  rownames(m2) <- rownames(m1)
  colnames(m2) <- c("GCAGCTCTGTTTCT", "GATATAACACGCAT", "AATGTTGACAGTCA",
                    "AGGTCATGAGTGTC", "AGAGATGATCTCGC", "GGGTAACTCTAGTG",
                    "CATGAGACACGGGA", "TACGCCACTCCGAA")
  m1[1, 3] <- 4.753095
  m1[3, 4] <- 4.378773
  m1[3, 7] <- 5.057837
  m2[2, 4] <- 5.089387
  m2[3, 3] <- 4.615121
  m2[1, 5] <- 3.976987
  m2[2, 6] <- 7.303943
  m2[3, 5] <- 7.326021
  m1 <- Matrix::Matrix(as.matrix(m1))
  m2 <- Matrix::Matrix(as.matrix(m2))
  matrix_list <- list(g1 = m1, g2 = m2)
  DEG_df <- data.frame(IGLL5 = c(0.2, 0.0), PPBP = c(0, 1), VDAC3 = c(1.0, 0.2))
  rownames(DEG_df) <- c("g1", "g2")
  expected_matrices <- list(DEG_df = DEG_df, matrices = matrix_list)
  matrices <- .get_matrices(seu = test_pbmc, meta = test_pbmc$cell_type,
      layer = "data", DEG_list = test_DEGs, genes = c("PPBP", "IGLL5", "VDAC3"))
  expect_equal(matrices, expected_matrices, tolerance = 6e-8)
})

test_that("test .get_matrices when only one cell is present in one of the
           idents", {
  test_pbmc <- pbmc_tiny
  test_pbmc$cell_type <- "g1"
  test_pbmc$cell_type[15] <- "g2"
  test_DEGs <- DEG_list
  test_DEGs$global <- NULL
  matrices <- suppressWarnings(.get_matrices(seu = test_pbmc, layer = "data",
              meta = test_pbmc$cell_type, DEG_list = test_DEGs,
              genes = c("PPBP", "IGLL5", "VDAC3")))
  expect_null(matrices$matrices$g2)
})

test_that("test .get_union_FCs works with default call", {
  expected <- data.frame(g1 = c(1, 0), g2 = c(0, 1))
  rownames(expected) <- c("VDAC3", "PPBP")
  expected <- as.matrix(expected)
  output <- .get_union_FCs(DEG_list = DEG_list)
  expect_equal(output, expected)
})

test_that("test .get_union_FCs properly applies FC filter", {
  expected <- data.frame(g1 = c(1, 0), g2 = c(0, 1))
  rownames(expected) <- c("VDAC3", "PPBP")
  expected <- as.matrix(expected)
  output <- .get_union_FCs(DEG_list, min_log2FC = 1, p_val_cutoff = 1)
  expect_equal(output, expected)
})

test_that("test .get_union_FCs works when specifying top argument", {
  expected <- data.frame(g1 = c(1, 0), g2 = c(0.2, 1))
  rownames(expected) <- c("VDAC3", "PPBP")
  expected <- as.matrix(expected)
  output <- .get_union_FCs(DEG_list, top = 1, min_log2FC = 0, p_val_cutoff = 1)
  expect_equal(output, expected)
})

test_that("test .get_union_FCs can handle FC and top arguments at once", {
  expected <- data.frame(g1 = c(1, 0), g2 = c(0, 1))
  rownames(expected) <- c("VDAC3", "PPBP")
  expected <- as.matrix(expected)
  output <- .get_union_FCs(DEG_list, min_log2FC = 1, top = 1, p_val_cutoff = 1)
  expect_equal(output, expected)
})

test_that("test .get_union_FCs applies p-value filter", {
  expected <- data.frame(g1 = c(0.0, 0.2, 1.0, -0.1),
                         g2 = c(1.0, 0.0, 0.2, 0.0))
  rownames(expected) <- c("PPBP", "IGLL5", "VDAC3", "GNLY")
  expected <- as.matrix(expected)
  output <- .get_union_FCs(DEG_list = DEG_list, p_val_cutoff = 0,min_log2FC = 0)
  expect_equal(output, expected)
})

test_that("subset_seurat works as intended", {
  output <- subset_seurat(seu = pbmc_tiny, column = "groups", value = "g2")
  expected <- c("ATGCCAGAACGACT", "GAACCTGATGAACC", "TGACTGGATTCTCA",
               "AGTCAGACTGCACA", "AGGTCATGAGTGTC", "GGGTAACTCTAGTG",
               "CATGAGACACGGGA", "TACGCCACTCCGAA")
  expect_equal(colnames(output), expected)
})

test_that("subset_seurat works for alternate operator", {
  output <- subset_seurat(seu = pbmc_tiny, column = "groups", value = "g1",
                             operator = "!=")
  expected <- c("ATGCCAGAACGACT", "GAACCTGATGAACC", "TGACTGGATTCTCA",
               "AGTCAGACTGCACA", "AGGTCATGAGTGTC", "GGGTAACTCTAGTG",
               "CATGAGACACGGGA", "TACGCCACTCCGAA")
  expect_equal(colnames(output), expected)
})

test_that("filter_sc_counts works as intended", {
  expected <- matrix(0, 10, 15)
  expected[2, 13] <- 15
  expected[3, 12] <- 29
  expected <- Matrix::Matrix(expected)
  rownames(expected) <- c("PPBP", "IGLL5", "VDAC3", "CD1C", "AKR1C3", "PF4",
                          "MYL9", "GNLY", "TREML1", "CA2")
  colnames(expected) <- c("ATGCCAGAACGACT", "CATGGCCTGTGCAT", "GAACCTGATGAACC",
                          "TGACTGGATTCTCA", "AGTCAGACTGCACA", "TCTGATACACGTGT",
                          "TGGTATCTAAACAG", "GCAGCTCTGTTTCT", "GATATAACACGCAT",
                          "AATGTTGACAGTCA", "AGGTCATGAGTGTC", "AGAGATGATCTCGC",
                          "GGGTAACTCTAGTG", "CATGAGACACGGGA", "TACGCCACTCCGAA")
  output <- suppressMessages(filter_sc_counts(object = pbmc_tiny,
                             min_counts = 2, layer = "counts"))
  expect_equal(output$RNA$counts, expected)
})

DEG_tags <- DEG_g1
DEG_tags$pct.1 <- c(0, 0, 0.5, 0.5)
DEG_tags$pct.2 <- c(0, 0.5, 0, 0.5)
DEG_tags$p_val_adj <- c(1, 0, 0.5, 0.1)
DEG_tags$avg_log2FC <- c(10, 5, 0, 0.5)
rownames(DEG_tags) <- DEG_tags$gene

test_that("tag_DEGs works with default values", {
  expected <- c("High_P-val,Low_proportion", "Low_proportion",
                   "Low_FC,High_P-val,Low_proportion", "None")
  output <- tag_DEGs(DEG_df = DEG_tags)
  expect_equal(output$cause_for_rejection, expected)
  expect_equal(output$DEG, c(rep(FALSE, 3), TRUE))
})

test_that("tag_DEGs works with lax values", {
  expected <- rep("None", 4)
  output <- tag_DEGs(DEG_df = DEG_tags, p_val_cutoff = 1, min_avg_log2FC = 0,
                     min_cell_proportion = 0)
  expect_equal(output$cause_for_rejection, expected)
  expect_true(all(output$DEG == TRUE))
})

test_that("tag_DEGs works with strict values", {
  expected <- rep("Low_FC,High_P-val,Low_proportion", 4)
  expected[2] <- "Low_FC,Low_proportion"
  output <- tag_DEGs(DEG_df = DEG_tags, p_val_cutoff = 0, min_avg_log2FC = 99,
                     min_cell_proportion = 100)
  expect_equal(output$cause_for_rejection, expected)
  expect_true(all(output$DEG == FALSE))
})

test_that("tag_DEGs works with NULL or FALSE data frames", {
  expect_no_error(tag_DEGs(DEG_df = NULL))
  expect_no_error(tag_DEGs(DEG_df = data.frame(FALSE)))
})

test_that("get_DEG_table works as intended", {
  test_table <- data.frame(matrix(ncol = 5, nrow = 2))
  colnames(test_table) <- c("irrelevant", "irrelevant2", "irrelevant3",
                            "relevant", "relevant2")
  test_table$irrelevant <- c("no", "thanks")
  test_table$irrelevant2 <- c("really", "no")
  test_table$irrelevant3 <- c("don't", "showthis")
  test_table$relevant <- c("yes", "please")
  test_table$relevant2 <- c("thank", "you")
  output <- get_DEG_table(test_table, c("relevant", "relevant2"))
  expect_equal(output, test_table[, c(4, 5)])
})

