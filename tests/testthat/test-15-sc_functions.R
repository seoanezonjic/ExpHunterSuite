
markers_df <- data.frame(samples = 1:4, cluster = 1:4,
                           gene = toupper(letters[1:4]))
markers_df$p_val_adj <- rep(1e-5, nrow(markers_df))
markers_df$avg_log2FC <- rep(1, nrow(markers_df))

test_that("match_cell_types, simple case", {
  cell_annotation <- data.frame(markers = c("A", "B", "C"),
               type = c("Type1", "Type2", "Type3"))
  expected_df <- markers_df
  expected_df$cell_type <- c("1. Type1", "2. Type2", "3. Type3", "4. Unknown")
  expected_df <- expected_df[order(expected_df$cluster), ]
  output_df <- match_cell_types(markers_df = markers_df,
                  cell_annotation = cell_annotation)$stats_table
  expect_equal(output_df, expected_df)
})

test_that("match_cell_types, unsorted simple case", {
  test_markers_df <- markers_df
  test_markers_df$samples <- 4:1
  test_markers_df$cluster <- c(2, 3, 1, 4)
  test_markers_df$gene <- c("A", "C", "B", "D")
  cell_annotation <- data.frame(markers = c("C", "B", "A"),
                                type = c("Type3", "Type2", "Type1"))
  expected_df <- test_markers_df
  expected_df$cell_type <- c("2. Type1", "3. Type3", "1. Type2", "4. Unknown")
  expected_df <- expected_df[order(expected_df$cluster), ]
  output_df <- match_cell_types(test_markers_df, cell_annotation)$stats_table
  expect_equal(output_df, expected_df)
})

test_that("match_cell_types, complex case", {
  test_markers_df <- data.frame(samples = seq(1, 50),
               cluster = c(rep(1, 20), rep(2, 20), rep (3, 10)))
  test_markers_df$p_val_adj <- rep(1e-5, nrow(test_markers_df))
  test_markers_df$avg_log2FC <- rep(1, nrow(test_markers_df))
  genes <- paste0("gene", seq(1, 50))
  test_markers_df$gene <- genes
  types <- c(rep("type1", 20), rep("type2", 13), rep("type3", 17))
  cell_annotation <- data.frame(markers = genes, type = types)
  expected_df <- test_markers_df
  expected_df$cell_type <-  c(rep("1. type1", 20), rep("2. type2", 20),
                rep("3. type3", 10))
  expected_df <- expected_df[order(expected_df$cluster), ]
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
  test_markers_df$cluster <- 4:1
  genes <- toupper(letters[1:3])
  types <- c("type1", "type2", "type3")
  cell_annotation <- data.frame(markers = genes, type = types)
  expected_df <- test_markers_df[order(test_markers_df$cluster), ]
  expected_df$gene[4] <- "None"
  expected_df$cell_type <- c("1. Unknown (a)", "2. type3",
                             "3. type2", "4. Unknown (b)")
  expect_warning(match_cell_types(test_markers_df, cell_annotation),
                 "WARNING: cluster 4 contains no significant markers")
  output_df <- suppressWarnings(match_cell_types(test_markers_df,
                              cell_annotation))$stats_table
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
  expect_false(suppressWarnings(.has_exclusive_idents(seu = test_pbmc,
                                  idents = "seurat_clusters", cond = "groups")))
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
  output_pct <- get_clusters_distribution(test_pbmc, 2)
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
  output_exp <- get_query_distribution(test_pbmc, query, 3)
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
  output_df <- suppressMessages(get_query_pct(test_pbmc, query, "sample"))
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
  warnings <- capture_warnings(suppressMessages(get_query_pct(test_pbmc,
                               missing_query, "sample")))
  expect_match(warnings, "NOEXPA, NOEXPB")
  output_df <- suppressWarnings(suppressMessages(get_query_pct(test_pbmc,
                               missing_query, "sample")))
  expect_equal(output_df, expected_df)
})

test_that("get_query_pct works with query of length one", {
  single_query <- "PPBP"
  expected_df <- matrix(nrow = 3, ncol = 1)
  rownames(expected_df) <- c("A", "B", "C")
  expected_df[, 1] <- c(20, 0, 20)
  colnames(expected_df) <- single_query
  output_df <- suppressMessages(get_query_pct(test_pbmc, single_query,
                                "sample"))
  expect_equal(output_df, expected_df)
})

test_that("get_query_pct works with alternate 'by' arguments", {
  single_query <- "PPBP"
  expected_df <- matrix(nrow = 5, ncol = 1)
  rownames(expected_df) <- c(paste0(0:4, ". type", toupper(letters[1:5])))
  expected_df[, 1] <- c(0, 33, 33, 0, 0)
  colnames(expected_df) <- single_query
  output_df <- suppressMessages(get_query_pct(test_pbmc, single_query,
                                "cell_type"))
  testthat::expect_equal(output_df, expected_df)
})

# Test disabled until I figure out a way to add a "counts" element to assays
# in a way that lets me access it in the same way that get_query_pct does.
# That method breaks in test dataset, I have not figured out what makes it
# different. Does not have to do with different seurat version.
# test_that("get_query_pct works with 'by' argument of length 2", {
#   pbmc_updated <- Seurat::CreateSeuratObject(counts = pbmc_tiny$RNA$counts)
#   pbmc_updated@meta.data$sample <- c(rep("A", 5), rep("B", 5), rep("C", 5))
#   pbmc_updated@meta.data$seurat_clusters <- rep(0:4, 3)
#   types <- paste0("type", toupper(letters[1:5]))
#   matA <- matrix(data = 0, nrow = 5, ncol = 3)
#   matB = matC <- matA
#   matB[, 3] <- c(0, 100, rep(0, 3))
#   matC[, 1] <- rep(100, 5)
#   matC[, 2] <- c(0, rep(100, 4))
#   matC[, 3] <- c(rep(100, 3), 0, 100)
#   colnames(matA) = colnames(matB) = colnames(matC) <- query
#   rownames(matA) = rownames(matB) = rownames(matC) <- types
#   expected_list <- list(A = matA, B = matB, C = matC)
#   output_list <- get_query_pct(pbmc_updated, query,
#                                by = c("sample", "seurat_clusters"))
#   expect_equal(output_list, expected_list)
# })

test_pbmc <- pbmc_tiny
test_pbmc@meta.data$groups <- "g1"
test_pbmc@meta.data$groups[7:15] <- "g2"
test_pbmc@meta.data$seurat_clusters <- 0
test_pbmc@meta.data$seurat_clusters[3:5] <- 1
test_pbmc@meta.data$seurat_clusters[15] <- 1
test_pbmc@meta.data$cell_types <- test_pbmc@meta.data$seurat_clusters

test_that("get_sc_markers skips exclusive clusters in DEG analysis", {
  expected_warning <- paste0("Cluster 2 contains less than three cells for ",
                "condition 'g2'")
  expect_warning(suppressMessages(get_sc_markers(seu = test_pbmc,
                 cond = "groups", DEG = TRUE, verbose = FALSE,
                 subset_by = "seurat_clusters")), expected_warning)
  expect_false(suppressWarnings(
                suppressMessages(get_sc_markers(seu = test_pbmc,
                cond = "groups", DEG = TRUE, verbose = FALSE,
                subset_by = "seurat_clusters")$markers[[2]][[1]])
                ))
})

test_that("get_sc_markers skips exclusive clusters in DEG analysis, alternate
           idents", {
  expected_warning <- paste0("Cluster 2 contains less than three cells for ",
                "condition 'g2'")
  expect_warning(suppressMessages(get_sc_markers(seu = test_pbmc,
                 cond = "groups", DEG = TRUE, verbose = FALSE,
                 subset_by = "cell_types")), expected_warning)
  expect_false(suppressWarnings(
                suppressMessages(get_sc_markers(seu = test_pbmc,
                cond = "groups", DEG = TRUE, verbose = FALSE,
                subset_by = "cell_types")$markers[[2]][[1]])
                ))
})

test_that("annotate_clusters simply assigns names to clusters", {
  test_pbmc$seurat_clusters <- 1:15
  Seurat::Idents(test_pbmc) <- test_pbmc$seurat_clusters
  new_clusters <- paste0("Type", toupper(letters[1:15]))
  annotated <- suppressWarnings(annotate_clusters(test_pbmc, new_clusters))
  output <- annotated@meta.data$cell_type
  expect_equal(as.character(output), new_clusters)
})

test_pbmc <- pbmc_tiny
test_pbmc@meta.data$sample <- c(rep("A", 5), rep("B", 5), rep("C", 5))
test_pbmc@meta.data$seurat_clusters <- c(rep(0, 5), rep (1, 5), rep (2, 5))

test_that("breakdown_query works in simple case", {
  expected_df <- c(0.130, 0.130, 0.270)
  names(expected_df) <- query
  output_df <- signif(breakdown_query(test_pbmc, query), 2)
  expect_equal(output_df, expected_df)
})

test_that("breakdown_query gives warning if any query genes are not found", {
  library(Seurat)
  missing_query <- c(query, "NOEXPA", "NOEXPB")
  expected_df <- c(0.130, 0.130, 0.270)
  names(expected_df) <- query
  expect_warning(breakdown_query(test_pbmc, missing_query), "NOEXPA, NOEXPB")
  output_df <- suppressWarnings(signif(breakdown_query(test_pbmc,
                                                       missing_query), 2))
  expect_equal(output_df, expected_df)
})

test_that("breakdown_query works with query of length one", {
  single_query <- c("PPBP")
  expected_df <- 0.13
  names(expected_df) <- single_query
  output_df <- signif(breakdown_query(test_pbmc, single_query), 2)
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
  expect_equal(get_top_genes(test_pbmc, top = 1), expected_vector)
  expect_equal(get_top_genes(test_pbmc, top = 10), expected_big)
  expect_error(get_top_genes(test_pbmc, top = 0), "greater than 1, was 0")
})

test_that("get_top_genes works even if N is greater than number of expressed
       genes", {
  expect_equal(get_top_genes(test_pbmc, top = 10e99999), expected_big)
})
