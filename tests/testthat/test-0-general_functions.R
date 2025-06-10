
test_that("parse_filter with a seurat object", {
	data(pbmc_tiny)
	output <- parse_filter(object = "pbmc_tiny@meta.data",
						   expression = "nCount_RNA == 2")
	expected <- data.frame(nCount_RNA = rep(FALSE, 15))
	expected[4, ] <- TRUE
	rownames(expected) <- rownames(pbmc_tiny@meta.data)
	expect_equal(as.data.frame(output), expected)
})
