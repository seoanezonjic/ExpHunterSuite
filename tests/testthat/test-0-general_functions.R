
test_that("parse_filter with numeric filter", {
	fellowship <<- data.frame(height = c(5.6, 6, 6.4, 6.6, 4.1),
						     race = c("Maia", "Sindar", "Man", "Man", "Hobbit"))
	rownames(fellowship) <- c("Gandalf", "Legolas", "Boromir", "Aragorn",
							  "Frodo")
	expected <- data.frame(height = rep(FALSE, 5))
	expected[5, ] <- TRUE
	rownames(expected) <- rownames(fellowship)
	output <- eval(parse_filter(object = "fellowship",
								expression = "height < 5"))
    expect_equal(as.data.frame(output), expected)
})

test_that("parse_filter with character filter", {
	fellowship <<- data.frame(height = c(5.6, 6, 6.4, 6.6, 4.1),
						     race = c("Maia", "Sindar", "Man", "Man", "Hobbit"))
	rownames(fellowship) <- c("Gandalf", "Legolas", "Boromir", "Aragorn",
							  "Frodo")
	expected <- data.frame(race = rep(FALSE, 5))
	expected[3:4, ] <- TRUE
	rownames(expected) <- rownames(fellowship)
	output <- eval(parse_filter(object = "fellowship",
						   expression = "race == \"Man\""))
	expect_equal(as.data.frame(output), expected)
})

test_that("parse_filter with a seurat object", {
	data(pbmc_tiny)
	output <- eval(parse_filter(object = "pbmc_tiny@meta.data",
						   expression = "nCount_RNA == 2"))
	expected <- data.frame(nCount_RNA = rep(FALSE, 15))
	expected[4, ] <- TRUE
	rownames(expected) <- rownames(pbmc_tiny@meta.data)
	expect_equal(as.data.frame(output), expected)
})

test_that("make_top_n_expression_table", {
	count_data <- data.frame(matrix(0, ncol = 10, nrow = 6))
	rownames(count_data) <- c(paste0("ctrl", 1:3), paste0("treat", 1:3))
	colnames(count_data) <- paste0("gene", 1:10)
	count_data[, paste0("gene", 1:4)] <- 100
	count_data$gene3 <- 300
	output <- make_top_n_expression_table(count_data, 5)
	expected <- count_data[, c(3, 1, 2, 4, 5)]
	expected$gene3 <- 50
	expected[, 2:4] <- 16.667
	expected$gene5 <- 0
	expect_equal(output, as.matrix(expected))
})
