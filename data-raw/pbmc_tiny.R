## R code to create `pbmc_tiny` dataset

library(SeuratObject)
pbmc_tiny <- pbmc_small[, 1:15]

usethis::use_data(pbmc_tiny, overwrite = TRUE)
