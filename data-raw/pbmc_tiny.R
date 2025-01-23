## R code to create `pbmc_tiny` dataset

library(SeuratObject)
pbmc_tiny <- pbmc_small[, 1:15]
gene_list <- c("PPBP", "IGLL5", "VDAC3", "CD1C", "AKR1C3", "PF4", "MYL9", "GNLY", "TREML1", "CA2")
pbmc_tiny <- downsample_seurat(pbmc_tiny, features = gene_list)

usethis::use_data(pbmc_tiny, overwrite = TRUE)
