## R code to create `pbmc_tiny` dataset

library(SeuratObject)
pbmc_tiny <- pbmc_small[, 1:15]
gene_list <- c("PPBP", "IGLL5", "VDAC3", "CD1C", "AKR1C3", "PF4", "MYL9", "GNLY", "TREML1", "CA2")
cell_list <- c("ATGCCAGAACGACT", "CATGGCCTGTGCAT", "GAACCTGATGAACC", "TGACTGGATTCTCA",
			   "AGTCAGACTGCACA", "TCTGATACACGTGT", "TGGTATCTAAACAG", "GCAGCTCTGTTTCT",
			   "GATATAACACGCAT", "AATGTTGACAGTCA", "AGGTCATGAGTGTC", "AGAGATGATCTCGC",
			   "GGGTAACTCTAGTG", "CATGAGACACGGGA", "TACGCCACTCCGAA")
pbmc_tiny <- ExpHunterSuite::downsample_seurat(pbmc_tiny, features = gene_list, cells = cell_list)
counts <- pbmc_tiny$RNA@counts
data <- pbmc_tiny$RNA@data
scale.data <- pbmc_tiny$RNA@scale.data
pbmc_tiny <- CreateSeuratObject(counts = counts, assay = "RNA", meta.data = pbmc_tiny@meta.data)
pbmc_tiny$RNA$data <- data
pbmc_tiny$RNA$scale.data <- scale.data
usethis::use_data(pbmc_tiny, overwrite = TRUE)
