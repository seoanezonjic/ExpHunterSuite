#! /usr/bin/env Rscript

######################################################################################################
####################### LIBRARY INSTALLING SCRIPT for DEgenes Hunter #################################
######################################################################################################

install.packages("BiocManager", repos='https://cloud.r-project.org')
print("Installing libraries from Bioconductor")
packages_list_biocond <- c("KEGG.db", "ReactomePA","preprocessCore","impute","limma", "edgeR", "DESeq2", "NOISeq", "biomaRt", "topGO", "KEGGREST", "clusterProfiler", "Rgraphviz", "org.Hs.eg.db", "diffcoexp")
BiocManager::install()
installed <- library()$results[,1]
packages_list_biocond <- setdiff(packages_list_biocond, installed)
BiocManager::install(packages_list_biocond)


print("Installing libraries from CRAN")
packages_list <-c("DT","FSA","ggplot2","VennDiagram","gplots","optparse","stringr","plyr","reshape2","PerformanceAnalytics", "PCIT")
installed <- library()$results[,1]
packages_list <- setdiff(packages_list, installed)
install.packages(packages_list, repos='https://cloud.r-project.org')

