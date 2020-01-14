#! /usr/bin/env Rscript

######################################################################################################
####################### LIBRARY INSTALLING SCRIPT for DEgenes Hunter #################################
######################################################################################################

############## Installing latest version of Bioconductor #####################
install.packages("BiocManager")
BiocManager::install(version = "3.10")
library('BiocManager')

####################  LIBRARIES FOR DIFFERENTIAL EXPRESSION ANALYSIS SCRIPT (degenes_Hunter.R) ############
print("Installing libraries for differential expression analysis script")

## PACKAGES IN CRAN
install.packages(c("FSA","ggplot2","VennDiagram","gplots","optparse","stringr","plyr","reshape2","PerformanceAnalytics"))

## PACKAGES IN BIOCONDUCTOR
BiocManager::install(c("limma", "edgeR", "DESeq2", "NOISeq"))

####################  LIBRARIES required specifically by the FUNCTIONAL ANALYSIS SCRIPT (functional_Hunter.R) ############
print("Installing additional libraries for functional analysis script")

## PACKAGES IN BIOCONDUCTOR
BiocManager::install(c("biomaRt", "topGO", "KEGGREST", "clusterProfiler", "Rgraphviz"))

################## Annotation Databases for functional analysis
## PACKAGES IN BIOCONDUCTOR
BiocManager::install(c("org.Hs.eg.db")
