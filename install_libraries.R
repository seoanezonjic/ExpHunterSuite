#! /usr/bin/env Rscript


######################################################################################################
####################### LIBRARY INSTALLING SCRIPT for DEgenes Hunter #################################
######################################################################################################

############### 
###############
#For using type source("install_libraries.R") in your R console (version R-3.2.3 or higher)


############## Installing latest version of Bioconductor #####################
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite()


####################  LIBRARIES FOR DIFFERENTIAL EXPRESSION ANALYSIS SCRIPT (degenes_Hunter.R) ############

print("Installing libraries for differential expression analysis script")

## PACKAGES IN CRAN
install.packages(c("ggplot2","VennDiagram","gplots","optparse","stringr","plyr"))

## PACKAGES IN BIOCONDUCTOR
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("limma")
biocLite("edgeR")
biocLite("DESeq2")
biocLite("NOISeq")



####################  LIBRARIES required specifically by the FUNCTIONAL ANALYSIS SCRIPT (functional_Hunter.R) ############

print("Installing additional libraries for functional analysis script")

## PACKAGES IN BIOCONDUCTOR
## try http:// if https:// URLs are not supported
## Functional analysis script requires "optparse", "stringr" and "plyr" as well.

source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
biocLite("topGO")
biocLite("KEGGREST")
