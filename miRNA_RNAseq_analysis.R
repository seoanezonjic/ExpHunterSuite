#!/usr/bin/env Rscript

################### INITILAIZE
options(scipen = 1) 
options(digits = 2)
suppressPackageStartupMessages(library(optparse)) 

 
 #### NOTAS PARA EL SCRIPT

full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile),  # works when using source
               error=function(e) # works when using R CMD
              normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', commandArgs())], '='))[2]))
main_path_script <- dirname(full.fpath)

################### OPTIONS
option_list <- list(
	make_option(c("-r", "--RNAseq_folder"), type="character", default=NULL,
		help="DEgenesHunter RNAseq execution folder"),
	make_option(c("-m", "--miRNAseq_folder"), type="character", default=NULL,
		help="DEgenesHunter miRNAseq execution folder"),
	make_option(c("-d", "--deg_tag"), type="character", default="PREVALENT_DEG",
		help="DEG tag. Set the DEG tag of DEGenesHunter to filter results. 'PREVALENT_DEG' or 'POSSIBLE_DEG'. Default : %default"),
	make_option(c("-o", "--output_files"), type="character", default=".", 
		help = "Output folder"),
	make_option(c("-a", "--aproaches"), type="character", default="EE,Eh,Ed,hd,hE",
		help = "EE = Anticorrelation between RNAseq modules Eigengene and miRNAseq modules Eigengene, Eh = Anticorrelation between RNAseq modules Eigengene and miRNAseq hub genes, Ed = Anticorrelation between RNAseq modules Eigengene and miRNAseq DEG expression profiles, hd = Anticorrelation between RNAseq hub genes and miRNAseq DEG expression profiles, hE = Anticorrelation between RNAseq hub genes and miRNAseq modules Eigengene. Default : %default"),
	make_option(c("--organism"), type ="character", default="hsa",
		help= "Reference organism to use. Available 'hsa' for human and 'mmu' for mouse."),
	make_option(c("-M", "--multimir_db"), type ="character", default=NULL,
		help= "Indicate multimir db download path which include 'organism'.RData file."),
	make_option(c("--only_known"), type="logical", default=FALSE, action = "store_true",
		help= "Perform analysis only using known miRNAs. Faster but you lose information of novel miRNAs. Only available when a translation file is given"),
	make_option(c("-p", "--p_val_cutoff"), type="double", default=0.05,
    	help="Correlation P value threshold . Default=%default"),
	make_option(c("-c", "--corr_cutoff"), type="double", default=-0.7,
    	help="Correlation threshold . Default=%default"),
	make_option(c("-R", "--report"), ,type = "character", default="miRNA_RNA_comparison.html",
	    help="Name of the html file. Default : %default"),
	make_option(c("-t", "--translation_file"), type = "character", default = NULL,
		help = 'Two columns (\t) file with miRNA names translation: Same IDs as imput on first column, and miRNA mature ID on second column (no need to header)')
)

suppressPackageStartupMessages(library(pryr))
print(mem_used())
opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(require(rmarkdown))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(WGCNA))
suppressPackageStartupMessages(require(ggplot2))  
suppressPackageStartupMessages(require(stringr))  
suppressPackageStartupMessages(require(data.table))  
suppressPackageStartupMessages(require(tidyverse))  
suppressPackageStartupMessages(require(DT))  
suppressPackageStartupMessages(require(miRBaseConverter))  
suppressPackageStartupMessages(require(gridExtra))

source(file.path(main_path_script, 'lib', 'plotting_functions.R')) 
source(file.path(main_path_script, 'lib', 'miRNA_RNA_functions.R'))
source(file.path(main_path_script, 'lib', 'functional_analysis_library.R'))

dir.create(opt$output_files, recursive = T)

time_control <- Sys.time()
opt$aproaches <- unlist(strsplit(opt$aproaches, ","))
print("packages")
print(mem_used())
translation_miRNA <- opt$translation_file
if (!is.null(translation_miRNA)) {
	translation_miRNA <- read.table(translation_miRNA, sep = "\t")
	colnames(translation_miRNA) <- c("mature_ID", "input_ID")
}

RNAseq <- load_DEGH_information(opt$RNAseq_folder)
miRNAseq <- load_DEGH_information(opt$miRNAseq_folder, translation_table =  translation_miRNA)
organism_info <- read.table(file.path(main_path_script, "lib", "organism_table.txt"), header = TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE, fill = NA)
organism_info <- subset(organism_info, organism_info$KeggCode %in% opt$organism)  




print("files")
print(mem_used())

if (!is.null(opt$multimir_db)){
	multimir_path <- file.path(opt$multimir_db, paste0("parsed_", opt$organism, ".RData"))
	load(multimir_path)
	multimir_summary <- summarize_multimir(multimir_summary)
	message("multiMiR database has been parsed and summarized")
	print(gc())
print(mem_used())
}

strategies <- list()
strategies[["zero"]] <- correlate_profiles(RNAseq[["normalized_counts"]], miRNAseq[["normalized_counts"]])
strategies[["zero"]]$plot_obj <- corr2plot(strategies[["zero"]])
strategies[["zero"]][["name"]] <- "counts_RNA_vs_counts_miRNA" 
strategies[["zero"]]$plot_obj <- add_multimir_info(strategies[["zero"]]$plot_obj, expand_miRNA = NULL , expand_RNA = NULL, multimir_summary = multimir_summary)
print(paste0("zero"," strategy has finished"))

print(gc())
print(mem_used())

if("EE" %in% opt$aproaches){
# Aproach 1: Anticorrelation between RNAseq modules Eigengene and miRNAseq modules Eigengene
# ###### RNAseq modules as rows and miRNA modules as columns
	strategies[["EE"]] <- correlate_profiles(RNAseq[["Eigengene"]], miRNAseq[["Eigengene"]])
	strategies[["EE"]]$plot_obj <- corr2plot(strategies[["EE"]])
	strategies[["EE"]][["name"]] <- "Eigen_RNA_v_Eigen_miRNA" 
	strategies[["EE"]]$plot_obj <- add_multimir_info(strategies[["EE"]]$plot_obj, expand_miRNA = miRNAseq[['DH_results']] , expand_RNA = RNAseq[['DH_results']], multimir_summary = multimir_summary)
	print(paste0("EE"," strategy has finished"))

	print(gc())
	print(mem_used()) 
}

if("Eh" %in% opt$aproaches){
# # Aproach 2: Anticorrelation between RNAseq modules Eigengene and miRNAseq hub genes 
# ###### RNAseq modules as rows and miRNA hub genes as columns
	miRNAseq$hub_1 <- get_hub_genes_by_MM(miRNAseq[["normalized_counts"]], miRNAseq[["DH_results"]], top = 1)
	strategies[["Eh"]] <- correlate_profiles(RNAseq[["Eigengene"]], miRNAseq[["hub_1"]])
	strategies[["Eh"]]$plot_obj <- corr2plot(strategies[["Eh"]])
	strategies[["Eh"]][["name"]] <- "Eigen_RNA_v_Hub1_miRNA" 
	strategies[["Eh"]]$plot_obj <- add_multimir_info(strategies[["Eh"]]$plot_obj, expand_miRNA = miRNAseq[['DH_results']], expand_RNA = RNAseq[['DH_results']], multimir_summary = multimir_summary) 
	print(paste0("Eh"," strategy has finished"))

	print(gc())
	print(mem_used()) 
}

if("Ed" %in% opt$aproaches){
# Aproach 3: Anticorrelation between RNAseq modules Eigengene and miRNAseq expression profiles
###### RNAseq modules as rows and miRNA genes as columns
	strategies[["Ed"]] <- correlate_profiles(RNAseq[["Eigengene"]], miRNAseq[["normalized_counts"]])
	strategies[["Ed"]]$plot_obj <- corr2plot(strategies[["Ed"]])
	strategies[["Ed"]][["name"]] <- "Eigen_RNA_v_counts_miRNA" 
	strategies[["Ed"]]$plot_obj <- add_multimir_info(strategies[["Ed"]]$plot_obj, expand_miRNA = NULL , expand_RNA = RNAseq[['DH_results']], multimir_summary = multimir_summary) 
	print(paste0("Ed"," strategy has finished"))

	print(gc())
	print(mem_used()) 

}

if("hd" %in% opt$aproaches){
# Aproach 4: Anticorrelation between RNAseq hub genes and miRNAseq expression profiles
	RNAseq$hub_1 <- get_hub_genes_by_MM(RNAseq[["normalized_counts"]], RNAseq[["DH_results"]], top = 1)
	strategies[["hd"]] <- correlate_profiles(RNAseq[["hub_1"]], miRNAseq[["normalized_counts"]])
	strategies[["hd"]]$plot_obj <- corr2plot(strategies[["hd"]])
	strategies[["hd"]][["name"]] <- "Hub1_RNA_v_counts_miRNA" 
	strategies[["hd"]]$plot_obj <- add_multimir_info(strategies[["hd"]]$plot_obj, expand_miRNA = NULL , expand_RNA = RNAseq[['DH_results']], multimir_summary = multimir_summary) 
	print(paste0("hd"," strategy has finished"))

	print(gc())
	print(mem_used()) 
}

if("hE" %in% opt$aproaches){
# Aproach 5: Correlation between RNAseq hub genes and miRNAseq modules Eigengene
	RNAseq$hub_1 <- get_hub_genes_by_MM(RNAseq[["normalized_counts"]], RNAseq[["DH_results"]], top = 1)
	strategies[["hE"]] <- correlate_profiles(RNAseq[["hub_1"]], miRNAseq[["Eigengene"]])
	strategies[["hE"]]$plot_obj <- corr2plot(strategies[["hE"]])
	strategies[["hE"]][["name"]] <- "Hub1_RNA_v_Eigen_miRNA" 
	strategies[["hE"]]$plot_obj <- add_multimir_info(strategies[["hE"]]$plot_obj, expand_miRNA = miRNAseq[['DH_results']] , expand_RNA = RNAseq[['DH_results']], multimir_summary = multimir_summary) 
	print(paste0("hE"," strategy has finished"))
	
	print(gc())
	print(mem_used()) 
}

all_strategies <- list()
for (strategy in names(strategies)) {
	strategies[[strategy]]$plot_obj[is.na(strategies[[strategy]]$plot_obj$predicted_c), "predicted_c"] <- 0
	strategies[[strategy]]$plot_obj[is.na(strategies[[strategy]]$plot_obj$validated_c), "validated_c"] <- 0
	all_strategies[[strategies[[strategy]][["name"]]]] <- as.data.table(strategies[[strategy]]$plot_obj)[correlation <= opt$corr_cutoff & pval <= opt$p_val_cutoff]
}

all_strategies <- as.data.frame(rbindlist(all_strategies, use.names=TRUE, idcol="strategy"))

filters_summary <- all_strategies %>%
					group_by(strategy) %>%
					summarize(sig_correlated = length(RNAseq),
								predicted = sum(predicted_c > 0),
								validated = sum(validated_c > 0), 
								both = sum(predicted_c > 0 & validated_c > 0))
filters_summary <- as.data.frame(filters_summary)
filters_summary <- reshape2::melt(filters_summary)
names(filters_summary) <- c("strategy", "type", "pairs")

################# GET IDS TARNSALTIONS
mirna_names <- unique(all_strategies$miRNAseq)
mirna_names <- mirna_names[grepl("MIMAT", mirna_names)]
mirna_names <- miRNA_AccessionToName(mirna_names, targetVersion = "v22")

gene_id_translation <- NULL
if (! organism_info$Bioconductor_VarName_SYMBOL[1] == "") {

	input_to_entrezgene <- id_translation_orgdb(input_ids = unique(all_strategies$RNAseq),
												 organism_db = organism_info$Bioconductor_DB[1],
												 org_var_name = organism_info$Bioconductor_VarName[1])

	colnames(input_to_entrezgene) <- c("ensembl_gene_id", "entrezgene") # Fix names


	input_to_symbol <- id_translation_orgdb(input_ids = input_to_entrezgene$entrezgene,
										 		organism_db = organism_info$Bioconductor_DB[1],
												org_var_name = organism_info$Bioconductor_VarName_SYMBOL[1])
	colnames(input_to_symbol) <- c("entrezgene", "Symbol")
	input_to_entrezgene <- as.data.table(input_to_entrezgene)
	input_to_symbol <- as.data.table(input_to_symbol)
	gene_id_translation <- merge(input_to_entrezgene, input_to_symbol, by = "entrezgene", all.x = TRUE)

}

rmarkdown::render(file.path(main_path_script, 'templates', 'miRNA_RNA.Rmd'), 
                  output_file = file.path(opt$output_files, opt$report), intermediates_dir = opt$output_files)