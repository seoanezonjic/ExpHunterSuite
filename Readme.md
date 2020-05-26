## DOCUMENTATION - README


# DEgenes Hunter(R) Version 1.0 
> RNA-seq Analysis Pipeline for customized differential gene expression and functional analysis

[![NPM Version][npm-image]][npm-url]
[![Build Status][travis-image]][travis-url]
[![Downloads Stats][npm-downloads]][npm-url]


## REQUIREMENTS

* R version R-3.3.1 and Bioconductor 3.4 or higher
* LaTeX 	


Included:
	
	* Main program scripts (degenes_Hunter.R and functional_Hunter.R)
	* All functional libraries:
		* general_functions.R
		* dif_expression_packages.R
		* qc_and_benchmarking_functions.R
		* functional_analysis_library.R
	* Table with all organisms provided: organism_table.txt
	* This README
	* Example data set (example_count_data.txt)
	* Example output of a differential analysis report dumped out by DEgenes Hunter (example_DE_report.pdf)
	* Example output of a functional analysis report dumped out by DEgenes Hunter (example_functional_report.pdf)



## Installation


Install latest R-Version

Go to page https://cloud.r-project.org/ and install the latest R version on your computer
Install also the latest Bioconductor version in http://bioconductor.org/install/


## Install DEgenes Hunter and required R-packages


To download DEgenes Hunter: https://github.com/Isabelggayte/DEgenesHunter

For installing the latest versions of all R-packages required to run DEgenes Hunter, use the install_libraries.R script which is contained in the main DEgenes_Hunter folder. 



## Make DEgenesHUnter path-accesible

Modify your .bashrc or .profile files to:

```sh
PATH=~path/to/DEgenesHunter:$PATH
export PATH

Then reload your terminal session or execute:

source ~/.bashrc

OR

source ~/.profile
```



## RUN DEgenes Hunter


DIFFERENTIAL EXPRESSION ANALYSIS with the degenes_Hunter.R script



## General Usage Notes

R script to perform differential expression analysis on RNA-seq count data. 



## Command Line Arguments


Launching example: 

```sh
	degenes_Hunter.R -i path/to/mapping_table -C G1_rep1,G1_rep2,G1_rep3 -T G2_rep1,G2_rep2,G2_rep3 -o path/to/output
```


	-i | -C | -T | -o 
	(required) Specify the path to the input counts/mapping table, names of control 
	and treatment columns and the path to the output folder

	-i - Input file with read counts
	-C - Columns considered as control samples in the count table provided with -i. 
		 Please indicate column names of control samples separated by commas
	-T - Columns considered as treatment samples in the count table provided with -i. 
		 Please indicate column names of treatment samples separated by commas

	-o - Output path
	  (optional) Output folder. Default = "hunter_DE_results"
	-r 0 | any whole number
	  (optional) Number of minimum mapped reads required in order to not be filtered out. Lesser number of reads are discarded. -r 0 = No filtering. 
	  By default, reads less than 2 are discarded.
	-l any whole number <= samples provided in the experiment.
	  (optional) Minimum number of mapped reads that must have a transcript in order to not to be filtered 
	  By default, minimum libraries required are 2.
	-p value between 0.01 and 0.1
	  (optional) Adjusted p-value for the differential expression analysis. Default = 0.05
	-f value between 1.5 and 2
	  (optional) Fold Change Value threshold. Default = 1.5
	-q value between 0.95 and 0.99
	  (optional) q value threshold for NOISeqBIO analysis. Default = 0.95 (recommended)
	-a "BH" | "bonferroni" | "holm" | "hochberg" | "hommel" | "BY"
	  (optional) adjust method for the combined nominal p-values. By default the BH method is performed.
	-n name of your experiment
	  (optional) Your experiment name. Default = Experiment1
	-m D | E | L | N
	  (optional) Differential expression packages to analyse data with.
	  D = DESeq2, E = edgeR, L = limma, N = NOISeq (NOISeqBIO function within NOISeq package is used)
	  Default = DELN.




## Output Files


Output folders tree structure:
	
	Main Folder (the folders' name is set with option -o)
		* boxplot_before_normalization.pdf
		* boxplot_normalized_data.pdf
		* group_dendogram.pdf
		* group_dendrogram_normalized.pdf
		* filtered_count_data.txt
		* genenumbers.pdf
		* statistics_report.txt
		* DE_report.pdf
		* Functional_analysis_report.pdf
	Subfolders:
		* Results_DESeq2
			* MA_plot_DESeq2.pdf
			* Normalized_counts_DESeq2.txt
			* DEgenes_DESeq2.txt
			* allgenes_DESeq2.txt
			* PCAplot.pdf
		* Results_edgeR
			* MA_plot_edgeR.pdf
			* Normalized_counts_edgeR.txt
			* DEgenes_edgeR.txt
			* allgenes_edgeR.txt
			* MDSplot.pdf
			* MDSplot_norm.pdf
		* Results_limma
			* Volcanoplot_limma.pdf
			* Normalized_counts_limma.txt
			* DEgenes_limma.txt
			* allgenes_limma.txt
		* Results_NOISeq
			* Expressionplot_NOISeq.pdf
			* Normalized_counts_NOISeq.txt
			* DEgenes_NOISeq.txt
			* allgenes_NOISeq.txt
		* Common_Results
			* VennDiagram.pdf
			* hunter_results_file.txt
			* Prevalent_geneIDs.txt
			* pos_prevalentDEGs_logFCs.txt
			* neg_prevalentDEGs_logFCs.txt
			* padj_possible_DEGs.pdf
			* padj_prevalent_DEGs.pdf
			* padj_all_genes.pdf
			* top20_genes.txt





## FUNCTIONAL ANALYSIS with the functional_Hunter.R script


## General Usage Notes

R script to perform functional analysis on the degenes_hunter.R output file. 




## Command Line Arguments


Launching example: 

```sh
	functional_Hunter.R -i path/to/complete_genes_statistics.txt -m Grapevine -t E -o path/to/output
```

	-i | -m | -t | -o 
	(required) Specify the path to the degenes_Hunter.R output file "complete_genes_statistics.txt", the model organism, the type of gene identifier and the path to the output folder

	-i - Path to the DEgenes Hunter's differential expression 
		 analysis output file "hunter_results_table.txt"
	-m - Ortologue species to be used as model organism to perform the functional analysis with.
	-t E | R
		 Gene ID provided. E = ENSEMBLE gene ID, R = REFSEQ peptide. Default = E.
	-o - Output path

	-L (optional) List all organisms provided.
	-a (optional) Path to file for providing own annotations for functional analysis.
	-f G | K 
	   (optional) Functional analysis choice.
	   G = Gene Ontology Enrichment (GOs), K = Pathway enrichment (KEGG)
	   Default = GK.
	-G M | B | C
	  (optional) Kinds of gene enrichment analysis to perform.
	  M = Molecular Function (MF), B = Biological Process (BP), C = Celular Components (CC)
	  Default = MBC.
	-K (optional) Ortologue species to be used to perform the pathway enrichment analysis in case
	  the model organism indicated in -m is not provided in the KEGG database.
	-q (optional) If indicated, biomaRt query is saved in an .RDS file.




## Output Files



Output folders tree structure:
	
	Main output Folder (the folders' name is set with option -o)
	* Functional_Hunter_Report.pdf

	Subfolders:
		* topGO_maps
			* GOgraph_allpos_(+).pdf
			* GOgraph_allpos_overex_(+).pdf
			* GOgraph_allpos_underex_(+).pdf
			* GOgraph_preval_(+).pdf
			* GOgraph_preval_overex_(+).pdf
			* GOgraph_preval_underex_(+).pdf

			(+) Type of enrichment analysis (MF, BP, or CC)

		* KEGG_pathways
			* KEGG_paths.html


## META

[https://github.com/Isabelggayte/github-link](https://github.com/DEgenes_Hunter/)


[npm-image]: https://img.shields.io/npm/v/datadog-metrics.svg?style=flat-square
[npm-url]: https://npmjs.org/package/datadog-metrics
[npm-downloads]: https://img.shields.io/npm/dm/datadog-metrics.svg?style=flat-square
[travis-image]: https://img.shields.io/travis/dbader/node-datadog-metrics/master.svg?style=flat-square
[travis-url]: https://travis-ci.org/dbader/node-datadog-metrics



## Publication

DEgenes Hunter was published in:

González Gayte, I., Bautista Moreno, R., Seoane Zonjic, P., & Claros, M. (2017). DEgenes Hunter - A Flexible R Pipeline for Automated RNA-seq Studies in Organisms without Reference Genome. Genomics And Computational Biology, 3(3), e31. doi:http://dx.doi.org/10.18547/gcb.2017.vol3.iss3.e31



Note to DEgenes Hunter publication:

Bootstraps with fold changes 1.5 and 2.5 were made with a sample size of 50.000 genes. The sample size employed does not affect in any way the results obtained in the paper.