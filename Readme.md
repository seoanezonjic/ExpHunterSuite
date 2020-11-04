## DOCUMENTATION - README

# DEgenes Hunter
> Protocol to undertake the key steps in the analysis of transcriptomics experiments, including differentially expressed gene detection, functional enrichment and the co-expression analysis.

[![NPM Version][npm-image]][npm-url]
[![Build Status][travis-image]][travis-url]
[![Downloads Stats][npm-downloads]][npm-url]

## REQUIREMENTS

* R version 4.

## INSTALLATION

### Install latest R-Version

Go to page https://cran.r-project.org and install the latest R version on your computer.

### Install DEgenes Hunter from console (using devtools)
To install DEgenes Hunter from console, we use the devtools utility to install R packages from GitHub. Please type the following commands: 

``` bash
    Rscript -e 'install.packages("devtools")'
    Rscript -e 'devtools::install_github("seoanezonjic/DEgenesHunter")'
```

Then, please create the folder where you want to install DEgenes Hunter command line scripts, copy the scripts there and make them command-line accesible using these commands: 

```bash
    mkdir install_folder
    Rscript -e "DEgenesHunter::install_DEgenes_hunter('install_folder')"
    export PATH=path_to_install_folder:$PATH
```
This export PATH can also be added to the .bash_rc or .bash_profile files.

## Usage

There are two ways for using DEgenes Hunter: as command line scripts or as an R package. We will first explain the procedure to perform the differential expression and functional analysis using the command line, and then how to use this tool as and R package.

### 1.1 Command line scripts for differential expression analysis.

Once installed, DEgenes Hunter performs the expression analysis from a raw count table. For this, the user must first create a targets file, including for each sample its name in the count table, it condition (if it is a control, mutant, disease, knock-out, etc), and to which group it belongs. This file must contain this information separated by tabs. 

    Note: we recommend to use ENSEMBL identifiers for the functional analysis.

Here we include an example in which the targets file must include the samples (CTL, TreatA and TreatB), the samples condition (Ctrl or Treat) and to which group they belong (ctrl, treatA, treatB):

| sample                                                                | treat               | group                    |
|---------------------------------------------------------------------|-------------------|--------------------------------|
|     CTL_1                                                     |     Ctrl    |     ctrl           |
|     CTL_2                                                     |     Ctrl    |     ctrl           |
|     TreatA_1                                                     |     Treat    |     treatA           |
|     TreatA_2                                                     |     Treat    |     treatA           |
|     TreatB_1                                                     |     Treat    |     treatB           |
|     TreatB_2                                                    |     Treat    |     treatB           |

Once generated, the expression analysis can be performed using degenes_Hunter.R script. For this, we must call degenes_Hunter.R and give it the following input arguments.

Here we show an example of usage:

```bash
    degenes_Hunter.R -p 0.05 -m DELN -c 4 -f 1 -t path_to_targets_file -S group -w FALSE -i path_to_counts_table -o path_to_DEgenesHunter_results
```

Mandatory arguments:

	-i | -t | -S
	(mandatory) Specify the path to the input counts/mapping table and to the targets file.

	-i Input read counts file.
	-t Targets file.
	-S Columns from the targets file to be used as categorical factors for the correlation analysis. If more than one to be used, should be comma separated.

[//]: <> (ERR: Pepe, please confirma si son todos estos los obligatorios)

Optional arguments:

    -o Output path.
	  (optional) Output folder. 
	  Default = "hunter_DE_results"
    -r any integer.
	  (optional) Number of minimum mapped reads required in order to not be filtered out. Lesser number of reads are discarded. 
	  0 = No filtering. 
	  By default, reads less than 2 are discarded.
	-l any integer <= samples provided in the experiment.
	  (optional) Minimum number of mapped reads that must have a transcript in order to not to be filtered. 
	  By default, minimum libraries required are 2.
	-p value between 0.01 and 0.1
	  (optional) Adjusted p-value for the differential expression analysis. 
	  Default = 0.05
	-f float number.
	  (optional) Minimum log2 fold change in expression. Please, consider this is on a log2 scale, so a value of 1 would mean a 2 fold change.
	  Default = 1.
	-q value between 0.95 and 0.99
	  (optional) q value threshold for NOISeqBIO analysis. 
	  Default = 0.95 (recommended)
	-a "BH" | "bonferroni" | "holm" | "hochberg" | "hommel" | "BY"
	  (optional) adjust method for the combined nominal p-values. 
	  By default the BH method is performed.
	-n name of your experiment
	  (optional) Your experiment name. 
	  Default = Experiment1
	-m D | E | L | N | W
	  (optional) Differential expression packages to analyse data with.
	  D = DESeq2, E = edgeR, L = limma, N = NOISeq (NOISeqBIO function within NOISeq package is used), W = WGCNA (for co-expression analysis).
	  Default = DELN.
	 -c Any integer.
	 (optional) Minimum number of packages to consider a gene as 'Prevalent' DEG.
	 Default = 4.
	 -e External DEG data file.
	 (optional) External file with pre-analysed DE data. Must consist of three columns containing p-value, logFC and FDR/p-adjust. Please, respect the columns order.
	 Default = NULL.

[//]: <> (ERR: Me deben faltar unos cuantos parámetros, ¿podrías echarle un vistazo y corregirlos?)

Results files will be included in the output_path:
    
    * DEG_report.html: file that encompass and summarizes all the information provided by the analysis. 
    * control_treatment.txt: file that includes information about the samples classification as determined in the targets file.
    * filtered_count_data.txt: filtered counts table resulting from the differential expression analysis.
    * opt_input_values.txt: summary of the parameters used for the differential expression analysis with DEgenes Hunter.

This folder will also include a Common_results folder with a file (table) with all methods used for the differential expression analysis and their logFC, FDR and p-value calculated, the number of DEGs and values for combined_FDR, FDR_labeling, mean_logFCs and genes_tag, and results for the WGCNA analysis (if established): Cluster_ID and Cluster_MM (MM: module membership).

In addition, the results folder will include subfolders generated in accordance to the methods used for the differential expression analysis (results_DESeq2, Results_edgeR, Results_limma, Results_NOISeq). All these folder include two files, one with the normalized counts for all samples and another one with the results given by each package.

In the case of performing the co-expression analysis with WGCNA, it will be created a Results_WGCNA folder including [...]

[//]: <> (ERR: Pepe, añade aquí lo que falta con respecto al WGCNA, please)


#### Non-canonical usage scenarios: 

##### A. Genes co-expression analysis with WGCNA

Co-expression analysis has been included in DEgenes Hunter to detect gene modules with related biological functions. Here we show an example of use with WGCNA:

```bash
    degenes_Hunter.R -m WDELN -c 4 -f 1 --WGCNA_mergecutHeight 0.1 --WGCNA_min_genes_cluster 15 --WGCNA_detectcutHeight 0.995 -S group -t path_to_targets_file -w FALSE -i path_to_counts_table -o path_to_DEgenesHunter_results
```

WGCNA arguments provided consists of:

    --WGCNA_memory integer.
    Maximum block size value to be passed to the blockwiseModules function of WGCNA as the maxBlockSize argument
    Default=5000.
    --WGCNA_norm_method Normalize counts for WGCNA method.
    Method used to normalized the table of counts for WGCNA. Must also run this method in the --modules argument. 
    Default = "DESeq2".
    --WGCNA_deepsplit integer.
    This option controls the module building process and is defined as 1, 2, 3 and 4 values. 1 for rough clustering and 4 for accurate clustering.
    Default = 2.
    --WGCNA_min_genes_cluster integer
    Minimum number of genes to keep a cluster.
    Default = 20
    --WGCNA_detectcutHeight float number.
    Cut height to split modules.
    Default = 0.995
    --WGCNA_mergecutHeight float number.
    Value to merge two similar modules: Maximum dissimilarity (i.e., 1-correlation).
    Default = 0.25.
    --WGCNA_all
    Run WGCNA for treated only, control only, and both as 3 separate runs. Needed if using PCIT. If false, WGCNA runs once, on the table including treament and control.
    Default = FALSE,
    --WGCNA_blockwiseNetworkType
    NetworkType option to be passed to blockwiseModules function.
    Default = "signed".
    --WGCNA_blockwiseTOMType
    TOMType option to be passed to blockwiseModules function.
    Default = "signed".

##### B. Analysing pre-normalized data with WGCNA

DEgenes Hunter requires a table of counts with integers. However, in some situtations, the user may wish to reanalyse a dataset consisting of non-integers, such as microarray data or pre-normalized data.

In this situation, the user can run WGCNA using the data in the table of counts directly, without performing normalization. To do this, they must run degenes_hunter.R with the argument --WGCNA_norm_method equal to "none" and the argument --modules must include "wl", i.e. specify limmaas it will accept normalised values. However the DE results will likely not make much sense.

[//]: <> (ERR: Esto estaba en el readme previo, pero no sé si se sigue haciendo: me suena "raro" eso del 'limmaas'. Os lo dejo para que lo corrijáis vosotros)

### 1.2 Command line scripts for functional enrichment analysis.

To perform the functional enrichment analysis we will use the functional_Hunter.R script. This tool will use the 

[//]: <> (ERR: Pepe please, te dejo para que escribas esto porque seguro que lo sabes 1000 veces mejor que yo xD)


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


[//]: <> (ERR: Te dejo igualmente que escribas esta parte)

### 2.1 R packages for differential expression analysis.

### 2.2 R packages for functional enrichment analysis.


## REFERENCES AND CITATION

Please cite DEgenes Hunter as:

González Gayte, I., Bautista Moreno, R., Seoane Zonjic, P., & Claros, M. (2017). DEgenes Hunter - A Flexible R Pipeline for Automated RNA-seq Studies in Organisms without Reference Genome. Genomics And Computational Biology, 3(3), e31. doi:http://dx.doi.org/10.18547/gcb.2017.vol3.iss3.e31

Fernando M. Jabato, Jose Cordoba-Caballero, Elena Rojano, Carlos Romá-Mateo, Pascual Sanz, Belén Pérez, Diana Gallego, Pedro Seoane, Juan A.G. Ranea y and James R. Perkins. Differential expression, co-expression and functional analysis of RNA-seq data using the DEgenes Hunter suite and its applicability to rare disease (in preparation). Briefings in Bioinformatics.



[npm-image]: https://img.shields.io/npm/v/datadog-metrics.svg?style=flat-square
[npm-url]: https://npmjs.org/package/datadog-metrics
[npm-downloads]: https://img.shields.io/npm/dm/datadog-metrics.svg?style=flat-square
[travis-image]: https://img.shields.io/travis/dbader/node-datadog-metrics/master.svg?style=flat-square
[travis-url]: https://travis-ci.org/dbader/node-datadog-metrics