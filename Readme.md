## DOCUMENTATION - README

# DEgenes Hunter
> Protocol to undertake the key steps in the analysis of transcriptomics experiments, including differentially expressed gene detection, a co-expression analysis and the functional enrichment.

## REQUIREMENTS

* R version 4

## INSTALLATION

### Install latest R-Version

Go to page https://cran.r-project.org and install the latest R version on your computer.    

### Install DEgenes Hunter from console (using devtools)
To install DEgenes Hunter from console, we use the devtools utility to install R packages from GitHub. Please type the following commands:

``` bash
    Rscript -e 'install.packages("devtools", repos="http://cran.us.r-project.org")'
    Rscript -e 'devtools::install_github("seoanezonjic/ExpHunterSuite", dependencies=TRUE)'
```

Sometimes reactome.db download fails because this package is too large (>400 Mb) and R has a download timeout of 60 seconds by default. 

``` bash
    Rscript -e "getOption('timeout')"
        [1] 60
```

This error can be solved by installing reactome.db from source or extending the timeout threshold before installing the ExpHunterSuite package with the command:

``` bash
    Rscript -e 'options(timeout=1000); devtools::install_github("seoanezonjic/ExpHunterSuite")'
```

Then, please create the folder where you want to install the command line scripts, copy the scripts there and make them command line accesible using these commands:

```bash
    mkdir install_folder
    Rscript -e "ExpHunterSuite::install_DEgenes_hunter('install_folder')"
    export PATH=path_to_install_folder:$PATH
```
This export PATH can also be added to the .bashrc or .bash_profile files.

## Usage

There are two ways for using DEgenes Hunter: as command line scripts or as an R package. We will explain the procedure to perform both 1. differential expression and 2. functional analysis using the command line.

### 1. Command line scripts for differential expression analysis.

Once installed, DEgenes Hunter performs the expression analysis from a raw count table. For this, the user must first create a targets file, including for each *sample* its name in the count table, it *treat* condition (Treat or Ctrl).This file must contain this information separated by tabs. 

    Note: we recommend to use ENSEMBL identifiers for the functional analysis.

Here we include an example in which the targets file must include the samples (CTL, TreatA and TreatB), the samples condition (Ctrl or Treat) and to which age_group they belong (adult or child). 
The correction including additional factors (-v and -M) or co-expression analysis using extra measures (-S and -C) require additional information that must be included in targets file. Extra measures are named as _traits_. These options use the _traits_ column names as arguments.

| sample                                | treat        |  age_group |
|---------------------------------------|--------------|------------|    
|     CTL_1                             |     Ctrl     |    adult   | 
|     CTL_2                             |     Ctrl     |    child   | 
|     TreatA_1                          |     Treat    |    adult   |
|     TreatA_2                          |     Treat    |    child   |
|     TreatB_1                          |     Treat    |    adult   |
|     TreatB_2                          |     Treat    |    adult   | 

Once generated, the expression analysis can be performed using degenes_Hunter.R script. For this, we must call degenes_Hunter.R and give it the following input arguments.

Here we show an example of basic usage:

```bash
    degenes_Hunter.R -t path_to_target_file -i path_to_counts_table -o path_to_results
```

Mandatory arguments:

    -i | -t 
    (mandatory) Specify the path to the input counts/mapping table and to the targets file.

    -i Input read counts file.
    -t Targets file.

Differential expresion analysis arguments:

    -o Output path.
      (optional) Output folder. 
      Default = "./hunter_DE_results"
    -r any integer.
      (optional) Number of minimum mapped reads required in order to not be filtered out. Lesser number of reads are discarded. 
      0 = No filtering. 
      By default, reads less than 2 are discarded.
    -l any integer <= samples provided in the experiment.
      (optional) Minimum number of libraries that must have reads of a transcript in order to not to be filtered. 
      By default, minimum libraries required are 2.
    -p value between 0.01 and 0.1
      (optional) Adjusted p-value for the differential expression analysis. 
      Default = 0.05
    -f float
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
      D = DESeq2, E = edgeR, L = limma, N = NOISeq (NOISeqBIO function within NOISeq package is used), W = WGCNA (this activates the co-expression analysis).
      Default = DELN.
    -c 1-4
      (optional) Minimum number of packages to consider a gene as 'Prevalent' DEG.
      Default = 4.
    -e External DEG data file.
      (optional) External file with pre-analysed DE data. Must consist of three columns containing p-value, logFC and FDR/p-adjust. Please, respect the columns order.
      Default = NULL.
    -v model variables
      (optional) Variables to include in the model. Must be comma separated and each variable must be a column in the target file.
      Default = NULL.
    -M model_text
      Text for this variable will be given directly to the model construction, overwriting the previous configuration.

Co-expression analysis arguments:
    
    -b Any integer.
      Maximum block size value to be given to the WGCNA blockwiseModules function as the maxBlockSize argument.
      Default = 5000.
    --WGCNA_norm_method NOISeq | DESeq2 | edgeR | limma
      Method used to normalized the table of counts for WGCNA. Must also run this method in the --modules argument. Raw counts are used if an empty string is given.
      Default=DESeq2
    --WGCNA_deepsplit 1-4
      This option controls the module building process and is defined as 1,2,3 and 4 values. 1 for rough clustering and 4 for accurate clustering.
      Default = 2.
    --WGCNA_min_genes_cluster integer
      Minimum number of genes to keep a cluster.
      Default = 20.
    --WGCNA_detectcutHeight 0 - 1 float
      Cut height to split modules.
      Default = 0.995.
    --WGCNA_mergecutHeight 0 - 1 float
      Value to merge two similar modules: Maximum dissimilarity (i.e., 1-correlation).
      Default = 0.25.
    -w
      Run WGCNA for treated only, control only, and both as 3 separate runs. Needed if using PCIT. If false, WGCNA runs once, on the table including treament and control.
    --WGCNA_blockwiseNetworkType unsigned | signed | signed hybrid
      NetworkType option to be passed to blockwiseModules function
      Default = signed.
    --WGCNA_blockwiseTOMType none | unsigned | signed | signed Nowick | unsigned 2 | signed 2 | signed Nowick 2
      TOMType option to be passed to blockwiseModules function. 
      Default = signed.
    -S comma sepparated text
      Columns in the target file to be used as categorical factors for the correlation analysis. If more than one to be used, should be comma separated
    -N comma sepparated text
      Columns in the target file to be used as numeric (continuous) factors for the correlation analysis. If more than one is specified, they must separated by commas.

Results files will be included in the output_path:
    
    * DEG_report.html: file that encompass and summarizes all the information provided by the expression analysis. 
    * control_treatment.txt: file that includes information about the samples classification as determined in the targets file.
    * filtered_count_data.txt: filtered counts table used the differential expression analysis. Filtering has been performed according to -r, -l and -F options.
    * opt_input_values.txt: summary of the parameters used for the differential expression analysis

This folder will also include a Common_results folder with a file (table) with all methods used for the differential expression analysis and their logFC, FDR and p-value calculated, the number of DEGs and values for combined_FDR, FDR_labeling, mean_logFCs and genes_tag, and results for the WGCNA analysis (if established): Cluster_ID and Cluster_MM (MM: module membership).

In addition, the results folder will include subfolders generated in accordance to the methods used for the differential expression analysis (results_DESeq2, Results_edgeR, Results_limma, Results_NOISeq). All these folder include two files, one with the normalized counts for all samples and another one with the results given by each package.

In the case of performing the co-expression analysis with WGCNA, it will be created a Results_WGCNA folder including tables with correlations between modules, genes and traits.

#### Non-canonical usage scenarios: 

##### A. Differential expression corrected by extra factors

RNA-seq samples can have several additional attributes in addition to the treatment/control status required to detect differential expression. In the previous example, one of these variables is defined in the targets file as age_group. It is possible to include these attribute in the model design as additional factors to control for. In these cases, the differential expression model can be completed by adding the additional argument -v <factors>, e.g. -v age_group

An example of code for this analysis is the follows:

```bash
    degenes_Hunter.R -t path_to_target_file -i path_to_counts_table -v age_group -o path_to_results
```

##### B. Genes co-expression analysis with WGCNA

Co-expression analysis has been included in DEgenes Hunter to detect gene modules with related biological functions. WGCNA can be activated using -m "W" option. Additional traits can be correlated with module mean profile (use -S for discrete columns and -N for continuous column). Here we show an example of using WGCNA with restrictive options:

```bash
    degenes_Hunter.R -m WDELN -c 4 -f 1 --WGCNA_mergecutHeight 0.1 --WGCNA_min_genes_cluster 15 --WGCNA_detectcutHeight 0.995 -S age_group -t path_to_targets_file -i path_to_counts_table -o path_to_results
```

##### C. Analysing pre-normalized data with WGCNA

DEgenes Hunter requires a table of counts with integers. However, in some situtations, the user may wish to reanalyse a dataset consisting of non-integers, such as microarray data or pre-normalized data.

In this situation, the user can run WGCNA using the data in the table of counts directly, without performing normalization. To do this, they must run degenes_hunter.R with the argument --WGCNA_norm_method equal to "none" and the argument --modules must include "WL", i.e. specify limma is the only algorithm that will accept normalised values. However the DE results will likely not make much sense.

```bash
    degenes_Hunter.R --WGCNA_norm_method none -m WL -c 1 -f 1 -S age_group -t path_to_targets_file -i path_to_normalized_table -o path_to_results
```

##### D. Using a pre-calculated list of DEA genes

In some cases, as well as using pre-normalized count data, we wish to use a precalculated list of DE genes. This can be useful if we want to run functional enrichment but not the DE analysis modules.

To do so, the user must provide a preeanalyzed list of DE analysis results. This should consist of four columns, with the following names and corresponding information: Entrez (or other gene id supported by functional hunter), P.Value, logFC and adj.P.Val, in that order.

In such a scenario, the user can provide a target and counts file, in which case the DE output reports will be generated using this information. They can also choose not to provide them, in which case the DE output report will be rather limited.

To run DEgenes Hunter using a pre-calculated gene list, the following command can be used:

```bash
    degenes_Hunter.R  -m "F" -t path_to_targets_file -i path_to_normalized_table -e path_to_precalculated_deg_file -o path_to_results
```

##### E. Multifactorial (2x2) analysis to look for interactions between factors and effects in distinct groups

Currently only a 2x2 factorial design is possible for interactions, and 2xn for group effects. 

In the case of a 2x2 design, if we consider the following experimental design, similar to the one shown above:

| sample    | treat |  age_group  |
|-----------|-------|-------------|
| ad_CTL_1  | Ctrl  |  ad         |
| ad_CTL_2  | Ctrl  |  ad         |
| ad_CTL_3  | Ctrl  |  ad         |
| ch_CTL_1  | Ctrl  |  ch         |
| ch_CTL_2  | Ctrl  |  ch         |
| ch_CTL_3  | Ctrl  |  ch         |
| ad_TRT_1  | Treat |  ad         |
| ad_TRT_2  | Treat |  ad         |
| ad_TRT_3  | Treat |  ad         |
| ch_TRT_1  | Treat |  ch         |
| ch_TRT_2  | Treat |  ch         |
| ch_TRT_3  | Treat |  ch         |

We can look for an interaction between treatment and age_group, the effects of treatment in a specific age group, or the differences between adults and children among untreated or treated samples.

Interaction can be thought of seeing whether the effect of treatment is different between age groups.

The required contrast (i.e. interaction or effect) must be specified in the following manner, via the flag --multifactorial:

"FactorA,FactorB:contrast"

In the case of an interaction between the factors, the contrast should be specified as "interaction,baseA,baseB", where baseA and baseB should be the base levels for each factor. The resulting logFC values detected by this contrast would represent [numA_numB - baseA_numB] - [numA_baseB - baseA_baseB] with numA/B representing the non-base levels for the factorA.

So, for our example, if one wished to see the interaction between treatment and age_group they should use the following:

```bash
degenes_Hunter.R -m "DEL" -i path_to_normalized_table -t path_to_targets_file -o results_2x2_interaction --multifactorial "treat,age_group:interaction,Ctrl,ch"
```

In the case of the effects of one factor in a group of samples specified by another factor, the contrast should be specified in the form "effect,baseA,groupB, where the baseA should be the level in FactorA that should be used as the base for FC calculation, and groupB represents the level in Factor B that is the group we are looking for the change in.

So, for our example if one wished to see the effect of treatment in children only, they should use the following:

```bash
degenes_Hunter.R -m "DEL" -i path_to_normalized_table -t path_to_targets_file -o results_2x2_effect_treat_ch --multifactorial "treat,age_group:effect,Ctrl,ch"
```

Similarly, if they wished to see the difference between age groups in the ctrl samples, they should use the following:

```bash
degenes_Hunter.R -m "DEL" -i path_to_normalized_table -t path_to_targets_file -o results_2x2_effect_age_ctrl --multifactorial "age_group,treat:effect,ch,Ctrl"
```

Note FactorB in the effects contrast can have more than 2 groups.


### 2. Command line scripts for functional enrichment analysis.

To perform the functional enrichment analysis we will use the functional_Hunter.R script. This tool will use the hypergeometric test to enrich genes in functions and pathways from GO, KEGG and Reactome. Depending on the expression analysis performed, the DEgenes Hunter functional analysis tool, functional_Hunter.R, will execute different enrichments:

    * When differential expression analysis is launched, all prevalent DEGs will be used to perform the functional enrichment. An html summary will be returned.
    * If co-expression analysis is set up, genes from each WGCNA independent module will be used to perform the functional enrichment. An html summary for each module and a global module enrichments summary will be returned.

Here we show an example of basic usage:

```sh
    functional_Hunter.R -i path_to_results -m Organism -o path_to_output
```

Mandatory arguments: 

    -i | -m | -o 
    (required) Specify the path to the degenes_Hunter.R output folder, the model organism to use and the path to the output folder.

    -i path
      Path to the ExpHunterSuite differential expression results.
    -m organism
      Ortologue species to be used as model organism to perform the functional analysis with. Run 'functional_Hunter.R' -L to display all available organisms.
    -o Output path.

Optional input arguments:

    -t input_ID
      Input gene IDs of counts table. Available IDs are: ENSEMBL (E), entrezgene (e), TAIR/Arabidopsis (T), Gene Names (G). 
      Default = E.
    -L 
      (optional) List all organisms provided.
    -a tab_file
      (optional) Path to file with own annotations for functional analysis.
    -f G | g | K | R
      Nomenclature and enrichment method(s) to use (topGO: G = GO | clusterProfiler: K = KEGG, g = GO, R = Reactome).
      Default = gKR.
    -G M | B | C
     Gene Ontology sub-classification to perform functional enrichment.
      M = Molecular Function (MF), B = Biological Process (BP), C = Celular Components (CC)
      Default = MBC.
    -A analysis_type
      Analysis performance (g = Gene Set Enrichment Analysis, o = Over Representation Analysis). 
      Default = go.
    -P float
      Enrichment p-value threshold. 
      Default = 0.1.
    -Q float
      Enrichment q-value threshold. 
      Default = 0.2.
    -c integer
      Cores to be used to parallelize clusters enrichments. 
      Default = 1
    -C files
      Files with custom functional annotation database (in GMT format) separated by commas (,)
    -r mode
       Flags to activate remote query from enrichments and genes translation. Use (b) to launch biomaRt translation; (k) to use KEGG remote database. Requires internet connection.
       Default = NULL
    -q 
    (optional) If indicated, biomaRt query is saved in an .RDS file.

#### DEgenes Hunter functional enrichment examples of use 

Here we show an example of use for DEgenes Hunter functional enrichment, changing some input parameters. 

*Functional enrichment in GO biological processes (-G B) using topGO (-f G) for H. sapiens (-m Human), using a overrepresentation analysis (-A o). P-value threshold set to 0.1 (-P 0.1). ctrl_vs_mut is the input folder with data from the functional expression analysis performed with degenes_Hunter.R (-i). Gene identifiers provided as entrez codes (-t E). Execution parallelized using 6 cores (-c 6).*
```
functional_Hunter.R -f G -G B -A o -P 0.1 -m Human -i ctrl_vs_mut -t E -c 6 -o functional_enrichment
```

### 3. R console pipeline.
Clone this repository and in R console move to the root folder. Then, you can execute the basic pipeline as follows:
```R
library('ExpHunterSuite')

target <- target_generation(from_file='inst/example/target.txt') # Read experiment design which describes the sample groups
raw_count_table <- read.table('inst/example/counts.txt', header=TRUE, row.names=1, sep="\t") # Read the table counts with the number of reads per gene.

final_results <- main_degenes_Hunter( # Perform the expresion analysis with default parameters
  target=target,
  raw=raw_count_table,
  modules='DE' # Use DEseq2 and EdgeR as expresion packages and other complementary analysis as WGCNA (coexpresion) are disabled. 
)

write.table(final_results[['raw_filter']], "filtered_count_data.txt", quote=FALSE, col.names=NA, sep="\t") # Raw table filtered by minreads parameter
write_df_list_as_tables(final_results[['all_data_normalized']], prefix = 'Normalized_counts_') #Normalized table by each expresion package used in modules argument.
write_df_list_as_tables(final_results[['all_counts_for_plotting']], prefix = 'allgenes_')
write.table(final_results[['DE_all_genes']], "hunter_results_table.txt", quote=FALSE, row.names=TRUE, sep="\t") # Table with all the expresion packages and the integrated results (gene_tag, combined p-value and combined log2FC)
write_expression_report(final_results) # Generate friendly html report with expresion data.


func_results <- functional_hunter( #Perform enrichment analisys
        final_results,
        'Mouse', #Use specified organism database
        func_annot_db = "R", #Perform enrichment analysis using only Reactome nomenclature
        analysis_type= "o" #Perform enrichment  using only Overepresentation analysis (Not GSEA)
)

write_enrich_files(func_results) #Write enrichment tables
write_functional_report(final_results, func_results) # Generate friendly html report with functional data.
```


## REFERENCES AND CITATION

Please cite DEgenes Hunter as:

Isabel González Gayte, Rocío Bautista Moreno, Pedro Seoane Zonjic & Manuel G. Claros (2017). DEgenes Hunter - A Flexible R Pipeline for Automated RNA-seq Studies in Organisms without Reference Genome. Genomics And Computational Biology, 3(3), e31. doi:http://dx.doi.org/10.18547/gcb.2017.vol3.iss3.e31

Fernando M. Jabato, José Córdoba-Caballero, Elena Rojano, Carlos Romá-Mateo, Pascual Sanz, Belén Pérez, Diana Gallego, Pedro Seoane, Juan A.G. Ranea y and James R. Perkins. Differential expression, co-expression and functional analysis of RNA-seq data using the DEgenes Hunter suite and its applicability to rare disease (in review). Briefings in Bioinformatics.

