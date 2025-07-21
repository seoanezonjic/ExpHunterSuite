#! /usr/bin/env Rscript

#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>

#############################################
### CONFIGURE 
#############################################

options(warn=1)
if( Sys.getenv('DEGHUNTER_MODE') == 'DEVELOPMENT' ){
  # Loading libraries
  # suppressPackageStartupMessages(require(optparse)) 
  # suppressPackageStartupMessages(require(TCC))

  # Obtain this script directory
  full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile),  
                 error=function(e) # works when using R CMD
                normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', 
                  commandArgs())], '='))[2]))
  main_path_script <- dirname(full.fpath)
  root_path <- file.path(main_path_script, '..', '..')
  # Load custom libraries
  devtools::load_all(root_path)
}else{
  require('ExpHunterSuite')
  root_path <- find.package('ExpHunterSuite')
}

#Loading libraries

option_list <- list(
  optparse::make_option(c("-m", "--mode"), type = "character", default = "compcodeR",
    help = "Function to use to simulate data counts. Possible values:
              \"Jabato\", \"compcodeR\"."),
  optparse::make_option(c("-r", "--replicates"), type="integer", default=3,
    help="Number of replicates for control/treatment group. Default: %default"),
  optparse::make_option(c("-n", "--ngenes"), type="integer", default=20000,
    help="Number of genes. Default : %default"),
  optparse::make_option(c("-d", "--DEGs_proportion"), type="double", 
    default=0.2,
    help=paste0("Proportion of differentially expressed genes (DEGs).",
      " Default: %default")),
  optparse::make_option(c("-f", "--FC_min"), type="double", default=1.3,
    help="Minimum Fold Change value. Default : %default"),
  optparse::make_option(c("-F", "--FC_max"), type="double", default=3.0,
    help="Maximum Fold Change value. Default : %default"),
  optparse::make_option(c("-T", "--P_up"), type="double", default=1,
    help=paste0("Proportion of DEGs that will be up-regulated in treatment",
      " group. Rest will be down-regulated. Default : %default")),
  optparse::make_option(c("-i", "--inputfile"), type="character", 
    default = NULL,
    help="[OPTIONAL] Input genes count matrix to be used for simulation"),
  optparse::make_option(c("-g", "--group"), type="character", 
    default = NULL,
    help=paste0("Columns from input file to be used as same group of ",
      "replicates. Specify indices separated by commas. Note: start at 1")),
  optparse::make_option(c("-o", "--outfile"), type="character",
    help=paste0("Output file basenames. A *_scount and *_predv files will",
      " be generated")),
  optparse::make_option(c("--dataset"), type = "character", default = "B_625_625",
    help = "A name or identifier for the data set/simulation settings."),
  optparse::make_option(c("--n.vars"), type = "integer", default = 12500,
    help = paste("The initial number of genes in the simulated data set. Based",
                 "on the filtering conditions (‘filter.threshold.total’ and",
                 "‘filter.threshold.mediancpm’), the number of genes in the",
                 "final data set may be lower than this number.", sep = "\n")),
  optparse::make_option("--samples.per.cond", type = "integer", default = 5,
    help = "The number of samples in each of the two conditions"),
  optparse::make_option("--n.diffexp", type = "integer", default = NULL,
    help = paste("Min number of cells for which a feature was recorded.",
                 "Overrides DEGs_proportion in compcodeR mode.", sep = "\n")),
  optparse::make_option("--repl.id", type = "integer", default = 1,
    help = paste("A replicate ID for the specific simulation instance. Useful",
                 "for example when generating multiple count matrices with the",
                 "same simulation settings.", sep = "\n")),
  optparse::make_option("--seqdepth", type = "integer", default = 1e+07,
    help = paste("The base sequencing depth (total number of mapped reads).",
                 "This number is multiplied by a value drawn uniformly between",
                 "‘minfact’ and ‘maxfact’ for each sample to generate data with",
                 "different actual sequencing depths.", sep = "\n")),
  optparse::make_option("--minfact", type = "numeric", default = 0.7,
    help = paste("The minimum for the uniform distribution",
                 "used to generate factors that are multiplied with ‘seqdepth’",
                 "to generate individual sequencing depths for the simulated",
                 "samples.", sep = "\n")),
  optparse::make_option("--maxfact", type = "numeric", default = 1.4,
    help = paste("The maximum for the uniform distribution",
                 "used to generate factors that are multiplied with ‘seqdepth’",
                 "to generate individual sequencing depths for the simulated",
                 "samples.", sep = "\n")),
  optparse::make_option("--relmeans", type = "character", default = "auto",
    help = paste("A vector of mean values to use in the simulation of data from",
                 "the Negative Binomial distribution, or ‘\"auto\"’. Note that",
                 "these values may be scaled in order to comply with the given",
                 "sequencing depth. With the default value (‘\"auto\"’), the mean",
                 "values are sampled from values estimated from the Pickrell",
                 "and Cheung data sets. If ‘relmeans’ is a vector, the provided",
                 "values will be used as mean values in the simulation for the",
                 "samples in the first condition. The mean values for the",
                 "samples in the second condition are generated by combining",
                 "the ‘relmeans’ and ‘effect.size’ arguments.", sep = "\n")),
  optparse::make_option("--dispersions", type = "character", default = "auto",
    help = paste("A vector or matrix of dispersions to use in the simulation",
                 "of data from the Negative Binomial distribution, or ‘\"auto\"’.",
                 "With the default value (‘\"auto\"’), the dispersion values are",
                 "sampled from values estimated from the Pickrell and Cheung",
                 "data sets. If both ‘relmeans’ and ‘dispersions’ are set to",
                 "‘\"auto\"’, the means and dispersion values are sampled in",
                 "pairs from the values in these data sets. If ‘dispersions’ is",
                 "a single vector, the provided dispersions will be used for",
                 "simulating data from both conditions. If it is a matrix with",
                 "two columns, the values in the first column are used for",
                 "condition 1, and the values in the second column are used for",
                 "condition 2.", sep = "\n")),
  optparse::make_option("--fraction.upregulated", type = "numeric", default = 1,
    help = paste("The fraction of the differentially expressed",
                   "genes that is upregulated in condition 2 compared to",
                   "condition 1", sep = "\n")),
  optparse::make_option("--between.group.diffdisp", type = "logical", default = FALSE, action = "store_true",
    help = paste("Whether or not the dispersion should be allowed",
                   "to be different between the conditions. Only applicable if",
                   "‘dispersions’ is ‘\"auto\"’.", sep = "\n")),
  optparse::make_option(c("--filter.threshold.total"), type = "numeric", default = 1,
    help = paste("The filter threshold on the total count for a",
                   "gene across all samples. All genes for which the total count",
                   "across all samples is less than the threshold will be",
                   "filtered out.", sep = "\n")),
  optparse::make_option("--filter.threshold.mediancpm", type = "numeric", default = 0,
    help = paste("The filter threshold on the median count",
                 "per million (cpm) for a gene across all samples. All genes",
                 "for which the median cpm across all samples is less than the",
                 "threshold will be filtered out.", sep = "\n")),
  optparse::make_option("--fraction.non.overdispersed", type = "numeric", default = 0,
    help = paste("The fraction of the genes that should be",
                 "simulated according to a Poisson distribution, without",
                 "overdispersion. The non-overdispersed genes will be divided",
                 "proportionally between the upregulated, downregulated and",
                 "non-differentially expressed genes.", sep = "\n")),
  optparse::make_option("--random.outlier.high.prob", type = "numeric", default = 0,
    help = "The fraction of 'random' outliers with unusually high counts."),
  optparse::make_option("--random.outlier.low.prob", type = "numeric", default = 0,
    help = "The fraction of 'random' outliers with unusually low counts."),
  optparse::make_option("--single.outlier.high.prob", type = "numeric", default = 0,
    help = "The fraction of 'single' outliers with unusually high counts."),
  optparse::make_option("--single.outlier.low.prob", type = "numeric", default = 0,
    help = "The fraction of 'single' outliers with unusually low counts."),
  optparse::make_option("--effect.size", type = "numeric", default = 1.5,
    help = paste("The strength of the differential expression, i.e., the",
                 "effect size, between the two conditions. If this is a single",
                 "number, the effect sizes will be obtained by simulating",
                 "numbers from an exponential distribution (with rate 1) and",
                 "adding the results to the ‘effect.size’. For genes that are",
                 "upregulated in the second condition, the mean in the first",
                 "condition is multiplied by the effect size. For genes that",
                 "are downregulated in the second condition, the mean in the",
                 "first condition is divided by the effect size. It is also",
                 "possible to provide a vector of effect sizes (one for each",
                 "gene), which will be used as provided. In this case, the",
                 "‘fraction.upregulated’ and ‘n.diffexp’ arguments will be",
                 "ignored and the values will be derived from the ‘effect.size’",
                 "vector.", sep = "\n")),
  optparse::make_option("--tree", type = "character", default = NULL,
            help = paste("A dated phylogenetic tree of class ‘phylo’ with",
                         "`samples.per.cond * 2` species.")),
  optparse::make_option("--prop.var.tree,", type = "numeric", default = 1,
            help = "Provided CPUs."),
  optparse::make_option("--model.process", type = "character", default = c("BM", "OU"),
            help = "Imported counts directory."),
  optparse::make_option("--selection.strength", type = "numeric", default = 0,
            help = "Cell types annotation file. Will be used to dynamically
                    annotate clusters."),
  optparse::make_option("--id.condition", type = "character", default = NULL,
            help = "Adjusted p-value cutoff."),
  optparse::make_option("--id.species", type = "character", default = NULL,
            help = "Verbosity of base Seurat and harmony function calls."),
  optparse::make_option("--check.id.species", type = "logical", default = TRUE, action = "store_true",
            help = "Randomly subset seurat object to 3000 cells, for quick testing."),
  optparse::make_option("--lengths.relmeans", type = "character", default = NULL,
            help = "Path to reference to use in SingleR annotation."),
  optparse::make_option("--lengths.dispersions", type = "character", default = NULL,
            help = "SingleR reference version."),
  optparse::make_option("--lengths.phylo", type = "logical", default = TRUE, action = "store_true", 
            help = "Column of reference metadata to use for annotation."),
  optparse::make_option("--output_dir", type = "character", default = NULL,
            help = "Directory where results will be written.")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

if(is.null(opt$n.diffexp)) {
  opt$n.diffexp <- opt$n.vars * opt$DEGs_proportion
}
if(is.null(opt$id.species)) {
  opt$id.species <- as.factor(rep(1, 2 * opt$samples.per.cond))
}

#############################################
### LOAD & PREPARE 
#############################################

if(opt$mode == "classic") {
  degsynth(inputfile = opt$inputfile, outfile = opt$outfile, group = opt$group,
  replicates = opt$replicates, ngenes = opt$ngenes, FC_min = opt$FC_min,
  DEGs_proportion = opt$DEGs_proportion, FC_max = opt$FC_max, P_up = opt$P_up)
}
saveRDS(opt, "opt.rds")
if(opt$mode == "compcodeR") {
  generate_synth_DEGs(dataset = opt$dataset, n.vars = opt$n.vars, samples.per.cond = opt$samples.per.cond,
        n.diffexp = opt$n.diffexp, repl.id = opt$repl.id, seqdepth = opt$seqdepth, fraction.upregulated = opt$fraction.upregulated,
        between.group.diffdisp = opt$between.group.diffdisp, filter.threshold.total = opt$filter.threshold.total,
        filter.threshold.mediancpm = opt$filter.threshold.mediancpm, fraction.non.overdispersed = opt$fraction.non.overdispersed,
        output_dir = opt$output_dir)
}

