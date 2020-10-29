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
  full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile),  # works when using source
                 error=function(e) # works when using R CMD
                normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', commandArgs())], '='))[2]))
  main_path_script <- dirname(full.fpath)
  root_path <- file.path(main_path_script, '..', '..')
  # Load custom libraries
  custom_libraries <- c('simulate_treatment_control.R')
  for (lib in custom_libraries){
    source(file.path(root_path, 'R', lib))
  }
}else{
  require('DEgenesHunter')
  root_path <- find.package('DEgenesHunter')
}




#Loading libraries

option_list <- list(
  optparse::make_option(c("-r", "--replicates"), type="integer", default=3,
    help="Number of replicates for control/treatment group. Default : %default"),
  optparse::make_option(c("-n", "--ngenes"), type="integer", default=20000,
    help="Number of genes. Default : %default"),
  optparse::make_option(c("-d", "--DEGs_proportion"), type="double", default=0.2,
    help="Proportion of differentially expressed genes (DEGs). Default : %default"),
  optparse::make_option(c("-f", "--FC_min"), type="double", default=1.3,
    help="Minimum Fold Change value. Default : %default"),
  optparse::make_option(c("-F", "--FC_max"), type="double", default=3.0,
    help="Maximum Fold Change value. Default : %default"),
  optparse::make_option(c("-T", "--P_up"), type="double", default=1,
    help="Proportion of DEGs that will be up-regulated in treatment group. Rest will be down-regulated. Default : %default"),
  optparse::make_option(c("-i", "--inputfile"), type="character", default = NULL,
    help="[OPTIONAL] Input genes count matrix to be used for simulation"),
  optparse::make_option(c("-g", "--group"), type="character", default = NULL,
    help="Columns from input file to be used as same group of replicates. Specify indices separated by commas. Note: start at 1"),
  optparse::make_option(c("-o", "--outfile"), type="character",
    help="Output file basenames. A *_scount and *_predv files will be generated")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))


#############################################
### LOAD & PREPARE 
#############################################

degsynth(
  inputfile = opt$inputfile,
  outfile = opt$outfile,
  replicates = opt$replicates,
  ngenes = opt$ngenes,
  DEGs_proportion = opt$DEGs_proportion,
  FC_min = opt$FC_min,
  FC_max = opt$FC_max,
  P_up = opt$P_up,
  group = opt$group)