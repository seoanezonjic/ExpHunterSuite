#! /usr/bin/env Rscript

#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>

#############################################
### CONFIGURE 
#############################################

options(warn=1)
if( Sys.getenv('DEGHUNTER_MODE') == 'DEVELOPMENT' ){
  # Loading libraries
  suppressPackageStartupMessages(require(optparse)) 
  suppressPackageStartupMessages(require(TCC))

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
### FUNCTIONS 
#############################################

#' Funtion used to scale a vector of numeric values from [x_min,x_max] range to [y_min,y_max] range.
#' 	f : R -> R
#'      [x_min, x_max] -> [y_min, y_max]
#'                         vect - x_min
#'  f(vect,y_min,y_max) =  ------------- x (y_max - y_min) + y_min
#'                         x_max - x_min
#'
#' @param vect vector to be transformated
#' @param nmin new minimum
#' @param nmax new maximum
#' @return transformated vector to new range
scale_range <- function(vect,nmin,nmax){
	((vect - min(vect))/(max(vect)-min(vect)))*(nmax-nmin)+nmin
}


#' Custom function to generate exponentially degradated foldchange range
#' @param means
#' @param fcmin minimum foldchange
#' @param fcmax maximum foldchange
#' @param meanLog param of rlnorm function
#' @param sdlog param of rlnorm function
#' @return a vector of foldchanges to be applied
fcfunc <- function(means,fcmin=1.4, fcmax=3, meanlog = 1, sdlog = 0.8){
	# Generate exponential distributio 
	xx <- rlnorm(length(means), meanlog = meanlog, sdlog = sdlog)
	xx <- scale_range(xx,fcmin, fcmax)
	ecdffun <- ecdf(xx)
	xx <- data.frame(FC = xx, Quant = unlist(lapply(xx,ecdffun)))
	xx <- xx[order(xx$Quant),]
	xx$Quant <- xx$Quant[seq(from = nrow(xx), to = 1)] # Apply inverse 
	# Prepare wuantiles of observed
	ecdffun <- ecdf(means)
	means <- data.frame(X = means, Quant = unlist(lapply(means,ecdffun)))
	# Merge values
	if(all(means$B %in% xx$B)){
		means <- merge(means,xx,by.y = "Quant",sort = FALSE) 
	}else{
		stop("Close search not implemented yet")
	}
	return(means)
}

#############################################
### LOAD & PREPARE 
#############################################

# Load count table if proceed
if(is.null(opt$inputfile)){
	bcount <- NULL
	group <- NULL
}else{
	bcount <- read.table(file = opt$inputfile, row.names = 1, header = TRUE, sep = "\t")
	if(is.null(opt$group)){
		group <- NULL
	}else{
		group <- as.numeric(unlist(strsplit(opt$group,",")))
	}
}

# Prepare DEGs proportion
posdeg <- opt$DEGs_proportion * opt$P_up
negdeg <- opt$DEGs_proportion - posdeg


#############################################
### SIMULATE 
#############################################

simul <- STC(Ngene      = opt$ngenes,
			 DEG.foldchange = function(means){fcfunc(means,fcmin = opt$FC_min,fcmax = opt$FC_max)}, 
			 replicates = opt$replicates, 
			 bcount     = bcount, 
			 group      = group, 
			 PDEG       = c(posdeg,negdeg))
prediction_vector <- simul$trueDEG
prediction_vector <- replace(prediction_vector, prediction_vector != 0, "TRUE")
prediction_vector <- replace(prediction_vector, prediction_vector == 0, "FALSE")
prediction_vector <- as.data.frame(prediction_vector)
prediction_vector <- cbind(rownames(prediction_vector),prediction_vector)
colnames(prediction_vector) <- c("Gene","Prediction")




#############################################
### EXPORT 
#############################################
write.table(simul$count, file=paste0(opt$outfile,"_scount"), quote = FALSE, col.names = TRUE, sep = "\t")
write.table(prediction_vector, file=paste0(opt$outfile,"_predv"), quote = FALSE, col.names = TRUE, sep = "\t", row.names = FALSE)
