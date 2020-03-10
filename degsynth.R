#! /usr/bin/env Rscript

#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>

#############################################
### CONFIGURE 
#############################################

#Loading libraries
suppressPackageStartupMessages(library(optparse)) 

option_list <- list(
  make_option(c("-r", "--replicates"), type="integer", default=3,
    help="Number of replicates for control/treatment group. Default : %default"),
  make_option(c("-n", "--ngenes"), type="integer", default=20000,
    help="Number of genes. Default : %default"),
  make_option(c("-d", "--DEGs_proportion"), type="double", default=0.2,
    help="Proportion of differentially expressed genes (DEGs). Default : %default"),
  make_option(c("-f", "--FC_min"), type="double", default=3.0,
    help="Minimum Fold Change value. Default : %default"),
  make_option(c("-F", "--FC_max"), type="double", default=1.3,
    help="Maximum Fold Change value. Default : %default"),
  make_option(c("-T", "--P_up"), type="double", default=1,
    help="Proportion of DEGs that will be up-regulated in treatment group. Rest will be down-regulated. Default : %default"),
  make_option(c("-i", "--inputfile"), type="character", default = NULL,
    help="[OPTIONAL] Input genes count matrix to be used for simulation"),
  make_option(c("-g", "--group"), type="character", default = NULL,
    help="Columns from input file to be used as same group of replicates. Specify indices separated by commas. Note: start at 1"),
  make_option(c("-o", "--outfile"), type="character",
    help="Output file basenames. A *_scount and *_predv files will be generated")
)

opt <- parse_args(OptionParser(option_list=option_list))

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


#' Function to Simulate Treatment-Control paired experiments for genes expression count tables. Pipeline inspired on TCC simulateReadCounts function
#' @param Ngene number of genes to be included in the experiment
#' @param PDEG Percentage of Differentially Expressed Genes. If one number is specified, only Up-regulated genes will be simulated. If a vector of two number is given, proportions will be applied for Up and Down regualted simulation respectively
#' @param DEG.foldchange foldchange to be applied. Can be a number or a function which must accept only a vector and return a vector of the same length with the fold-change values. Warning: NOT LOG FOLD CHANGE
#' @param replicates to be ge generated of each group (Treatment and Control)
#' @param bcount matrix of genes counts used to simulate. If NULL, arabidopsis dataset used (not from) in TCC package will be loaded
#' @param group columns to be used as group replicates. If NULL, all columns will be used.
#' @return a list with COUNT) a simulated gene counts matrix with simulated genes and sample names and TRUEDEG) a vector with real up/down regulated genes
#' @import TCC package(s)
#' @dependency scale_range function(s)   
STC <- function(Ngene = 1000, PDEG = 0.2, DEG.foldchange = 2, replicates = 3, bcount = NULL, group = NULL){
    # Load & prepare base data
	if(is.null(bcount)){
		require(TCC)
		data(arab)
		bcount <- arab
		group <- seq(3)
	}
	if(is.null(group)){
		group <- seq(ncol(bcount))
	}
	rpm.a <- sweep(bcount[, group], 2, median(colSums(bcount[, group]))/colSums(bcount[,group]), "*")
	rpm.a <- rpm.a[apply(rpm.a, 1, var) > 0, ]
	mean.a <- apply(rpm.a, 1, mean)
	var.a <- apply(rpm.a, 1, var)
	dispersion <- (var.a - mean.a)/(mean.a * mean.a)
	population <- data.frame(mean = mean.a, disp = dispersion)
	population <- population[population$disp > 0, ]
	population <- population[sample(1:nrow(population), Ngene, replace = TRUE), ]
	# Define which genes will be DEGs
	if(is.vector(PDEG)){
		ndeg_pos <- floor(Ngene*abs(PDEG[1]))
		ndeg_neg <- floor(Ngene*abs(PDEG[2]))
	}else{
		ndeg_pos <- floor(Ngene*abs(PDEG))
		ndeg_neg <- 0
	}
	ndeg <- ndeg_pos + ndeg_neg
	shuffle <- function(vector){vector[sample(seq_along(vector))]}
	DEGs <- c(rep(1.0,ndeg_pos),rep(-1.0,ndeg_neg),rep(0.0,Ngene - ndeg)) 
	DEGs <- shuffle(DEGs)
	# Check foldchange
	if(is.function(DEG.foldchange)){
		DEG.foldchange <- DEG.foldchange(population$mean)$FC
	}else{
		DEG.foldchange <- rep(DEG.foldchange,Ngene)
	}
	# Prepare foldchange
	experiments <- matrix(0,ncol = replicates * 2, nrow = Ngene)
	invisible(lapply(seq_along(DEGs),function(i){
		# Check
		if(DEGs[i] == 0){
			experiments[i,] <<- rep(1.0,ncol(experiments))
			return()
		} 
		# Generate foldchange
		experiments[i,] <<- c(rep(1.0,replicates),
						 DEGs[i] * runif(replicates,min = DEG.foldchange[i] - 0.2, max = DEG.foldchange[i] + 0.2))
	}))	
	# Generate
	experiments <- apply(experiments, 2, function(x, pp = population) {rnbinom(n = Ngene, mu = (abs(x)^(x/abs(x)))*pp$mean, size = 1/pp$disp)})
	# Change names
	colnames(experiments) <- c(paste0("Ctrl_rep",seq(replicates)),paste0("Treat_rep",seq(replicates)))
	rownames(experiments) <- paste0("gene_",seq(Ngene))
	names(DEGs) <- rownames(experiments)
	# END
	return(list(count = experiments, trueDEG = DEGs))
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
prediction_vector <- replace(prediction_vector, prediction_vector != 0, "DEG")
prediction_vector <- replace(prediction_vector, prediction_vector == 0, "NOTDEG")
prediction_vector <- as.data.frame(prediction_vector)
prediction_vector <- cbind(rownames(prediction_vector),prediction_vector)
colnames(prediction_vector) <- c("Gene","Prediction")




#############################################
### EXPORT 
#############################################
write.table(simul$count, file=paste0(opt$outfile,"_scount"), quote = FALSE, col.names = TRUE, sep = "\t")
write.table(prediction_vector, file=paste0(opt$outfile,"_predv"), quote = FALSE, col.names = TRUE, sep = "\t", row.names = FALSE)
