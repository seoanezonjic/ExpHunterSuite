#' Function to Simulate Treatment-Control paired experiments for genes
#' expression count tables. Pipeline inspired on TCC simulateReadCounts function
#' @param Ngene number of genes to be included in the experiment
#' @param PDEG Percentage of Differentially Expressed Genes. If one number
#' is specified, only Up-regulated genes will be simulated. If a vector of two
#' number is given, proportions will be applied for Up and Down regualted
#' simulation respectively
#' @param DEG.foldchange foldchange to be applied. Can be a number or a
#' function which must accept only a vector and return a vector of the same
#' length with the fold-change values. Warning: NOT LOG FOLD CHANGE
#' @param replicates to be ge generated of each group (Treatment and Control)
#' @param bcount matrix of genes counts used to simulate. If NULL,
#' arabidopsis dataset used (not from) in TCC package will be loaded
#' @param group columns to be used as group replicates. If NULL,
#' all columns will be used.
#' @return a list with COUNT) a simulated gene counts matrix with simulated
#' genes and sample names and TRUEDEG) a vector with real up/down regulated
#' genes
#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>
#' @importFrom utils data
#' @importFrom stats median var runif rnbinom
STC <- function(Ngene = 1000, 
    PDEG = 0.2, 
    DEG.foldchange = 2, 
    replicates = 3, 
    bcount = NULL, 
    group = NULL){
    # Load & prepare base data
    if(is.null(bcount)){
        arab <- NULL
        utils::data("arab", envir = environment())
        bcount <- arab
        group <- seq(3)
    }
    if(is.null(group)){
        group <- seq(ncol(bcount))
    }
    rpm.a <- sweep(bcount[, group], 2, stats::median(
        colSums(bcount[, group]))/colSums(bcount[,group]), "*")
    rpm.a <- rpm.a[apply(rpm.a, 1, stats::var) > 0, ]
    mean.a <- apply(rpm.a, 1, mean)
    var.a <- apply(rpm.a, 1, stats::var)
    dispersion <- (var.a - mean.a)/(mean.a * mean.a)
    population <- data.frame(mean = mean.a, disp = dispersion)
    population <- population[population$disp > 0, ]
    population <- population[sample(seq(nrow(population)), 
        Ngene, replace = TRUE), ]
    # Define which genes will be DEGs
    if(is.vector(PDEG) && length(PDEG) > 1){
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
                         DEGs[i] * stats::runif(replicates,
                             min = DEG.foldchange[i] - 0.2, 
                             max = DEG.foldchange[i] + 0.2))
    }))    
    # Generate
    experiments <- apply(experiments, 2, function(x, pp = population) {
        stats::rnbinom(n = Ngene, 
            mu = (abs(x)^(x/abs(x)))*pp$mean, size = 1/pp$disp)
    })
    # Change names
    colnames(experiments) <- c(paste0("Ctrl_rep",seq(replicates)),
        paste0("Treat_rep",seq(replicates)))
    rownames(experiments) <- paste0("gene_",seq(Ngene))
    names(DEGs) <- rownames(experiments)
    # END
    return(list(count = experiments, trueDEG = DEGs))
}
