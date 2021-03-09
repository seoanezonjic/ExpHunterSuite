#' Funtion used to scale a vector of numeric values from 
#' (x_min,x_max) range to 
#' (y_min,y_max) range.
#'     f : R -> R
#'      (x_min, x_max) -> (y_min, y_max)
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
#' @param means means vector
#' @param fcmin minimum foldchange
#' @param fcmax maximum foldchange
#' @param meanlog param of rlnorm function
#' @param sdlog param of rlnorm function
#' @return a vector of foldchanges to be applied
#' @importFrom stats ecdf rlnorm
fcfunc <- function(means,fcmin=1.4, fcmax=3, meanlog = 1, sdlog = 0.8){
    # Generate exponential distributio 
    xx <- stats::rlnorm(length(means), meanlog = meanlog, sdlog = sdlog)
    xx <- scale_range(xx,fcmin, fcmax)
    ecdffun <- stats::ecdf(xx)
    xx <- data.frame(FC = xx, Quant = unlist(lapply(xx,ecdffun)))
    xx <- xx[order(xx$Quant),]
    xx$Quant <- xx$Quant[seq(from = nrow(xx), to = 1)] # Apply inverse 
    # Prepare quantiles of observed
    ecdffun <- stats::ecdf(means)
    means <- data.frame(X = means, Quant = unlist(lapply(means,ecdffun)))
    # Merge values
    if(all(means$B %in% xx$B)){
        means <- merge(means,xx,by.y = "Quant",sort = FALSE) 
    }else{
        stop("Close search not implemented yet")
    }
    return(means)
}

#' Main function to generate synthetic data using an specific exponential 
#' distribution for logFC
#' @param outfile output file
#' @param inputfile input file
#' @param replicates number of replicates
#' @param ngenes number of genes
#' @param DEGs_proportion numeric (0,1) proportion of DEG genes to be simulated
#' @param FC_min minimum Fold-Change
#' @param FC_max maximum Fold-Change
#' @param P_up proportion og up-regulated genes
#' @param group optional group identifiers seprated by commas
#' @keywords synthetic
#' @export
#' @return synthetic data generated
#' @importFrom utils read.table write.table
#' @examples
#' synthetic_dataset <- degsynth()
#' # Returns simulated count dataset and indication of which are DE/not DE
degsynth <- function(
    outfile = NULL,
    inputfile = NULL,
    replicates = 3,
    ngenes = 20000,
    DEGs_proportion = 0.2,
    FC_min = 1.3,
    FC_max = 3.0,
    P_up = 1,
    group = NULL
    ){

    # Load count table if proceed
    if(is.null(inputfile)){
        bcount <- NULL
        group <- NULL
    }else{
        bcount <- utils::read.table(file = inputfile, row.names = 1, 
            header = TRUE, sep = "\t")
        if(is.null(group)){
            group <- NULL
        }else{
            group <- as.numeric(unlist(strsplit(group,",")))
        }
    }

    # Prepare DEGs proportion
    posdeg <- DEGs_proportion * P_up
    negdeg <- DEGs_proportion - posdeg


    #############################################
    ### SIMULATE 
    #############################################
    simul <- STC(Ngene      = ngenes,
                 DEG.foldchange = function(means){fcfunc(means,
                                                          fcmin = FC_min,
                                                          fcmax = FC_max)}, 
                 replicates = replicates, 
                 bcount     = bcount, 
                 group      = group, 
                 PDEG       = c(posdeg,negdeg))
    prediction_vector <- simul$trueDEG
    prediction_vector <- replace(prediction_vector, 
        prediction_vector != 0, "TRUE")
    prediction_vector <- replace(prediction_vector, 
        prediction_vector == 0, "FALSE")
    prediction_vector <- as.data.frame(prediction_vector)
    prediction_vector <- cbind(rownames(prediction_vector),prediction_vector)
    colnames(prediction_vector) <- c("Gene","Prediction")




    #############################################
    ### EXPORT OR RETURN
    #############################################
    if(is.null(outfile)) {
        return(list(simul_count=simul$count, 
                    prediction_vector=prediction_vector))
    } else {
        utils::write.table(simul$count, file=paste0(outfile,"_scount"), 
            quote = FALSE, col.names = TRUE, sep = "\t")
        utils::write.table(prediction_vector, file=paste0(outfile,"_predv"), 
            quote = FALSE, col.names = TRUE, sep = "\t", row.names = FALSE)
    }
}

