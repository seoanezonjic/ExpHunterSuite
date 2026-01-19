
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

#' main_compcodeR
#'
#' `main_compcodeR` Is a wrapper for the `generateSyntheticData` function
#' from package `compcodeR`.
#'
#' @inheritParams compcodeR::generateSyntheticData
#' @importFrom compcodeR generateSyntheticData
#' @importFrom utils write.table
#' @param output_dir Directory where counts matrix and DEG table will be saved.
#' @returns A list. Element `all` contains the entire results object. Element
#' `counts_table` contains just the counts matrix. Element `DEGs` contains
#' a data frame detailing which genes are differentially expressed in dataset.
#' @examples
#' B_625_625 <- main_compcodeR(dataset = "B_625_625", n.vars = 12500, 
#'                                  samples.per.cond = 5, n.diffexp = 1250, 
#'                                  repl.id = 1, seqdepth = 1e7, 
#'                                  fraction.upregulated = 0.5, 
#'                                  between.group.diffdisp = FALSE, 
#'                                  filter.threshold.total = 1, 
#'                                  filter.threshold.mediancpm = 0, 
#'                                  fraction.non.overdispersed = 0)
#' B_625_625
#' @export

main_compcodeR <- function(dataset, n.vars, samples.per.cond, n.diffexp,
    repl.id = 1, seqdepth = 1e+07, minfact = 0.7, maxfact = 1.4,
    relmeans = "auto", dispersions = "auto", fraction.upregulated = 1,
    between.group.diffdisp = FALSE, filter.threshold.total = 1,
    filter.threshold.mediancpm = 0, fraction.non.overdispersed = 0,
    random.outlier.high.prob = 0, random.outlier.low.prob = 0,
    single.outlier.high.prob = 0, single.outlier.low.prob = 0,
    effect_sizes = 1.5, tree = NULL, prop.var.tree = 1,
    nDEGs = as.character(nvars/2), model.process = c("BM", "OU"),
    selection.strength = 0, id.condition = NULL, overlap_size = 0,
    id.species = as.factor(rep(1, 2 * samples.per.cond)),
    check.id.species = TRUE, lengths.relmeans = NULL,
    lengths.dispersions = NULL, lengths.phylo = TRUE, output_dir = NULL,
    method = "vanilla", fixed_DEG_lists = "", fixed_upregulated_DEGs = "",
    condition_columns) {
    effect_sizes <- strsplit(effect_sizes, ",")[[1]]
    if(method != "vanilla") {
        synth_diffexp <- 0
        synth_diffdisp <- FALSE
    } else {
        if(length(effect_sizes) > 1) {
            stop("compcodeR DEG synthesis not compatible with multiple effect",
                 "sizes. Please provide only one.")
        }
        synth_diffexp <- n.diffexp
        synth_diffdisp <- between.group.diffdisp
    }
    synth_data <- compcodeR::generateSyntheticData(dataset = dataset,
        n.vars = n.vars, samples.per.cond = samples.per.cond,
        n.diffexp = synth_diffexp, repl.id = repl.id, seqdepth = seqdepth,
        minfact = minfact, maxfact = maxfact, relmeans = relmeans,
        dispersions = dispersions, effect.size = effect_sizes[1], tree = tree,
        fraction.upregulated = fraction.upregulated, 
        between.group.diffdisp = synth_diffdisp, 
        filter.threshold.total = filter.threshold.total,
        filter.threshold.mediancpm = filter.threshold.mediancpm,
        fraction.non.overdispersed = fraction.non.overdispersed,
        random.outlier.high.prob = random.outlier.high.prob,
        random.outlier.low.prob = random.outlier.low.prob,
        single.outlier.high.prob = single.outlier.high.prob,
        single.outlier.low.prob = single.outlier.low.prob,
        check.id.species = check.id.species,
        prop.var.tree = prop.var.tree, lengths.relmeans = lengths.relmeans,
        model.process = model.process, id.species = id.species,
        selection.strength = selection.strength,
        id.condition = id.condition, lengths.phylo = lengths.phylo,
        lengths.dispersions = lengths.dispersions)
    counts_table <- synth_data@count.matrix
    counts_table <- cbind(data.frame(gene = rownames(counts_table)),
                           counts_table)
    exp_design <- synth_data@sample.annotations["condition"]
    exp_design <- cbind(data.frame(sample = rownames(exp_design), exp_design))
    rownames(exp_design) <- NULL
    rownames(counts_table) <- NULL
    DEGs <- synth_data@variable.annotations["differential.expression"]
    DEGs <- cbind(data.frame(gene = rownames(DEGs), DEGs))
    rownames(DEGs) <- NULL
    if(method != "vanilla") {
        synth_data <- custom_synth(counts_table = counts_table, 
            exp_design = exp_design, nDEGs = nDEGs, effect_sizes = effect_sizes,
            columns = condition_columns, deg_method = method,
            fixed_DEG_lists = fixed_DEG_lists, overlap_size = overlap_size, 
            fixed_upregulated_DEGs = fixed_upregulated_DEGs, 
            fraction_upregulated = fraction.upregulated,
            samples_per_cond = samples.per.cond)
        counts_table <- synth_data$counts_table
        exp_design <- synth_data$exp_design
        DEGs <- synth_data$DEGs
    }
    res <- list(all = synth_data, counts_table = counts_table, DEGs = DEGs,
                exp_design = exp_design)
    if(!is.null(output_dir)) {
        utils::write.table(res$counts_table, file.path(output_dir,
            "synth_counts.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
        utils::write.table(res$DEGs, file.path(output_dir,
            "synth_DEGs.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
        utils::write.table(res$exp_design, file.path(output_dir,
            "synth_exp_design.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
    } else {
        return(res)
    }
}
