
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
#' `counts_matrix` contains just the counts matrix. Element `DEGs` contains
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
    effect.size = 1.5, tree = NULL, prop.var.tree = 1,
    model.process = c("BM", "OU"), selection.strength = 0, id.condition = NULL,
    id.species = as.factor(rep(1, 2 * samples.per.cond)),
    check.id.species = TRUE, lengths.relmeans = NULL,
    lengths.dispersions = NULL, lengths.phylo = TRUE, output_dir = NULL,
    method = "vanilla", fixed_DEG_lists = "", fixed_upregulated_DEGs = "") {
    effect_sizes <- strsplit(effect.size, ",")[[1]]
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
        between.group.diffdisp = synth_diffidsp, 
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
    counts_matrix <- synth_data@count.matrix
    counts_matrix <- cbind(data.frame(gene = rownames(counts_matrix)),
                           counts_matrix)
    exp_design <- synth_data@sample.annotations["condition"]
    exp_design <- cbind(data.frame(sample = rownames(exp_design), exp_design))
    rownames(exp_design) <- NULL
    rownames(counts_matrix) <- NULL
    DEGs <- synth_data@variable.annotations["differential.expression"]
    DEGs <- cbind(data.frame(gene = rownames(DEGs), DEGs))
    rownames(DEGs) <- NULL
    if(method != "vanilla") {
        synth_data <- main_custom_synth(counts_table = counts_table, 
            experimental_design = experimental_design, nDEGs = nDEGs,
            fixed_DEG_lists = fixed_DEG_lists, effect_sizes = effect_sizes,
            fixed_upregulated_DEGs = fixed_upregulated_DEGs, 
            fraction_upregulated = fraction_upregulated, overlap_overlap_size, deg_method, samples_per_cond)
    }
    res <- list(all = synth_data, counts_matrix = counts_matrix, DEGs = DEGs,
                exp_design = exp_design)
    if(!is.null(output_dir)) {
        utils::write.table(res$counts_matrix, file.path(output_dir,
            "synth_counts.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
        utils::write.table(res$DEGs, file.path(output_dir,
            "synth_DEGs.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
        utils::write.table(res$exp_design, file.path(output_dir,
            "synth_exp_design.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
    } else {
        return(res)
    }
}


#' main_custom_synth
#'
#' `main_custom_synth` Is a wrapper for the `generateSyntheticData` function
#' from package `compcodeR`.
#'
#' @inheritParams compcodeR::generateSyntheticData
#' @importFrom compcodeR generateSyntheticData
#' @importFrom utils write.table
#' @param output_dir Directory where counts matrix and DEG table will be saved.
#' @returns A list. Element `all` contains the entire results object. Element
#' `counts_matrix` contains just the counts matrix. Element `DEGs` contains
#' a data frame detailing which genes are differentially expressed in dataset.
#' @examples
#' B_625_625 <- main_custom_synth(dataset = "B_625_625", n.vars = 12500, 
#'                                  samples.per.cond = 5, n.diffexp = 1250, 
#'                                  repl.id = 1, seqdepth = 1e7, 
#'                                  fraction.upregulated = 0.5, 
#'                                  between.group.diffdisp = FALSE, 
#'                                  filter.threshold.total = 1, 
#'                                  filter.threshold.mediancpm = 0, 
#'                                  fraction.non.overdispersed = 0)
#' B_625_625
#' @export

main_custom_synth <- function(counts_table, experimental_design, nDEGs,
    fixed_DEG_lists, fixed_upregulated_DEGs, effect_sizes,
    fraction_upregulated, overlap_size, deg_method, samples_per_cond) {
    if(opt$fixed_DEG_lists == "") {
        DEG_lists <- create_DEG_lists(universe = rownames(counts_table),
                                      nDEGs = nDEGs, overlap = overlap_size)    
    } else {
        DEG_lists <- read_DEG_lists(paths_string = fixed_DEG_lists)
    }
    new_tables <- rename_samples(exp_design = experimental_design,
                    counts_table = counts_table, conditions_column = c(2, 3))
    experimental_design <- new_tables$exp_design
    counts_table <- new_tables$counts_table
    message("Experimental design is")
    print(experimental_design)
    message("Count table head is ")
    head(counts_table)
    effect_sizes <- c(opt$effect_size_1, opt$effect_size_2)
    factor_columns <- c("condition", "condition_2")
    for(i in seq(factor_columns)) {
        if(opt$fixed_upregulated_DEGs == "") {
            n_ups <- floor(opt$fraction_upregulated * length(DEG_lists[[i]]))
            up_DEGs <- sample(DEG_lists[[i]], n_ups, replace = FALSE)
        } else {
            up_DEGs <- read_DEG_lists(opt$fixed_upregulated_DEGs)[[i]]
        }
        counts_table <- make_DEG_table(table = counts_table, up_DEGs = up_DEGs,
            experimental_design = experimental_design, DEGs = DEG_lists[[i]],
            effect_size = effect_sizes[i], factor_column = factor_columns[i],
            samples_per_group = opt$samples_per_group, method = opt$method)
    }
}
