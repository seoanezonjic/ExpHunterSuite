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
    return(((vect - min(vect))/(max(vect)-min(vect)))*(nmax-nmin)+nmin)
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

#' Add tag to duplicate values in vector

#' `mark_duplicates` finds duplicate values inside a vector and appends
#' a tag to each of them. This tag can either be an index or a letter
#'
#' @param vector Input vector
#' @param index_as_letter A boolean
#'   * `TRUE`:  Duplicates will be marked as value_A, value_B
#'   * `FALSE` (the default): Duplicates will be marked as value_1, value_2
#' @returns A vector with marked duplicates

mark_duplicates <- function(string, index_as_letter = FALSE) {
	if(index_as_letter) {
		mark_fun <- function(i) LETTERS[seq_along(i)]
	} else {
		mark_fun <- seq_along
	}
	res <- ave(string, string, FUN = function(i) paste0(i, "_", mark_fun(i)))
	return(res)
}

#' Use a dictionary to rename a vector

#' `rename_with_dict` renames a vector by consulting a list that links
#' the original names with the desired names
#'
#' @param vector Vector to rename
#' @param dict Dictionary containing name-new_name pairs
#' @returns Renamed vector

rename_with_dict <- function(vector, dict) {
    vector_IDs <- match(names(dict), vector)
    dict_IDs <- match(vector, names(dict))
    dict_IDs <- dict_IDs[!is.na(dict_IDs)]
    vector <- dict[dict_IDs]
    return(vector)
}

#' Give samples a more descriptive name
#'
#' `rename_samples` takes an experimental design table and renames samples
#' with their experimental conditions, making the design more descriptive at
#' a first glance and easier to interpret
#'
#' @param exp_design A data frame describing the experimental design
#' @param condition_columns An integer vector. Indices of columns to use as
#' conditions
#' @param counts_table Counts table where renaming will be propagated
#' @returns A list. Exp design is the renamed experimental design, counts_table
#' is the renamed counts table

rename_samples <- function(exp_design, condition_columns, counts_table) {
    new_counts_table <- counts_table[, !colnames(counts_table) == "gene"]
	subs_design <- exp_design[, c(1, condition_columns), drop = FALSE]
	new_names <- apply(exp_design[, -1, drop = FALSE], 1, paste, collapse = "_")
	new_names <- mark_duplicates(new_names)
	dict <- new_names
	names(dict) <- exp_design[, 1]
	colnames(new_counts_table) <- rename_with_dict(colnames(new_counts_table),
                                                   dict)
    subs_design$sample <- rename_with_dict(subs_design$sample, dict)
    subs_design <- subs_design[order(subs_design$sample), ]
    new_counts_table <- match_counts_to_design(exp_design = subs_design,
                                               counts_table = new_counts_table)
    new_counts_table <- cbind(gene = counts_table$gene, new_counts_table)
    new_counts_table <- htmlreportR::col_to_rownames(new_counts_table)
	return(list(exp_design = subs_design, counts_table = new_counts_table))
}

#' Generate synthetic DEGs following a random distribution.
#'
#' @description
#' `synth_normal_DEGs` generates two random vectors by which to multiply
#' expression values of samples in a counts table. Randomisation is taken from
#' the normal distribution, with a default standard deviation of two and an
#' effect size supplied by the user
#'
#' `synth_exponential_DEGs` follows an exponential distribution instead
#'
#' `synth_deterministic_DEGs` is a simple multiplication, with no noise
#'
#' @inheritParams synth_DEGs_by_condition
#' @returns Two vectors of values by which to multiply expression levels in
#' counts table. Upregulated vector includes effect size, static vector only
#' includes noise from the normal distribution.

synth_normal_DEGs <- function(up_vector, stat_vector, effect_size,
                              stdev = 0.2) {
    all_samples <- c(up_vector, stat_vector)
    n_samples <- length(all_samples)
    nscale <- length(up_vector)
    error_size <- stdev
    noise <- rnorm(n = n_samples, mean = 0, sd = error_size)
    scale_factor <- effect_size + 1 + abs(noise[1:nscale])
    distr_factor <- stat_vector * 2 * noise[(nscale + 1):n_samples] # Emulate
    # distribution with no effect size
    up_vector <- up_vector * scale_factor
    return(list(up_vector = round(up_vector), stat_vector = round(stat_vector)))
}

#' @rdname synth_normal_DEGs

synth_exponential_DEGs <- function(up_vector, stat_vector, effect_size,
                                   stdev = NULL) {
    all_samples <- c(up_vector, stat_vector)
    n_samples <- length(all_samples)
    nscale <- length(up_vector)
    ncontrol <- length(stat_vector)
    noise <- rexp(n = n_samples, rate = 1)
    samples <- up_vector * (effect_size + 1 + noise[1:nscale] ** 2)
    stat_vector <- 1 + noise[(nscale + 1):n_samples]
    return(list(up_vector = round(up_vector), stat_vector = round(stat_vector)))
}

#' @inheritParams synth_normal_DEGs

synth_deterministic_DEGs <- function(up_vector, stat_vector, effect_size,
                                     stdev = NULL){
    up_vector <- up_vector * (effect_size + 1)
    return(list(up_vector = round(up_vector), stat_vector = round(stat_vector)))
}

#' Generate DEGs using desired method in specified vectors.
#'
#' `synth_DEGs_by_condition` is a wrapper for the synthetic DEG generation
#' functions. Provide sample vectors and desired method, and this function will
#' dispatch the corresponding synthetic generator
#'
#' @param up_vector Samples that will be upregulated
#' @param stat_vector Samples that will not be upregulated, but still need
#' to be randomised in order to not introduce a bias
#' @param effect_size Effect size to emulate in expression values
#' @param method Method to use. Values: "norm" (`synth_normal_DEGs`), det
#' (`synth_deterministic_DEGs`) and "exp" (`synth_exponential_DEGs`)
#' @param stdev Standard deviation for distribution. Used in `synth_normal_DEGs`
#' @returns Two vectors of values by which to multiply expression levels in
#' counts table. Upregulated vector includes effect size, static vector only
#' includes noise from the normal distribution.

synth_DEGs_by_condition <- function(up_vector, stat_vector, effect_size,
                                    method = "norm", stdev = 0.2) {
    function_list <- c(norm = synth_normal_DEGs, det = synth_deterministic_DEGs,
    				   exp = synth_exponential_DEGs)
    synth_function <- function_list[[method]]
    res <- synth_function(up_vector, stat_vector, effect_size, stdev)
    return(list(up_vector = res$up_vector, stat_vector = res$stat_vector))
}

#' Randomly select which genes will be differentially expressed 
#' `create_DEG_lists` takes a list of genes and chooses which among them will
#' be differentially expressed.
#'
#' @param universe List of genes from which DEGs will be chosen
#' @param nDEGs Vector of number of DEGs to select. One number for every list
#' of DEGs to generate
#' @param overlap_size Proportion of DEGs that multiple experimental conditions
#' will share
#'
#' @returns A list containing a vector of DEGs for every experimental condition

create_DEG_lists <- function(universe, nDEGs, overlap_size) {
	names(nDEGs) <- seq(nDEGs)
	nDEGs <- sort(nDEGs)
    DEG_lists <- vector(mode = "list", length = length(nDEGs))
    names(DEG_lists) <- names(nDEGs)
    nsample <- nDEGs[[1]]
    overlap <- NULL
    for(i in seq(nDEGs)) {
        overlap_abs <- floor(overlap_size * min(nDEGs[[i]])) # Simpson overlap
    	if(i != 1) {
    		universe <- universe[!universe %in% DEG_lists[[i - 1]]]
    		nsample <- nDEGs[[i]] - overlap_abs
    	}
    	DEG_lists[[i]] <- sample(universe, nsample, replace = FALSE)
    	DEG_lists[[i]] <- c(DEG_lists[[i]], overlap)
    	if(i == 1) {
    		overlap <- sample(DEG_lists[[i]], overlap_abs, replace = FALSE)
    	}
    }
    DEG_lists <- DEG_lists[order(names(DEG_lists))]
    return(DEG_lists)
}

#' Calculate the coefficient of variation in a vector
#' `CV` performs a simple calculation to retrieve the coefficient of variation
#' of a vector
#'
#' @param vector Vector whose CV to calculate
#' @returns The coefficient of variation of the vector

CV <- function(vector) {
    return(var(vector) / mean(vector))
}

#' Extract samples by their value of requested experimental condition.
#' `define_conditions` takes an experimental design and the index of a factor
#' column and creates two vectors, one for each value of the supplied condition
#'
#' @param exp_design Experimental design
#' @param factor_column Index of column whose factor to consider

define_conditions <- function(exp_design, factor_column) {
    factors <- unique(exp_design[,factor_column])
    if(length(factors) > 2) {
        stop("Factor column has more than two unique levels. DEG synthesis impossible.")
    }
    cond_1 <- exp_design$sample[which(exp_design[factor_column] == factors[1])]
    cond_2 <- exp_design$sample[which(exp_design[factor_column] == factors[2])]
    return(list(cond_1 = cond_1, cond_2 = cond_2))
}

#' Sort counts table by samples in experimental design
#' `match_counts_to_design` makes sure that samples are sorted in the same order
#' in counts table and experimental design.
#'
#' @param counts_table Counts table
#' @param exp_design Experimental design
#' @returns Reordered counts table

match_counts_to_design <- function(counts_table, exp_design) {
    counts_table <- counts_table[, colnames(counts_table) != "gene"]
    new_count_order <- match(exp_design$sample, colnames(counts_table))
    res <- counts_table[, new_count_order]
    return(res)
}

#' Create synthetic DEG table
#' `make_DEG_table` adds a synthetic differential expression effect to a
#' counts table. It controls effect size, number of samples, standard deviation
#' and condition to consider.
#'
#' @inheritParams define_conditions
#' @inheritParams synth_DEGs_by_condition
#' @param table Counts table to modify
#' @param up_DEGs Genes to upregulate
#' @param samples_per_group Number of samples that each group should have
#' @returns Table with differential expression emulation

make_DEG_table <- function(table, exp_design, factor_column, up_DEGs,
                           effect_size, method, DEGs, samples_per_group,
                           stdev = 0.2) {
    conds <- define_conditions(exp_design, factor_column)
    res <- table
    for(gene in DEGs) {
        if(gene %in% up_DEGs) {
            up_cond <- conds$cond_2
            stat_cond <- conds$cond_1
        } else {
            up_cond <- conds$cond_1
            stat_cond <- conds$cond_2
        }
        up_samples <- exp_design$sample %in% up_cond
        stat_samples <- exp_design$sample %in% stat_cond
        up_vector <- table[gene, up_samples]
        stat_vector <- table[gene, stat_samples]
        synths <- synth_DEGs_by_condition(up_vector = up_vector,
            stat_vector = stat_vector, effect_size = effect_size,
            method = method, stdev = stdev)
        res[gene, up_samples] <- synths$up_vector
        res[gene, stat_samples] <- synths$stat_vector
    }
    return(res)
}

#' Read list of DEGs from disk
#' `read_DEG_lists` was developed for testing purposes. It loads a list of DEGs
#' from disk and returns it in a format that makes it usable in our DEG
#' synthesis library.

read_DEG_lists <- function(paths_string) {
    paths <- strsplit(paths_string, ",")[[1]]
    res <- lapply(paths, function(x) readLines(x))
    return(res)
}

#' custom_synth
#'
#' `custom_synth` Is a wrapper for the `generateSyntheticData` function
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

custom_synth <- function(counts_table, exp_design, nDEGs, columns,
    fixed_DEG_lists, fixed_upregulated_DEGs, effect_sizes,
    fraction_upregulated, overlap_size, deg_method, samples_per_cond) {
    nDEGs <- as.integer(strsplit(nDEGs, ",")[[1]])
    if(fixed_DEG_lists == "") {
        DEG_lists <- create_DEG_lists(universe = counts_table$gene,
                                nDEGs = nDEGs, overlap = overlap_size)
    } else {
        DEG_lists <- read_DEG_lists(paths_string = fixed_DEG_lists)
    }
    columns <- as.integer(strsplit(columns, ",")[[1]])
    if(any(duplicated(columns))) {
        stop("Supplied duplicate column indexes. Please fix value of columns",
        " argument")
    }
    if(length(columns) != length(DEG_lists)) {
        stop("Factor columns and DEG conditions do not match. Please specify
              as many columns as there are conditions")
    }
    new_tables <- rename_samples(exp_design = exp_design,
                    counts_table = counts_table, condition_columns = columns)
    exp_design <- new_tables$exp_design
    counts_table <- new_tables$counts_table
    factor_columns <- colnames(new_tables$exp_design)[columns]
    for(i in seq(factor_columns)) {
        if(fixed_upregulated_DEGs == "") {
            n_ups <- floor(fraction_upregulated * length(DEG_lists[[i]]))
            up_DEGs <- sample(DEG_lists[[i]], n_ups, replace = FALSE)
        } else {
            up_DEGs <- read_DEG_lists(fixed_upregulated_DEGs)[[i]]
        }
        counts_table <- make_DEG_table(table = counts_table, up_DEGs = up_DEGs,
            exp_design = exp_design, DEGs = DEG_lists[[i]],
            effect_size = effect_sizes[i], factor_column = factor_columns[i],
            samples_per_group = samples_per_cond, method = deg_method)
    }
    DEGs <- get_diffexp_info(counts_table, unique(unlist(DEG_lists)))
    return(list(exp_design = exp_design, counts_table = counts_table,
                DEGs = DEGs))
}

#' Add a column to a counts table specifying DEG status
#'
#' `get_diffexp_info` adds a column to a counts table that specifies whether
#' or not the gene is differentially expressed, imitating the output of counts
#' tables generated by the compcodeR library.
#'
#' @param counts_table Counts table
#' @param DEGs List of genes to mark as DEGs
#' @returns A counts table with a new column with value `1` for DEGs and value
#' `0` for control genes.

get_diffexp_info <- function(counts_table, DEGs) {
    res <- data.frame(gene = rownames(counts_table),differential.expression = 0)
    res[res$gene %in% DEGs, ]$differential.expression <- 1
    return(res)
}
