#! /usr/bin/env Rscript


## Function has been deleted in a FRASER update. Will redefine it here
## until we have time to update workflow to new FRASER version
## https://github.com/gagneurlab/drop/tree/master/drop/modules

resultsByGenes <- function(res, geneColumn="hgncSymbol", method="BY"){
    # sort by pvalue
    res <- res[order(res$pValue)]

    # extract subset
    if(is(res, "GRanges")){
        ans <- data.table::as.data.table(S4Vectors::mcols(res)[,c(geneColumn, "pValue", "sampleID")])
        colnames(ans) <- c("features", "pval", "sampleID")
    } else {
        ans <- featureNames <- res[,.(
                features=get(geneColumn), pval=pvalue, sampleID=sampleID)]
    }

    # remove NAs
    naIdx <- ans[,is.na(features)]
    ansNoNA <- ans[!is.na(features)]

    # compute pvalues by gene
    ansNoNA[,pByFeature:=min(p.adjust(pval, method="holm")),
            by="sampleID,features"]

    # subset to lowest pvalue by gene
    dupIdx <- duplicated(ansNoNA[,.(features,sampleID)])
    ansGenes <- ansNoNA[!dupIdx]

    # compute FDR
    ansGenes[,fdrByFeature:=p.adjust(pByFeature, method=method),
            by="sampleID"]

    # get final result table
    finalAns <- res[!naIdx][!dupIdx]
    finalAns$pValueGene  <- ansGenes$pByFeature
    finalAns$padjustGene <- ansGenes$fdrByFeature
    finalAns
}

# Results formatting functions

get_aberrants <- function(df, pval_cutoff = 0.01) {
  if(is.null(df)) {
    return(NULL)
  }
  colnames(df)[colnames(df)=="hgncSymbol"] <- "geneID"
  df$aberrant = FALSE
  df$aberrant[df$padjust<=pval_cutoff] <- TRUE
  aberrants <- df$aberrant==TRUE
  genes <- df[aberrants, ]$geneID
  samples <- df[aberrants, ]$sampleID
  res <- df[df$geneID %in% genes & df$sampleID %in% samples, ]
  if((any(is.na(res)))) {
    warning("NAs in aberrant genes, removing")
    res <- res[!is.na(res$geneID), ]
  }
  if(nrow(res)==0) {
    warning("No aberrants found in input")
    return(NULL)
  }
  return(res)
}

format_aberrants <- function(df) {
  if(is.null(df)) {
    warning("NULL input. This can be due to missing module results or no
      aberrants found.")
    return(NULL)
  }
  names <- c("sampleID", "geneID", "padjust", "type", "zScore", "altRatio")
  matches <- colnames(df) %in% names
  res <- df[, matches]
  res$padjust <- -log(res$padjust)
  colnames(res)[colnames(res)=="padjust"] <- "p_padjust"
  return(res)
}

subset_types <- function(splicing_df) {
  if(is.null(splicing_df)) {
    return(NULL)
  }
  working_df <- splicing_df[, colnames(splicing_df)!="zScore"]
  pval_column <- grep("p_padjust", colnames(working_df))
  p_padjust <- working_df$p_padjust[!duplicated(working_df[, -pval_column])]
  working_df <- working_df[, -pval_column]
  working_df <- plyr::ddply(working_df, colnames(working_df), nrow)
  working_df <- cbind(working_df, p_padjust)
  colnames(working_df)[colnames(working_df)=="V1"] <- "spliceSites"
  jaccard <- working_df[working_df$type == "jaccard", ]
  psi3 <- working_df[working_df$type == "psi3", ]
  psi5 <- working_df[working_df$type == "psi5", ]
  theta <- working_df[working_df$type == "theta", ]
  res <- list(jaccard = jaccard, psi3 = psi3, psi5 = psi5, theta = theta)
  return(res)
}

option_list <- list(
  optparse::make_option(c("-a", "--sample_annotation"), type="character", default=NULL,
    help="Sample annotation table in tsv format."),
  optparse::make_option(c("-d", "--dataset"), type="character", default=NULL,
    help="Name of the dataset."),
  optparse::make_option(c("-g", "--genome_UCSC"), type="character", default=NULL,
    help="Genome version in UCSC format."),
   optparse::make_option(c("-n", "--cpu"), type="integer", default=1,
    help="Number of CPUs provided to job."),
  optparse::make_option(c("-r", "--recount"), type="logical", default=NULL, action = "store_true",
    help="Recount existing split reads."),
  optparse::make_option(c("-k", "--keep_non_standard_chrs"), type="logical", default=NULL, action = "store_true",
    help="Whether to keep non standard chromosomes."),
  optparse::make_option(c("-m", "--min_expression_in_one_sample"), type="integer", default=NULL,
    help="Minimum expression filter."),
  optparse::make_option(c("-l", "--long_read"), type="logical", default=NULL, action = "store_true",
    help="Boolean."),
  optparse::make_option(c("-q", "--quantile"), type = "integer", default = NULL,
    help="Quantile cutoff."),
  optparse::make_option(c("--quantile_min_expression"), type = "integer", default = NULL,
    help="Min expression cutoff in selected quantile."),
  optparse::make_option(c("--min_delta_psi"), type = "integer", default = NULL,
    help="Only introns for which the maximal difference in the psi value of a sample to the mean psi of the intron is larger than this value pass the filter."),
  optparse::make_option(c("-f" ,"--filter"), type = "logical", default = FALSE, action = "store_true",
    help="Discard junctions that do not pass the filter."),
  optparse::make_option(c("--anno_database"), type="character", default=NULL,
    help="TxDb annotation database in db format."),
   optparse::make_option(c("-s", "--seed"), type="character", default="",
    help="Use a random or set seed."),
  optparse::make_option(c("-i", "--implementation"), type="character", default=NULL,
    help="Implementation to use in fitting."),
  optparse::make_option(c("--max_tested_dimension_proportion"), type="integer", default=NULL,
    help="Number of CPUs provided to job."),
  optparse::make_option(c("-o", "--add_hpo_cols"), type="character", default=NULL,
    help="Path to add_HPO_cols.R script. We use a modified version which does not require internet connection,
      	  looking for the hpo file locally."),
  optparse::make_option(c("--gene_mapping_file"), type="character", default=NULL,
    help="Gene mapping file in tsv format."),
  optparse::make_option(c("--p_adj_cutoff"), type="integer", default=0.05,
    help="Adjusted p-value cutoff to consider a splicing metric as aberrant."),
  optparse::make_option(c("--delta_psi_cutoff"), type="integer", default=0.05,
    help="Minimum delta psi value to consider a splicing metric as aberrant."),
  optparse::make_option(c("--output_dir"), type="character", default=NULL,
    help="Directory were output will be written."),
  optparse::make_option(c("--hpo_file"), type="character", default=NULL,
    help="Genes and associated HPO terms in compressed tsv format."),
  optparse::make_option(c("-z", "--z_score_cutoff"), type="integer", default=0.05,
    help="Minimum delta psi value to consider a splicing metric as aberrant."),
  optparse::make_option(c("-t", "--top_N"), type="integer", default=10,
    help="Top N genes by adjusted p-value to be selected.")
  )

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
options("FRASER.maxSamplesNoHDF5"=0)
options("FRASER.maxJunctionsNoHDF5"=-1)

BiocParallel::register(BiocParallel::MulticoreParam(opt$cpu))
# Limit number of threads for DelayedArray operations
DelayedArray::setAutoBPPARAM(BiocParallel::MulticoreParam(opt$cpu))

sample_annotation <- data.table::fread(opt$sample_annotation)

annoGroups <- stringr::str_split(sample_annotation$EXP_GROUP,",")
subset_index <- unlist(lapply(annoGroups, function(x) any(x==opt$dataset)))
subset_ids <- sample_annotation$RNA_ID[subset_index]
annoSub <- sample_annotation[RNA_ID %in% subset_ids]
data.table::setnames(annoSub, "RNA_ID", "sampleID")
data.table::setnames(annoSub, "RNA_BAM_FILE", "bamFile")
data.table::setcolorder(annoSub, unique(c("sampleID", "STRAND", "PAIRED_END", "bamFile"), colnames(annoSub)))

RNA_BAM_FILE <- annoSub$bamFile
col_data_file<- tibble::add_column(annoSub, RNA_BAM_FILE, .after = "bamFile")

# Create initial FRASER object

fds <- FRASER::FraserDataSet(colData = col_data_file, workingDir = getwd(),
                             name = paste0("raw-local-", opt$dataset))

# Add paired end and strand specificity to the fds
FRASER::pairedEnd(fds) <- SummarizedExperiment::colData(fds)$PAIRED_END
FRASER::strandSpecific(fds) <- 'no'
if(data.table::uniqueN(SummarizedExperiment::colData(fds)$STRAND) == 1){
  FRASER::strandSpecific(fds) <- unique(SummarizedExperiment::colData(fds)$STRAND)
}
message(date(), ": FRASER object initialized for ", opt$dataset)
sample_ids <- colnames(fds) # Dumb adjustment I had to make due to their insistence in mixing external and internal data @alvaro
if(all(FRASER::strandSpecific(fds) != 0)){
  genome <- NULL
} else {
  genome <- BSgenome::getBSgenome(opt$genome_UCSC)
}
sample_results <- vector(mode='list', length=length(sample_ids))
for (i in seq(length(sample_ids))) { # This needs parallelisation. Furthermore, right now every sample is counted once for every group to which it belongs, whereas they should be counted only once AND THEN be divided by groups
                                     # FRASER has a function called getSplitReadCountsForAllSamples. Might be useful
  sample_results[[i]] <- FRASER::countSplitReads(sampleID = sample_ids[i], fds = fds, NcpuPerSample = opt$cpu,
                                 recount = opt$recount, keepNonStandardChromosomes = opt$keep_non_standard_chrs, genome = genome)
  message(date(), ": ", opt$dataset, ", ", sample_ids[i],
        " no. splice junctions (split counts) = ", length(sample_results[[i]]))
}
opt$min_expression_in_one_sample <- opt$min_expression_in_one_sample

# If samples are recounted, remove the merged ones
splitCountsDir <- file.path(getwd(), "savedObjects", 
                          paste0("raw-local-", opt$dataset), 'splitCounts')

if(opt$recount == TRUE & dir.exists(splitCountsDir)){
  unlink(splitCountsDir, recursive = TRUE)
}
# Why is this not use to count on the first place??? Can be run in parallel, would remove need for the other loop
splitCounts <- FRASER::getSplitReadCountsForAllSamples(fds=fds, recount=FALSE)
# Extract, annotate and save granges
splitCountRanges <- SummarizedExperiment::rowRanges(splitCounts)
# Annotate granges from the split counts
splitCountRanges <- FRASER:::annotateSpliceSite(splitCountRanges)
# Create ranges for non split counts
# Subset by minExpression
maxCount <- DelayedMatrixStats::rowMaxs(SummarizedExperiment::assay(splitCounts, "rawCountsJ"))
passed <- maxCount >= opt$min_expression_in_one_sample
# extract granges after filtering
nonSplitCountRanges <- splitCountRanges[passed, ]
# Extract splitSiteCoodinates: extract donor and acceptor sites
# take either filtered or full fds
spliceSiteCoords <- FRASER:::extractSpliceSiteCoordinates(splitCountRanges)
message(date(), ": ", opt$dataset, " total no. splice junctions = ", 
        length(splitCounts))
# Get sample id from wildcard
# Read splice site coordinates from RDS
# Count nonSplitReads for all samples
sample_results <- vector(mode = 'list', length = length(sample_ids))
for (i in seq(length(sample_ids))) {
        sample_results[[i]] <- FRASER::countNonSplicedReads(sampleID = sample_ids[i], ## Rsubread is incapable of following the symbolik link created by AutoFlow, used in tmpdir argument. I'll have to look for a way to pass a custom tmpdir. Idea: setting workingDir(fds) to getwd()
                                              splitCountRanges = NULL, fds = fds,
                                              NcpuPerSample = opt$cpu, minAnchor=5,
                                              recount=opt$recount,
                                              spliceSiteCoords=spliceSiteCoords,
                                              longRead=opt$long_read)
        message(date(), ": ", opt$dataset, ", ", sample_ids[i],
                " no. splice junctions (non split counts) = ", length(sample_results[[i]]))
}
# If samples are recounted, remove the merged ones
nonSplitCountsDir <- file.path(getwd(), "savedObjects", 
                            paste0("raw-local-", opt$dataset), 'nonSplitCounts')
if(dir.exists(nonSplitCountsDir)){
  unlink(nonSplitCountsDir, recursive = TRUE)
}
# Get and merge nonSplitReads for all sample ids
nonSplitCounts <- FRASER::getNonSplitReadCountsForAllSamples(fds = fds,
						  minAnchor = 5, splitCountRanges = nonSplitCountRanges, 
                          recount = FALSE, longRead = opt$long_read)                                                                                                                                                                                 
message(date(), ":", opt$dataset, " nonSplit counts done")

# Get splitReads and nonSplitRead counts in order to store them in FRASER object
splitCounts_h5 <- HDF5Array::HDF5Array(file.path(getwd(), paste0("savedObjects/raw-local-", opt$dataset, "/rawCountsJ.h5")), "rawCountsJ")
splitCounts_se <- SummarizedExperiment::SummarizedExperiment(colData = SummarizedExperiment::colData(fds),
      rowRanges = splitCountRanges, assays = list(rawCountsJ=splitCounts_h5))
nonSplitCounts_h5 <- HDF5Array::HDF5Array(file.path(getwd(), paste0("savedObjects/raw-local-", opt$dataset, "/rawCountsSS.h5")), "rawCountsSS")
nonSplitCounts_se <- SummarizedExperiment::SummarizedExperiment(colData = SummarizedExperiment::colData(fds),
  rowRanges = spliceSiteCoords, assays = list(rawCountsSS=nonSplitCounts_h5)
)
fds <- FRASER::addCountsToFraserDataSet(fds=fds, splitCounts=splitCounts_se,
                                nonSplitCounts=nonSplitCounts_se)
fds <- FRASER::calculatePSIValues(fds, overwriteCts = TRUE)

# External recognition block @alvaro. Not yet in working order
if (any(fds$EXTERNAL == "external")) {
    annoGroups <- stringr::str_split(sample_annotation$EXP_GROUP,",")
    exCountIndex <- which(annoGroups==opt$dataset) # annoGroups will be a list of vectors, this line will recognize those which ONLY contain fraser_external @alvaro
    exCountFiles <- sample_annotation$SPLICE_COUNTS_DIR[exCountIndex]
    exCountIDs <- sample_annotation$RNA_ID[exCountIndex]
    # Add external data if provided by dataset
	# Block explanation (@alvaro): if exCountIDs is greater than zero, it means the dataset contains external data.
	# If that is the case, the fds object is saved in a new directory, named "raw". Then, external information is
	# processed, and the updated object overwrites the one that has just been saved. Seems redundant. In any case,
	# now an fds object exists inside a "raw" folder, and it contains information on external counts.
    message("create new merged fraser object")
    for(resource in unique(exCountFiles)){
        exSampleIDs <- exCountIDs[exCountFiles == resource]
        exAnno <- fread(opt$sample_annotation, key="RNA_ID")[J(exSampleIDs)] # Reads the sample annotation table. Does not read it for group recognition, but does here.
        data.table::setnames(exAnno, "RNA_ID", "sampleID")
        ctsNames <- c("k_j", "k_theta", "n_psi3", "n_psi5", "n_theta")
        ctsFiles <- paste0(resource, "/", ctsNames, "_counts.tsv.gz") # Edited this line. Was previously dirname(resource), which deleted the name of the directory.
        # Merging external counts restricts the junctions to those that 
        # are only present in both the counted (fromBam) junctions AND the 
        # junctions from the external counts.
        fds <- FRASER::mergeExternalData(fds = fds, countFiles = ctsFiles,
                sampleIDs = exSampleIDs, annotation = exAnno)
        fds@colData$isExternal <- as.factor(!is.na(fds@colData$SPLICE_COUNTS_DIR))
    }
} else {
    fds@colData$isExternal <- as.factor(FALSE)
    FRASER::name(fds) <- paste0("raw-", opt$dataset)
}

fds <- FRASER::filterExpressionAndVariability(fds, 
        minExpressionInOneSample = opt$min_expression_in_one_sample,
        quantile = opt$quantile,
        quantileMinExpression = opt$quantile_min_expression,
        minDeltaPsi = opt$min_delta_psi,
        filterOnJaccard = FALSE,
        filter=FALSE)

# Keep junctions that pass filter
FRASER::name(fds) <- opt$dataset
if (opt$filter == TRUE) {
    filtered <- S4Vectors::mcols(fds, type="j")[,"passed"]
    fds <- fds[filtered,]
    message(paste("filtered to", nrow(fds), "junctions"))
}

out_k_files <- c("k_j_counts.tsv.gz", "k_theta_counts.tsv.gz")
out_n_files <- c("n_psi5_counts.tsv.gz", "n_psi3_counts.tsv.gz", "n_theta_counts.tsv.gz")

txdb <- AnnotationDbi::loadDb(opt$anno_database)
introns <- unique(unlist(GenomicFeatures::intronsByTranscript(txdb)))
introns <- GenomeInfoDb::keepStandardChromosomes(introns, pruning.mode = 'coarse')

GenomeInfoDb::seqlevels(fds) <- GenomeInfoDb::seqlevelsInUse(fds)
SummarizedExperiment::colData(fds)$sampleID <- as.character(SummarizedExperiment::colData(fds)$sampleID)
GenomeInfoDb::seqlevelsStyle(fds) <- GenomeInfoDb::seqlevelsStyle(introns)[1]
fds_known <- fds[unique(to(IRanges::findOverlaps(introns, MatrixGenerics::rowRanges(fds, type="j"), type="equal"))), ]

## Unavoidable due to the call of BiocGenerics::get, I have not been able to specify namespace of function to get so we need the library
library(FRASER)

# save k/n counts
sapply(c(out_k_files, out_n_files), function(i){
  ctsType <- toupper(strsplit(basename(i), "_")[[1]][1])
  psiType <- strsplit(basename(i), "_")[[1]][2]
  cts <- data.table::as.data.table(BiocGenerics::get(ctsType)(fds_known, type = psiType))
  grAnno <- MatrixGenerics::rowRanges(fds_known, type=psiType)
  anno <- data.table::as.data.table(grAnno)
  anno <- anno[,.(seqnames, start, end, strand)]
  
  data.table::fwrite(cbind(anno, cts), file=paste(opt$dataset, i, sep="_"), quote=FALSE, row.names=FALSE, sep="\t", compress="gzip")
}) |> invisible()

has_external <- (all(is.null(fds@colData$SPLICE_COUNTS_DIR)) || is.null(fds@colData$SPLICE_COUNTS_DIR))
if(has_external){
    fds@colData$isExternal <- as.factor(!is.na(fds@colData$SPLICE_COUNTS_DIR))
}else{
    fds@colData$isExternal <- as.factor(FALSE)
}

if(has_external){
    externalCountIDs <- SummarizedExperiment::colData(fds)[as.logical(SummarizedExperiment::colData(fds)[, "isExternal"]), "sampleID"]
    localCountIDs <- SummarizedExperiment::colData(fds)[!as.logical(SummarizedExperiment::colData(fds)[, "isExternal"]), "sampleID"]

    cts <- FRASER::K(fds,"psi5")
    ctsLocal<- cts[, localCountIDs]
    ctsExt<- cts[, externalCountIDs]

    rowMeanLocal <- rowMeans(ctsLocal)
    rowMeanExt <- rowMeans(ctsExt)

    dt <- data.table("Mean counts of local samples" = rowMeanLocal,
                     "Mean counts of external samples" = rowMeanExt)
                 
    ggplot2::ggplot(dt, ggplot2::aes(x = `Mean counts of local samples`, y= `Mean counts of external samples`)) +
       ggplot2::geom_hex() + cowplot::theme_cowplot(font_size = 16) +
	     ggplot2::theme_bw() + ggplot2::scale_x_log10() + ggplot2::scale_y_log10() + 
       ggplot2::geom_abline(slope = 1, intercept =0) +
       ggplot2::scale_color_brewer(palette="Dark2") + 
       ggplot2::ggtitle(opt$dataset) # Added ggtitles to help identify output @alvaro
}else{
	print("No external counts, comparison is ommitted")
}

if(!is.null(opt$seed)) {
  if(opt$seed != "") {
    set.seed(as.numeric(opt$seed))
  }
}

mp <- opt$max_tested_dimension_proportion

# Get range for latent space dimension
a <- 2 
b <- min(ncol(fds), nrow(fds)) / mp   # N/mp

maxSteps <- 12
if(mp < 6){
  maxSteps <- 15
}

Nsteps <- min(maxSteps, b)
pars_q <- round(exp(seq(log(a), log(b), length.out = Nsteps))) |> unique()


if(opt$genome_UCSC == "hg38") {
  GRCh <- 38
} else if(opt$genome_UCSC == "hg19") {
  GRCh <- 37
} else {
  stop("Invalid genome version. Please select hg38 or hg19. Selected ", opt$genome_UCSC)
}

for(type in c("jaccard")){
    message(date(), ": ", type)
    fds <- FRASER::optimHyperParams(fds, type = type, implementation = opt$implementation,
                                    q_param = pars_q, plot = FALSE)
    FRASER::currentType(fds) <- type
    q <- FRASER::bestQ(fds, type)
    FRASER::verbose(fds) <- 3   # Add verbosity to the FRASER object
    fds <- FRASER::fit(fds, q = q, type = type, iterations = 15, implementation = opt$implementation)
    fds <- FRASER::calculateZscore(fds, type=type)
    # Pvalues
    fds <- FRASER::calculatePvalues(fds, type=type)
    # Adjust Pvalues
    fds <- FRASER::calculatePadjValues(fds, type=type)
}

source(opt$add_hpo_cols)

orgdb <- data.table::fread(opt$gene_mapping_file)
seqlevels_fds <- GenomeInfoDb::seqlevelsStyle(fds)[1]
GenomeInfoDb::seqlevelsStyle(orgdb$seqnames) <- seqlevels_fds
GenomeInfoDb::seqlevelsStyle(txdb) <- seqlevels_fds

fds <- FRASER::annotateRangesWithTxDb(fds, txdb = txdb, orgDb = orgdb, 
			   feature = 'gene_name', featureName = 'hgnc_symbol', keytype = 'gene_id')

res_junc <- FRASER::results(fds, psiType = c("j"), ## This is wrong, but quickest way to get it to work before we can update to new workflow compatible with new FRASER version
                            padjCutoff = opt$p_adj_cutoff, deltaPsiCutoff = opt$delta_psi_cutoff)
res_junc_dt   <- data.table::as.data.table(res_junc)
print('Results per junction extracted')

if(nrow(res_junc_dt) > 0){

  # number of samples per gene and variant
  res_junc_dt[, numSamplesPerGene := data.table::uniqueN(sampleID), by = hgncSymbol]
  res_junc_dt[, numEventsPerGene := .N, by = "hgncSymbol,sampleID"]
  res_junc_dt[, numSamplesPerJunc := data.table::uniqueN(sampleID), by = "seqnames,start,end,strand"]

  # add colData to the results
  res_junc_dt <- merge(res_junc_dt, data.table::as.data.table(SummarizedExperiment::colData(fds)), by = "sampleID")
  res_junc_dt[, c("bamFile", "pairedEnd", "STRAND") := NULL]
} else{
  warning("The aberrant splicing pipeline gave 0 results for the ", opt$dataset, " dataset.")
}

# Aggregate results by gene
if(length(res_junc) > 0){
  res_genes_dt <- resultsByGenes(res_junc) |> data.table::as.data.table()
  res_genes_dt <- merge(res_genes_dt, data.table::as.data.table(SummarizedExperiment::colData(fds)), by = "sampleID")
  res_genes_dt[, c("bamFile", "pairedEnd", "STRAND") := NULL]
  # add HPO overlap information
  sa <- data.table::fread(opt$sample_annotation, 
              colClasses = c(RNA_ID = 'character', DNA_ID = 'character'))
  if(!is.null(sa$HPO_TERMS)){
    if(!all(is.na(sa$HPO_TERMS)) & ! all(sa$HPO_TERMS == '')){
      res_genes_dt <- add_HPO_cols(res_genes_dt, hpo_file = opt$hpoFile) # Changed this so the file path is looked up on config file @alvaro
    }
  }
} else{
  res_genes_dt <- data.table::data.table()
  warning("The aberrant splicing pipeline gave 0 results for the ", opt$dataset, " dataset.")
}

# Results
write.table(res_junc_dt, file=paste0(opt$dataset,'_results_per_junction.tsv'), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(res_genes_dt, file=paste0(opt$dataset,'_results.tsv'), quote = FALSE, sep = "\t", row.names = FALSE)

hasExternal <- length(levels(SummarizedExperiment::colData(fds)$isExternal) > 1)
dataset_title <- paste0("Dataset: ", opt$dataset, "--", opt$genome_UCSC)
psiTypes <- "jaccard"

for(type in psiTypes){
  g <- OUTRIDER::plotEncDimSearch(fds, type=type) 
  if (!is.null(g)) {
    g <- g + cowplot::theme_cowplot(font_size = 16) + 
      ggplot2::ggtitle(paste(opt$dataset, "Q estimation", type, sep= ", ")) + ggplot2::theme(legend.position = "none")
    print(g)
  }
}

## We need to recompute for some reason

fds <- FRASER::calculateZscore(fds, type="jaccard")
# Pvalues
fds <- FRASER::calculatePvalues(fds, type="jaccard")
# Adjust Pvalues
fds <- FRASER::calculatePadjValues(fds, type="jaccard")

#' ## Aberrantly spliced genes per sample
OUTRIDER::plotAberrantPerSample(fds, type=psiTypes, padjCutoff = opt$p_adj_cutoff,
                      zScoreCutoff = opt$z_score_cutoff,
                      deltaPsiCutoff = opt$delta_psi_cutoff,
                      aggregate = TRUE, main = dataset_title) + 
cowplot::theme_cowplot(font_size = 16) +
ggplot2::theme(legend.position = "top")

#' ## Batch Correlation: samples x samples
topN <- 30000
topJ <- 10000
anno_color_scheme <- RColorBrewer::brewer.pal(n = 3, name = 'Dark2')[1:2]

for(type in psiTypes){
  ## WAS FAILING WITH NORMALIZED = TRUE
  ## ORIGINAL LINE WAS THE FOLLOWING:
  ## for(normalized in c(F, T)) {
  for(normalized in c(F)){
    hm <- OUTRIDER::plotCountCorHeatmap(
      object=fds,
      type = type,
      logit = TRUE,
      topN = topN,
      topJ = topJ,
      plotType = "sampleCorrelation",
      normalized = normalized,
      genome_version_col = "isExternal",
      genome_version_row = NA,
      sampleCluster = NA,
      minDeltaPsi = minDeltaPsi,
      plotMeanPsi = FALSE,
      plotCov = FALSE,
      genome_version_legend = TRUE,
      genome_version_colors = list(isExternal = c("FALSE" = anno_color_scheme[1],
                                              "TRUE" = anno_color_scheme[2]))
    )
    hm
  }
}

#' ## Results
res <- res_genes_dt
file <- paste0('FRASER_results_', opt$dataset,'--', opt$genome_UCSC,'.tsv')
write.table(res, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
#+ echo=FALSE, results='asis'
cat(paste0("<a href='./", basename(file), "'>Download FRASER results table</a>"))

# round numbers
if(nrow(res_genes_dt) > 0){
  res_genes_dt[, pValue := signif(pValue, 3)]
  res_genes_dt[, padjust := signif(padjust, 3)]
  res_genes_dt[, deltaPsi := signif(deltaPsi, 2)]
  res_genes_dt[, psiValue := signif(psiValue, 2)]
  res_genes_dt[, pValueGene := signif(pValueGene, 2)]
  res_genes_dt[, padjustGene := signif(padjustGene, 2)]
  padjGene_cols <- grep("padjustGene", colnames(res_genes_dt), value=TRUE)
  for(padj_col in padjGene_cols){
      res_genes_dt[, c(padj_col) := signif(get(padj_col), 2)]
  }
}

DT::datatable(
  head(res_genes_dt, 1000),
  caption = 'FRASER results (up to 1,000 rows shown)',
  options=list(scrollX=TRUE),
  escape=FALSE,
  filter = 'top'
)

samples <- unique(res_genes_dt$sampleID)

topSpliceSites <- function(fds, metric, top_N) { # help @alvaro
  # We want the maximum value for that metric
  spliceData <- S4Vectors::mcols(fds)[, metric]
  # save original splice site row
  spliceData <- data.frame(spliceData, seq(1, length(spliceData)))
  colnames(spliceData) <- c(metric, "spliceSiteID")
  # Sort splice sites decreasingly by maximum value for input metric
  res <- spliceData[order(spliceData[, metric], decreasing=TRUE), ]
  return(head(res$spliceSiteID, top_N))
}

jaccardSigSites <- topSpliceSites(fds, "b_jaccard", opt$top_N)
psi3SigSites <- topSpliceSites(fds, "maxDPsi3", opt$top_N)
psi5SigSites <- topSpliceSites(fds, "maxDPsi5", opt$top_N)
## But what about acceptor?? Hope new workflow addresses this (else I'll just to put both, but I liked it better when they were integrated)
thetaSigSites <- topSpliceSites(fds, "maxDThetaDonor", opt$top_N)

#' ### Volcano plot
# set basePlot to FALSE to create an interactive plot

AS_Sample_Overview <- function(sample){
  pdf(file=paste(sample, opt$dataset, "FRASER.pdf",sep="_"))
  VolcanoJ <- FRASER::plotVolcano(fds, sample, type = 'jaccard', basePlot = TRUE,
                      deltaPsiCutoff = opt$delta_psi_cutoff,
                      padjCutoff = opt$p_adj_cutoff)
  # New FRASER doesn't like other splicing metrics, I guess
  # Volcano3 <- FRASER::plotVolcano(fds, sample, type = 'psi3', basePlot = TRUE,
  #                     deltaPsiCutoff = opt$delta_psi_cutoff,
  #                     padjCutoff = opt$p_adj_cutoff)
  # Volcano5 <- FRASER::plotVolcano(fds, sample, type = 'psi5', basePlot = TRUE,
  #                     deltaPsiCutoff = opt$delta_psi_cutoff,
  #                     padjCutoff = opt$p_adj_cutoff)
  # VolcanoZ <- FRASER::plotVolcano(fds, sample, type = 'theta', basePlot = TRUE,
  #                     deltaPsiCutoff = opt$delta_psi_cutoff,
  #                     padjCutoff = opt$p_adj_cutoff)
  print(VolcanoJ)
  graphics.off()
}

#' ### Expression plot
#' ### Expected vs observed PSI (or theta)

AS_SpliceSite_Overview <- function(site){
  pdf(file=paste("splice_site", site, opt$dataset,"FRASER.pdf",sep="_"))
  # exPlot not working
  # expPlotJ <- FRASER::plotExpression(fds, type = 'jaccard', site = site, basePlot = TRUE)
  # expPlot3 <- FRASER::plotExpression(fds, type = 'psi3', site = psi3SigSites, basePlot = TRUE)
  # expPlot5 <- FRASER::plotExpression(fds, type = 'psi5', site = psi5SigSites, basePlot = TRUE)
  # expPlotZ <- FRASER::plotExpression(fds, type = 'theta', site = thetaSigSites, basePlot = TRUE)
  # vsPlot3 <- FRASER::plotExpectedVsObservedPsi(fds, type = 'psi3',
  #                                 idx = psi3SigSites, basePlot = TRUE)
  # vsPlot5 <- FRASER::plotExpectedVsObservedPsi(fds, type = 'psi5',
  #                                 idx = psi5SigSites, basePlot = TRUE)
  vsPlotJ <- FRASER::plotExpectedVsObservedPsi(fds, type = "jaccard", idx = site, basePlot = TRUE)
  # vsPlotZ <- FRASER::plotExpectedVsObservedPsi(fds, type = 'theta',
  #                                 idx = thetaSigSites, basePlot = TRUE)
  # print(expPlotJ)
  print(vsPlotJ)
  graphics.off()
}

BiocParallel::register(BiocParallel::MulticoreParam(opt$cpu))

BiocParallel::bplapply(samples, AS_Sample_Overview)
BiocParallel::bplapply(jaccardSigSites, AS_SpliceSite_Overview)

## Results formatting
aberrants <- get_aberrants(as.data.frame(res_junc_dt), pval_cutoff = opt$p_adj_cutoff)
formatted <- format_aberrants(aberrants)
subset <- subset_types(formatted)
write.table(subset$jaccard, "AS_jaccard_results.tsv", quote = FALSE, row.names = FALSE)
write.table(subset$psi3, "AS_psi3_results.tsv", quote = FALSE, row.names = FALSE)
write.table(subset$psi5, "AS_psi5_results.tsv", quote = FALSE, row.names = FALSE)
write.table(subset$theta, "AS_theta_results.tsv", quote = FALSE, row.names = FALSE)
