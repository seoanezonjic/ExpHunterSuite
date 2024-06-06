#' Preprocess GTF for DROP analysis
#' 
#' `preprocess_gtf` takes a GTF file and a directory path. It outputs a txdb,
#'count ranges, and gene_mapping files, which it writes to provided directory.
#' @importFrom GenomicFeatures makeTxDbFromGFF exonsBy
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom AnnotationDbi saveDb
#' @param gtf GTF file to process.
#' @param outdir Processed files directory.
#' @returns A list containing all three generated objects.
#' @examples
#' gtf <- system.file("extData/testdata", "gencode.v45.toy.annotation.gtf",
#'                     package = "ExpHunterSuite")
#' @export

preprocess_gtf <- function(gtf) {
  message("Generating TxDb from provided GTF file")
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf) # Adjusted (came from snakemake)
  txdb <- GenomeInfoDb::keepStandardChromosomes(txdb)
  message("Calculating count ranges")
  count_ranges <- GenomicFeatures::exonsBy(txdb, by = "gene")
  message("Creating gene-name mapping file")
  gene_name_mapping <- .map_genes(gtf)
  preprocessing_results <- list(txdb = txdb,
                             count_ranges = count_ranges,
                             gene_name_mapping = gene_name_mapping)
  return(preprocessing_results)
}

#' Create gene name mapping file
#' 
#' `map_genes` takes a GTF file and builds a gene-name mapping data frame out
#' of it.
#' @inheritParams preprocess_gtf
#' @returns A data frame containing GTF file gene annotation.
#' @importFrom AnnotationDbi saveDb
#' @importFrom tools file_path_sans_ext
#' @importFrom rtracklayer import

.map_genes <- function(gtf) {
  gtf_name <- basename(tools::file_path_sans_ext(gtf))
  gtf_df <- as.data.frame(rtracklayer::import(gtf))
    if (!"gene_name" %in% colnames(gtf_df)) {
      gtf_df$gene_name <- gtf_df$gene_id
    }
    gtf_df <- gtf_df[gtf_df$type == 'gene',]
    if('gene_biotype' %in% colnames(gtf_df)) {
      colnames(gtf_df)[colnames(gtf_df) == "gene_biotype"] <- "gene_type"
    }
    # Subset to the following columns only
    columns <- c('seqnames', 'start', 'end', 'strand', 'gene_id',
                 'gene_name', 'gene_type', 'gene_status')
    columns <- intersect(columns, colnames(gtf_df))
    gtf_df <- gtf_df[, columns]
    return(gtf_df) 
}

#' Add HPO terms to table
#' 
#' `.add_HPO_cols` takes a data table and adds a hgncSymbol column.
#' @param RES Data table to process.
#' @param sample_id_col Name of sample id column.
#' @param gene_name_col Name of column where hgnc symbols will be included.
#' @param hpo_file File containing HPO terms associated to genes.
#' @importFrom data.table fread setnames
#' @returns Data frame with new hgnc symbol column and (optionally) HPO terms
#' associated to symbols.

.add_HPO_cols <- function(RES, sample_id_col = 'sampleID', 
                         gene_name_col = 'hgncSymbol', hpo_file = NULL){
  
  filename <- file.path(hpo_file) # Adjusted. Instead of trying to download it
                                  # it finds it locally. Same change as in original
                                  # DROP (thanks Jim!)

  hpo_dt <- data.table::fread(filename)
  
  # change column names
  data.table::setnames(RES, old = c(sample_id_col, gene_name_col), new = c('sampleID', 'hgncSymbol'))
  
  # Get HPO terms of all genes
  f2 <- merge(unique(RES[, .(sampleID, hgncSymbol)]), 
              hpo_dt[,.(hgncSymbol, HPO_id, HPO_label, Ngenes_per_pheno, Nphenos_gene)], 
              by = 'hgncSymbol', allow.cartesian = TRUE)
  sa[, HPO_TERMS := gsub(', ', ',', HPO_TERMS)]
  
  if(nrow(f2) > 0){
    f3 <- merge(f2, sa[,.(RNA_ID, HPO_TERMS)], by.x = 'sampleID', by.y = 'RNA_ID')
    f3 <- f3[HPO_TERMS != ''] # remove cases with no HPO terms to speed up
    if(nrow(f3) > 0){
      f3[, HPO_match := any(grepl(HPO_id, HPO_TERMS)), by = 1:nrow(f3)]
      f3 <- f3[HPO_match == TRUE] #only take those that match HPOs
      if(nrow(f3) > 0){
        f4 <- f3[, .(HPO_label_overlap = paste(HPO_label, collapse = ', '),
                     HPO_id_overlap = paste(HPO_id, collapse = ', '),
                     Ngenes_per_pheno = paste(Ngenes_per_pheno, collapse = ', '),
                     Nphenos_gene = unique(Nphenos_gene)), 
                 by = .(sampleID, hgncSymbol)]
        RES <- merge(RES, f4, by = c('sampleID', 'hgncSymbol'), all.x = TRUE)
      }
    }
  }
  
  data.table::setnames(RES, old = c('sampleID', 'hgncSymbol'), new = c(sample_id_col, gene_name_col))
  return(RES)
}

#' Get counts table from rds file or table
#' 
#' `.get_file_counts` takes a file from which to extract a counts table. File may
#' be in RDS or table format.
#' @param file Counts table file.
#' @returns Loaded counts table.
#' @importFrom data.table fread

.get_file_counts <- function(file) {
  print(file)
    if(grepl('rds$', file))
      counts <- assay(readRDS(file))
    else {
      counts <- as.matrix(data.table::fread(file), rownames = "geneID")
    }
  return(counts)
}

#' Get unique rownames of a data frame
#' 
#' `.get_unique_rownames` takes a data frame and extracts its unique row names.
#' @param df A data frame.
#' @returns Unique rownames of data frame.

.get_unique_rownames <- function(df) {
  res <- sort(unique(rownames(df)))
  return(res)
}

#' Get base sample ID from RNA_ID string generated for sample annotation file
#' 
#' `.get_base_ID` takes an RNA_ID string from the DROP sample annotation file.
#' String generated in this way contain not only sample ID, but the name of
#' the program that generated them. This function extracts the sample ID from
#' that string.
#' @param RNA_ID A string containing the sample ID and the program used for
#' read counting, separated by underscores.
#' @returns Sample ID.

.get_base_ID <- function(RNA_ID) {
  base_ID <- .split_string_by_char(RNA_ID, "_", 1)
  return(base_ID)
}

#' Add base ID column to sample annotation data frame
#' 
#' `.add_base_IDs` takes a sample annotation data frame and calls .get_base_ID
#' on it, then assigns the results to new BASE_ID column.
#' @param col_data Data frame containing RNA_ID column.
#' @returns The same dataframe, with a new BASE_ID column containing sample ID.
#' If RNA_ID column only contains sample ID, new column will be an exact copy.

.add_base_IDs <- function(col_data) {
  base_IDs <- sapply(col_data$RNA_ID, .get_base_ID)
  col_data$BASE_ID <- base_IDs
  return(col_data)
}

#' Estimate BCV (Biological Coefficient of Variation) from Outrider Dataset.
#' 
#' `.estimateThetaWithoutAutoCorrect` takes an ODS (Outrider DataSet) object
#' and calculates its biological coefficient of variation, adding it to the ODS.
#' @param ods Outrider DataSet object.
#' @returns An updated ODS containing estimated BCV.
#' @importFrom OUTRIDER OutriderDataSet normalizationFactors counts fit theta
#' @importFrom SummarizedExperiment colData

.estimateThetaWithoutAutoCorrect <- function(ods){
  
  ods1 <- OUTRIDER::OutriderDataSet(
                                   countData = OUTRIDER::counts(ods),
                                   colData = SummarizedExperiment::colData(ods)
                                   )
  # use rowMeans as expected means
  OUTRIDER::normalizationFactors(ods1) <- matrix(
                                            rowMeans(OUTRIDER::counts(ods1)), 
                                            ncol=ncol(ods1), nrow=nrow(ods1)
                                            )
  ods1 <- OUTRIDER::fit(ods1)
  OUTRIDER::theta(ods1)
  
  return(OUTRIDER::theta(ods1))
}

#' Simple wrapper for OUTRIDER::plotVolcano function. Calls it for specified
#' sample in input Outrider DataSet (ods)
#' 
#' `AE_Sample_Overview` takes a sample and an Outrider DataSet (ods), and
#' prints a volcano plot of aberrantly expressed genes in said sample.
#' @param ods Outrider DataSet object.
#' @param sample Sample to plot.
#' @returns Volcano plot of aberrantly expressed genes in specified sample
#' in provided ods.
#' @importFrom OUTRIDER plotVolcano
#' @export

AE_Sample_Overview <- function(sample, ods, dataset, cfg, z_score_cutoff,
                               p_adj_cutoff){
  pdf(file=paste(sample, dataset,"OUTRIDER.pdf",sep="_"))
  plot <- OUTRIDER::plotVolcano(ods, sample, basePlot = TRUE,
                      zScoreCutoff = z_score_cutoff,
                      padjCutoff = p_adj_cutoff)
  print(plot)
  graphics.off()
}

#' Simple wrapper for OUTRIDER::plotExpressionRank and
#' plotExpectedVsObservedCounts functions. Calls them for specified
#' genes in input Outrider DataSet (ods)
#' 
#' `AE_Gene_Overview` takes a gene and an Outrider DataSet (ODS), and
#' prints two plots. The first one is the Expression Rank, which represents
#' normalized counts of specified gene by sample. Samples are sorted in
#' ascending order by normalized counts (i.e expression rank). Outliers are
#' coloured red. ExpectedVsObserved plot creates a graph which represents
#' raw (observed) versus expected counts. Expected counts are an estimation
#' based on OUTRIDER normalization and total expression levels in sample.
#' @param gene Gene to plot.
#' @param ods Outrider DataSet object.
#' @returns Volcano plots of provided gene expression rank and expected vs
#' observed counts across all samples.
#' @importFrom OUTRIDER plotExpressionRank plotExpectedVsObservedCounts
#' @export

AE_Gene_Overview <- function(gene, ods, dataset, cfg){
  pdf(file=paste(gene, dataset,"OUTRIDER.pdf",sep="_"))
  expPlot <- OUTRIDER::plotExpressionRank(ods, gene, basePlot = TRUE,
                      zScoreCutoff = z_score_utoff,
                      padjCutoff = p_adj_cutoff)
  vsPlot <- OUTRIDER::plotExpectedVsObservedCounts(ods, gene, basePlot = TRUE)
  print(expPlot)
  print(vsPlot)
  graphics.off()
}

#' Function to separate locally processed and GTEx samples from OUTRIDER
#' results.
#' `.processed_vs_imported` takes a data frame and splits it into two data frames
#' contained in a list. The two elements are "processed" and "imported".
#' All samples whose name contain "GTEX" are considered imported, the rest
#' are considered exported.
#' @param df A data frame containing read counts by sample.
#' @returns A list with two elements, one for processed counts and another
#' for imported counts.

.processed_vs_imported <- function(df) {
  extIndex <- grep("GTEX", df$sampleID)
  if(any(extIndex)) {
    processed <- df[-grep("GTEX", df$sampleID), ]
    imported <- df[grep("GTEX", df$sampleID), ]
    return(list(processed = processed, imported = imported))
  } else {
    return(list(processed = df))
  }
}

#' Function to identify aberrantly expressed genes in an outrider results
#' data frame.
#' `get_aberrants` takes an outrider results data frame and creates a column
#' named "aberrant", set to FALSE by default. It then compares each gene's
#' zScore and adjusted p-value to the cutoffs provided and sets that gene's
#' aberrant value to TRUE if it passes the check (zScore greater than cutoff and
#' adjusted p-value lower than cutoff)'. It then subsets the input data frame
#' to samples with at least one aberrantly expressed gene and genes aberrantly
#' expressed in at least one sample.
#' @inheritParams main_abgenes_Hunter
#' @param df An OUTRIDER results data frame.
#' @returns A data frame containing only samples with at least one aberrantly
#' expressed genes and genes aberrantly expressed in at least one sample.
#' @export

get_aberrants <- function(df, z_score_cutoff, p_adj_cutoff) {
  if(is.null(df)) {
    return(NULL)
  }
  df$aberrant <- FALSE
  df$aberrant[abs(df$zScore) > z_score_cutoff & df$padjust <= p_adj_cutoff] <- TRUE
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

#' Function to format DROP results
#' `format_aberrants` takes a data frame and subsets it to columns named
#' "sampleID", "geneID", "padjust", "type", "zScore" and "altRatio". It then
#' updates padjust column and takes its negative decimal logarithm, and renames
#' it to p_padjust (as in -log(padjust)).
#' @param input_table A DROP results data frame or data table.
#' @returns A subset of the data frame containing only "sampleID", "geneID",
#' "p_padjust" (calculated from "padjust" column), "type", "zScore" and
#' "altRatio" columns, if they existed in the input data frame.
#' @import data.table
#' @export

format_aberrants <- function(input_table) {
  if(is.null(input_table)) {
    warning("NULL input. This can be due to missing module results or no
      aberrants found.")
    return(NULL)
  }
  names <- c("sampleID", "geneID", "padjust", "type", "zScore", "altRatio")
  matches <- colnames(input_table) %in% names
  if("data.table" %in% class(input_table)){
    res <- input_table[, matches, with = FALSE]
  } else {
    res <- input_table[, matches]
  }
  res$padjust <- -log(res$padjust)
  colnames(res)[colnames(res)=="padjust"] <- "p_padjust"
  return(res)
}

#' Function to parse imported data information in merge_counts function.
#' `.parse_externals` subsets the sample annotation table to entries with
#' "EXTERNAL" column set to "external", and sets up a custom object to allow
#' the incorporation of imported data into OUTRIDER analysis.
#' @inheritParams .get_metadata
#' @param local_counts Counts table to merge with external table
#' @returns A list containing three objects. "run" is a logical that is TRUE
#' if the is an imported counts table. Counts is the counts table extracted
#' from the file. IDs is a vectors containing imported sample IDs.

.parse_externals <- function(sample_anno, local_counts) {
  external_anno <- sample_anno[sample_anno$EXTERNAL=="external", ]
  external_file <- unique(external_anno$GENE_COUNTS_FILE)
  ex_counts <- data.frame(fread(external_file),
                          check.names = FALSE, row.names = "geneID")
  ex_counts <- ex_counts[, colnames(ex_counts)!="Description"]
  exCountIDs <- external_anno$RNA_ID
  if(exCountIDs=="all") {
    # Replace string "all" with sample names from external counts
    exCountIDs <- colnames(ex_counts)
  } else {
    ex_counts <- subset(ex_counts, select = exCountIDs)
  }
  if(!identical(.get_unique_rownames(local_counts),
                .get_unique_rownames(ex_counts))){
                            stop('The rows (genes) of the count matrices to be
                                merged are not the same.')
  }
  return(list(run = nrow(external_anno) > 0, counts = ex_counts,
              IDs = exCountIDs))
}

#' Function to add rowRanges to merged counts table, and imported counts
#' sample information if imported counts table exists in dataset.
#' `.get_metadata` takes a merged counts table, a sample annotation table,
#' a path to a count_ranges file and a logical, and adds some metadata to
#' merged counts table.
#' @inheritParams main_abgenes_Hunter
#' @param total_counts A merged counts table.
#' @param sample_anno A sample annotation table.
#' @param external A boolean.
#'   * `TRUE` : function will check for imported table sample information, and
#'              update the merged_table.
#'   * `FALSE`: function will only set new rowRanges.
#' @returns A merged counts table with updated rowRanges and imported sample
#'          metadata.

.get_metadata <- function(total_counts, sample_anno, count_ranges, external) {
  SummarizedExperiment::rowRanges(total_counts) <- readRDS(count_ranges)
  if(external) {
    col_data[col_data$RNA_ID%in%exCountIDs, ]$EXTERNAL <- "external"
  }
  col_data <- data.table::data.table(RNA_ID = as.character(
                                            colnames(total_counts)))
  col_data <- .add_base_IDs(col_data)
  colnames(sample_anno)[colnames(sample_anno) == "RNA_ID"] <- "BASE_ID"
  col_data <- dplyr::left_join(col_data, sample_anno, by = "BASE_ID")
  rownames(col_data) <- col_data$RNA_ID
  col_data <- as(col_data, "DataFrame")
  rownames(col_data) <- col_data$RNA_ID
  SummarizedExperiment::colData(total_counts) <- col_data
  return(total_counts)
}

#' Function to extract count tables and merge them into one.
#' `merge_counts` takes multiple paths to count tables and reads them, merging
#' them into one.
#' @inheritParams .get_metadata
#' @inheritParams main_abgenes_Hunter
#' @returns A merged counts table.
#' @export

merge_counts <- function(cpu, sample_anno, count_files, count_ranges) {
  BiocParallel::register(BiocParallel::MulticoreParam(cpu))
  count_files <- strsplit(count_files, " ")[[1]]
  local_counts_list <- BiocParallel::bplapply(count_files, .get_file_counts)
  merged_assays <- do.call(cbind, local_counts_list)
  message(paste("read", length(count_files), 'files'))
  external_anno <- sample_anno[sample_anno$EXTERNAL=="external", ]
  if(nrow(external_anno) > 0) {
    external <- .parse_externals(sample_anno, merged_assays)
    merged_assays <- cbind(merged_assays, external$counts)
  } else {
    external <- list(run = FALSE)
  }
  merged_counts_table <- cbind(rownames(merged_assays), merged_assays)
  rownames(merged_counts_table) <- NULL
  colnames(merged_counts_table)[1] <- "gene_ID"
  total_counts <- SummarizedExperiment::SummarizedExperiment(assays=list(
                         counts=as.matrix(merged_assays)))
  total_counts <- .get_metadata(total_counts = total_counts,
                               sample_anno = sample_anno,
                               count_ranges = count_ranges,
                               external = external$run)
  return(total_counts)
}

#' Wrapper for DROP OUTRIDER filtering script.
#' `filter_counts` takes a counts table, a txdb object and a fpkm_cutoff to
#' build an OutriderDataSet object and filter it by fpkm.
#' @inheritParams main_abgenes_Hunter
#' @param counts A counts table.
#' @returns A filtered OutriderDataSet built from the merged table.
#' @export

filter_counts <- function(counts, txdb, fpkm_cutoff) {
  ods <- OUTRIDER::OutriderDataSet(counts)
  col_data <- SummarizedExperiment::colData(ods)
  row_data <- SummarizedExperiment::rowData(ods)
  col_data$EXTERNAL[is.na(col_data$EXTERNAL)] <- "no"
  ods <- OUTRIDER::filterExpression(ods, gtfFile=txdb, filter=FALSE,
                          fpkm_cutoff=fpkm_cutoff, addExpressedGenes=TRUE)

  # add column for genes with at least 1 gene
  row_data$counted1sample = rowSums(assay(ods)) > 0

  col_data$isExternal <- col_data$EXTERNAL=="external"
  return(ods)
}

#' Wrapper for DROP main OUTRIDER script.
#' `run_outrider` takes a filtered Outrider dataset and runs the OUTRIDER
#' algorithm on it.
#' @inheritParams main_abgenes_Hunter
#' @param ods_unfitted The filtered Outrider dataset on which to run OUTRIDER.
#' @returns A fitted outrider dataset.
#' @export

run_outrider <- function(ods_unfitted, implementation, max_dim_proportion) {
  ods_unfitted <- ods_unfitted[
           SummarizedExperiment::mcols(ods_unfitted)$passedFilter, ]

  gr <- unlist(S4Vectors::endoapply(SummarizedExperiment::rowRanges(
                             ods_unfitted), range))
  if(length(gr) > 0){
      rd <- SummarizedExperiment::rowData(ods_unfitted)
      SummarizedExperiment::rowRanges(ods_unfitted) <- gr
      SummarizedExperiment::rowData(ods_unfitted) <- rd
  }

  ods_unfitted <- OUTRIDER::estimateSizeFactors(ods_unfitted)

  b <- min(ncol(ods_unfitted), nrow(ods_unfitted)) / max_dim_proportion
  
  if(max_dim_proportion < 4){
      maxSteps <- 20
  } else {
    maxSteps <- 15
  }

  Nsteps <- min(maxSteps, b)
  pars_q <- unique(round(exp(seq(log(5), log(b), length.out = Nsteps))))
  ods_unfitted <- OUTRIDER::findEncodingDim(ods_unfitted, params = pars_q,
                      implementation = implementation)

  ods <- OUTRIDER::OUTRIDER(ods_unfitted, implementation = implementation)
  message("outrider fitting finished")
  return(ods)   
}

#' Wrapper for DROP OUTRIDER results script.
#' `get_ods_results` extracts the results from an ods object.
#' @inheritParams main_abgenes_Hunter
#' @param ods A fitted outrider dataset.
#' @returns A list with two items. The first one, named "all", contains all
#' results. The second one, "table", only contains significant results.
#' @export

get_ods_results <- function(ods, p_adj_cutoff, z_score_cutoff,
                 gene_mapping_file, sample_anno) {
    OUTRIDER_results_all <- OUTRIDER::results(ods, padjCutoff = p_adj_cutoff,
                                      zScoreCutoff = z_score_cutoff, all = TRUE)
    OUTRIDER_results_all[, foldChange := round(2^l2fc, 2)]
    res <- OUTRIDER_results_all[padjust <= p_adj_cutoff &
                   abs(zScore) > z_score_cutoff]
    gene_annot_dt <- data.table::fread(gene_mapping_file)
    if(!is.null(gene_annot_dt$gene_name)){
      if(grepl('ENSG00', res[1,geneID]) & grepl('ENSG00', gene_annot_dt[1,gene_id])){
        res <- merge(res, gene_annot_dt[, .(gene_id, gene_name)],
                     by.x = 'geneID', by.y = 'gene_id', sort = FALSE, all.x = TRUE)
        data.table::setnames(res, 'gene_name', 'hgncSymbol')
        res <- cbind(res[, .(hgncSymbol)], res[, - 'hgncSymbol'])
      }
    }
    if(!is.null(sample_anno$HPO_TERMS) & nrow(res) > 0){
      if(!all(is.na(sample_anno$HPO_TERMS)) & ! all(sa$HPO_TERMS == '')){
        res <- .add_HPO_cols(res, hpo_file = hpo_file)
      }
    }
    return(list(all = OUTRIDER_results_all,
          table = res))
  }

#' Function to get the Biological Coefficient Variation of an OutriderDataSet
#' object.
#' `get_bcv` extracts the biological coefficient variation of an outrider
#' dataset.
#' @inheritParams get_ods_results
#' @returns A data table containing the biological coefficient variation
#' of the outrider dataset before and after normalisation.
#' @export

get_bcv <- function(ods) {
    before <- data.table::data.table(when = "Before",
                 BCV = 1/sqrt(.estimateThetaWithoutAutoCorrect(ods)))
    after <- data.table::data.table(when = "After",
                                    BCV = 1/sqrt(OUTRIDER::theta(ods)))
    bcv_dt <- rbind(before, after)
    return(bcv_dt)
}

#' Function to merge bam stats.
#' `merge_bam_stats` Merges the bam stat files contained in the specified path.
#' @inheritParams main_abgenes_Hunter
#' @returns A list containing two merged counts matrices, one for local and
#'  one for xternal counts, and a merged bam stats table.
#' @export

merge_bam_stats <- function(ods, stats_path) {
  stats <- list.files(stats_path, pattern=".txt", full.names = TRUE)
  bam_coverage <- lapply(stats, read.table)
  bam_coverage <- do.call(rbind, bam_coverage)
  colnames(bam_coverage) <- c("sampleID", "record_count")
  merged_bam_stats <- .extract_coverage_info(ods = ods,
                                             bam_coverage = bam_coverage)
  return(merged_bam_stats)
}

#' Function to convert aberrant expression results into a format that allows
#' a heatmap report of aberrants found.
#' `format_for_report` Takes an OUTRIDER results table and reformats it,
#' returning only information pertaining genes aberrantly expressed in at least
#' one sample, and only samples containing at least one aberrantly expressed
#' gene. It keeps p-value and z-score information.
#' @inheritParams main_abgenes_Hunter
#' @returns A table containing the merged bam stats.
#' @export

format_for_report <- function(results, z_score_cutoff, p_adj_cutoff) {
  data <- .processed_vs_imported(results)
  aberrants <- lapply(data, get_aberrants, z_score_cutoff, p_adj_cutoff)
  formatted <- lapply(aberrants, format_aberrants)
  return(formatted)
}

#' Function to split a string by a certain character and return only the desired
#' element.
#' `.split_string_by_char` Splits a string by the provided character, and
#' keeps the result specified by the provided index.
#' @param string String to split.
#' @param char Char by which the string will be split.
#' @param index Index of the element after splitting that will be kept
#' @returns The element of the split string specified by provided index.

.split_string_by_char <- function(string, char, index) {
  if (!char %in% strsplit(string, "")[[1]]) {
    warning("WARNING: Attempted to split string by character not present")
  }
  split <- strsplit(x = string, split = char, fixed = TRUE)[[1]]
  if(index > length(split) || index < 1) {
    warning("WARNING: Attempted to split string by out-of-bounds index.")
  }
  return(split[index])
}

#' Extract correlation between samples according to gene counts.
#' `get_counts_correlation` Calculates the correlation between all samples of an
#' ods object according to gene count levels.
#' @inheritParams .estimateThetaWithoutAutoCorrect
#' @param normalized A boolean.
#'   * `TRUE` : raw counts.
#'   * `FALSE`: normalized counts.
#' @returns The correlation data frame, calculated with the spearman method.
#' It constains a column of clusters (calculated with euclidean distance) and a
#' row of EXTERNAL metadata ('yes' for imported counts and 'no' for locally
#' processed counts).
#' @export
get_counts_correlation <- function(ods, normalized) {
  counts <- OUTRIDER::counts(ods, normalized = normalized)
  fc_matrix <- as.matrix(log2(counts + 1))
  counts_corr <- cor(fc_matrix, method = "spearman")
  clust_annotation <- .get_clusters_df(ods = ods, matrix = counts_corr,
                                       groups = "EXTERNAL", nClust = 4)
  external_row <- clust_annotation["EXTERNAL"]
  external_row <- t(rbind(NA, external_row))
  colnames(external_row)[1] <- "nClust"
  res <- cbind(clust_annotation["nClust"], counts_corr)
  res <- rbind(external_row, res)
  rownames(res)[1] <- "EXTERNAL"
  return(res)
}

#' Clusterize correlation data.
#' `.get_clusters_df` Clusterizes a count correlation matrix, adding metadata
#' specified in the ods object from which it was calculated.
#' @inheritParams .estimateThetaWithoutAutoCorrect
#' @param matrix log2 counts correlation matrix.
#' @param groups A string. Column of ods colData to use as secondary grouping
#' atribute. Default: "EXTERNAL".
#' @param nClust Number of clusters in which data should be divided.
#' @returns A data frame. Rows are named after samples. Columns are nClust
#' (cluster where sample belongs) and specified groups column.

.get_clusters_df <- function(ods, matrix, groups = character(),
                             nClust = 4, type = "sample") {
  if(type == "gene") {
    data <- SummarizedExperiment::rowData(ods)
  } else {
    data <- SummarizedExperiment::colData(ods)
  }
  ans <- as.data.frame(data[, groups])
  colnames(ans) <- groups
  clusters <- cutree(hclust(dist(matrix)), nClust)
  ans[["nClust"]] <- as.character(clusters)
  rownames(ans) <- rownames(data)
  return(ans)
}

#' Get gene counts correlations across all samples of an ods object.
#' samples from an ods object. Genes are then clusterized.
#' `get_gene_sample_correlations` Extracts the correlation between counts across
#' all samples and clusterizes it.
#' @inheritParams .estimateThetaWithoutAutoCorrect
#' @param normalized A boolean.
#'   * `TRUE` : raw counts.
#'   * `FALSE`: normalized counts.
#' @returns A data frame containing the correlations between the log2 of counts
#' calculated with the spearman method. It constains a column of clusters
#' (calculated with euclidean distance) and a row of EXTERNAL metadata ('yes'
#' for imported counts and 'no' for locally processed counts).
#' @export

get_gene_sample_correlations <- function(ods, normalized = TRUE, nGenes = 500,
                                       rowCentered = TRUE, bcvQuantile = 0.9) {
  bcv <- 1/sqrt(OUTRIDER::theta(ods))
  SummarizedExperiment::rowData(ods)$BCV <- bcv
  ods_sub <- ods[!is.na(bcv) & bcv > stats::quantile(bcv, probs=bcvQuantile,
                                                     na.rm=TRUE),]
  if(!is.null(nGenes)){
    ods_sub <- ods_sub[BiocGenerics::rank(
                       SummarizedExperiment::rowData(ods_sub)$BCV) <= nGenes,]
  }
  if(normalized){
    fc_mat <- as.matrix(log2(OUTRIDER::counts(
                                                 ods_sub, normalized=TRUE) + 1))
  } else {
    fc_mat <- as.matrix(log2(OUTRIDER::counts(
                                  ods_sub, normalized=FALSE) + 1))
    fc_mat <- t(t(fc_mat)/
                       DESeq2::estimateSizeFactorsForMatrix(fc_mat))
  }
  if(rowCentered) {
    fc_mat <- fc_mat - rowMeans(fc_mat)
  }
  cluster_col <- .get_clusters_df(ods = ods_sub, matrix = fc_mat, 
                                     nClust = 4, type = "gene")
  external_row <- .get_clusters_df(ods = ods_sub, matrix = t(fc_mat), 
                                     groups = "EXTERNAL", nClust = 4)["EXTERNAL"]
  external_row <- t(rbind(NA, external_row))
  colnames(external_row)[1] <- "nClust"
  res <- cbind(cluster_col, fc_mat)
  res <- rbind(external_row, res)
  return(res)
}

.extract_coverage_info <- function(bam_coverage, ods) {
  cnts_mtx <- OUTRIDER::counts(ods, normalized = F)
  local_columns <- SummarizedExperiment::colData(ods)$EXTERNAL == "no"
  cnts_mtx_local <- cnts_mtx[, local_columns]
  rownames(bam_coverage) <- bam_coverage$sampleID
  coverage_df <- data.frame(sampleID = colnames(ods),
                            read_count = colSums(cnts_mtx))
  coverage_df <- merge(bam_coverage, coverage_df, by = "sampleID", sort = FALSE)
  coverage_dt <- data.table::data.table(coverage_df)
  data.table::setorder(coverage_dt, read_count)
  coverage_dt[, count_rank := .I]
  coverage_dt[, counted_frac := read_count/record_count]
  data.table::setorder(coverage_dt, counted_frac)
  coverage_dt[, frac_rank := .I]
  size_factors <- OUTRIDER::sizeFactors(ods)
  locals <- names(size_factors) %in%  colnames(cnts_mtx_local)
  local_size_factors <- size_factors[locals]
  coverage_dt[, size_factors := local_size_factors]
  data.table::setorder(coverage_dt, size_factors)
  coverage_dt[, sf_rank := 1:.N]
  coverage_df <- as.data.frame(coverage_dt)
  return(list(coverage_df = coverage_df, counts_matrix = cnts_mtx,
              local_counts_matrix = cnts_mtx_local, size_factors = size_factors,
              local_size_factors = local_size_factors))
}

filter_matrix <- function(ods, cnts_mtx, cnts_mtx_local, has_external) {
  if(has_external){
    filter_mtx <- list(
      local = cnts_mtx_local,
      all = cnts_mtx,
      `passed FPKM` = cnts_mtx[SummarizedExperiment::rowData(ods)$passedFilter,],
      `min 1 read` = cnts_mtx[MatrixGenerics::rowQuantiles(cnts_mtx, probs = 0.95) > 1, ],
      `min 10 reads` = cnts_mtx[MatrixGenerics::rowQuantiles(cnts_mtx, probs = 0.95) > 10, ]
    )
    filter_dt <- lapply(names(filter_mtx), function(filter_name) {
      mtx <- filter_mtx[[filter_name]]
      data.table::data.table(gene_ID = rownames(mtx), median_counts = rowMeans(mtx), filter = filter_name)
    }) |> data.table::rbindlist()
    filter_dt[, filter := factor(filter, levels = c('local', 'all', 'passed FPKM', 'min 1 read', 'min 10 reads'))]
  } else {
    filter_mtx <- list(
      all = cnts_mtx,
      `passed FPKM` = cnts_mtx[SummarizedExperiment::rowData(ods)$passedFilter,],
      `min 1 read` = cnts_mtx[MatrixGenerics::rowQuantiles(cnts_mtx, probs = 0.95) > 1, ],
      `min 10 reads` = cnts_mtx[MatrixGenerics::rowQuantiles(cnts_mtx, probs = 0.95) > 10, ]
    )
    filter_dt <- lapply(names(filter_mtx), function(filter_name) {
      mtx <- filter_mtx[[filter_name]]
      data.table::data.table(gene_ID = rownames(mtx), median_counts = rowMeans(mtx), filter = filter_name)
    }) |> data.table::rbindlist()
    filter_dt[, filter := factor(filter, levels = c('all', 'passed FPKM', 'min 1 read', 'min 10 reads'))]
  }
  return(as.data.frame(filter_dt))
}


# Data manipulation extracted from OUTRIDER::plotExpressedGenes, WITHOUT
# the plotting part.
get_expressed_genes <- function(ods) {
  exp_genes_cols <- c(sampleID = "sampleID", Rank = "expressedGenesRank",
                      Expressed = "expressedGenes",
                      Union = "unionExpressedGenes",
                      Intersection = "intersectionExpressedGenes", 
                      Passed = "passedFilterGenes")
  col_data <- SummarizedExperiment::colData(ods)
  if(!all(exp_genes_cols %in% names(col_data))) {
      stop("Compute expressed genes first by executing \n\tods <- ",
          "filterExpression(ods, addExpressedGenes=TRUE)")
    }
  dt <- data.table::as.data.table(col_data[, exp_genes_cols])
  colnames(dt) <- names(exp_genes_cols)
  df <- as.data.frame(dt)
  return(df[order(df$Rank), ])
}
