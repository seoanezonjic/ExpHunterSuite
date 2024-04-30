#' Preprocess GTF for DROP analysis
#' 
#' `preprocess_gtf` takes a GTF file and a directory path. It outputs a txdb,
#'count ranges, and gene_mapping files, which it writes to provided directory.
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom AnnotationDbi saveDb
#' @param gtf GTF file to process.
#' @param outdir Processed files directory.
#' @examples
#' gtf <- system.file("extData/testdata", "gencode.v45.toy.annotation.gtf",
#'                     package = "ExpHunterSuite")
#' @export
preprocess_gtf <- function(gtf, outdir) {
  gtf_name <- basename(tools::file_path_sans_ext(gtf))
  db_file <- file.path(outdir, paste0(gtf_name, "_txdb.db"))
  message("Generating TxDb from provided GTF file")
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf) # Adjusted (came from snakemake)
  txdb <- GenomeInfoDb::keepStandardChromosomes(txdb)
  AnnotationDbi::saveDb(txdb, db_file)
  message("Calculating count ranges")
  cr_file <- file.path(outdir, paste0(gtf_name, "_count_ranges.rds"))
  count_ranges <- GenomicFeatures::exonsBy(txdb, by = "gene")
  saveRDS(count_ranges, cr_file)
  message("Creating gene-name mapping file")
  map_file <- file.path(outdir, paste0(gtf_name, "_gene_name_mapping.tsv"))
  gene_name_mapping <- map_genes(gtf)
  write.table(gene_name_mapping, file = map_file, sep = "\t", quote = FALSE,
              row.names = FALSE)
  return("GTF processed successfully!")
}

#' Create gene name mapping file
#' 
#' `map_genes` takes a GTF file and builds a gene-name mapping data frame out
#' of it.
#' @param gtf GTF file to process.
#' @returns A data frame containing GTF file gene annotation.
#' @export

map_genes <- function(gtf) {
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
#' `add_HPO_cols` takes a data table and adds a hgncSymbol column.
#' @param RES Data table to process.
#' @param sample_id_col Name of sample id column.
#' @param gene_name_col Name of column where hgnc symbols will be included.
#' @param hpo_file File containing HPO terms associated to genes.
#' @returns Data frame with new hgnc symbol column and (optionally) HPO terms
#' associated to symbols.
#' @export

# Function to add HPO terms to results table

add_HPO_cols <- function(RES, sample_id_col = 'sampleID', 
                         gene_name_col = 'hgncSymbol', hpo_file = NULL){
  require(data.table)
  
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
#' `get_counts` takes a file from which to extract a counts table. File may
#' be in RDS or table format.
#' @param file Counts table file.
#' @returns Loaded counts table.
#' @export

get_counts <- function(file) {
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
#' `get_unique_rownames` takes a data frame and extracts its unique row names.
#' @param df A data frame.
#' @returns Unique rownames of data frame.
#' @export

get_unique_rownames <- function(df) {
  res <- sort(unique(rownames(df)))
  return(res)
}

#' Get base sample ID from RNA_ID string generated for sample annotation file
#' 
#' `get_base_ID` takes an RNA_ID string from the DROP sample annotation file.
#' String generated in this way contain not only sample ID, but the name of
#' the program that generated them. This function extracts the sample ID from
#' that string.
#' @param RNA_ID A string containing the sample ID and the program used for
#' read counting, separated by underscores.
#' @returns Sample ID.
#' @export

get_base_ID <- function(RNA_ID) {
  base_ID <- strsplit(RNA_ID, "_")[[1]][1]
  return(base_ID)
}

#' Add base ID column to sample annotation data frame
#' 
#' `add_base_IDs` takes a sample annotation data frame and calls get_base_ID
#' on it, then assigns the results to new BASE_ID column.
#' @param col_data Data frame containing RNA_ID column.
#' @returns The same dataframe, with a new BASE_ID column containing sample ID.
#' If RNA_ID column only contains sample ID, new column will be an exact copy.
#' @export

add_base_IDs <- function(col_data) {
  base_IDs <- sapply(col_data$RNA_ID, get_base_ID)
  col_data$BASE_ID <- base_IDs
  return(col_data)
}

#' Estimate BCV (Biological Coefficient of Variation) from Outrider Dataset.
#' 
#' `estimateThetaWithoutAutoCorrect` takes an ODS (Outrider DataSet) object
#' and calculates its biological coefficient of variation, adding it to the ODS.
#' @param ods Outrider DataSet object.
#' @returns An updated ODS containing estimated BCV.
#' @export

estimateThetaWithoutAutoCorrect <- function(ods){
  
  ods1 <- OutriderDataSet(countData=counts(ods), colData=SummarizedExperiment::colData(ods))
  # use rowMeans as expected means
  normalizationFactors(ods1) <- matrix(rowMeans(counts(ods1)), 
                                       ncol=ncol(ods1), nrow=nrow(ods1))
  ods1 <- fit(ods1)
  theta(ods1)
  
  return(theta(ods1))
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
#' @export

AE_Sample_Overview <- function(ods, sample){
  pdf(file=paste(sample, opt$dataset,"OUTRIDER.pdf",sep="_"))
  plot <- OUTRIDER::plotVolcano(ods, sample, basePlot = TRUE,
                      zScoreCutoff = cfg$aberrantExpression$zScoreCutoff,
                      padjCutoff = cfg$aberrantExpression$padjCutoff)
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
#' @returns Volcano plots of provided gene expression rank and expected vs
#' observed counts across all samples.
#' @export

AE_Gene_Overview <- function(gene){
  pdf(file=paste(gene, opt$dataset,"OUTRIDER.pdf",sep="_"))
  expPlot <- OUTRIDER::plotExpressionRank(ods, gene, basePlot = TRUE,
                      zScoreCutoff = cfg$aberrantExpression$zScoreCutoff,
                      padjCutoff = cfg$aberrantExpression$padjCutoff)
  vsPlot <- OUTRIDER::plotExpectedVsObservedCounts(ods, gene, basePlot = TRUE)
  print(expPlot)
  print(vsPlot)
  graphics.off()
}

#' Function to separate locally processed and GTEx samples from OUTRIDER
#' results.
#' `processed_vs_imported` takes a data frame and splits it into two data frames
#' contained in a list. The two elements are "processed" and "imported".
#' All samples whose name contain "GTEX" are considered imported, the rest
#' are considered exported.
#' @param df A data frame containing read counts by sample.
#' @returns A list with two elements, one for processed counts and another
#' for imported counts.

processed_vs_imported <- function(df) {
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
#' @param df An OUTRIDER results data frame.
#' @param zScoreCutoff Cutoff that gene zScore must surpass in order to be
#' considered aberrantly expressed.
#' @param padjCutoff Adjusted p-value cutoff. Gene adjusted p-value must be
#' lower than this number in order to be considered aberrantly expressed.
#' @returns A data frame containing only samples with at least one aberrantly
#' expressed genes and genes aberrantly expressed in at least one sample.

get_aberrants <- function(df, zScoreCutoff, padjCutoff) {
  if(is.null(df)) {
    return(NULL)
  }
  df$aberrant <- FALSE
  df$aberrant[abs(df$zScore) > zScoreCutoff & df$padjust <= padjCutoff] <- TRUE
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

#' Function to format aberrant
#' `format_aberrants` takes a data frame and subsets it to columns named
#' "sampleID", "geneID", "padjust", "type", "zScore" and "altRatio". It then
#' updates padjust column and takes its negative decimal logarithm, and renames
#' it to p_padjust (as in -log(padjust), as in pH meaning -log[H]).
#' @param df A DROP results data frame.
#' @returns A subset of the data frame containing only "sampleID", "geneID",
#' "p_padjust" (calculated from "padjust" column), "type", "zScore" and
#' "altRatio" columns, if they existed in the input data frame.

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
