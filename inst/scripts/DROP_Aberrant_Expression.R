#! /usr/bin/env Rscript


### ORIGIN:  MergeCounts.R
option_list <- list(
  optparse::make_option(c("-s", "--sample_annotation"), type="character", default=NULL,
    help="Sample annotation table in tsv format."),
  optparse::make_option(c("-c", "--count_ranges"), type="character", default=NULL,
    help="Count ranges R object in rds format."),
  optparse::make_option(c("-g", "--gene_mapping_file"), type="character", default=NULL,
    help="Gene mapping file in tsv format."),
  optparse::make_option(c("-r", "--count_files"), type="character", default=NULL,
    help="Count tables, can be rds objects or text files."),
  optparse::make_option(c("-d", "--dataset"), type="character", default=NULL,
    help="Dataset name."),
  optparse::make_option(c("-n", "--cpu"), type="integer", default=NULL,
    help="Number of CPUs provided to job."),
  optparse::make_option(c("-a", "--anno_database"), type="character", default=NULL,
    help="TxDb annotation database in db format."),
  optparse::make_option(c("-c", "--config_file"), type="character", default=NULL,
    help="Configuration file in yaml format. Will be loaded as a list."),
  optparse::make_option(c("-f", "--hpo_file"), type="character", default=NULL,
    help="Genes and associated HPO terms in compressed tsv format."),
  optparse::make_option(c("-r", "--results_outrider"), type="character", default=NULL,
    help="OUTRIDER analysis results in tsv format."),
  optparse::make_option(c("-i", "--input_bams"), type="character", default=NULL,
    help="Bam stats files in plain text format."),
  optparse::make_option(c("-b", "--bam_stats"), type="character", default=NULL,
    help="Merged Bam Stats in tsv format."),
  optparse::make_option(c("-t", "--top_N"), type="integer", default=10,
    help="Top N genes by adjusted p-value to be selected.")
  )

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

register(MulticoreParam(opt$cpu)) # From snakemake@threads, was 30.

# Read counts

sample_anno <- fread(opt$sample_annotation, # First argument from snakemake@config$sampleAnnotation, the Sample Annotation Table. Needed earlier than
                    colClasses = c(RNA_ID = 'character', DNA_ID = 'character')) # original in order to recognize groups.
count_files <- strsplit(opt$count_files, " ")[[1]]

### Local counts

local_counts_list <- bplapply(count_files, get_counts)

message(paste("read", length(count_files), 'files'))

merged_locals <- do.call(cbind, local_counts_list)

### External counts

external_anno <- sample_anno[sample_anno$EXTERNAL=="external", ]
if(nrow(external_anno) > 0) {
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
  if(!identical(get_unique_rownames(merged_locals), get_unique_rownames(ex_counts))){
    stop('The rows (genes) of the count matrices to be merged are not the same.')
  }
}

# merge counts

if(nrow(external_anno) > 0) {
  merged_assays <- cbind(merged_locals, ex_counts)
} else {
  merged_assays <- merged_locals
}

assays_table <- cbind(rownames(merged_assays), merged_assays)
rownames(assays_table) <- NULL
colnames(assays_table)[1] <- "gene_ID"
write.table(assays_table, 'total_counts.txt', sep = '\t', quote = FALSE, row.names = FALSE)

## merge local and external

total_counts <- SummarizedExperiment(assays=list(counts=as.matrix(merged_assays)))

# assign ranges
rowRanges(total_counts) <- readRDS(opt$count_ranges)

# Add sample annotation data (colData)

col_data <- data.table(RNA_ID = as.character(colnames(total_counts)))
col_data <- add_base_IDs(col_data)
colnames(sample_anno)[colnames(sample_anno) == "RNA_ID"] <- "BASE_ID"
col_data <- left_join(col_data, sample_anno, by = "BASE_ID")
if(nrow(external_anno) > 0)
{
  col_data[col_data$RNA_ID%in%exCountIDs, ]$EXTERNAL <- "external"
}

rownames(col_data) <- col_data$RNA_ID
SummarizedExperiment::colData(total_counts) <- as(col_data, "DataFrame")
rownames(SummarizedExperiment::colData(total_counts)) <- SummarizedExperiment::colData(total_counts)$RNA_ID
# save in RDS format
saveRDS(total_counts, 'total_counts.rds')

# ORIGIN: exportCounts.R

data.table::fwrite(as.data.table(assay(total_counts), keep.rownames = 'geneID'),
       file = paste0(opt$dataset, "_geneCounts.tsv.gz"),
       quote = FALSE, row.names = FALSE, sep = '\t', compress = 'gzip')

# ORIGIN: filterCounts.R

counts <- readRDS(total_counts) # snakemake@input$counts, total_counts.rds.
ods <- OUTRIDER::OutriderDataSet(counts)
SummarizedExperiment::colData(ods)$EXTERNAL[is.na(SummarizedExperiment::colData(ods)$EXTERNAL)] <- "no"
txdb <- AnnotationDbi::loadDb(opt$anno_database)

# filter not expressed genes
fpkmCutoff <- yaml::read_yaml(file=opt$config_file)$aberrantExpression$fpkmCutoff
ods <- OUTRIDER::filterExpression(ods, gtfFile=txdb, filter=FALSE,
                        fpkmCutoff=fpkmCutoff, addExpressedGenes=TRUE)

# add column for genes with at least 1 gene
rowData(ods)$counted1sample = rowSums(assay(ods)) > 0

# External data check
for (i in seq(1, length(SummarizedExperiment::colData(ods)$EXTERNAL))) {
  message(i)
  if (OUTRIDER::SummarizedExperiment::colData(ods)$EXTERNAL[i]=="no"){
    OUTRIDER::SummarizedExperiment::colData(ods)$isExternal[i] <- FALSE
  } else {
    OUTRIDER::SummarizedExperiment::colData(ods)$isExternal[i] <- TRUE
  }
}

# Save the ods before filtering to preserve the original number of genes
ods_unfitted <- ods
saveRDS(ods_unfitted, 'ods_unfitted.rds')

# ORIGIN: runOutrider.R

cfg <- yaml::read_yaml(opt$config_file)
implementation <- cfg$aberrantExpression$implementation # From snakemake@config$aberrantExpression$implementation ("autoencoder")
mp <- cfg$aberrantExpression$maxTestedDimensionProportion # From snakemake@config$aberrantExpression$maxTestedDimensionProportion
ods <- readRDS(opt$ods_unfitted) # From snakemake@input$ods (ods_unfitted.Rds)
register(MulticoreParam(opt$cpu)) # From snakemake@threads, is 30.

## subset filtered
ods_unfitted <- ods_unfitted[mcols(ods_unfitted)$passedFilter, ]

# add gene ranges to rowData
gr <- unlist(S4Vectors::endoapply(rowRanges(ods_unfitted), range))
if(length(gr) > 0){
    rd <- rowData(ods_unfitted)
    rowRanges(ods_unfitted) <- gr
    rowData(ods_unfitted) <- rd
}

ods_unfitted <- OUTRIDER::estimateSizeFactors(ods_unfitted)

## find optimal encoding dimension
a <- 5 
b <- min(ncol(ods_unfitted), nrow(ods_unfitted)) / mp   # N/3

maxSteps <- 15
if(mp < 4){
    maxSteps <- 20
}

Nsteps <- min(maxSteps, b)   # Do at most 20 steps or N/3
# Do unique in case 2 were repeated
pars_q <- round(exp(seq(log(a),log(b),length.out = Nsteps))) %>% unique
ods_unfitted <- OUTRIDER::findEncodingDim(ods_unfitted, params = pars_q, implementation = implementation)

## fit OUTRIDER
ods <- OUTRIDER::OUTRIDER(ods, implementation = implementation)
message("outrider fitting finished")

saveRDS(ods, 'ods_fitted.rds') # snakemake@output$ods, 'ods.Rds'.

# ORIGIN: OUTRIDER_Results.R

cfg <- yaml::read_yaml(opt$config_file)
res <- OUTRIDER::results(ods, padjCutoff = cfg$aberrantExpression$padjCutoff, # snakemake@params$padjCutoff, was '`sm cfg.AE.get("padjCutoff")`'
			   zScoreCutoff = cfg$aberrantExpression$zScoreCutoff, all = TRUE) # snakemake@params$pzScoreCutoff, zScoreCutoff: '`sm cfg.AE.get("zScoreCutoff")`'

# Add fold change
res[, foldChange := round(2^l2fc, 2)] 

# Save all the results and significant ones
OUTRIDER_results_all <- res
saveRDS(OUTRIDER_results_all, paste0(opt$dataset,'_OUTRIDER_results_all.rds')) # snakemake@output$results_all

# Subset to significant results
res <- res[padjust <= cfg$aberrantExpression$padjCutoff &
               abs(zScore) > cfg$aberrantExpression$zScoreCutoff]

gene_annot_dt <- data.table::fread(opt$gene_mapping_file) # snakemake@input$gene_name_mapping
if(!is.null(gene_annot_dt$gene_name)){
  if(grepl('ENSG00', res[1,geneID]) & grepl('ENSG00', gene_annot_dt[1,gene_id])){
    res <- merge(res, gene_annot_dt[, .(gene_id, gene_name)],
                 by.x = 'geneID', by.y = 'gene_id', sort = FALSE, all.x = TRUE)
    data.table::setnames(res, 'gene_name', 'hgncSymbol')
    res <- cbind(res[, .(hgncSymbol)], res[, - 'hgncSymbol'])
  }
}

# Add HPO terms, requires online connection (not anymore! @alvaro) and for there to be annotated HPO terms.
sa <- data.table::fread(opt$sample_annotation, colClasses = c(RNA_ID = 'character', DNA_ID = 'character')) # From snakemake@config$sampleAnnotation.
if(!is.null(sa$HPO_TERMS) & nrow(res) > 0){
  if(!all(is.na(sa$HPO_TERMS)) & ! all(sa$HPO_TERMS == '')){
    res <- add_HPO_cols(res, hpo_file = opt$hpo_file) # From snakemake@params$hpoFile.
  }
}

# Save results
OUTRIDER_results_table <- res
data.table::fwrite(OUTRIDER_results_table, paste0(opt$dataset, '_OUTRIDER_results.tsv'), sep = "\t", quote = F) # From snakemake@output$results

# ORIGIN: OUTRIDER_Summary.R

cfg <- yaml::read_yaml(opt$config_file)

# used for most plots
dataset_title <- paste("Dataset:", paste(opt$dataset, names(cfg$geneAnnotation), sep = '--')) # snakemake@wildcards$dataset es "outrider" y "outrider_external". 'v29' viene de snakemake@wildcards$annotation.

#' Number of samples: `r ncol(ods)`  
#' Number of expressed genes: `r nrow(ods)`  

#'
#' ## Visualize
#' ### Encoding dimension
OUTRIDER::plotEncDimSearch(ods) +
          ggplot2::labs(title = dataset_title) +
          cowplot::theme_cowplot() +
          cowplot::background_grid() +
          ggplot2::scale_color_brewer(palette="Dark2")

#' ### Aberrantly expressed genes per sample
OUTRIDER::plotAberrantPerSample(ods, main = dataset_title, 
                                padjCutoff = cfg$aberrantExpression$padjCutoff,
                                zScoreCutoff = cfg$aberrantExpression$zScoreCutoff)

#' ### Batch correction
#+ countCorHeatmap, fig.height=8, fig.width=8
OUTRIDER::plotCountCorHeatmap(ods, normalized = FALSE, colGroups = "isExternal", colColSet = "Dark2",
                              main = paste0('Raw Counts (', dataset_title, ')'))
OUTRIDER::plotCountCorHeatmap(ods, normalized = TRUE, colGroups = "isExternal", colColSet = "Dark2",
                              main = paste0('Normalized Counts (', dataset_title, ')'))

#' ### Expression by gene per sample
#+ geneSampleHeatmap, fig.height=12, fig.width=8
OUTRIDER::plotCountGeneSampleHeatmap(ods, normalized = FALSE, nGenes = 50, colGroups = "isExternal", colColSet = "Dark2",
                                     main = paste0('Raw Counts (', dataset_title, ')'),
                                     bcvQuantile = .95, show_names = 'row')
OUTRIDER::plotCountGeneSampleHeatmap(ods, normalized = TRUE, nGenes = 50, colGroups = "isExternal", colColSet = "Dark2",
                                     main = paste0('Normalized Counts (', dataset_title,')'),
                                     bcvQuantile = .95, show_names = 'row')

before <- data.table::data.table(when = "Before",
                     BCV = 1/sqrt(estimateThetaWithoutAutoCorrect(ods)))
after <- data.table::data.table(when = "After", BCV = 1/sqrt( theta( ods )))
bcv_dt <- rbind(before, after)

# boxplot of BCV Before and After Autoencoder
#+ BCV, fig.height=5, fig.width=6
ggplot2::ggplot(bcv_dt, ggplot2::aes(when, BCV)) +
                ggplot2::geom_boxplot() +
                ggplot2::theme_bw(base_size = 14) +
                ggplot2::labs(x = "Autoencoder correction",
                              y = "Biological coefficient \nof variation",
                              title = dataset_title)

#' ## Results
res <- OUTRIDER_results_table # Venía de snakemake@input$results

## This block looks more of a report than anything
#' Total number of expression outliers: `r nrow(res)`
#' Samples with at least one outlier gene: `r res[, uniqueN(sampleID)]`  

#'
#' ### Aberrant samples
#' An aberrant sample is one that has more than 0.1% of the expressed genes called as outliers.
if (nrow(res) > 0) {
  ab_table <- res[AberrantBySample > nrow(ods)/1000, .("Outlier genes" = .N), by = .(sampleID)] %>% unique
  if (nrow(ab_table) > 0) {
    data.table::setorder(ab_table, "Outlier genes") 
    DT::datatable(ab_table)
  } else {
    print("no aberrant samples")
  }
} else {
  print('no results')
}

# #' ### Results table (report @alvaro)

# ## Save the results table in the html folder and provide link to download
# file <- paste0(opt$dataset, '_OUTRIDER_results.tsv') # De snakemake@output$res_html
# #+ echo=FALSE, results='asis'
# cat(paste0("<a href='./", basename(file), "'>Download OUTRIDER results table</a>"))

# res[, pValue := format(pValue, scientific = T, digits = 3)]
# res[, padjust := format(padjust, scientific = T, digits = 3)]

# # Move this here, was previously written before the final steps @alvaro
# fwrite(res, file, sep = '\t', quote = F)

# # This is only useful for the report @alvaro

# DT::datatable(
#   head(res, 1000),
#   caption = 'OUTRIDER results (up to 1,000 rows shown)',
#   options=list(scrollX=TRUE),
#   filter = 'top'
# )

## End of report block

# ORIGIN: mergeBamStats.R. Depends on bamStats, bash script.

sa <- read.table(file = opt$sample_annotation, sep = '\t', header = TRUE) # es un data frame
bams_folder <- dirname(opt$input_bams)
bams <- list.files(bams_folder, pattern=".txt", full.names = TRUE)
expTable <- data.frame("sampleID" = rep(NA, length(bams)), "record_count" = NA)

for (i in 1:length(bams)) {
    expTable[i,] <- read.table(bams[i])
}

expTable$record_count <- as.numeric(expTable$record_count)
write.table(expTable, file = paste0(opt$dataset,"_outrider.tsv"), row.names = FALSE, sep = "\t")

# ORIGIN: counting_Summary.R

# Old, was useful when multiple bam_stats files were loaded. Leaving it here just in case.
# bam_folder <- stringr::str_sub(opt$bam_stats,1,-6)
# bam_stats_files <- list.files(bam_folder, pattern=stringr::str_sub(opt$bam_stats,-6,-1))
# parsed_bam_stats <- sapply(bam_stats_files, function(x) paste0(bam_folder,x))
parsed_bam_stats <- opt$bam_stats

ods <- readRDS(opt$input_ods)
has_external <- any(as.logical(SummarizedExperiment::colData(ods)$isExternal))
cnts_mtx_local <- OUTRIDER::counts(ods, normalized = F)[, !as.logical(ods@colData$isExternal)]
cnts_mtx <- OUTRIDER::counts(ods, normalized = F)


#' ## Number of samples:  
#' Local: `r sum(!as.logical(ods@colData$isExternal))`  
#' External: `r sum(as.logical(ods@colData$isExternal))`  
#' 
#' # Count Quality Control
#' 
#' Compare number of records vs. read counts  
#' The `Obtained Read Count Ratio` plot does not include external counts
#' because there are no raw reads to be counted.
#' 
  bam_coverage <- data.frame(data.table::fread(parsed_bam_stats))
  rownames(bam_coverage) <- bam_coverage$sampleID
  coverage_df <- data.frame(sampleID = colnames(ods),
                            read_count = colSums(cnts_mtx))
  coverage_df <- merge(bam_coverage, coverage_df, by = "sampleID", sort = FALSE)
  # read count
  coverage_dt <- data.table::data.table(coverage_df)
  data.table::setorder(coverage_dt, read_count)
  coverage_dt[, count_rank := .I]
  # ratio
  coverage_dt[, counted_frac := read_count/record_count]
  data.table::setorder(coverage_dt, counted_frac)
  coverage_dt[, frac_rank := .I]

  # size factors 
  ods <- OUTRIDER::estimateSizeFactors(ods)
  local_size_factors <- OUTRIDER::sizeFactors(ods)[names(OUTRIDER::sizeFactors(ods)) %in% rownames(bam_coverage)]
  coverage_dt[, size_factors := local_size_factors]
  data.bale::setorder(coverage_dt, size_factors)
  coverage_dt[, sf_rank := 1:.N]

  ### @alvaro This block introduces the concept of sample rank. It seems to be
  ### the order by which the samples are sorted according to read counts,
  ### read count ratio and size factors. Need to properly understand the concept
  ### of size factor. Plots should be improved to include sample IDs.

  p_depth <- ggplot2::ggplot(coverage_dt, ggplot2::aes(x = count_rank, y = read_count)) +
    ggplot2::geom_point(size = 3, show.legend = has_external) +
    cowplot::theme_cowplot() +
    cowplot::background_grid() +
    ggplot2::labs(title = "Read Counts", x="Sample Rank", y = "Reads Counted") +
    ggplot2::ylim(c(0,NA)) +
    ggplot2::scale_color_brewer(palette="Dark2")

  p_frac <- ggplot2::ggplot(coverage_dt, ggplot2::aes(x = frac_rank, y = counted_frac)) +
    ggplot2::geom_point(size = 3, show.legend = has_external) +
    cowplot::theme_cowplot() +
    cowplot::background_grid() +
    ggplot2::labs(title = "Read Count Ratio", x = "Sample Rank", y = "Percent Reads Counted") +
    ggplot2::ylim(c(0,NA)) +
    ggplot2::scale_color_brewer(palette="Dark2")

  #+ QC, fig.height=6, fig.width=12
  cowplot::plot_grid(p_depth, p_frac) 

  p_sf <- ggplot2::ggplot(coverage_dt, ggplot2::aes(sf_rank, local_size_factors)) +
    ggplot2::geom_point(size = 3, show.legend = has_external) +
    ggplot2::ylim(c(0,NA)) +
    cowplot::theme_cowplot() +
    cowplot::background_grid() +
    ggplot2::labs(title = 'Size Factors', x = 'Sample Rank', y = 'Size Factors') +
    ggplot2::scale_color_brewer(palette="Dark2")

  p_sf_cov <- ggplot2::ggplot(coverage_dt, ggplot2::aes(read_count, size_factors)) +
    ggplot2::geom_point(size = 3, show.legend = has_external) +
    ggplot2::ylim(c(0,NA)) +
    cowplot::theme_cowplot() +
    cowplot::background_grid() +
    ggplot2::labs(title = 'Size Factors vs. Read Counts',
         x = 'Reads Counted', y = 'Size Factors') +
    ggplot2::scale_color_brewer(palette="Dark2")

  #+ sizeFactors, fig.height=6, fig.width=12
  cowplot::plot_grid(p_sf, p_sf_cov)

  #' # Filtering
  #' **local**: A pre-filtered summary of counts using only the local (from BAM) counts. Omitted if no external counts  
  #' **all**: A pre-filtered summary of counts using only the merged local (from BAM) and external counts  
  #' **passed FPKM**: Passes the user defined FPKM cutoff in at least 5% of genes  
  #' **min 1 read**: minimum of 1 read expressed in 5% of genes  
  #' **min 10 reads**: minimum of 10 reads expressed in 5% of genes  

  quant <- .95

  if(has_external){
      filter_mtx <- list(
        local = cnts_mtx_local,
        all = cnts_mtx,
        `passed FPKM` = cnts_mtx[SummarizedExperiment::rowData(ods)$passedFilter,],
        `min 1 read` = cnts_mtx[MatrixGenerics::rowQuantiles(cnts_mtx, probs = quant) > 1, ],
        `min 10 reads` = cnts_mtx[MatrixGenerics::rowQuantiles(cnts_mtx, probs = quant) > 10, ]
      )
      filter_dt <- lapply(names(filter_mtx), function(filter_name) {
        mtx <- filter_mtx[[filter_name]]
        data.table::data.table(gene_ID = rownames(mtx), median_counts = rowMeans(mtx), filter = filter_name)
      }) %>% rbindlist
      filter_dt[, filter := factor(filter, levels = c('local', 'all', 'passed FPKM', 'min 1 read', 'min 10 reads'))]
  } else {
      filter_mtx <- list(
        all = cnts_mtx,
        `passed FPKM` = cnts_mtx[SummarizedExperiment::rowData(ods)$passedFilter,],
        `min 1 read` = cnts_mtx[MatrixGenerics::rowQuantiles(cnts_mtx, probs = quant) > 1, ],
        `min 10 reads` = cnts_mtx[MatrixGenerics::rowQuantiles(cnts_mtx, probs = quant) > 10, ]
      )
      filter_dt <- lapply(names(filter_mtx), function(filter_name) {
        mtx <- filter_mtx[[filter_name]]
        data.table::data.table(gene_ID = rownames(mtx), median_counts = rowMeans(mtx), filter = filter_name)
      }) %>% rbindlist
      filter_dt[, filter := factor(filter, levels = c('all', 'passed FPKM', 'min 1 read', 'min 10 reads'))]
  }

  binwidth <- .2
  p_hist <- ggplot2::ggplot(filter_dt, ggplot2::aes(x = median_counts, fill = filter)) +
    ggplot2::geom_histogram(binwidth = binwidth) +
    ggplot2::scale_x_log10() +
    ggplot2::facet_wrap(.~filter) +
    ggplot2::labs(x = "Mean counts per gene", y = "Frequency", title = 'Mean Count Distribution') +
    ggplot2::guides(col = ggplot2::guide_legend(title = NULL)) +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    cowplot::theme_cowplot() +
    ggplot2::theme(legend.position = "none")

  p_dens <- ggplot2::ggplot(filter_dt, ggplot2::aes(x = median_counts, col = filter)) +
    ggplot2::geom_density(ggplot2::aes(y=binwidth * ..count..), size = 1.2) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "Mean counts per gene", y = "Frequency") +
    ggplot2::guides(col = ggplot2::guide_legend(title = NULL)) +
    ggplot2::scale_color_brewer(palette = "Paired") +
    cowplot::theme_cowplot() +
    ggplot2::theme(legend.position = "top",
          legend.justification="center",
          legend.background = element_rect(color = NA))

  #+ meanCounts, fig.height=6, fig.width=12
  cowplot::plot_grid(p_hist, p_dens)

  #' ### Expressed Genes
  exp_genes_cols <- c(Rank = "expressedGenesRank",`Expressed\ngenes` = "expressedGenes", 
                      `Union of\nexpressed genes` = "unionExpressedGenes", 
                      `Intersection of\nexpressed genes` = "intersectionExpressedGenes", 
                      `Genes passed\nfiltering` = "passedFilterGenes")

  expressed_genes <- data.table::as.data.table(SummarizedExperiment::colData(ods)[, exp_genes_cols])
  colnames(expressed_genes) <- names(exp_genes_cols)

  #+ expressedGenes, fig.height=6, fig.width=8
  OUTRIDER::plotExpressedGenes(ods) + 
    cowplot::theme_cowplot() +
    cowplot::background_grid(major = "y") +
    ggplot2::geom_point(data = data.table::melt(expressed_genes, id.vars = c("Rank")),
               ggplot2::aes(x = Rank, y = value, col = variable), show.legend = has_external)

  if(has_external){
      DT::datatable(expressed_genes[order(Rank)], rownames = F)
  } else{
      DT::datatable(expressed_genes[order(Rank), -"Is External"], rownames = F)
  }

  write.table(expressed_genes, paste0(opt$dataset, "_expressed_genes.tsv"), sep="\t", row.names= F)

# ORIGIN: OUTRIDER_Overview.R

#+ include=FALSE
# source(snakemake@input$functions) # html, not needed @alvaro

#+ eval=TRUE, echo=FALSE
# get parameters
cfg <- yaml::read_yaml(opt$config_file)
annotations <- names(cfg$geneAnnotation)
#htmlDir <- snakemake@params$htmlDir

# count_links <- sapply(
#   annotations, function(x) build_link_list(
#     file_paths = file.path(htmlDir, "Counting", x, paste0('Summary_', datasets, '.html')),
#     captions = datasets)
# )

# results_links <- sapply(
#  annotations, function(x) build_link_list(
#     file_paths = file.path(htmlDir, "Outrider", x, paste0('Summary_', datasets, '.html')),
#     captions = datasets)
# )

# ods_links <- build_link_list(snakemake@input$odsFiles)
# results_tables <- build_link_list(snakemake@input$resultTables)

## start html

#'
#' **Datasets:** `r paste(datasets, collapse = ', ')`
#'
#' **Gene annotations:** `r paste(annotations, collapse = ', ')`
#'
#' ## Summaries
#'
#' ### Counts summary
#'
#' `r display_text(caption = 'Gene annotation version ', links = count_links)`
#'
#' ### OUTRIDER summary
#'
#' `r display_text(caption = 'Gene annotation version ', links = results_links)`
#'
#' ## Files
#' `r display_text(caption = 'OUTRIDER datasets (ods)', links = ods_links)`
#' `r display_text(caption = 'Results tables', links = results_tables)`
#'

#' ## Analyze Individual Results
#+ echo=FALSE

#' Display the results table of the first dataset
#+ echo=FALSE
DT::datatable(OUTRIDER_results_table, filter = 'top')

#' Choose genes and samples to plot from significant table.

sortedRes <- OUTRIDER_results_table[order(OUTRIDER_results_table$padj_rank), ]
sigSamples <- unique(sortedRes$sampleID)
sigGenes <- head(unique(sortedRes$geneID), opt$top_N)

BiocParallel::register(BiocParallel::MulticoreParam(opt$cpu))

#' ### Volcano plot
#' setting basePlot = FALSE creates an interactive plot
#' that allows finding the gene(s) of interest

BiocParallel::bplapply(sigSamples,AE_Sample_Overview)

#' ### Gene expression plot (normalized counts)
#' ### Expected vs observed counts
BiocParallel::bplapply(sigGenes,AE_Gene_Overview)

# ORIGIN: format_for_report.R, custom script

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

cfg <- yaml::read_yaml(opt$config_file)
zScoreCutoff <- cfg$aberrantExpression$zScoreCutoff
padjCutoff <- cfg$aberrantExpression$padjCutoff

data <- processed_vs_imported(OUTRIDER_results_all)
aberrants <- lapply(data, get_aberrants)
formatted <- lapply(aberrants, format_aberrants)

write.table(formatted$processed, "processed_AE_results.tsv", quote = FALSE, row.names = FALSE)
write.table(formatted$imported, "imported_AE_results.tsv", quote = FALSE, row.names = FALSE)
