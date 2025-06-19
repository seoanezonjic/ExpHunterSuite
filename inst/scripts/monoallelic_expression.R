#! /usr/bin/env Rscript


##########################################
## LOAD LIBRARIES
##########################################

if( Sys.getenv('DEGHUNTER_MODE') == 'DEVELOPMENT' ){
  # Obtain this script directory
  full.fpath <- normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', 
                  commandArgs())], '='))[2])
  main_path_script <- dirname(full.fpath)
  root_path <- file.path(main_path_script, '..', '..')
  # Load custom libraries
  custom_libraries <- c('main_abgenes_Hunter.R', 'DROP_functions.R')
  for (lib in custom_libraries){
    source(file.path(root_path, 'R', lib))
  }
  template_folder <- file.path(root_path, 'inst', 'templates')
} else {
  require('ExpHunterSuite')
  root_path <- find.package('ExpHunterSuite')
  template_folder <- file.path(root_path, 'templates')
}

## Function definitions

run_deseq <- function(file) {
	mae_counts <- data.table::fread(file, fill=TRUE)
	mae_counts <- mae_counts[contig != '']
	mae_counts[, position := as.numeric(position)]
	# Sort by chr
	mae_counts <- mae_counts[!grep("Default|opcode", contig)]
	mae_counts[,contig := factor(contig, 
	                levels = unique(stringr::str_sort(mae_counts$contig, numeric = TRUE)))]
	mae_counts <- mae_counts[!is.na(position)]
	# Function from tMAE pkg
	rmae <- DESeq4MAE(mae_counts) ## negative binomial test for allelic counts

	### Add AF information from gnomAD
	if (opt$add_af == TRUE) {
	    # obtain the assembly from the config
	    rmae <- add_gnomAD_AF(rmae, genome_assembly = opt$genome_version,
	        max_af_cutoff = opt$max_af, populations = c("AF", "AF_afr", "AF_amr", "AF_eas", "AF_nfe"))
	} else {
	    rmae[, rare := NA]
	}
	return(rmae)
}

run_deseq_qc <- function(file) {
	qc_counts <- data.table::fread(file, fill=TRUE)
	qc_counts <- qc_counts[!is.na(position)]
	# Run DESeq
	rmae <- DESeq4MAE(qc_counts, minCoverage = 10)
	rmae[, RNA_GT := '0/1']
	rmae[altRatio < .2, RNA_GT := '0/0']
	rmae[altRatio > .8, RNA_GT := '1/1']
	rmae[, position := as.numeric(position)]

	# Convert to granges
	qc_gr <- GenomicRanges::GRanges(seqnames = rmae$contig, 
	                 ranges = IRanges::IRanges(start = rmae$position, end = rmae$position), 
	                 strand = '*')
	S4Vectors::mcols(qc_gr) = S4Vectors::DataFrame(RNA_GT = rmae$RNA_GT)
	return(qc_gr)
}

option_list <- list(
  optparse::make_option(c("-m", "--mae_counts"), type="character", default=NULL,
    help="MAE counts file."),
  optparse::make_option(c("-q", "--qc_counts"), type="character", default=NULL,
    help="QC counts file."),
  optparse::make_option(c("-g", "--genome_version"), type="character", default=NULL,
    help="Genome release in UCSC format."),
  optparse::make_option(c("-a", "--add_af"), type="logical", default = FALSE, action = "store_true",
    help="Whether or not to add allelic frequency info sourced from gnomAD."),
  optparse::make_option(c("-x", "--max_af"), type="integer", default = 0.001,
    help="Max allelic frequency cutoff to consider as a rare variant."),
  optparse::make_option(c("-f", "--gene_mapping_file"), type="character", default=NULL,
    help="Gene mapping file in tsv format."),
  optparse::make_option(c("-d", "--dataset"), type="character", default=NULL,
    help="Dataset name."),
  optparse::make_option(c("--allelic_ratio_cutoff"), type="integer", default=0.8,
    help="Min allelic ratio cutoff to consider as allelic imbalace."),
  optparse::make_option(c("-p", "--p_adj_cutoff"), type="integer", default=0.05,
    help="Min adjusted P value cutoff to consider an observation significant."),
  optparse::make_option(c("-c", "--max_var_freq_cohort"), type="integer", default = 0.04,
    help="Max allelic frequency in cohort cutoff to consider as a rare variant."),
  optparse::make_option(c("-s", "--sample_annotation"), type="character", default=NULL,
    help="Sample annotation table in tsv format."),
  optparse::make_option(c("-v", "--qc_vcf"), type="character", default=NULL,
    help="High-quality VCF file to use as golden standard"),
  optparse::make_option(c("-n", "--cpu"), type="integer", default=1,
    help="Number of CPUs provided to job.")
  )

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

## Configuration

HDF5Array::setHDF5DumpDir(file.path(getwd(), "tmp"))
input_files <- strsplit(opt$mae_counts, " ")[[1]]
qc_input_files <- strsplit(opt$qc_counts, " ")[[1]]

BiocParallel::register(BiocParallel::MulticoreParam(opt$cpu))

deseq_list <- BiocParallel::bplapply(input_files, run_deseq)

rmae <- lapply(deseq_list, function(rt){
  # force consistant UCSC chromosome style
  rt <- rt[!grepl("chr",contig),contig:= paste0("chr",contig)]
  return(rt)
}) |> data.table::rbindlist()

# re-factor contig
rmae$contig <- factor(rmae$contig)

# Convert results into GRanges
rmae_ranges <- GenomicRanges::GRanges(seqnames = rmae$contig, 
                       IRanges::IRanges(start = rmae$position, end = rmae$position), strand = '*')

# Read annotation and convert into GRanges
gene_annot_dt <- data.table::fread(opt$gene_mapping_file)
gene_annot_ranges <-  GenomicRanges::GRanges(seqnames = gene_annot_dt$seqnames, 
                             IRanges::IRanges(start = gene_annot_dt$start, end = gene_annot_dt$end), 
                             strand = gene_annot_dt$strand)
gene_annot_ranges <- GenomeInfoDb::keepStandardChromosomes(gene_annot_ranges, pruning.mode = 'coarse')

# Keep the chr style of the annotation in case the results contain different styles
GenomeInfoDb::seqlevelsStyle(rmae_ranges) <- GenomeInfoDb::seqlevelsStyle(gene_annot_ranges)

# Overlap results and annotation
fo <- IRanges::findOverlaps(rmae_ranges, gene_annot_ranges)

# Add the gene names
MAE_annot <- cbind(rmae[S4Vectors::from(fo), ],  gene_annot_dt[S4Vectors::to(fo), .(gene_name, gene_type)])

# Prioritize protein coding genes
MAE_annot <- rbind(MAE_annot[gene_type == 'protein_coding'], 
                   MAE_annot[gene_type != 'protein_coding'])

# Write all the other genes in another column
MAE_annot[, aux := paste(contig, position, sep = "-")]
rvar <- unique(MAE_annot[, .(aux, gene_name)])
rvar[, N := 1:.N, by = aux]

r_other <- rvar[N > 1, .(other_names = paste(gene_name, collapse = ',')), by = aux]
mae_res <- merge(MAE_annot, r_other, by = 'aux', sort = FALSE, all.x = TRUE) 
mae_res[, c('aux') := NULL]
mae_res <- mae_res[, .SD[1], by = .(ID, contig, position)]

# Bring gene_name column front
mae_res <- cbind(mae_res[, .(gene_name)], mae_res[, -"gene_name"])

# Calculate variant frequency within cohort 
maxCohortFreq <- opt$max_var_freq_cohort
mae_res[, N_var := .N, by = .(gene_name, contig, position)]
mae_res[, cohort_freq := round(N_var / data.table::uniqueN(ID), 3)]
mae_res[, rare := (rare | is.na(rare)) & cohort_freq <= maxCohortFreq] 

# Add significance columns
allelicRatioCutoff <- opt$allelic_ratio_cutoff
mae_res[, MAE := padj <= opt$p_adj_cutoff &
      (altRatio >= allelicRatioCutoff | altRatio <= (1-allelicRatioCutoff)) 
    ] 
mae_res[, MAE_ALT := MAE == TRUE & altRatio >= allelicRatioCutoff]
#'
#' Number of samples: `r uniqueN(mae_res$ID)`
#'
#' Number of genes: `r uniqueN(mae_res$gene_name)`
#'
#' Number of samples with significant MAE for alternative events: `r uniqueN(mae_res[MAE_ALT == TRUE, ID])`

#+echo=F

# Save full mae_results zipped
mae_res[, altRatio := round(altRatio, 3)]
data.table::fwrite(mae_res, paste('MAE_results', opt$genome_version, opt$dataset, "all.tsv.gz", sep="_"), sep = '\t', 
       row.names = FALSE, quote = FALSE, compress = 'gzip')

# Save significant mae_results
data.table::fwrite(mae_res[MAE_ALT == TRUE], paste('MAE_results', opt$genome_version, opt$dataset, "alt.tsv", sep="_"), 
       sep = '\t', row.names = FALSE, quote = FALSE)

# Save significant mae_results
data.table::fwrite(mae_res[MAE_ALT == TRUE & rare == TRUE], paste('MAE_results', opt$genome_version, opt$dataset, 'alt_rare.tsv', sep="_"), 
       sep = '\t', row.names = FALSE, quote = FALSE)


# Add columns for plot
mae_res[, N := .N, by = ID]
plot_mae_res <- mae_res[,.(N = .N,
             N_MAE = sum(MAE==TRUE),
             N_MAE_REF=sum(MAE==TRUE & MAE_ALT == FALSE),
             N_MAE_ALT=sum(MAE_ALT == TRUE),
             N_MAE_REF_RARE = sum(MAE ==TRUE & MAE_ALT==FALSE & rare == TRUE),
             N_MAE_ALT_RARE = sum(MAE_ALT ==TRUE & rare ==TRUE)
			 ),by = ID]


# Made changes here so the output table is not messed up. Original format worked great for wBuild, but returned
# janky tables when writing to disk @alvaro

melt_dt <- data.table::melt(plot_mae_res, id.vars = 'ID')
melt_dt[variable == 'N', variable := '>10 counts']
melt_dt[variable == 'N_MAE', variable := 'MAE']
melt_dt[variable == 'N_MAE_REF', variable := 'MAE for REF']
melt_dt[variable == 'N_MAE_ALT', variable := 'MAE for ALT']
melt_dt[variable == 'N_MAE_REF_RARE', variable := 'MAE for REF & rare']
melt_dt[variable == 'N_MAE_ALT_RARE', variable := 'MAE for ALT & rare']

#' 
#' ## Cascade plot 
#' a cascade plot that shows a progmae_ression of added filters  
#'   - >10 counts: only variants supported by more than 10 counts
#'   - +MAE: and shows mono allelic expmae_ression
#'   - +MAE for REF : the monoallelic expmae_ression favors the reference allele
#'   - +MAE for ALT : the monoallelic expmae_ression favors the alternative allele
#'   - rare: 
#'     - if `add_AF` is set to true in config file must meet minimum AF set by the config value `max_AF`
#'     - must meet the inner-cohort frequency `maxVarFreqCohort` cutoff

ggplot2::ggplot(melt_dt, ggplot2::aes(variable, value)) + ggplot2::geom_boxplot() +
  ggplot2::scale_y_log10(limits = c(1,NA)) + ggplot2::theme_bw(base_size = 14) +
  ggplot2::labs(y = 'Heterozygous SNVs per patient', x = '') + 
    ggplot2::annotation_logticks(sides = "l")

#'
#' ## Variant Frequency within Cohort
ggplot2::ggplot(unique(mae_res[,cohort_freq,by =.(gene_name, contig, position)]), ggplot2::aes(x = cohort_freq)) + ggplot2::geom_histogram( binwidth = 0.02)  +
  ggplot2::geom_vline(xintercept = maxCohortFreq, col = "red",linetype="dashed") + ggplot2::theme_bw(base_size = 14) +
  ggplot2::xlim(0,NA) + ggplot2::xlab("Variant frequency in cohort") + ggplot2::ylab("Variants")

#' Median of each category
DT::datatable(melt_dt[, .(median = stats::median(value, na.rm = TRUE)), by = variable])
data.table::fwrite(melt_dt, paste0(opt$dataset, '_medians.tsv'), sep = '\t', 
       row.names = FALSE, quote = FALSE)


# round numbers
if(nrow(mae_res) > 0){
  mae_res[, pvalue := signif(pvalue, 3)]
  mae_res[, padj := signif(padj, 3)]
  mae_res[, log2FC := signif(log2FC, 3)]
}
#' 
#' ## MAE Results table
DT::datatable(
  head(mae_res[MAE_ALT == TRUE], 1000),
  caption = 'MAE mae_results (up to 1,000 rows shown)',
  options=list(scrollX=TRUE),
  filter = 'top'
)

annotations <- opt$genome_version
#htmlDir <- snakemake@params$htmlDir
#resultsDir <- snakemake@params$resultsDir # Not needed, task directory is results directory

# results_links <- sapply(
#   annotations, function(x) build_link_list(
#     file_paths = file.path(htmlDir, paste0(datasets, '--', x, '_results.html')),
#     captions = datasets
#   )
# )

# table_links <- sapply(
#   annotations, function(v) build_link_list(
#     file_paths = file.path(resultsDir, paste0(datasets, '/MAE_results_', v, '.tsv')),
#     captions = paste0(datasets)
#   )
# )

#'
#' **Datasets:** `r paste(datasets, collapse = ', ')`
#'
#' **Gene annotations:** `r paste(annotations, collapse = ', ')`
#'
#' ## MAE results
#' `r display_text(caption = 'Gene annotation version ', links = results_links)`
#'
#' ## Files
#' * [Allelic counts](`r file.path(cfg$root, 'processed_data/mae/allelic_counts/')`)
#' * [Results data tables of each sample (.Rds)](`r file.path(cfg$root, 'processed_results/mae/samples/')`)  

#'
#' `r display_text(caption = 'Significant MAE results tables ', links = table_links)`


#' ## Quality Control: VCF-BAM Matching
#+ eval=TRUE, echo=FALSE
# qc_groups <- sort(snakemake@params$qc_groups)
# qc_links <- build_link_list(
#     file_paths = file.path(htmlDir, paste0('QC/', qc_groups, '.html')),
#     captions = qc_groups
# )

# qc_matrix_links <- build_link_list(
#     file_paths = file.path(snakemake@input$qc_matrix),
#     captions = qc_groups
# )

#' `r display_text(caption = 'QC Overview ', links = qc_links)`
#' `r display_text(caption = 'DNA-RNA matrix ', links = qc_matrix_links)`
#'

#+ eval=TRUE, echo=TRUE
#' ## Analyze Individual Results
# Read the first results table

MAE_Overview <- function(mae_result) {
  print(unique(mae_result$ID))
  if(all(is.na(mae_result$rare))){
    g1 <- plotMA4MAE(mae_result,
                    title=mae_result$ID,
                    padjCutoff = opt$p_adj_cutoff,
                    allelicRatioCutoff = opt$allelic_ratio_cutoff )
    g2 <- plotAllelicCounts(mae_result,
                    padjCutoff = opt$p_adj_cutoff,
                    allelicRatioCutoff = opt$allelic_ratio_cutoff )
  } else {
    g1 <- plotMA4MAE(mae_result, rare_column = 'rare',
                    title=mae_result$ID,
                    padjCutoff = opt$p_adj_cutoff,
                    allelicRatioCutoff = opt$allelic_ratio_cutoff )
    g2 <- plotAllelicCounts(mae_result, rare_column = 'rare',
                    title=mae_result$ID,
                    padjCutoff = opt$p_adj_cutoff,
                    allelicRatioCutoff = opt$allelic_ratio_cutoff )
  }
  pdf(file=paste(mae_result$ID,opt$dataset,"MAE.pdf",sep="_"))
  #' ### MA plot: fold change vs RNA coverage
  #+echo=F
  print(g1)
  #' ### Alternative vs Reference plot
  #+echo=F
  print(g2)
  graphics.off()
}

mae_list <- vector(mode = "list", length = length(unique(mae_res$ID)))
names(mae_list) <- unique(mae_res$ID)
for(sample in unique(mae_res$ID)) {
	mae_list[[sample]] <- subset(mae_res, ID == sample)
}

BiocParallel::register(BiocParallel::MulticoreParam(opt$cpu))
BiocParallel::bplapply(mae_list, MAE_Overview)

deseq_list_qc <- BiocParallel::bplapply(qc_input_files, run_deseq_qc)

qc_list <- vector(mode = "list", length = length(unique(deseq_list_qc$ID)))
names(qc_list) <- unique(deseq_list_qc$ID)
for(sample in unique(deseq_list_qc$ID)) {
	qc_list[[sample]] <- subset(deseq_list_qc, ID == sample)
}

sa <- data.table::fread('filtered_table.tsv', 
            colClasses = c(RNA_ID = 'character', DNA_ID = 'character'))
rna_samples <- sa$RNA_ID
dna_samples <- sa$DNA_ID
vcf_files <- sa$DNA_VCF_FILE

# Read the test vcf as GRanges
gr_test <- GenomicRanges::granges(VariantAnnotation::readVcf(opt$qc_vcf))
S4Vectors::mcols(gr_test)$GT <- "0/0"

# For every DNA, find the overlapping variants with each RNA
N <- length(vcf_files)
lp <- lapply(1:N, function(i){
  # Read sample vcf file
  sample <- as.character(dna_samples[i])
  param <-  VariantAnnotation::ScanVcfParam(fixed=NA, info='NT', geno='GT', samples=sample, trimEmpty=TRUE) 
  vcf_sample <- VariantAnnotation::readVcf(vcf_files[i], param = param, row.names = FALSE)
  # Get GRanges and add Genotype
  gr_sample <- GenomicRanges::granges(vcf_sample)
  
  if(!is.null(VariantAnnotation::geno(vcf_sample)$GT)){
    gt <- VariantAnnotation::geno(vcf_sample)$GT
    gt <- gsub('0|0', '0/0', gt, fixed = TRUE)
    gt <- gsub('0|1', '0/1', gt, fixed = TRUE)
    gt <- gsub('1|0', '0/1', gt, fixed = TRUE)
    gt <- gsub('1|1', '1/1', gt, fixed = TRUE)
  } else if(!is.null(VariantAnnotation::info(vcf_sample)$NT)){
    gt <- VariantAnnotation::info(vcf_sample)$NT
    gt <- gsub('ref', '0/0', gt)
    gt <- gsub('het', '0/1', gt)
    gt <- gsub('hom', '1/1', gt)
 }
  
  S4Vectors::mcols(gr_sample)$GT <- gt
  
  # Find overlaps between test and sample
  gr_res <- data.table::copy(gr_test)
  GenomeInfoDb::seqlevelsStyle(gr_res) <- GenomeInfoDb::seqlevelsStyle(GenomeInfoDb::seqlevelsInUse(gr_sample))[1]
  ov <- IRanges::findOverlaps(gr_res, gr_sample, type = 'equal')
  S4Vectors::mcols(gr_res)[S4Vectors::from(ov),]$GT <- S4Vectors::mcols(gr_sample)[S4Vectors::to(ov),]$GT
  
  # Find similarity between DNA sample and RNA sample
  x <- vapply(deseq_list_qc, function(gr_rna){
    # @alvaro
    # Following line throws:
    # Error in (function (classes, fdef, mtable)  :
    # unable to find an inherited method for function ‘seqinfo’ for signature ‘"data.table"’
    GenomeInfoDb::seqlevelsStyle(gr_rna) <- GenomeInfoDb::seqlevelsStyle(gr_res)[1]

    ov <- IRanges::findOverlaps(gr_res, gr_rna, type = 'equal')
    gt_dna <- gr_res[S4Vectors::from(ov)]$GT
    gt_rna <- gr_rna[S4Vectors::to(ov)]$RNA_GT
    sum(gt_dna == gt_rna) / length(gt_dna)
  }, 1.0)
  return(x)
})

# Create a matrix
mat <- do.call(rbind, lp)
row.names(mat) <- dna_samples
colnames(mat) <- rna_samples
save.image('Testing.RData')
melt_mat <- data.table::as.data.table(reshape2::melt(mat))

identityCutoff <- .85

ggplot2::ggplot(melt_mat, ggplot2::aes(value)) + ggplot2::geom_histogram(fill = 'cadetblue4', binwidth = 0.05, center = .025) + 
  ggplot2::theme_bw(base_size = 14) + 
  ggplot2::labs(x = 'Proportion of matching DNA-RNA variants', y = 'DNA-RNA combinations', title = opt$dataset) + 
  ggplot2::scale_y_log10() + ggplot2::annotation_logticks(sides = "l") + 
  ggplot2::expand_limits(x=c(0,1)) +
  ggplot2::geom_vline(xintercept=identityCutoff, linetype='dashed', color = 'firebrick')


#' ## Identify matching samples

#' Number of RNA samples: `r ncol(qc_mat)`
#'
#' Number of DNA samples: `r nrow(qc_mat)`
#' 
#' Number of samples that match RNA and DNA: `r length(qc_mat[qc_mat > identityCutoff])`
#'
#' Median of proportion of matching variants in matching samples: `r round(median(qc_mat[qc_mat > identityCutoff]), 2)`
#'
#' Median of proportion of matching variants in not matching samples: `r round(median(qc_mat[qc_mat < identityCutoff]), 2)`
#'
#' **Considerations:**
#' On our experience, the median of the proportion of matching variants in matching samples is around 0.95,
#' and the median of the proportion of matching variants in not matching samples is around 0.58.
#' Sometimes we do see some values between 0.7 - 0.85. That could mean that the DNA-RNA combination is 
#' not from the same person, but from a relative. It could also be due to a technical error. For those cases, 
#' check the following:
#' 
#' * RNA sequencing depth (low seq depth that can lead to variants not to be found in the RNA)
#' * Number of variants (too many variants called due to sequencing errors)
#' * Ratio of heterozygous/homozygous variants (usually too many called variants means too many heterozygous ones)
#' * Is the sample a relative of the other?
#' 

sa <- data.table::fread('filtered_table.tsv', 
            colClasses = c(RNA_ID = 'character', DNA_ID = 'character'))[, .(DNA_ID, RNA_ID)]
sa[, ANNOTATED_MATCH := TRUE]
colnames(melt_mat)[1:2] <- c('DNA_ID', 'RNA_ID')

#' ### Samples that were annotated to match but do not 
false_matches <- merge(sa, melt_mat, by = c('DNA_ID', 'RNA_ID'), 
                       sort = FALSE, all.x = TRUE)
false_matches <- false_matches[value < identityCutoff]

#' ### Samples that were not annotated to match but actually do
unexpected_matches <- merge(melt_mat, sa, by = c('DNA_ID', 'RNA_ID'), 
                          sort = FALSE, all.x = TRUE)
unexpected_matches <- unexpected_matches[is.na(ANNOTATED_MATCH) & value > identityCutoff]
unexpected_matches$ANNOTATED_MATCH <- NULL

data.table::fwrite(false_matches, paste0(opt$dataset,'_false_matches.tsv'), sep = '\t', 
       row.names = FALSE, quote = FALSE)
data.table::fwrite(unexpected_matches, paste0(opt$dataset,'_unexpected_matches.tsv'), sep = '\t', 
       row.names = FALSE, quote = FALSE)

get_aberrants <- function(df) {
	if(is.null(df)) {
		return(NULL)
	}
	colnames(df)[colnames(df)=="MAE"] <- "aberrant"
	colnames(df)[colnames(df)=="gene_name"] <- "geneID"	
	colnames(df)[colnames(df)=="ID"] <- "sampleID"
	colnames(df)[colnames(df)=="padj"] <- "padjust"
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

data <- as.data.frame(mae_res)
aberrants <- get_aberrants(data)
formatted <- format_aberrants(aberrants)

write.table(formatted, "MAE_results.tsv", quote = FALSE, row.names = FALSE)
