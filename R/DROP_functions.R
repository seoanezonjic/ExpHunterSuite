#' Preprocess GTF for DROP analysis
#' 
#' `process_gtf` takes a GTF file and a directory path. It outputs a txdb,
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
  txdb <- GenomicFeatures::makeTxDbFromGFF(opt$database_gtf) # Adjusted (came from snakemake)
  txdb <- GenomeInfoDb::keepStandardChromosomes(txdb)
  AnnotationDbi::saveDb(txdb, db_file)
  cr_file <- file.path(outdir, paste0(gtf_name, "_count_ranges.rds"))
  count_ranges <- exonsBy(txdb, by = "gene")
  saveRDS(count_ranges, cr_file)
  map_file <- file.path(outdir, paste0(gtf_name, "_gene_name_mapping.tsv"))
  gene_name_mapping <- map_genes()
  write_table(gene_name_mapping, file = map_file, sep = "\t", quote = FALSE,
              row.names = FALSE)
  return("GTF processed successfully!")
}

#' Create gene name mapping file
#' 
#' `map_gnes` takes a GTF file and builds a gene-name mapping data frame out
#' of it.
#' @param gtf GTF file to process.
#' @returns A data frame containing GTF file gene annotation.
#' @examples
#' None yet
#' @export

map_genes <- function(gtf) {
  gtf_name <- basename(tools::file_path_sans_ext(gtf))
  gm_file <- file.path(outdir, paste0(gtf_name, "_gene_name_mapping.tsv"))
  if(file.exists(gm_file)) {
    message("Gene name mapping already exist, skipping step")
  } else {
    gtf_df <- as.data.frame(tracklayer::import(gtf))
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
}
