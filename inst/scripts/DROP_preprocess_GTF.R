#! /usr/bin/env Rscript

devtools::load_all("~aestebanm/dev_R/ExpHunterSuite")
option_list <- list(
  optparse::make_option(c("-g", "--gtf"), type="character", default=NULL,
    help="GTF file for txdb."),
  optparse::make_option(c("-o", "--output_dir"), type="character", default=NULL,
    help="Path to output directory.")
  )
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

preprocessing_results <- preprocess_gtf(opt$gtf)

if(is.null(opt$output_dir)) {
  ### I want this to point do inst/DROP_data, or root/DROP_data in installed package
  output_dir <- system.file("DROP_data", package= "ExpHunterSuite")
} else {
  output_dir <- opt$output_dir
}

gtf_name <- basename(tools::file_path_sans_ext(opt$gtf))
db_file <- file.path(outdir, paste0(gtf_name, "_txdb.db"))
cr_file <- file.path(outdir, paste0(gtf_name, "_count_ranges.rds"))
map_file <- file.path(outdir, paste0(gtf_name, "_gene_name_mapping.tsv"))

AnnotationDbi::saveDb(preprocessing_results$txdb, db_file)
saveRDS(preprocessing_results$count_ranges, cr_file)
write.table(preprocessing_results$gene_name_mapping, file = map_file, sep = "\t", quote = FALSE,
            row.names = FALSE)
