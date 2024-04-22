#! /usr/bin/env Rscript


option_list <- list(
  optparse::make_option(c("-g", "--gtf"), type="character", default=NULL,
    help="GTF file for txdb."),
  optparse::make_option(c("-o", "--output_dir"), type="character", default=NULL,
    help="Path to output directory.")
  )
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

if(is.null(opt$output_dir)) {
  ### I want this to point do inst/DROP_data, or root/DROP_data in installed package
  output_dir <- system.file("DROP_data", package= "ExpHunterSuite")
} else {
  output_dir <- opt$output_dir
}

process_gtf(opt$gtf, output_dir)
