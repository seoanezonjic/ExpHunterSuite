#!/usr/bin/env Rscript


render_pca_report <- function(multivar_res, output_files,template_folder, multivar_type) {
  source_folder <- find.package("htmlreportR")
  if( Sys.getenv('HTMLREPORTER_MODE') == 'DEVELOPMENT' )
    source_folder <- file.path(source_folder, "inst")

  plotter <- htmlreportR:::htmlReport$new(title_doc = "PCA report", 
                          container = c(multivar_res, list(multivar_type = multivar_type)), 
                          tmp_folder = file.path(normalizePath(output_files), "tmp"),
                          src = source_folder,
                          compress_obj = FALSE,
                          type_index = "menu")
  
  plotter$build(file.path(template_folder, 'main_PCA.txt'))
  plotter$write_report(file.path(opt$output_files, "PCA_report.html"))
}




############################################################
##                      SETUP PROGRAM                     ##
############################################################
options(warn=1)
if( Sys.getenv('DEGHUNTER_MODE') == 'DEVELOPMENT' ){
  # Obtain this script directory
  full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile), 
                 error=function(e) # works when using R CMD
                normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', 
                  commandArgs())], '='))[2]))
  main_path_script <- dirname(full.fpath)
  root_path <- file.path(main_path_script, '..', '..')
  # Load custom libraries
  custom_libraries <- list.files(file.path(root_path, "R"), full.names = TRUE)
    # template_folder <- file.path(root_path, 'inst/templates')
  template_folder <- file.path(root_path, 'inst/templates')

  invisible(sapply(custom_libraries, source)) # source(file.path(root_path, 'R', lib))
  
}else{
  require('ExpHunterSuite')
  root_path <- find.package('ExpHunterSuite')
  # template_folder <- file.path(root_path, 'templates')
  template_folder <- file.path(root_path, 'templates')
}




# Prepare command line input 
option_list <- list(
  optparse::make_option(c("-i", "--input_file"), type="character", default=NULL,
    help="Input file where rows are samples and columns are variables. Samples and variable names must be given as first row and column respectively. To change this behaviour, please see -t option"),
  optparse::make_option(c("-o", "--output_files"), type="character", 
    default=file.path(getwd(), "results"),
    help="Output path. Default=%default"),
   optparse::make_option(c("--supp_file"), type="character", default=NULL,
    help="Additional file with supplementary variables"),
  optparse::make_option(c("-a", "--analysis_type"), type="character", 
    default="pca",
    help="Indicate which type of multivariate analysis must be performed. Default=%default"),
   optparse::make_option(c("-t", "--transpose"), type="logical", action = "store_true", 
    default=FALSE,
    help="If activated, rows are treated as variables and columns as samples."),
   optparse::make_option(c("-S", "--add_qualitative_vars"), type="character", 
    default=NULL,
    help="Variables in the input to be treated as supplementary QUALITATIVE variables"),
  optparse::make_option(c("-N", "--add_quantitative_vars"), type="character", 
    default=NULL,
    help="Variables in the input to be treated as supplementary QUANTITATIVE variables"),
  optparse::make_option(c("-s", "--add_samples"), type="character", 
    default=NULL,
    help="Comma seppatared list of samples to be used as supplementary samples/individuals")
 )
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))


input_file <- read.table(file.path(opt$input_file), header = TRUE, row.names = 1, check.names = F)
target <- NULL 
if (!is.null(opt$supp_file)) {
  target <- read.table(opt$supp_file, header = TRUE)
}
if (!is.null(opt$add_quantitative_vars)) {
	opt$add_quantitative_vars <- split_str(opt$add_quantitative_vars, ",")
}

if (!is.null(opt$add_qualitative_vars)) {
	opt$add_qualitative_vars <- split_str(opt$add_qualitative_vars, ",")
}

if (!is.null(opt$add_samples)) {
	opt$add_samples <- split_str(opt$add_samples, ",")
}

if(!file.exists(opt$output_files))
	dir.create(opt$output_files, recursive = TRUE)

if (opt$analysis_type == "pca") {
	pca_res <- compute_pca(pca_data = input_file,
							transpose = opt$transpose,
							string_factors = opt$add_qualitative_vars, 
							numeric_factors = opt$add_quantitative_vars,
							add_samples = opt$add_samples,
              target = target)

} else if (opt$analysis_type == "mca") {
  pca_res <- compute_mca(mca_data = input_file,
              transpose = opt$transpose,
              string_factors = opt$add_qualitative_vars, 
              numeric_factors = opt$add_quantitative_vars,
              add_samples = opt$add_samples, 
              target = target)
}
	pca_output <- file.path(opt$output_files, "PCA_results")
	dir.create(pca_output)
  save(pca_res, file = file.path(pca_output, "pca_tmp.rdata"))
	write_general_pca(pca_res, pca_output)
  write.table(pca_res$pca_data$var$eta2, file = file.path(pca_output, "pca_vars.txt"), quote = F, sep = "\t")
  pca_res$string_factors <-opt$add_qualitative_vars
  pca_res$numeric_factors <- opt$add_quantitative_vars

  render_pca_report(multivar_res = pca_res, output_files = opt$output_files, 
                    template_folder = template_folder, multivar_type = opt$analysis_type)

