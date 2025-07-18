#!/usr/bin/env Rscript


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
  devtools::load_all(root_path)
  template_folder <- file.path(root_path, 'inst/templates')
}else{
  require('ExpHunterSuite')
  root_path <- find.package('ExpHunterSuite')
  template_folder <- file.path(root_path, 'templates')
}

option_list <- list(
  optparse::make_option(c("-i", "--input_files"), type="character", default=NULL,
    help="Comma separated input files. Indicate active and supplementary tables throught -A and -S respectively. All files are active by default."),
  optparse::make_option(c("-o", "--output_files"), type="character", 
    default=file.path(getwd(), "results"),
    help="Output path. Default=%default"),
  optparse::make_option(c("-A", "--act_des"), type="character", 
    default=NULL,
    help="Indicate which files are ACTIVE using a comma seppatared string. 
    The file name and the type of data must be sepparated by ':'. 
    The data types are 'c' (quantitative), 's' (quantitative that must be scaled) and 'n' (qualitative) Example: 'file1.txt:n,file2.txt:s,file3.txt:c'"),
  optparse::make_option(c("-S", "--supp_desc"), type="character", 
    default=NULL,
    help="Indicate which files are SUPPLEMENTARY using a comma seppatared string. 
    The file name and the type of data must be sepparated by ':'. 
      The data types are 'c' (quantitative) and 'n' (qualitative) Example: 'file1.txt:n,file2.txt:n,file3.txt:c'"),
  optparse::make_option(c("-s", "--supp_samples"), type="character", 
    default=NULL,
    help="Comma seppatared list of samples to be used as supplementary samples/individuals"),
  optparse::make_option(c("-c", "--hcpc_consol"), type="logical", action = "store_false",
    default=TRUE,
    help="Deactivate HCPC consolidation through k-means."),
  optparse::make_option(c("-n", "--n_clusters"), type="integer",
    default=-1,
    help="Number of HCPC clusters."),
  optparse::make_option(c("--seed"), type="character",
    default=NULL, help="Seed to define in degenes_Hunter.R script. Will affect PCA results. Leave empty (\"\") to use a random seed, else provide an integer.")
 )
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

if(!is.null(opt$seed)) set.seed(opt$seed)

if(!file.exists(opt$output_files))
  dir.create(opt$output_files, recursive = TRUE)

input_files <- split_str(opt$input_files, ",")
input_tables <- lapply(input_files, read.table, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
names(input_tables) <- gsub("\\..*", "", basename(input_files))

act_des <- parse_multivar_input(opt$act_des)
colnames(act_des) <- act_des[1,]

numeric_factors <- NULL
string_factors <- NULL
merged_supp_tables <- NULL
supp_desc <- NULL
if (!is.null(opt$supp_samples)) 
  opt$supp_samples <- split_str(opt$supp_samples, ",")

if (!is.null(opt$supp_desc)){
  supp_desc <- parse_multivar_input(opt$supp_desc)
  supp_tables <- process_supp_files_ind(input_tables, supp_desc)
  qual_supp_tables <- input_tables[names(supp_tables$supp_str_files)]
  for(table in names(qual_supp_tables)) {
    input_tables[[table]] <- data.frame(lapply(qual_supp_tables[[table]],
                    as.factor), row.names = rownames(qual_supp_tables[[table]]))
  }
  numeric_factors <- unlist(lapply(supp_tables[["supp_num_files"]], colnames))
  string_factors <- unlist(lapply(supp_tables[["supp_str_files"]], colnames))
  merged_supp_tables <- merge_all_df(unlist(supp_tables, recursive = FALSE))
  merged_supp_tables$sample <- rownames(merged_supp_tables)
}

pca_res <- lapply(act_des, perform_individual_analysis,
                          all_files = input_tables, 
                          numeric_factors = numeric_factors,
                          string_factors = string_factors, 
                          target = merged_supp_tables,
                          hcpc_consol = opt$hcpc_consol,
                          n_clusters = opt$n_clusters)

pca_res$ind_analysis <- names(pca_res)

if(length(act_des) > 1) {
  pca_res$mfa <- compute_mfa(act_des,supp_desc,input_tables, 
                          hcpc_consol = opt$hcpc_consol,
                          n_clusters = opt$n_clusters)
}

pca_output <- file.path(opt$output_files, "PCA_results")
pca_res$target <- merged_supp_tables
pca_res$act_groups <- colnames(act_des)
pca_res$numeric_factors <- numeric_factors
pca_res$string_factors <- string_factors
dir.create(pca_output)
save(pca_res,act_des, file = file.path(pca_output, "multivar_tmp.rdata"))

for(ind_an in pca_res$ind_analysis){
  write_general_pca(pca_res[[ind_an]], pca_output, tag = paste0(ind_an, "_"))
}

if (!is.null(pca_res$mfa)) {
  mfa_res <-pca_res$mfa
  mfa_res$pca_data <- mfa_res$pca_data$global.pca 
  write_general_pca(mfa_res, pca_output, tag = "mfa_")
}

names(pca_res$act_groups) <- ifelse(act_des["type",] == "n", "mca","pca")
render_multivar_report(multivar_res = pca_res, output_files = opt$output_files, opt = opt,
                    template_folder = template_folder,  string_factors = string_factors, numeric_factors = numeric_factors)
