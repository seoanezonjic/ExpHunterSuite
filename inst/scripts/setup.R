#! /usr/bin/env Rscript


##########################################
## OPTION PARSER
##########################################

option_list <- list(
  optparse::make_option(c("-f", "--force"), type = "character", default = FALSE, action = "store_true",
            help = "Force reinstallation of already-installed packages"),
  optparse::make_option(c("-d", "--check_dependencies"), type = "character", default = FALSE, action = "store_true",
            help = "Update package dependencies")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

if(!"loupeR" %in% .packages(all =TRUE) | isTRUE(opt$force)) {
  remotes::install_github("10XGenomics/loupeR", dependencies = opt$check_dependencies)
  loupeR::setup()
}
