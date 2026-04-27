#! /usr/bin/env Rscript

#################################################################################################
## FUNCTIONS
#################################################################################################

load_file <- function(path, stranded) {
	file <- read.table(path, sep = "\t", header = FALSE, row.names = 1)
	if(stranded == "no") {
		res <- data.frame(file[, 2, drop = FALSE])
	} else {
		res <- as.data.frame(apply(file[, -1], 1, max))
	}
	colnames(res) <- NULL
	return(res)
}

#################################################################################################
## INPUT PARSING
#################################################################################################

option_list <- list(
  optparse::make_option(c("-i", "--input_file PATH"), type="character", default=NULL,
                        help="File to process"),
  optparse::make_option(c("-s", "--stranded"), type="character", default='no',
                        help="Strand attribute to select column counts. Default \'no\'")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

selected_counts <- load_file(opt$input, opt$stranded)
write.table(selected_counts, col.names = FALSE, quote = FALSE, sep = "\t", file = stdout())
