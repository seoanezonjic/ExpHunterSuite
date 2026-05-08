#! /usr/bin/env Rscript


#################################################################################################
## FUNCTIONS
#################################################################################################

load_all_tables <- function(all_files) {
	all_tables <- lapply(all_files, load_table)
	return(all_tables)
}

load_table <- function(path) {
  file <- read.table(path, sep = "\t", header = FALSE, row.names = 1)
	rownames(file) <- lapply(strsplit(rownames(file), "\\."), `[[`, 1) # Split rownames by . and retrieve first element of every resulting sublist (gene ID without version identifier)
  return(file)
}

merge_all_tables <- function(all_tables, tags) {
	merged_tables <- do.call(cbind, all_tables)
	colnames(merged_tables) <- unlist(strsplit(tags, ","))
	return(merged_tables)
}

#################################################################################################
## INPUT PARSING
#################################################################################################

option_list <- list(
  optparse::make_option(c("-i", "--input"), type="character", default=NULL,
                        help="File to process"),
  optparse::make_option(c("-t", "--tags"), type="character", default=NULL,
                        help="")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

opt$input <- unlist(strsplit(opt$input, ","))

#############################
## MAIN
#############################
all_counts_tables <- load_all_tables(opt$input)
merged_tables <- merge_all_tables(all_counts_tables, opt$tags)
merged_tables <- rbind(colnames(merged_tables), merged_tables)
merged_tables <- cbind(rownames(merged_tables), merged_tables)
merged_tables[1, 1] <- "ID"
write.table(merged_tables, file = "", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

