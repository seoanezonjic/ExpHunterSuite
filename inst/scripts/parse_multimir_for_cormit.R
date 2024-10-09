#!/usr/bin/env Rscript

library(dplyr)
library(stringr)
library(data.table)
library(optparse)


parse_validated_db <- function(multimir_summary, org_mirna_targets){
	for (db in c("mirecords", "mirtarbase", "tarbase")) {
		database_pairs <- as.data.table(org_mirna_targets[org_mirna_targets$database == db, c("target_ensembl", "mature_mirna_acc")])
		database_pairs <- unique(database_pairs)
		database_pairs <- paste0(database_pairs$target_ensembl, "_AND_", database_pairs$mature_mirna_acc)
		db_pairs <- multimir_summary$pairs %in% database_pairs
		multimir_summary[[db]] <- db_pairs
	}
	return(multimir_summary)

}

parse_predicted_db <- function(multimir_summary, org_mirna_targets, scale){
	prediction_databases <- c("diana_microt", "elmmo", "microcosm", "miranda","mirdb", "pictar", "pita", "targetscan")
	prediction_databases <- prediction_databases[prediction_databases %in% unique(org_mirna_targets$database)]
	org_mirna_targets_split <- split(org_mirna_targets, org_mirna_targets$database)

	
	parsed_db_pairs <- lapply(prediction_databases, function(db){
		db_info <- org_mirna_targets_split[[db]][, c("target_ensembl", "mature_mirna_acc", "score")]
		db_info$score <- as.numeric(db_info$score)

		if (scale == "db_quantile"){
			db_info <- as.data.frame(scale_quantile(db_info = db_info, database = db))
		} else if (scale == "pairs_quantile"){
			db_info <- scale_r_score(db_info = db_info, database = db)
		} else {
			db_info$raw_score <- db_info$score
		}
		
		database_pairs <- paste0(db_info$target_ensembl, "_AND_", db_info$mature_mirna_acc)
		parsed_multimir_summary <- db_info[match(multimir_summary$pairs, database_pairs), c("raw_score","score")]
		return(parsed_multimir_summary)
	})

	names(parsed_db_pairs) <- prediction_databases 
	for (database in names(parsed_db_pairs)) {
		multimir_summary[, database] <- parsed_db_pairs[[database]]$score
		multimir_summary[, paste0("raw_",database)] <- parsed_db_pairs[[database]]$raw_score
	}
	
	return(multimir_summary)
}

scale_quantile <- function(db_info, database){
	db_info <- as.data.frame(db_info)
	db_info$raw_score <- db_info$score
	if (database %in% c("targetscan", "pita","miranda")){
		db_info$score <- -1 * db_info$score
	}
	db_info$score <- rank(db_info$score) / nrow(db_info)
	return(db_info)
}


scale_r_score <- function(db_info, database){
	
	mirnas_pairs <- split(db_info, db_info$mature_mirna_acc)

	parsed_mirna_pairs <- lapply(mirnas_pairs, function(mirna_pairs){
		mirna_pairs$raw_score <- mirna_pairs$score
		if (database %in% c("targetscan", "pita", "miranda")){
			mirna_pairs$score <- mirna_pairs$score * -1
		}
		mirna_pairs$score <- rank(mirna_pairs$score) / nrow(mirna_pairs)		
		return(mirna_pairs)
	})

	parsed_mirna_pairs <- data.table::rbindlist(parsed_mirna_pairs)

	return(as.data.frame(parsed_mirna_pairs))
}


option_list <- list(
	optparse::make_option(c("-i", "--input"), type= "character", default = ".", 
		help = "Set input folder"),
	optparse::make_option(c("-r", "--report_mode"), type= "logical", action = "store_true", default = FALSE, 
		help = "Activate report mode, load input/parsed_[org].RData and print report."),
	optparse::make_option(c("--organism"), type = "character", default = NULL, 
		help = "Set the model organisms available on multimiR (hsa, mmu or rno)"),
	optparse::make_option(c("-s", "--scale"), type = "character", default = "raw_score", 
		help = "Scaling method can be set 'raw_score', 'pairs_quantile' (specific miRNA quantile) or 'db_quantile' (whole database quantiles). Default=%default"),
	optparse::make_option(c("-o","--output"), type = "character", default = ".", 
		help = "Set the output path. Parsed multimir will be saved in parsed_[organism]_[scale].RData") 
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
load(file.path(opt$input, paste0(opt$organism, ".RData")))

message("Filtering starts")
org_mirna_targets <- org_mirna_targets[org_mirna_targets$target_ensembl != ""   &
 									   org_mirna_targets$mature_mirna_acc != "" & 
 									   !is.na(org_mirna_targets$database),]
databases_names <- unique(org_mirna_targets$database)

all_pairs <- as.data.table(org_mirna_targets[, c("target_ensembl", "mature_mirna_acc")])
all_pairs <- unique(all_pairs)
unique_pairs <- paste0(all_pairs$target_ensembl, "_AND_", all_pairs$mature_mirna_acc)

message("Unique starts")
multimir_summary <- data.table(pairs = unique_pairs)

message("Parsing validated databases")
multimir_summary <- parse_validated_db(multimir_summary = multimir_summary,
										org_mirna_targets = org_mirna_targets)

message("Parsing prediction databases")

multimir_summary <- parse_predicted_db(multimir_summary= multimir_summary, org_mirna_targets = org_mirna_targets, scale = opt$scale)

new_columns <- str_split_fixed(multimir_summary$pairs, "_AND_", n = 2)
multimir_summary$target_ensembl <- new_columns[,1]
multimir_summary$mature_mirna_acc <- new_columns[,2]
multimir_summary$pairs <- NULL
multimir_summary <- as.data.frame(multimir_summary)

# }
if (opt$scale != "raw_score"){
	message("Rendering report")
	rmarkdown::render(file.path(".", 'multiMiR_parsing.Rmd'), 
	                  output_file = file.path(opt$output, paste0("multiMiR_",opt$organism,"_", opt$scale, "_test_report.html")), intermediates_dir = opt$output)
} else {
	message("Scores will not be scaled")
}

message("Saving data")
multimir_summary[,grepl("raw", colnames(multimir_summary))] <- NULL
save(multimir_summary, file = file.path(opt$output, paste0("parsed_", opt$scale, "_", opt$organism, ".RData")))
