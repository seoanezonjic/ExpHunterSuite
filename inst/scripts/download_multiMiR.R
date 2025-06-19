#!/usr/bin/env Rscript

# library(multiMiR)
# library(optparse)
library(data.table)
library(dplyr)
library(stringr)

option_list <- list(
    optparse::make_option(c("-s", "--chunk_size"), type="integer",default=100,
        help="Database chunks size"),
    optparse::make_option(c("-o", "--output"), type="character",
        help="Tabulated file with information about each sample"),
    optparse::make_option(c("-O", "--organism"), type="character",
        help="Tabulated file with information about each sample"),
    optparse::make_option(c("-c", "--cache_mode"), type="logical", 
        default=FALSE, action = "store_true",
        help=paste0("In this mode each chunk wil be launched in a diferent",
            " job. You must launch this script many times as chunks you have",
            " divided the data. You can incorporate this mode inside a loop",
            " and check if output/temp/finished exist"))
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))


temp_dir <- file.path(opt$output, "temp", opt$organism)
remaining_miRNAs_file <- file.path(temp_dir, "remaining_mirnas.txt")
log_file <- file.path(opt$output, paste0(opt$organism, "_multimir_log"))
# print(opt)
# print(all_organism)
##### Prepare data
if(opt$cache_mode && file.exists(remaining_miRNAs_file)){
    remaining_mirnas <- readLines(remaining_miRNAs_file)
}else{
    miRNA_info <- multiMiR::list_multimir("mirna")
  
    allmiRNAs <- unique(miRNA_info[miRNA_info$org == opt$organism, 
             "mature_mirna_id"])
    remaining_mirnas <- allmiRNAs[allmiRNAs != ""]
        
    if(opt$cache_mode){
        dir.create(temp_dir, recursive=TRUE)
        writeLines(remaining_mirnas, con = remaining_miRNAs_file)
    }
}


### Download data 
org_path <- file.path(temp_dir, opt$organism)
if(!dir.exists(org_path)){
    dir.create(org_path, recursive = TRUE)
}

while(length(remaining_mirnas) > 0){
    mirna_chunk <- unique(sample(remaining_mirnas, opt$chunk_size, 
        replace = TRUE))
    multimir <- multiMiR::get_multimir(mirna = mirna_chunk, table = "all", 
        org = opt$organism, predicted.cutoff.type="n", 
        predicted.cutoff = 1000000000, predicted.site = "all", limit = NULL)
    multimir <- multimir@data
    #str(multimir)
    if(nrow(multimir) == 0) {
        no_downloaded_mirnas <- mirna_chunk
        multimir <- data.frame(mature_mirna_id = no_downloaded_mirnas)
        remaining_mirnas <- remaining_mirnas[! remaining_mirnas %in% 
             no_downloaded_mirnas]
        cat(paste0(paste0(no_downloaded_mirnas, sep = ","),
           " miRNAs could not been found on multiMiR. ", length(remaining_mirnas), 
           " miRNAs remaining..."), file = log_file, append=TRUE, sep = "\n")
    } else {
        downloaded_mirnas <- multimir[,"mature_mirna_id"] %>% unique
        if (length(remaining_mirnas[remaining_mirnas %in%
            downloaded_mirnas]) == 0) {
            #para cuando una query devuelve mirnas con los identificadores
            # cambiados (actualizados o obsoletos)
            remaining_mirnas <- remaining_mirnas[! remaining_mirnas %in% 
            mirna_chunk]   
        } else {
            remaining_mirnas <- remaining_mirnas[! remaining_mirnas %in% 
               downloaded_mirnas]
        }
        cat(paste0(length(downloaded_mirnas), " of ", length(downloaded_mirnas) + length(remaining_mirnas), 
            " miRNAs has been downloaded. ", length(remaining_mirnas), 
            " miRNAs remaining..."), file = log_file, append=TRUE, sep = "\n")
    }

    save(multimir, file = file.path(org_path, paste0("data_", 
        length(list.files(org_path, pattern = ".RData")) + 1, ".RData")))
    
    if(opt$cache_mode){
        writeLines(remaining_mirnas, con = remaining_miRNAs_file)
        stop("Program stopped in a controlled way. Launch it again.")
    }
}

print("Merging results...")    
org_mirna_targets <- lapply(list.files(org_path), function(cache_file){
    load(file.path(org_path, cache_file))
    return(as.data.table(multimir))
})

org_mirna_targets <- as.data.frame(unique(rbindlist(org_mirna_targets, 
            fill = TRUE))) # This unique is recommended by the author

print("Filtering starts")

org_mirna_targets <- org_mirna_targets[,c("target_ensembl", 
    "mature_mirna_acc", "database", "score")] 
## Remove empty strings for miRNAs included in databases but with no data
org_mirna_targets <- org_mirna_targets[org_mirna_targets$target_ensembl != "" & 
           org_mirna_targets$mature_mirna_acc != "" & 
           !is.na(org_mirna_targets$database),]

databases_names <- unique(org_mirna_targets$database)
print("Paste starts")
print(databases_names)

# all_pairs <- as.data.table(org_mirna_targets[, c("target_ensembl", 
#            "mature_mirna_acc")])
# all_pairs <- unique(all_pairs)
# unique_pairs <- paste0(all_pairs$target_ensembl, "_AND_", 
#           all_pairs$mature_mirna_acc)

# print("Unique starts")
# multimir_summary <- data.table(pairs = unique_pairs)

# for (db in databases_names) { #This block is for chnagin data structure
#     print(db)
#     database_pairs <- as.data.table(org_mirna_targets[
#               org_mirna_targets$database == db, c("target_ensembl", 
#                  "mature_mirna_acc")])
#     database_pairs <- unique(database_pairs)
#     database_pairs <- paste0(database_pairs$target_ensembl, "_AND_", 
#               database_pairs$mature_mirna_acc)
#     db_pairs <- multimir_summary$pairs %in% database_pairs
#     print(str(db_pairs))
#     multimir_summary[[db]] <- db_pairs
# } 
# print(str(multimir_summary))

# new_columns <- str_split_fixed(multimir_summary$pairs, "_AND_", n = 2)
# multimir_summary$target_ensembl <- new_columns[,1]
# multimir_summary$mature_mirna_acc <- new_columns[,2]
# all_pair2 <- multimir_summary$pairs
# multimir_summary$pairs <- NULL
# multimir_summary <- as.data.frame(multimir_summary)

save(org_mirna_targets, file = file.path(opt$output, 
    paste0(opt$organism, ".RData")))
# save(multimir_summary, file = file.path(opt$output, paste0("parsed_", 
#                     opt$organism, ".RData")))
#suppressPackageStartupMessages(require(mirbase.db)) 







cat("FINISHED", file=file.path(opt$output, paste0(opt$organism, "_finished")), 
    sep="\n")
# unlink(temp_dir, recursive = TRUE, force = TRUE)
