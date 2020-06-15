#!/usr/bin/env Rscript

library(multiMiR)
library(optparse)
library(data.table)

option_list <- list(
	make_option(c("-s", "--chunk_size"), type="integer",default=100,
		help="Database chunks size"),
	make_option(c("-o", "--output"), type="character",
		help="Tabulated file with information about each sample"),
	make_option(c("-O", "--organism"), type="character",
		help="Tabulated file with information about each sample"),
	make_option(c("-c", "--cache_mode"), type="logical", default=FALSE, action = "store_true",
		help="In this mode each chunk wil be launched in a diferent job. You must launch this script many times as chunks you have divided the data. You can incorporate this mode inside a loop and check if output/temp/finished exist")
)

opt <- parse_args(OptionParser(option_list=option_list))


temp_dir <- file.path(opt$output, "temp")
cache <- file.path(temp_dir, "cache.RData")
# print(opt)
# print(all_organism)
##### Prepare data
if(file.exists(cache) && opt$cache_mode){
	load(cache)
}else{
	multimir_info <- list_multimir()
	# list_organism <- list(rep(0, length(all_organism)))
 # 	names(list_organism) <- all_organism
	org_multimir <- multimir_info[multimir_info$org == opt$organism, ]
	org_multimir_list <- unique(org_multimir$mature_mirna_id[org_multimir$mature_mirna_id != ""])
	splitted_list <- split(org_multimir_list, ceiling(seq_along(org_multimir_list)/opt$chunk))
	all_multimir_info <- list(
		"splitted_list" = splitted_list,
		"original_chunks" = length(splitted_list),
		"remaining_chunks" = length(splitted_list)
	)
	
	if(opt$cache_mode){
		dir.create(temp_dir, recursive=TRUE)
		save(opt, all_multimir_info, file = cache)
	}
}


### Download data 
org_path <- file.path(temp_dir, opt$organism)
if(!dir.exists(org_path)){
	dir.create(org_path, recursive = T)
}
remaining_chunks <- all_multimir_info[["remaining_chunks"]] 
while(remaining_chunks != 0){
	mirna_chunk <- all_multimir_info[["splitted_list"]][[remaining_chunks]]
	multimir <- get_multimir(mirna = mirna_chunk, table = "all")
	multimir <- multimir@data
	save(multimir, file = file.path(org_path, paste0("data_", remaining_chunks, ".RData")))
	all_multimir_info[["splitted_list"]][remaining_chunks] <- NULL
	 remaining_chunks <- remaining_chunks - 1
	if(opt$cache_mode){
		all_multimir_info[["remaining_chunks"]] <- remaining_chunks 
		save(opt, all_multimir_info, file = cache)
		q()
	}
}
	

last_chunk <- all_multimir_info[["original_chunks"]]

org_mirna_targets <- lapply(seq(1, last_chunk), function(i_chunk){
	load(file.path(org_path, paste0("data_", i_chunk, ".RData")))
	return(multimir)
})

org_mirna_targets <- rbindlist(org_mirna_targets, fill = TRUE)
save(org_mirna_targets, file = file.path(opt$output, paste0(opt$organism, ".RData")))
system(paste0("rm -r ", temp_dir))
cat("FINISHED", file=file.path(opt$output, "finished") ,sep="\n")
