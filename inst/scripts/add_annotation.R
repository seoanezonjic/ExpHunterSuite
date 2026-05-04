#!/usr/bin/env Rscript


option_list <- list(
  optparse::make_option(c("-i", "--input_file"), type="character", default=NULL,
                        help="Table to annotate"),
  optparse::make_option(c("-o", "--output_file"), type="character", default='Annotated_table.txt',
                        help="Define the output path. Default = %default"),
  optparse::make_option(c("-c", "--column"), type="character", default=1,
                        help="Column name or index or 'rownames' with IDs to translate. Default = %default"),
  optparse::make_option(c("-m", "--mirna"), type= "logical", default = FALSE, action ="store_true",
                        help= "Indicate if the ids to translate are miRNA (from miRBase to mature ID). If active, -O is assumed to be \"human\". Combine with -l to see valid keytypes. Set flag -I to \"NAME\" to retrieve accession IDs for input_keytype"),
  optparse::make_option(c("-I", "--input_keytype"), type="character", default=NULL,
                        help="Set the input keytype (use flag -l to list valid keytypes). Default=%default: All valid keytypes will be printed"),
  optparse::make_option(c("-K", "--output_keytype"), type="character", default=NULL,
                        help="Set the output keytype (use flag -l to list valid keytypes)"),
  optparse::make_option(c("-O", "--organism"), type="character", default="Human",
                        help="Set the model organism. Default = %default"),
  optparse::make_option(c("-l", "--list_keytypes"), type="logical", default=FALSE,
                        action = "store_true", help="List valid keytypes and exit. Shows mRNA keytypes by default, combine with flag --mirna to see valid miRNA keytypes")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

if(opt$list_keytypes) {
	options(show.error.messages = FALSE)
	if(!opt$mirna) {
		valid_keytypes <- paste(AnnotationDbi::columns(org.Hs.eg.db::org.Hs.eg.db), collapse = ", ")
	} else {
		valid_keytypes <- paste(miRBaseVersions.db::keytypes(miRBaseVersions.db::miRBaseVersions.db), collapse = ", ")
	}
	message("Valid keytypes: ", valid_keytypes)
	stop()
}

options(warn=1)
if( Sys.getenv('DEGHUNTER_MODE') == 'DEVELOPMENT' ){
  # Obtain this script directory
  full.fpath <- normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', 
                  commandArgs())], '='))[2])

  main_path_script <- dirname(full.fpath)
  root_path <- file.path(main_path_script, '..', '..')
  organisms_table_file <- file.path(root_path, "inst","external_data", 
        "organism_table.txt")
  # Load custom libraries
  devtools::load_all(root_path)
  source_folder <- file.path(root_path, 'inst')
}else{
  require('ExpHunterSuite')
  root_path <- find.package('ExpHunterSuite')
  source_folder <- file.path(root_path)
  organisms_table_file <- file.path(root_path, "external_data", 
        "organism_table.txt")
}

organisms_table <- get_organism_table(organisms_table_file)
current_organism_info <- organisms_table[rownames(organisms_table) %in% opt$organism,]
org_db <- get_org_db(current_organism_info)

if (grepl("[0-9]+", opt$column)) opt$column <- as.numeric(opt$column)

table_to_annot <- read.table(opt$input_file, sep = "\t", header = TRUE)

if (opt$column == "rownames") {
  ids_to_translate <- rownames(table_to_annot)
} else {
  ids_to_translate <- table_to_annot[,opt$column]
}

if (opt$mirna) {
  output_col <- "miRNA_name"
  if(opt$input_keytype == "NAME") output_col <- "MIMAT"
  table_to_annot[, output_col] <- translate_miRNA_ids(ids_to_translate, input_keytype = opt$input_keytype, output_keytype = opt$output_keytype)

} else {
  translated_keytypes <- translate_ids_orgdb(ids = ids_to_translate, 
                      input_id=opt$input_keytype, output_id = opt$output_keytype, org_db=org_db)
  translated_ids <- translated_keytypes[match(ids_to_translate, 
                                          translated_keytypes[,opt$input_keytype]) ,opt$output_keytype]
  table_to_annot[,opt$output_keytype] <- translated_ids
  if(any(is.na(table_to_annot$SYMBOL))) {
    NAs <- which(is.na(table_to_annot$SYMBOL))
    id_warn <- paste(length(NAs), opt$input_keytype, "IDs could not be translated")
    warning(id_warn)
    table_to_annot$SYMBOL[NAs] <- ids_to_translate[NAs]
  }
}

write.table(table_to_annot, sep = "\t", quote = FALSE, row.names = (opt$column == "rownames"), file = opt$output_file)
