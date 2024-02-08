#'Main function to parse pseudocounts
#' @param input_folder folder with pseudocounts files
#' @param mapper string pseudomapping program
#' @param output_type string indicate i fthe output must be [G] genes, [T] transcripts or both
#' @param annotation_file gtf or gff path, requested for parsin salmon results
#' @examples
#'	\dontrun{
#'		parse_pseudocounts(input_folder = "Path_to_folder",
#'							mapper = "salmon",
#'							annotation_file = "path_to_gtf")	
#'}
#' @export
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom AnnotationDbi keys select
#' @importFrom tximport tximport
parse_pseudocounts <- function(input_folder, mapper, output_type, annotation_file){
	results <- list()
	file_type <- list("salmon" = "*.sf", 
					  "kallisto" = "*abundance.h5",
					  "stringtie" = "*t_data.ctab",
					  "rsem" = "*.genes.results.gz"
					  )


	if (mapper == "stringtie") {
		stop("StringTie input has not yet been included.")
		#https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#StringTie
	}
	
	tx2gene <- NULL
	if (mapper == "salmon") {
		annot_txdb <- GenomicFeatures::makeTxDbFromGFF(annotation_file)
		db_keys <- AnnotationDbi::keys(annot_txdb,  keytype = "TXNAME")
		tx2gene <- AnnotationDbi::select(annot_txdb, db_keys, "GENEID", "TXNAME")
		tx2gene <- remove_pattern_from_df(tx2gene, "[.][0-9]+")
	}
	input_files <- Sys.glob(file.path(input_folder, "*", file_type[[mapper]]))
	names(input_files) <- basename(Sys.glob(file.path(input_folder, "*")))

	if (grepl("G", output_type)) {
		txi_genes <- tximport::tximport(input_files, type = mapper, tx2gene = tx2gene, ignoreTxVersion = TRUE)

		results[["genes"]] <- txi_genes$counts
	}
	if (grepl("T", output_type)) {
		txi_transcripts <- tximport::tximport(input_files, type = mapper, tx2gene = tx2gene, txOut = TRUE, ignoreTxVersion = TRUE)
		results[["transcripts"]] <- txi_transcripts$counts
	}
	return(results)
}



