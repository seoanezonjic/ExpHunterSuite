#' Function to generate target
#'
#' Load target file if it exists, otherwise use the -C and -T flags. Note target takes precedence over target.
#' @param from_file 
#' @param ctrl_samples
#' @param treat_samples
#' @keywords design
#' @export
#' @examples
#' target_generation()

target_generation <- function(from_file=NULL, ctrl_samples=NULL, treat_samples=NULL){
	target <- NULL
	if(! is.null(from_file)) {
	  target <- read.table(from_file, header=TRUE, sep="\t", check.names = FALSE)
	  # Check there is a column named treat
	  if(! "treat" %in% colnames(target)) {
	    stop(cat("No column named treat in the target file.\nPlease resubmit"))
	  }
	} else {
	  index_control_cols <- unlist(strsplit(ctrl_samples, ",")) 
	  index_treatmn_cols <- unlist(strsplit(treat_samples, ","))
	  # Create the target data frame needed in the calls to the DE detection methods
	  target <- data.frame(sample=c(index_control_cols, index_treatmn_cols), 
	    treat=c(rep("Ctrl", length(index_control_cols)), rep("Treat", length(index_treatmn_cols))))
	}	
	return(target)
}


#' Function to write expression package results
#'
#' Write to disk the results of each diff expression package
#' @param df_list 
#' @param prefix
#' @param root
#' @keywords output
#' @export
#' @examples
#' write_df_list_as_tables()

write_df_list_as_tables <- function(df_list, prefix, root) {
  invisible(
    lapply(1:length(df_list), function(i) {
      pack <- names(df_list)[i]
      folder <- file.path(root, paste0("Results_", pack))
      if(!dir.exists(folder)){dir.create(folder)}
      write.table(df_list[[pack]],
      file=file.path(folder, paste0(prefix, pack, '.txt')), quote=FALSE, col.names = NA, sep="\t")
    })
  )
}