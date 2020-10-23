#' Hunter installation Function
#'
#' This function allows you to install scripts.
#' @param to Path to folder in which copy Hunter scripts
#' @keywords installation
#' @export
#' @examples
#' install_DEgenes_hunter()
install_DEgenes_hunter <- function(to){
	scripts_path <- file.path(find.package('DEgenesHunter'), 'scripts')
	hunter_scripts <- list.files(scripts_path, pattern = ".R", recursive = FALSE, include.dirs = FALSE)
	for(hunter_script in hunter_scripts){
		file.copy(file.path(scripts_path, hunter_script), to, overwrite = TRUE)
	}
}
