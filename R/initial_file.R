#' Initial file structure
#'
#' @param GEO_name GSE number
#' @param dir Dest dir
#'
#' @return GSE download File structer
#' @export
#'
#' @examples
#' library(Dbio)
#' initial_file("GSE54566", "./")
initial_file <- function(GEO_name = GEO_name, dir = dir) {
  des_dir <- paste0(dir, GEO_name)
  if (!dir.exists(des_dir)) {
    dir.create(des_dir)
  }
  if (!dir.exists(paste0(des_dir, "/1.rawdata"))) {
    dir.create(paste0(des_dir, "/1.rawdata"))
  }
  if (!dir.exists(paste0(des_dir, "/2.output"))) {
    dir.create(paste0(des_dir, "/2.output"))
  }
}
