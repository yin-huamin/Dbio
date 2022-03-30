#' Get matrixs data from GEO
#'
#' @param GEO_name GSE Number
#' @param dir Dest dir
#'
#' @return
#' @export
#'
#' @examples
#' library(Dbio)
#' get_matrixs("GSE54566", "./")
get_matrixs <- function(GEO_name = GEO_name, dir = dir) {

  Dbio::initial_file(GEO_name = GEO_name, dir = dir)

  if (requireNamespace("GEOquery", quietly = TRUE)) {
    des_dir <- paste0(dir, GEO_name)
    gset <- GEOquery::getGEO(
      GEO = GEO_name,
      destdir = paste0(des_dir, "/1.rawdata"),
      AnnotGPL = F,
      getGPL = F
    )
    expmat <- Biobase::exprs(gset[[1]])
  } else {
    stop(
      "Package \"GEOquery\" must be installed to use this function.",
      call. = FALSE
    )
  }
  return(expmat)
}
