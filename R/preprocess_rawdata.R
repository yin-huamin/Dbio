#' Process raw data from GEO
#'
#' @param GEO_name GSE number
#' @param annotation Probe anotation platform
#' @param arrayType one of "AFFY", "AGI4-44K"
#' @param method.probe one of "medain", "mean", "max", "min", default is "median"
#' @param dir Des dir
#'
#' @return
#' @export
#'
#' @examples
#' library(Dbio)
#' # make file structure
#' initial_file("GSE54566", "./")
#' # get probe annotation data
#' annotation <- get_probe2symbol("GPL570", "GSE54566", "./")
#' preprocess_rawdata(
#'   GEO_name = "GSE54566",
#'   annotation = annotation,
#'   arrayType = "AFFY",
#'   method.probe = "median",
#'   dir = "./"
#' )
preprocess_rawdata <- function(GEO_name = GEO_name,
                               GPL_name = GPL,
                               arrayType = Type,
                               method.probe = "median",
                               dir = dir) {
  if (requireNamespace("GEOquery", quietly = TRUE)) {
    #######################################################
    #               1. download raw data
    #######################################################
    des_dir <- paste0(dir, GEO_name)
    GEOquery::getGEOSuppFiles(GEO_name, filter_regex = "RAW", baseDir = paste0(des_dir, "/1.rawdata"))
    #######################################################
    #               2. preprocess raw data
    #######################################################
    if (arrayType == "AFFY") {
      exp<-preprocess_rawdata.affy(GEO_name = GEO_name,
                                   GPL_name = GPL,
                                   method.probe = "median",
                                   dir = dir)
    }
    if (arrayType == "AGI4-44K") {
      exp<-preprocess_rawdata.agi4_44k(GEO_name = GEO_name,
                                                 GPL_name = GPL,
                                                 method.probe = "median",
                                                 dir = dir)
    }
    if (arrayType == "ILMN"){
      exp<-preprocess_rawdata.illumina(GEO_name = GEO_name,
                                                 GPL_name = GPL,
                                                 method.probe = "median",
                                                 dir = dir)
    }
  } else {
    stop(
      "Package \"GEOquery\" must be installed to use this function.",
      call. = FALSE
    )
  }
}
