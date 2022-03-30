#' Pipeline of processing GEO data
#'
#' @param GEO_name GSE number
#' @param GPL_number Probe annotation platform
#' @param dir Des dir
#' @param quiet logical vector, default is FALSE
#' @param rawdata logical vector, default is FALSE
#' @param Type Type one of "AFFY", "AGI4-44K" when \code{rawdata=TRUE}
#' @param multiprobe multiprobe one of "median", "mean", "max", "min", default is median
#' @param removeRawdata logical vector, default is TRUE
#'
#' @return
#' @export
#'
#' @examples
#' library(Dbio)
#' # Download matrixs and preprocess
#' deal_geo("GSE54566", "GPL570", "./")
#' # Download, preprocess matrixs and keep rawdata
#' deal_geo("GSE54566", "GPL570", "./", removeRawdata = FALSE)
#' # Download, preprocess rawdata of affymetrixs
#' deal_geo("GSE54566", "GPL570", "./", rawdata = TRUE, Type = "AFFY")
deal_geo <- function(GEO_name = GEO_name,
                     GPL_number = GPL_number,
                     dir = dir,
                     quiet = TRUE,
                     rawdata = FALSE,
                     Type = "AFFY",
                     multiprobe = "median",
                     removeRawdata = TRUE) {
  # 1. --- make dir
  Dbio::initial_file(dir, GEO_name)
  if (!quiet) {
    message("1. Make File Structure")
  }


  # 2. --- clinical and clinical process
  sample_post <- Dbio::deal_sample(GEO_name, dir)

  # 3 ------ GPL processing
  gpl <- Dbio::get_probe2symbol(GPL_number, GEO_name, dir)

  if (!quiet) {
    message("3. Get GPL")
  }

  # 4 ------ ExpMat processing
  if (rawdata) {
    exp <- Dbio::preprocess_rawdata(
      GEO_name = GEO_name,
      annotation = gpl,
      arrayType = Type,
      method.probe = multiprobe,
      dir = dir
    )
  } else {
    # 4.1 get expmat
    expmat <- Dbio::get_matrixs(GEO_name, dir)
    # 4.2 detect log2 condition
    expmat <- Dbio::detect_log2(expmat)
    # 4.3 gene expression profile
    exp <- Dbio::preprocess_matrixs(
      dataset = expmat,
      metadata = sample_post,
      annotation = gpl,
      method.probe = multiprobe
    )
  }

  if (!quiet) {
    message("4. Preprocess expressionMat")
  }

  des_dir <- paste0(dir, GEO_name)
  exp.process <- paste0(paste0(des_dir, "/2.output/"),
                        GEO_name,
                        "_",
                        nrow(exp),
                        ".txt")
  write.table(exp,
              exp.process,
              col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

  # 5 ------ Remove rawdata
  if (removeRawdata) {
    files <- paste0(paste0(des_dir, "/1.rawdata/"), list.files(paste0(des_dir, "/1.rawdata/"), recursive = T))

    for (i in seq_len(length(files))) {
      file.remove(files, showWarnings = FALSE)
    }
  } else {
    message(paste0("rawdata have kept in", paste0(des_dir, "/1.rawdata/")))
  }

  return(message("Finished"))
}
