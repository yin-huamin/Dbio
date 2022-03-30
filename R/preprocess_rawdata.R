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
                               annotation = probe2gene,
                               arrayType = Type,
                               method.probe = "median",
                               dir = dir) {
  if (requireNamespace("GEOquery", quietly = TRUE)) {
    # 1. download raw data
    des_dir <- paste0(dir, GEO_name)
    GEOquery::getGEOSuppFiles(GEO_name, filter_regex = "RAW", baseDir = paste0(des_dir, "/1.rawdata"))
    if (arrayType == "AFFY") {
      library(affy)
      # 2. find rawdata dir
      tar_dir <- paste0(des_dir, "/1.rawdata/", GEO_name, "/")
      untar(tarfile = paste0(tar_dir, list.files(tar_dir)), exdir = tar_dir)
      # 3. read .cel files
      mydata <- affy::ReadAffy(celfile.path = tar_dir)
      # 4. B bg correct N normalize S summarize
      eset <- affy::rma(mydata)
      # 5. get expression mat
      mat <- Biobase::exprs(eset)
      # 6. probe 2 symbol
      exp <- merge(annotation,
                   mat %>%
                     as.data.frame() %>%
                     rownames_to_column("ID")) %>%
        dplyr::mutate(gene = toupper(gene)) %>%
        Dbio::probe2multisymbol(data = ., method.probe)
      return(exp)
    }
    if (arrayType == "AGI4-44K") {
      library(limma)
      # reference : https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf
      # 2. find rawdata dir
      tar_dir <- paste0(des_dir, "/1.rawdata/", GEO_name, "/")
      untar(tarfile = paste0(tar_dir, list.files(tar_dir)), exdir = tar_dir)
      mydata <- limma::read.maimages(
        files = paste0(tar_dir, list.files(tar_dir, pattern = ".txt")),
        source = "agilent", green.only = TRUE
      ) # only green when 1 channel array
      mydata2 <- limma::backgroundCorrect(mydata, method = "normexp")
      mydata2 <- limma::normalizeBetweenArrays(mydata2, method = "quantile")

      Control <- mydata2$genes$ControlType == 1L
      isExpr <- rowSums(mydata2$E) >= 1
      mat <- cbind(gene = mydata2$genes$GeneName[!Control & isExpr], mydata2$E[!Control & isExpr, ]) %>% as.data.frame()
      mat[, -1] <- apply(mat[, -1], 2, as.numeric)

      exp <- mat %>%
        dplyr::mutate(gene = toupper(gene)) %>%
        dplyr::filter(!stringr::str_detect(gene, "^\\d")) %>%
        dplyr::filter(!stringr::str_detect(gene, "DARKC")) %>%
        Dbio::probe2multisymbol(data = ., method.probe)
      colnames(exp)[-1] <- str_extract(colnames(exp), "GSM\\d+")[-1]
      return(exp)
    }
  } else {
    stop(
      "Package \"GEOquery\" must be installed to use this function.",
      call. = FALSE
    )
  }
}
