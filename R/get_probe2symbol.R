#' Get probe annation data
#'
#' @param GPL_number Probe annotation platfrom number
#' @param GEO_name GSE number
#' @param dir Des dir
#'
#' @return
#' @export
#'
#' @examples
#' library(Dbio)
#' get_probe2symbol("GPL570", "GSE54566", "./")
get_probe2symbol <- function(GPL_number,
                             GEO_name = GEO_name,
                             dir = dir,
                             Array_type = "None") {

  des_dir <- paste0(dir, GEO_name)
  # gpl_dir <- "./data/"
  gpl.list<-data(package = "Dbio")[[3]][,3]

  if (any(str_detect(gpl.list, GPL_number))) { # if gpl exists

    gpl.name <- stringr::str_extract(
      gpl.list[stringr::str_detect(gpl.list,GPL_number)],
      "GPL.*\\d")

    data(list=gpl.name,package="Dbio",envir = environment())

  } else { # get gpl data from GEO
    if (requireNamespace("GEOquery", quietly = TRUE)) {
      gpl <- GEOquery::getGEO(
        GEO = GPL_number,
        destdir = paste0(des_dir, "/1.rawdata")
      )
    } else {
      stop(
        "Package \"GEOquery\" must be installed to use this function.",
        call. = FALSE
      )
    }
    if (Array_type == "ILMN") {
      # illumina
      gpl <- GEOquery::Table(gpl) %>%
        dplyr::select(any_of(c("ID", "ILMN_Gene", "ILMN_GENE"))) %>%
        dplyr::rename(ILMN_Gene = 2) %>%
        dplyr::mutate(gene = trimws(ILMN_Gene, which = "both")) %>%
        dplyr::select(-ILMN_Gene) %>%
        dplyr::filter(gene != "") %>%
        unique()
    } else if (Array_type == "AFFY") {
      # affy
      gpl <- GEOquery::Table(gpl) %>%
        dplyr::select(any_of(c("ID", "GeneSymbol", "Gene Symbol"))) %>%
        dplyr::rename(GeneSymbol = 2) %>%
        tidyr::separate_rows("GeneSymbol", sep = " /// ") %>%
        dplyr::mutate(gene = trimws(GeneSymbol, which = "both")) %>%
        dplyr::select(-GeneSymbol) %>%
        dplyr::filter(gene != "") %>%
        unique()
    } else if (Array_type == "PHLANX") {
      # Phalanx
      gpl <- GEOquery::Table(gpl) %>%
        dplyr::select(ID, GENE_SYMBOL, Gene_symbol, ORF) %>%
        dplyr::mutate(gene = trimws(GENE_SYMBOL, which = "both")) %>%
        dplyr::select(-GENE_SYMBOL) %>%
        dplyr::filter(gene != "") %>%
        unique()
    } else if (Array_type == "other") {
      message("Gpl hsa not been preprocessed. Welcome to issue the problem.
              https://github.com/yin-huamin/Dbio/issues")
    }
  }
  return(gpl)
}
