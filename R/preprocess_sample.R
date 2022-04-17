#' Preproces clinical data from GEO
#'
#' @param metadata Sample metadata get from GEO
#'
#' @return
#' @export
#'
#' @examples
#' library(Dbio)
#' library(GEOquery)
#' sample <- get_sample("GSE54566", "./")
#' sample_post <- preprocess_sample(sample)
preprocess_sample <- function(metadata = sample) {
  dat <- metadata[, c(1, grep("_ch1$", colnames(metadata)))]

  colnames(dat) <- gsub("ch1$", "", colnames(dat))

  res <- dat %>%
    rownames_to_column("GSM") %>%
    dplyr::distinct(GSM, .keep_all = TRUE)

  return(res)
}
