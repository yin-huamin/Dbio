#' Process matrixs from GEO
#'
#' @param dataset Expression matrixs
#' @param metadata sample metadata
#' @param annotation probe annotaion file
#' @param method.probe one of "medain","mean","max","min", default is "median"
#'
#' @return preprocessed expression matrixs
#' @export
#'
#' @examples
#' library(Dbio)
#' # make file structure
#' initial_file("GSE54566", "./")
#' # get sample data
#' sample <- get_sample("GSE54566", "./")
#' # process sample
#' sample_post <- preprocess_sample(sample)
#' # get expression data
#' expmat <- get_matrixs("GSE54566", "./")
#' # get probe annotation data
#' annotation <- get_probe2symbol("GPL570", "GSE54566", "./")
#' preprocess_matrixs(
#'   dataset = expmat,
#'   metadata = sample_post,
#'   annotation = annotation
#' )
preprocess_matrixs <- function(dataset = expmat,
                               metadata = sample_post,
                               annotation = probe2gene,
                               method.probe = "median") {
  if (is.null(metadata)) {
    expmat_cln <- data.frame(dataset) %>%
      tibble::rownames_to_column("ID") %>%
      dplyr::inner_join(annotation, by = "ID") %>% # merge with platform
      Dbio::probe2multisymbol(data = ., method.probe)
  } else {
    sid <- base::intersect(colnames(dataset), metadata$GSM)

    expmat_cln <- data.frame(dataset) %>%
      dplyr::select(all_of(sid)) %>% # same GSM ID with metadata
      tibble::rownames_to_column("ID") %>%
      dplyr::inner_join(annotation, by = "ID") %>% # merge with platform
      Dbio::probe2multisymbol(data = ., method.probe)
  }
  return(expmat_cln)
}
