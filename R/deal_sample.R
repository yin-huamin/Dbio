#' Download, Preprocess, and Save sample metadata
#'
#' @param GEO_name GSE number
#' @param dir Des dir
#'
#' @return
#' @export
#'
#' @examples
#' library(Dbio)
#' deal_sample(GEO_name = "GSE54566", dir = "./")
deal_sample <- function(GEO_name = GEO_name, dir = dir) {
  sample <- Dbio::get_sample(GEO_name, dir)
  sample_post <- Dbio::preprocess_sample(sample)
  Dbio::output_sample(sample, sample_post, GEO_name, dir)
  return(sample_post)
}
