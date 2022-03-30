#' Output clinical data from GEO to des dir
#'
#' @param sample Sample metadata get from get_sample function
#' @param sample_post Preprocessed sample metadata from preprocess_sample
#' @param GEO_name GSE Number
#' @param dir Des dir
#'
#' @return
#' @export
#'
#' @examples
#' library(Dbio)
#' sample <- get_sample("GSE54566", "./")
#' sample_post <- preprocess_sample(sample)
#' output_sample(sample, sample_post, "GSE54566", "./")
output_sample <- function(sample, sample_post, GEO_name = GEO_name, dir = dir) {
  des_dir <- paste0(dir, GEO_name)
  sample.origin <- paste0(paste0(des_dir, "/2.output/"), GEO_name, "_sample.txt")
  sample.process <- paste0(paste0(des_dir, "/2.output/"), GEO_name, "_sample-tiny.txt")
  write.table(sample, sample.origin, col.names = T, row.names = F, sep = "\t", quote = F)
  write.table(sample_post, sample.process, col.names = T, row.names = F, sep = "\t", quote = F)
}
