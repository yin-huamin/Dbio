#' Deal with Probe to Multi Symbol
#'
#' Calculate the unique symbol (expression) matrixs
#' The gene symbol column should name : "gene"
#' @param data expression matrixs
#' @param method.probe one of "medain","mean","max","min", default is "median"
#'
#' @return
#' @export
#'
#' @examples
#' library(Dbio)
#' # get median value of same gene symbol in matrix
#' exp1 <- probe2multisymbol(exp)
#' # get max  value of same gene symbol in matrix
#' exp2 <- probe2multisymbol(exp, "max")
probe2multisymbol <- function(data = exp, method.probe = "median") {
  # select the method, when occuring multi symbol
  method.probe.switch <- dplyr::case_when(
    method.probe == "median" ~ 1,
    method.probe == "mean" ~ 2,
    method.probe == "max" ~ 3,
    method.probe == "min" ~ 4
  )
  switch(method.probe.switch,
         data %>%
           dplyr::group_by(gene) %>%
           dplyr::summarise_if(is.numeric, median, na.rm = TRUE),
         data %>%
           dplyr::group_by(gene) %>%
           dplyr::summarise_if(is.numeric, mean, na.rm = TRUE),
         data %>%
           dplyr::group_by(gene) %>%
           dplyr::summarise_if(is.numeric, max, na.rm = TRUE),
         data %>%
           dplyr::group_by(gene) %>%
           dplyr::summarise_if(is.numeric, min, na.rm = TRUE)
  )
}
