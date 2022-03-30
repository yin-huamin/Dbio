#' Confirm expmat log2(+1) transfered by detection
#'
#' @param data expression matrixs
#'
#' @return log2(+1) transformed matrixs
#' @export
#'
#' @examples
#' library(Dbio)
#' expmat <- get_matrixs("GSE54566", "./")
#' expmat2 <- detect_log2(expmat)
detect_log2 <- function(data = expmat) {
  qx <- as.numeric(quantile(data, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
  LogC <- (qx[5] > 100) ||
    (qx[6] - qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

  if (LogC) {
    data[which(data <= 0)] <- NA
    data <- log2(data + 1)
    print("log2 transform finished")
  } else {
    print("log2 transform not needed")
  }
  return(data)
}
