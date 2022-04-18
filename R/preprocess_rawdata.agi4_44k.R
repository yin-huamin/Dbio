#' Process Agilent 4-44K raw data from GEO
#'
#' @param GEO_name GSE number
#' @param GPL_name Probe anotation platform
#' @param method.probe one of "medain", "mean", "max", "min", default is "median"
#' @param dir Des dir
#'
#' @return
#' @export
#'
#' @examples
preprocess_rawdata.agi4_44k<-function(GEO_name = GEO_name,
                                  GPL_name = GPL,
                                  method.probe = "median",
                                  dir = dir){
  if(!requireNamespace("limma")){
    BiocManager::install("limma")
  }
  library(limma)
  # reference : https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf
  # 2. find rawdata dir
  tar_dir <- paste0(des_dir, "/1.rawdata/", GEO_name, "/")
  if(length(tarfile)==1){
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
  }else{
    message("Multi platfrom study, will only preprocess selected GPL...")
    tarfile <- tarfile[stringr::str_detect(tarfile,GPL)]
    untar(tarfile = tarfile, exdir = tar_dir)
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
}
