#' Process AFFYMATRIX raw data from GEO
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
preprocess_rawdata.affy<-function(GEO_name = GEO_name,
                                  GPL_name = GPL,
                                  method.probe = "median",
                                  dir = dir){

  if(!requireNamespace("oligo")){
    BiocManager::install("oligo")
  }
  library(oligo)
  # 2.1 find rawdata dir
  des_dir <- paste0(dir, GEO_name)
  tar_dir <- paste0(des_dir, "/1.rawdata/", GEO_name, "/")

  if(length(tarfile)==1){
    untar(tarfile = tarfile, exdir = tar_dir)
    # 3. read .cel files
    celfiles <- list.files(tar_dir,pattern = "[CEL|cel]$")

    mydata <- oligo::read.celfiles(paste0(tar_dir,celfiles))
    # 4. B bg correct N normalize S summarize
    eset <- oligo::rma(mydata)
    # 5. get expression mat
    mat <- Biobase::exprs(eset)
    # 6. probe 2 symbol
    annotation <- Dbio::get_probe2symbol(GPL_name, GEO_name, dir,Array_type = Type)
    exp <- merge(annotation,
                 mat %>%
                   as.data.frame() %>%
                   rownames_to_column("ID")) %>%
      dplyr::mutate(gene = toupper(gene)) %>%
      Dbio::probe2multisymbol(data = ., method.probe)

  }else{
    message("Multi platfrom study, will only preprocess selected GPL...")
    tarfile <- tarfile[stringr::str_detect(tarfile,GPL)]
    untar(tarfile = tarfile, exdir = tar_dir)
    # 3. read .cel files
    celfiles <- list.files(tar_dir,pattern = "[CEL|cel]$")

    mydata <- oligo::read.celfiles(paste0(tar_dir,celfiles))
    # 4. B bg correct N normalize S summarize
    eset <- oligo::rma(mydata)
    # 5. get expression mat
    mat <- Biobase::exprs(eset)
    # 6. probe 2 symbol
    annotation <- Dbio::get_probe2symbol(GPL_name, GEO_name, dir,Array_type = Type)
    exp <- merge(annotation,
                 mat %>%
                   as.data.frame() %>%
                   rownames_to_column("ID")) %>%
      dplyr::mutate(gene = toupper(gene)) %>%
      Dbio::probe2multisymbol(data = ., method.probe)

  }
  return(exp)
}


