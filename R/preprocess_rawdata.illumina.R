#' Process ILLUMINA raw data from GEO
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
preprocess_rawdata.illumina<-function(GEO_name = GEO_name,
                                      GPL_name = GPL,
                                      method.probe = "median",
                                      dir = dir){

  if(!requireNamespace("limma")){
    BiocManager::install("limma")
  }
  library(limma)

  # 2.1 find rawdata dir
  des_dir <- paste0(dir, GEO_name)
  raw_dir <- paste0(des_dir, "/1.rawdata/")
  raw_file <- paste0(raw_dir,list.files(raw_dir,pattern = "non-normalized.*.txt$"))
  rawdt<-read_delim(raw_file,delim = "\t",col_names = T)

  rawdt.pval<-rawdt %>% dplyr::select(1,ends_with("Pval"))
  pval_detect<-apply(rawdt.pval[,-1], 1, function(x){length(x[x<0.05])/length(x)>0.5})

  rawdt.exp<-rawdt %>% dplyr::select(1,ends_with("Signal"))

  rawdt.exp<-rawdt.exp[pval_detect,] %>% column_to_rownames(var = "PROBE_ID")

  rawdt.exp <- Dbio::detect_log2(rawdt.exp)
  # background correction
  mat<-backgroundCorrect(rawdt.exp,method = "normexp")
  # quantile guiyihua
  mat<-normalizeBetweenArrays(mat,method = "quantile")
  # log2
  mat <- Dbio::detect_log2(mat)

  annotation <- Dbio::get_probe2symbol(GPL_name, GEO_name, dir,Array_type = Type)

  exp <- merge(annotation,
               mat %>%
                 as.data.frame() %>%
                 rownames_to_column("ID")) %>%
    dplyr::mutate(gene = toupper(gene)) %>%
    Dbio::probe2multisymbol(data = ., method.probe)
  return(exp)
}


