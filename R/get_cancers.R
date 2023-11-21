#' get_cancers
#'
#' @description Return the sample summary of 33 types of cancer in TCGA.
#'
#' @export
#' @examples get_cancers()
get_cancers=function(){
  table(tpm$Cancer,tpm$Group)
}

