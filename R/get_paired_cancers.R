#' get_paired_cancers
#'
#' @description Return the sample summary of 15 types of cancer containing more than 20 paired samples in TCGA.
#' @export
#' @examples get_paired_cancers()
get_paired_cancers=function(){
  table(paired_tpm$Cancer,paired_tpm$Group)
}

