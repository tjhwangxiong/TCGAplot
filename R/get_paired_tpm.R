#' get_paired_tpm
#'
#' @description Extract the TPM matrix of a specific type of cancer with paired samples (n>20) in TCGA.
#' @param cancer cancer name likes "COAD".
#' @export
#' @examples get_paired_tpm("COAD")
get_paired_tpm=function (cancer)
{
  exp = subset(paired_tpm, Cancer == cancer)
  return(exp)
}

