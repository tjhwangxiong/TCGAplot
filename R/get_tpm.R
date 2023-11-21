#' get_tpm
#'
#' @description Extract the TPM matrix of a specific type of cancer in TCGA.
#' @param cancer cancer name likes "COAD".
#' @export
#' @examples get_tpm("COAD")
get_tpm=function (cancer)
{
  data("tpm")
  exp = subset(tpm, Cancer == cancer)
  return(exp)
}

