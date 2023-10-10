#' get_meta
#'
#' @description Extract the clinical information of a specific type of cancer in TCGA.
#' @param cancer cancer name likes "COAD".
#' @export
#' @examples get_meta("COAD")
get_meta=function (cancer)
{
  cl = subset(meta, Cancer == cancer)
  return(cl)
}

