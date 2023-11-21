#' get_promoter_methy
#'
#' @description Extract promoter methylation information of a specific type of tumor.
#' @param cancer cancer name likes "COAD".
#' @export
#' @examples get_promoter_methy("COAD")
get_promoter_methy=function (cancer) {
  data("pro_methy","promoter_probe")
  df=rownames(subset(tpm,Cancer==cancer & Group=="Tumor"))

  exp=list()
  exp$probe = promoter_probe
  exp$methy=pro_methy[rownames(pro_methy) %in% df, ]
  return(exp)
}
