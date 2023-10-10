#' get_methy
#'
#' @description Extract the promoter methylation information of all samples in TCGA.
#' @export
#' @examples get_methy()
get_methy=function ()
{
  data("cg_probe","methy")
  exp=list()
  exp$probe = cg_probe
  exp$methy=methy
  return(exp)
}

