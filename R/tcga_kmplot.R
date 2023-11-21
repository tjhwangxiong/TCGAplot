#' tcga_kmplot
#'
#' @description K_M survival plot for a single gene in a specific type of cancer.
#' @param cancer cancer name likes "COAD".
#' @param gene gene name likes "KLF7".
#' @param palette the color palette to be used for coloring or filling by groups. Allowed values include scientific journal palettes from ggsci R package, e.g.: "npg", "aaas", "lancet", "jco".
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom tibble add_column
#' @importFrom tibble remove_rownames
#' @importFrom tibble column_to_rownames
#' @importFrom survival Surv
#' @importFrom survival survfit
#' @importFrom survminer ggsurvplot
#' @export
#' @examples tcga_kmplot("COAD","KLF7")
tcga_kmplot=function(cancer,gene,palette='jco'){
  data("tpm")
  exprSet=subset(tpm,Group=="Tumor" & Cancer==cancer)%>%
    dplyr::select(all_of(gene))%>%
    tibble::add_column(ID = stringr::str_sub(rownames(.),1,12)) %>%
    dplyr::filter(!duplicated(ID)) %>%
    tibble::remove_rownames(.) %>%
    tibble::column_to_rownames("ID")%>%
    dplyr::filter(rownames(.) %in% rownames(subset(meta,Cancer==cancer)))%>%
    t()%>%
    as.matrix()

  cl=meta[colnames(exprSet),]

  # plot
  cl$expression = ifelse(exprSet[gene,]> median(exprSet[gene,]),'high','low')
  sfit = survival::survfit(survival::Surv(time, event) ~expression, data = cl)
  survminer::ggsurvplot(sfit, pval = TRUE, palette = palette,
                        data = cl, legend = c(0.8, 0.8),
                        title =paste0('KMplot of ',gene, ' in ',cancer),risk.table = T)
}
