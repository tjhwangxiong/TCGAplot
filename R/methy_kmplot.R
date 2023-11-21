#' methy_kmplot
#'
#' @description Describes the K_M survival plot based on the promoter methylation of a single gene in a specific type of cancer. A pdf file will be generated in the working directory.
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
#' @examples methy_kmplot("COAD","KLF7")

methy_kmplot=function(cancer,gene,palette='jco'){

  if (gene %in% promoter_probe$gene) {
    print('A pdf file will be generated in the working directory.')

    cg_list=promoter_probe$probe[promoter_probe$gene==gene]
    cg_list = cg_list[cg_list %in% colnames(pro_methy)]

    exprSet=dplyr::select(pro_methy,all_of(cg_list))%>%
      tibble::add_column(ID = stringr::str_sub(rownames(.),1,12)) %>%
      dplyr::filter(!duplicated(ID)) %>%
      tibble::remove_rownames(.) %>%
      tibble::column_to_rownames("ID")%>%
      dplyr::filter(rownames(.) %in% rownames(subset(meta,Cancer==cancer)))%>%
      t()

    cl=meta[colnames(exprSet),]

    # plot
    pdf(paste0(gene,'_methylation_kmplot_',cancer,'.pdf'),width = 8,height = 6)
    for (cg in cg_list){
      cl1=cl%>%
        dplyr::mutate(methylation=ifelse(exprSet[cg,]> median(exprSet[cg,]),'high','low'))

      sfit = survival::survfit(survival::Surv(time, event) ~ methylation, data = cl1)
      p=survminer::ggsurvplot(sfit, pval = TRUE, palette = palette,
                              data = cl1, legend = c(0.8, 0.8),
                              title =paste0('KMplot of ',gene,' ',cg, ' in ',cancer),risk.table = T)
      print(p)
    }
    dev.off()
  } else {print('This gene is not included in the methylation matrix!')}
}
