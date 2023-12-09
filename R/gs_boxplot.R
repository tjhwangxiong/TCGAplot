#' gs_boxplot
#'
#' @description Create a tumor-normal box plot for a gene set in a specific type of cancer.
#'
#' @param cancer cancer name likes "BRCA".
#' @param geneset gene set name likes "HALLMARK_DNA_REPAIR".
#' @param add character vector for adding another plot element likes "none", "dotplot", "jitter".
#' @param palette the color palette to be used for coloring or filling by groups. Allowed values include scientific journal palettes from ggsci R package, e.g.: "npg", "aaas", "lancet", "jco".
#' @param legend legend position. Allowed values include "top","bottom","left","right" and "none".
#' @param label character string specifying label type. Allowed values include "p.signif" (shows the significance levels), "p.format" (shows the formatted p value). label="p.signif" is default.
#' @param method a character string indicating which method to be used for comparing means including "wilcox.text" and "t.test". method="wilcox.test" is default.
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom GSVA gsva
#' @import  ggpubr
#' @export
#'
#' @examples gs_boxplot("BRCA","HALLMARK_DNA_REPAIR")
gs_boxplot=function(cancer,geneset,add = "jitter",palette="jco",legend="none",label="p.signif",method="wilcox.test"){
  genelist=msigb_genesets[geneset]

  df=filter(tpm,Cancer==cancer)%>%
     select(-c(Cancer,Group))%>%
     t()

  #gsva
  df <- GSVA::gsva(expr=df,
                   gset.idx.list=genelist,
                   kcdf="Gaussian",
                   verbose=T,
                   parallel.sz = parallel::detectCores())
  df=t(df)

  df1=filter(tpm,Cancer==cancer)%>%
      select(Cancer,Group)

  identical(rownames(df),rownames(df1))

  dat=cbind(df1,df)

  p <- ggboxplot(dat, x = "Group", y = geneset,
                 color = "Group", add=add, palette = palette,
                 xlab=cancer,ylab=geneset)+
       theme(legend.title=element_text(size=14),
          legend.text=element_text(size=12),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size=14),
          axis.text.y = element_text(size = 12))+
    border("black")+
    theme(legend.position = legend)
  p + stat_compare_means(label = label,
                         method = method,
                         label.x.npc='center',
                         label.y.npc = 'top',
                         hide.ns = T)
}
