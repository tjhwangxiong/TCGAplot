#' gs_stage
#'
#' @description Create a box plot for a gene set grouped by stage in a specific type of cancer.
#' @param cancer cancer name likes "COAD".
#' @param geneset gene set name likes "HALLMARK_DNA_REPAIR".
#' @param add character vector for adding another plot element. likes "none", "dotplot", "jitter".
#' @param palette the color palette to be used for coloring or filling by groups. Allowed values include scientific journal palettes from ggsci R package, e.g.: "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty".
#' @param legend legend position. Allowed values include "top","bottom","left","right" and "none".
#' @param method a character string indicating which method to be used for comparing means including "wilcox.text" and "t.test". method="wilcox.test" is default.
#' @importFrom magrittr %>%
#' @import ggpubr
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom GSVA gsva
#' @importFrom tibble add_column
#' @importFrom stringr str_sub
#' @importFrom tibble remove_rownames
#' @importFrom tibble column_to_rownames
#' @export
#' @examples gs_stage("COAD","HALLMARK_DNA_REPAIR")
gs_stage=function(cancer,geneset,add = "jitter",palette="jco",legend="none",method="wilcox.test"){

  genelist=msigb_genesets[geneset]

  df=dplyr::filter(tpm,Group=="Tumor" & Cancer==cancer)%>%
    dplyr::select(-c(Cancer,Group))%>%
    tibble::add_column(ID = stringr::str_sub(rownames(.),1,12),
                       .before=1) %>%
    dplyr::filter(!duplicated(ID)) %>%
    tibble::remove_rownames(.) %>%
    tibble::column_to_rownames("ID")%>%
    dplyr::filter(rownames(.) %in% rownames(subset(meta,Cancer==cancer)))%>%
    t()

  #gsva
  df <- GSVA::gsva(expr=df,
                   gset.idx.list=genelist,
                   kcdf="Gaussian",
                   verbose=T,
                   parallel.sz = parallel::detectCores())
  df=t(df)

  cl=meta[rownames(df),]

  identical(rownames(df),rownames(cl))
  df=cbind(cl,df)
  df=df[,c(6:7)]
  df=na.omit(df)
  df$stage=factor(df$stage,levels = c('I','II','III','IV'))

  my_comparisons <- list(c("I","II"), c("I", "III"),c("I", "IV"),
                         c("II", "III"),c("II", "IV"),
                         c("III", "IV"))

  p <- ggboxplot(df, x = "stage", y = geneset,
                 color = "stage", add=add, palette = palette,
                 xlab=cancer,ylab=geneset)+
    theme(legend.title=element_text(size=14),
          legend.text=element_text(size=12),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size=14),
          axis.text.y = element_text(size = 12))+
    border("black")+
    theme(legend.position = legend)
  p + stat_compare_means(comparisons = my_comparisons,
                         label = "p.format",
                         method = method,
                         label.x.npc='center',
                         label.y.npc = 'top',
                         hide.ns = F)
}
