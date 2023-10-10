#' gene_gender
#'
#' @description Create a box plot for a single gene with symbols indicating statistical significance grouped by gender in a specific type of cancer.
#' @param cancer cancer name likes "BLCA".
#' @param gene gene name likes "KLF7".
#' @param add character vector for adding another plot element. likes "none", "dotplot", "jitter".
#' @param palette the color palette to be used for coloring or filling by groups. Allowed values include scientific journal palettes from ggsci R package, e.g.: "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty".
#' @param legend legend position. Allowed values include "top","bottom","left","right" and "none".
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom stringr str_sub
#' @importFrom tibble add_column
#' @importFrom tibble remove_rownames
#' @importFrom tibble column_to_rownames
#' @import  ggpubr
#' @export
#' @examples gene_gender("BLCA","KLF7")
gene_gender=function(cancer,gene,add = "jitter",palette="jco",legend="none"){
  df=dplyr::filter(tpm,Group=="Tumor" & Cancer==cancer)%>%
    dplyr::select(all_of(gene))%>%
    tibble::add_column(ID = stringr::str_sub(rownames(.),1,12),
                       .before=gene) %>%
    dplyr::filter(!duplicated(ID)) %>%
    tibble::remove_rownames(.) %>%
    tibble::column_to_rownames("ID")%>%
    dplyr::filter(rownames(.) %in% rownames(subset(meta,Cancer==cancer)))

  cl=meta[rownames(df),]

  identical(rownames(df),rownames(cl))
  df=cbind(cl,df)
  df=df[,c(5,7)]
  df=na.omit(df)

  p <- ggpubr::ggboxplot(df, x = "gender", y = gene,
                 color = "gender", add=add, palette = palette,
                 xlab=cancer,ylab=gene)+
    theme(legend.title=element_text(size=14),
          legend.text=element_text(size=12),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size=14),
          axis.text.y = element_text(size = 12))+
    border("black")+
    theme(legend.position = legend)
  p + stat_compare_means(label = "p.signif",
                         label.x.npc='center',
                         label.y.npc = 'top',
                         hide.ns = F)
}
