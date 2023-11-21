#' gs_3age
#'
#' @description Create a box plot for a gene set grouped by three age groups in a specific type of cancer.
#' @param cancer cancer name likes "ACC".
#' @param geneset gene set name likes "HALLMARK_DNA_REPAIR".
#' @param age1 numeric number of young age like 40.
#' @param age1 numeric number of old age like 60.
#' @param add character vector for adding another plot element. likes "none", "dotplot", "jitter".
#' @param palette the color palette to be used for coloring or filling by groups. Allowed values include scientific journal palettes from ggsci R package, e.g.: "npg", "aaas", "lancet", "jco".
#' @param legend legend position. Allowed values include "top","bottom","left","right" and "none".
#' @param method a character string indicating which method to be used for comparing means including "wilcox.text" and "t.test". method="wilcox.test" is default.
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom GSVA gsva
#' @importFrom stringr str_sub
#' @importFrom tibble add_column
#' @importFrom tibble remove_rownames
#' @importFrom tibble column_to_rownames
#' @import  ggpubr
#' @export
#' @examples gs_3age("COAD","HALLMARK_DNA_REPAIR")
gs_3age=function(cancer,geneset,age1=40,age2=60,add = "jitter",palette="jco",legend="none",method="wilcox.test"){
  data("tpm")
  genelist=msigb_genesets[geneset]

  df=dplyr::filter(tpm,Group=="Tumor" & Cancer==cancer)%>%
    dplyr::select(-c(Cancer,Group))%>%
    tibble::add_column(ID = stringr::str_sub(rownames(.),1,12),
                       .before=gene) %>%
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
  df$Group=ifelse(df$age<=age1,paste0('≤',age1),
                  ifelse(df$age<=age2,paste0(age1,'~',age2),paste0('>',age2)))

  df$Group=factor(df$Group,levels=c(paste0('≤',age1),paste0(age1,'~',age2),paste0('>',age2)))
  df=df[,7:8]
  df=na.omit(df)

  my_comparisons <- list(c(paste0('≤',age1),paste0(age1,'~',age2)),
                         c(paste0('≤',age1), paste0('>',age2)),
                         c(paste0(age1,'~',age2), paste0('>',age2)))

  p <- ggpubr::ggboxplot(df, x = "Group", y = geneset,
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
  p + stat_compare_means(comparisons = my_comparisons,
                         label = "p.format",
                         method = method,
                         label.x.npc='center',
                         label.y.npc = 'top',
                         hide.ns = F)
  }
