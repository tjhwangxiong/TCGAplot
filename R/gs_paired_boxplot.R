#' gs_paired_boxplot
#' @description Create a paired tumor-normal box plot for a gene set in a specific type of cancer. Only cancers with more than 20 paired samples could be analyzed, including "BLCA","BRCA","COAD","ESCA","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","PRAD","STAD","THCA", and "UCEC".
#' @param cancer name of cancer likes "BRCA"
#' @param geneset gene set name likes "HALLMARK_DNA_REPAIR".
#' @param palette the color palette to be used for coloring or filling by groups. Allowed values include scientific journal palettes from ggsci R package, e.g.: "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty".
#' @param legend legend position. Allowed values include "top","bottom","left","right" and "none".
#' @param method a character string indicating which method to be used for comparing means including "wilcox.text" and "t.test". method="wilcox.test" is default.
#' @import ggpubr
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom GSVA gsva
#' @importFrom dplyr filter
#' @export
#'
#' @examples gs_paired_boxplot("BRCA","HALLMARK_DNA_REPAIR")
gs_paired_boxplot=function(cancer,geneset,palette="jco",legend="none",method="wilcox.test"){
  pcancers=c("BLCA","BRCA","COAD","ESCA","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","PRAD","STAD","THCA","UCEC")

  if (cancer %in% pcancers) {
    genelist=msigb_genesets[geneset]

    df=filter(paired_tpm,Cancer==cancer)%>%
      select(-c(Cancer,Group,ID))%>%
      t()

    #gsva
    df <- GSVA::gsva(expr=df,
                     gset.idx.list=genelist,
                     kcdf="Gaussian",
                     verbose=T,
                     parallel.sz = parallel::detectCores())
    df=t(df)

    df1=filter(paired_tpm,Cancer==cancer)%>%
      select(Cancer,Group,ID)

    identical(rownames(df),rownames(df1))

    dat=cbind(df1,df)

    p=ggpaired(dat,x="Group", y=geneset,
               color = "Group",palette = palette,add="jitter",
               xlab=cancer,ylab=geneset,
               line.color = "gray",line.size = 0.4)+
      theme_classic()+
      theme(legend.title=element_text(size=14),
            legend.text=element_text(size=12),
            axis.title.x = element_text(size = 14),
            axis.text.x = element_text(size = 12,colour = "black"),
            axis.title.y = element_text(size=14),
            axis.text.y = element_text(size = 12,colour = "black"),
            strip.text.x = element_text(size = 12, color = "black", face = "bold"))+
      theme(legend.position = legend)+
      border("black")

    p + stat_compare_means(label="p.format",method = method,paired=T,label.x.npc='center',label.y.npc ='top')
  }
  else
    print(paste0("No sufficient paired samples (n>20) were found in ", cancer,"."))
}
