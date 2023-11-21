#' pan_paired_boxplot
#'
#' @description Create a pan-cancer paired box plot for a single gene with symbols indicating statistical significance.
#' @param gene gene name likes "KLF7"
#' @param palette the color palette to be used for filling by groups. Allowed values include scientific journal palettes from ggsci R package, e.g.: "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty".
#' @param legend legend position. Allowed values include "top","bottom","left","right" and "none".
#' @param method a character string indicating which method to be used for comparing means including "wilcox.text" and "t.test". method="wilcox.test" is default.
#' @import ggpubr
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @export
#' @examples pan_paired_boxplot("KLF7")
pan_paired_boxplot=function(gene,palette="jco",legend="right",method="wilcox.test"){
  data("paired_tpm")
  # select cancer with significant p value
  # pcancers=unique(paired_tpm$Cancer)
  pcancers=c("BLCA","BRCA","COAD","ESCA","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","PRAD","STAD","THCA","UCEC")
  # dat=list()
  # for (cancer in pcancers) {
  #   data=subset(paired_tpm,Cancer==cancer)%>%
  #     dplyr::select(Cancer,Group,ID,all_of(gene))
  #   fml = as.formula(paste0(gene,"~Group"))
  #   dat[[cancer]]=as.data.frame(ggpubr::compare_means(fml,data,paired=T))
  # }
  # dat=as.data.frame(do.call(rbind,dat))
  # selected=rownames(dat)[dat$p.signif!="ns"]

  # plot
  df =dplyr::select(paired_tpm,Cancer,Group,ID,all_of(gene))

  p=ggpaired(df,x ="Group",y =gene,id = "ID",
             color = "Group",palette = palette,add="jitter",
             xlab='',ylab=gene,
             line.color = "gray",line.size = 0.4,
             facet.by = "Cancer",nrow=1)+
    theme_classic()+
    rotate_x_text(angle = 90)+
    theme(legend.title=element_text(size=14),
          legend.text=element_text(size=12),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12,colour = "black"),
          axis.title.y = element_text(size=14),
          axis.text.y = element_text(size = 12,colour = "black"),
          strip.text.x = element_text(size = 12, color = "black"))+
    theme(legend.position = legend)

  p + stat_compare_means(label="p.signif",
                         method = method,
                         paired=T,
                         label.x.npc=0.4,
                         label.y.npc ='top')
}
