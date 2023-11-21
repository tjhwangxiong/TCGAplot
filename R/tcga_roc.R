#' tcga_roc
#' @description Diagnostic ROC curve of a single gene in a specific type of cancer.
#'
#' @param cancer cancer name likes "BRCA".
#' @param gene gene name likes "KLF7".
#' @import pROC
#' @import ggplot2
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @export
#' @examples tcga_roc("BRCA","KLF7")
tcga_roc=function(cancer,gene){

  df=dplyr::filter(tpm,Cancer==cancer)%>%
    dplyr::select(Group,all_of(gene))

  res<- pROC::roc_(df,"Group",gene,aur=TRUE,
            ci=TRUE, # 95%CI
            #percent=TRUE, # percentile
            smooth=TRUE,
            levels=c('Normal','Tumor'))

  p<- pROC::ggroc(res, color ="red",legacy.axes = TRUE)+
      geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype=4)+
      theme_bw() +
      ggtitle(paste0(cancer,' ROC Curve'))+
      theme(plot.title = element_text(hjust = 0.5,size = 14),
          legend.title=element_text(size=14,colour = "black"),
          legend.text=element_text(size=12,colour = "black"),
          axis.title.x = element_text(size = 14,colour = "black"),
          axis.text.x = element_text(size = 12,colour = "black"),
          axis.title.y = element_text(size=14,colour = "black"),
          axis.text.y = element_text(size = 12,colour = "black"))

 p + ggplot2::annotate("text",x=0.75,y=0.25,label=paste("AUC = ", round(res$auc,3)))+
     ggplot2::annotate("text",x=0.75,y=0.20,label=paste("95%CI: ", round(res$ci[1],3),'-',round(res$ci[3],3)))
}
