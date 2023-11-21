#' gene_immustimulator_heatmap
#'
#' @description Create a pan-cancer heatmap with symbols indicating statistical significance to reveal the correlation between the expression of a single gene and immune stimulators.
#' @param gene gene names likes "KLF7".
#' @param method method="pearson" is the default value. The alternatives to be passed to correlation are "spearman" and "kendall".
#' @param lowcol colour for low point
#' @param highcol colour for high point
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom psych corr.test
#' @importFrom reshape2 melt
#' @export
#' @examples gene_immustimulator_heatmap("KLF7")
gene_immustimulator_heatmap=function(gene,method = "pearson",lowcol="blue",highcol="red"){
  data("tpm")
  clist=list()
  plist=list()
  immustimulator=c("CD27","CD276","CD28","CD40","CD40LG","CD48","CD70","CD80","CD86","CXCL12","CXCR4","ENTPD1","HHLA2","ICOS","ICOSLG","IL2RA","IL6","IL6R","KLRC1","KLRK1","LTA","MICB","NT5E","PVR","RAET1E","TMIGD2","TNFRSF13B","TNFRSF13C","TNFRSF14","TNFRSF17","TNFRSF18","TNFRSF25","TNFRSF4","TNFRSF8","TNFRSF9","TNFSF13","TNFSF13B","TNFSF14","TNFSF15","TNFSF18","TNFSF4","TNFSF9","ULBP1")
  immustimulator=immustimulator[immustimulator %in% colnames(tpm)]
  for (cancer in cancers) {
    tg=subset(tpm,Group=="Tumor" & Cancer==cancer)%>%
      dplyr::select(all_of(gene))

    ic=subset(tpm,Group=="Tumor" & Cancer==cancer)%>%
      dplyr::select(all_of(immustimulator))
    identical(rownames(tg),rownames(ic))

    # Correlation
    cor <-psych::corr.test(tg, ic, method = method,adjust="none")
    # Cor, p value
    cmt <-cor$r
    pmt <- cor$p
    clist[[cancer]]=as.data.frame(cmt)
    plist[[cancer]]=as.data.frame(pmt)
  }
  clist=do.call(rbind,clist)
  plist=do.call(rbind,plist)
  clist[is.na(clist)]=0
  plist[is.na(plist)]=0.06
  plist=as.data.frame(ifelse(plist<0.01,'**',ifelse(plist<0.05,'*','')))

  clist$cancer=rownames(clist)
  plist$cancer=rownames(plist)
  clist=reshape2::melt(clist,id.vars = "cancer")
  plist=reshape2::melt(plist,id.vars = "cancer")

  df=cbind(clist,plist)
  df=df[,c(1:3,6)]
  colnames(df)=c("cancer","cell","cor","sig")
  df$cell=factor(df$cell,levels=rev(unique(df$cell)))

  ggplot(df, aes(cancer, cell, fill = cor)) +
    geom_tile(aes(width = 1, height = 1), size = 10) +
    geom_text(aes(cancer, cell, label = sig), color = "black") +
    scale_fill_gradient2(low = lowcol, mid = "white", high = highcol) +
    theme(panel.grid = element_blank(),
          panel.background = element_blank())+
    theme(axis.text.x = element_text(angle = 90,colour = "black",size=12),
          axis.text.y = element_text(colour = "black",size=12))+
    labs(x='',y='')+
    theme(axis.ticks = element_blank())
}
