#' gs_receptor_heatmap
#'
#' @description Create a pan-cancer heatmap with symbols indicating statistical significance to reveal the correlation between the expression of a gene set and chemokine receptors.
#' @param geneset gene set name likes "KEGG_APOPTOSIS".
#' @param method method="pearson" is the default value. The alternatives to be passed to correlation are "spearman" and "kendall".
#' @param lowcol colour for low point
#' @param highcol colour for high point
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom GSVA gsva
#' @importFrom psych corr.test
#' @importFrom reshape2 melt
#' @import  ggplot2
#' @export
#' @examples gs_receptor_heatmap("KEGG_APOPTOSIS")
gs_receptor_heatmap=function(geneset,method = "pearson",lowcol="blue",highcol="red"){
  data("tpm")
  clist=list()
  plist=list()
  receptor=c("CCR1","CCR2","CCR3","CCR4","CCR5","CCR6","CCR7","CCR8","CCR9","CCR10",
             "CXCR1","CXCR2","CXCR3","CXCR4","CXCR5","CXCR6","XCR1","CX3R1")
  receptor=receptor[receptor %in% colnames(tpm)]

  genelist=msigb_genesets[geneset]

  for (cancer in cancers) {
    tg=subset(tpm,Group=="Tumor" & Cancer==cancer)%>%
      dplyr::select(-c(Cancer,Group))%>%
      t()
    #gsva
    tg <- GSVA::gsva(expr=tg,
                     gset.idx.list=genelist,
                     kcdf="Gaussian",
                     verbose=T,
                     parallel.sz = parallel::detectCores()) # use all cores
    tg=t(tg)

    ic=subset(tpm,Group=="Tumor" & Cancer==cancer)%>%
      dplyr::select(all_of(receptor))
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
