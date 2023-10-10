#' gene_immunescore_heatmap
#'
#' @description Create a pan-cancer heatmap with symbols indicating statistical significance to reveal the correlation between the expression of a single gene and immune scores, including Stromal score, immune score, and ESTIMATE score.
#' @param gene gene names likes "KLF7".
#' @param method method="pearson" is the default value. The alternatives to be passed to correlation are "spearman" and "kendall".
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom psych corr.test
#' @importFrom reshape2 melt
#' @export
#' @examples gene_immunescore_heatmap("KLF7")
#'
gene_immunescore_heatmap=function(gene,method="pearson"){
  clist=list()
  plist=list()
  for (cancer in cancers) {
    tg=subset(tpm,Group=="Tumor" & Cancer==cancer)%>%
      dplyr::select(all_of(gene))%>%
      dplyr::filter(rownames(.) %in% rownames(immuscore))

    ic=immuscore[rownames(tg),]
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
  colnames(df)=c("cancer","score","cor","sig")
  df$score=factor(df$score,levels=rev(unique(df$score)))

  ggplot(df, aes(cancer, score, fill = cor)) +
    geom_tile(aes(width = 1, height = 1), linewidth = 10) +
    geom_text(aes(cancer, score, label = sig), color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    theme(panel.grid = element_blank(),
          panel.background = element_blank())+
    theme(axis.text.x = element_text(angle = 90,colour = "black",size=12),
          axis.text.y = element_text(colour = "black",size=12))+
    labs(x='',y='')+
    theme(axis.ticks = element_blank())+coord_equal()
}
