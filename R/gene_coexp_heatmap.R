#' gene_coexp_heatmap
#'
#' @description Heatmap and Go enrichment of the positive and negative co-expressed genes of a single gene in a specific type of cancer.
#' @param cancer cancer name likes "STAD".
#' @param gene gene name likes "KLF7".
#' @param top_n the number of co-expressed genes.
#' @param method method="pearson" is the default value. The alternatives to be passed to correlation were "spearman" and "kendall".
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom psych corr.test
#' @importFrom reshape2 melt
#' @importFrom enrichplot dotplot
#' @importFrom tinyarray ggheat
#' @import  ggplot2
#' @import clusterProfiler
#' @import org.Hs.eg.db
#' @import patchwork
#' @export
#' @examples gene_coexp_heatmap("STAD","KLF7")
gene_coexp_heatmap=function(cancer,gene,top_n=20, method="pearson"){

    tg=subset(tpm,Group=="Tumor" & Cancer==cancer)%>%
      dplyr::select(all_of(gene))

    ic=subset(tpm,Group=="Tumor" & Cancer==cancer)%>%
      dplyr::select(-c(1,2,all_of(gene)))
    identical(rownames(tg),rownames(ic))

    # Correlation
    cor <-psych::corr.test(tg, ic, method = method,adjust="none")
    # Cor, p value
    cmt <-cor$r
    pmt <- cor$p
    df=cbind(as.data.frame(t(cmt)),as.data.frame(t(pmt)))
    df=na.omit(df)
    colnames(df)=c("cor","pvalue")
    df=df[df$pvalue<0.05,]
    df=df[order(df$cor,decreasing = T),]
    gene_pos=rownames(df[1:top_n,])
    gene_neg=rownames(df[(nrow(df)-top_n+1):nrow(df),])

    exp = subset(tpm, Group == "Tumor" & Cancer == cancer)%>%
          dplyr::select(all_of(gene_pos),all_of(gene_neg),all_of(gene))

    ## code from tinyarray
    riskscore = exp[,ncol(exp)]
    cut = median(riskscore)
    exp$fp=riskscore

    exp$Risk = ifelse(exp$fp<median(exp$fp),"low","high")
    exp$Risk = factor(exp$Risk,levels = c("low","high"))

    fp_dat=data.frame(patientid=1:length(riskscore),
                      riskscore=as.numeric(sort(riskscore)),
                      group = exp$Risk[order(riskscore)])

    exp_dat = scale(exp[order(exp$fp),-c((ncol(exp)-2):ncol(exp))])


    p1=ggplot(fp_dat,aes(x=patientid,y=riskscore,color = group))+
      geom_point()+
      scale_color_manual(values = c("#2874C5","#f87669"))+
      geom_vline(xintercept = sum(riskscore<cut),lty = 2)+
      scale_x_continuous(expand=c(0,0))+
      ylab(paste0(gene))+
      theme_bw()

    p2 = tinyarray::ggheat(exp_dat,
                    exp$Risk[order(riskscore)],
                    show_rownames = F,
                    color = c("#2874C5","white","#f87669"),
                    legend_color = c("#2874C5","#f87669"))
    ## code from tinyarray

    go_pos <- enrichGO(gene = gene_pos,
                      OrgDb = org.Hs.eg.db,
                      keyType = "SYMBOL",
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 1,
                      qvalueCutoff = 1,
                      readable = T)

    go_neg <- enrichGO(gene = gene_neg,
                      OrgDb = org.Hs.eg.db,
                      keyType = "SYMBOL",
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 1,
                      qvalueCutoff = 1,
                      readable = T)

    p3 = enrichplot::dotplot(go_pos, showCategory = 5,
                             title=paste0("GO of ",gene," positively co-expressed genes"))
    p4 = enrichplot::dotplot(go_neg, showCategory = 5,
                             title=paste0("GO of ",gene," negatively co-expressed genes"))

    p1/p2+ plot_layout(design = "A
                         B
                         B
                         B") | p3/p4
}
