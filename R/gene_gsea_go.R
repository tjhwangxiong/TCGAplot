#' gene_gsea_go
#'
#' @description GSEA-GO analysis of DEGs grouped by the expression of a single gene in a specific type of cancer, and the top 5 GO BP pathways were shown.
#'
#' @param cancer cancer name likes "BLCA".
#' @param gene gene name likes "KLF7".
#' @param logFC_cutoff cutoff value of logFC, 2 was the default setting.
#' @param pvalue_cutoff cutoff value of pvalue, 0.05 was the default setting.
#' @import limma
#' @import edgeR
#' @import org.Hs.eg.db
#' @import clusterProfiler
#' @import dplyr
#' @importFrom enrichplot gseaplot2
#' @export
#' @examples gene_gsea_go("BLCA","KLF7")
gene_gsea_go=function(cancer,gene,logFC_cutoff=2,pvalue_cutoff = 0.05){
  data("tpm")
  exp=subset(tpm,Group=="Tumor" & Cancer==cancer)
  exp=exp[,-c(1:2)]
  exp=as.matrix(t(exp))

  Group = factor(ifelse(exp[gene,]> median(exp[gene,]),'high','low'),levels=c("low","high"))
  logFC_t = logFC_cutoff
  pvalue_t = pvalue_cutoff

  #DEG
  dat=normalizeBetweenArrays(exp)
  design=model.matrix(~Group)
  fit=lmFit(dat,design)
  fit=eBayes(fit)
  options(digits = 4)
  DEG=topTable(fit,coef=2,adjust='BH',n=Inf)
  DEG = na.omit(DEG)
  DEG$symbol=rownames(DEG)

  #GSEA
  s2e <- clusterProfiler::bitr(DEG$symbol,
              fromType = "SYMBOL",
              toType = "ENTREZID",
              OrgDb = org.Hs.eg.db)

  DEG <- dplyr::inner_join(DEG,s2e,by=c("symbol"="SYMBOL"))

  geneList=DEG$logFC
  names(geneList)=DEG$ENTREZID
  geneList=sort(geneList,decreasing = T)

  #GSEA-GO
  Go_gseresult <- clusterProfiler::gseGO(geneList, 'org.Hs.eg.db',
                                          keyType = "ENTREZID", ont="BP", pvalueCutoff=1)

  # plot
  enrichplot::gseaplot2(Go_gseresult, 1:5,
                        title = cancer)

}
