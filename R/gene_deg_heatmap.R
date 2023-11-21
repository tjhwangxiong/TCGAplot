#' gene_deg_heatmap
#'
#' @description Create a heatmap for differentially expressed genes grouped by the expression of a single gene in a specific type of cancer.
#'
#' @param cancer cancer name likes "BLCA".
#' @param gene gene name likes "KLF7".
#' @param top_n the number of top DEGS to be shown in the heatmap.
#' @import limma
#' @import edgeR
#' @import dplyr
#' @importFrom tinyarray draw_heatmap
#' @export
#' @examples gene_deg_heatmap("BLCA","KLF7")
gene_deg_heatmap=function (cancer, gene,top_n=20){
  exp = subset(tpm, Group == "Tumor" & Cancer == cancer)
  exp = exp[, -c(1:2)]
  exp = as.matrix(t(exp))
  Group = factor(ifelse(exp[gene, ] > median(exp[gene, ]),
                        "high", "low"), levels = c("low", "high"))

  dat = normalizeBetweenArrays(exp)
  design = model.matrix(~Group)
  fit = lmFit(dat, design)
  fit = eBayes(fit)
  options(digits = 4)
  DEG = topTable(fit, coef = 2, adjust = "BH", n = Inf)
  DEG = na.omit(DEG)
  write.csv(DEG,file="DEG.csv")
  DEG = subset(DEG, P.Value<0.05)
  DEG = DEG[order(DEG$logFC,decreasing = T),]
  DEG = DEG[c(1:top_n,(nrow(DEG)-top_n+1):nrow(DEG)),]

  df = exp[rownames(DEG),order(Group)]
  group_list=factor(c(rep("low",sum(Group=="low")),rep("high",sum(Group=="high"))), levels = c("low", "high"))

  draw_heatmap(df,group_list,n_cutoff = 2, legend = T, show_rownames = T, annotation_legend = T,
               cluster_cols=F,cluster_rows=F)
}

