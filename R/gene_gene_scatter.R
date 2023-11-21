#' gene_gene_scatter
#'
#' @description Scatter plot of gene and gene correlation in a specific type cancer.
#'
#' @param cancer cancer name likes "BLCA".
#' @param gene1 name of gene1 likes "CBX2".
#' @param gene2 name of gene2 likes "CBX3".
#' @param density whether density of gene expression was shown.
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom ggpubr ggscatter
#' @importFrom ggpubr ggscatterhist
#' @export
#' @examples gene_gene_scatter("BLCA","CBX2","CBX3"),
#'           gene_gene_scatter("BLCA","CBX2","CBX3",density="T")
gene_gene_scatter=function(cancer,gene1,gene2,density="F"){
   data("tpm")
  df=subset(tpm,Group=="Tumor" & Cancer==cancer)%>%
    select(all_of(gene1),all_of(gene2))

 if(density=="F"){
   p=ggscatter(df, x = gene1, y = gene2,
             add = "reg.line",
             add.params = list(color = "blue", fill = "lightgray"),
             conf.int = TRUE,
             cor.coef = TRUE, cor.method = "pearson",
             color = "black", size = 3, alpha = 0.6,
             xlab = gene1, ylab = gene2)
   print(p)
 }

 if(density=="T"){ggscatterhist(df, x = gene1, y = gene2,
               add = "reg.line",
               add.params = list(color = "blue", fill = "lightgray"),
               conf.int = TRUE,
               cor.coef = TRUE, cor.method = "pearson",
               xlab = gene1, ylab = gene2,
               color = "black", size = 3, alpha = 0.6,
               margin.params = list(fill = '#fccb7e', color = "#fccb7e", size = 0.3))
 }

}
