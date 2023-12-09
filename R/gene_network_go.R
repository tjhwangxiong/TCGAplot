#' gene_network_go
#'
#' @description Create a cnetplot to depict the linkages of gene(s) and GO terms as a network.
#' @param gene gene name likes "KLF7", or a vector of gene names like c("LAMA3","LAMC2","TNC","OSMR").
#' @import  clusterProfiler
#' @import enrichplot
#' @export
#' @examples gene_network_go(c("LAMA3","LAMC2","TNC","OSMR"))
gene_network_go=function(gene){
  print("This function is performed by the clusterProfiler and enrichplot R packages.")
  s2e <- clusterProfiler::bitr(gene,
                               fromType = "SYMBOL",
                               toType = "ENTREZID",
                               OrgDb = org.Hs.eg.db)
  ego <- enrichGO(s2e$ENTREZID, OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)

  ## remove redundent GO terms
  ego2 <- simplify(ego)

  cnetplot(ego2, foldChange=s2e$ENTREZID, circular = TRUE, colorEdge = TRUE)

}
