#' gene_network_kegg
#'
#' @description Create a cnetplot to depict the linkages of gene(s) and KEGG pathways as a network.
#' @param gene gene name likes "KLF7", or a vector of gene names like c("LAMA3","LAMC2","TNC","OSMR").
#' @import  clusterProfiler
#' @import enrichplot
#' @export
#' @examples gene_network_kegg(c("LAMA3","LAMC2","TNC","OSMR"))
gene_network_kegg=function(gene){
  print("This function is performed by the clusterProfiler and enrichplot R packages.")
  s2e <- clusterProfiler::bitr(gene,
                               fromType = "SYMBOL",
                               toType = "ENTREZID",
                               OrgDb = org.Hs.eg.db)

  ekk <- enrichKEGG(s2e$ENTREZID)

  #transform ENTREZID to gene symbol
  ekk = setReadable(ekk,
                               OrgDb = "org.Hs.eg.db",
                               keyType = "ENTREZID")

  cnetplot(ekk, foldChange=s2e$ENTREZID, circular = TRUE, colorEdge = TRUE)
}
