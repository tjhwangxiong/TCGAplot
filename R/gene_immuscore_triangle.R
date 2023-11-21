#' gene_immunescore_triangle
#'
#' @description Create a pan-cancer triangle reveals the correlation between the expression of a single gene and immune scores, including Stromal score, immune score, and ESTIMATE score.
#' @param gene gene name likes "KLF7".
#' @param method method="pearson" is the default value. The alternatives to be passed to cor are "spearman" and "kendall".
#' @import ComplexHeatmap
#' @import grid
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom circlize colorRamp2
#' @importFrom reshape2 melt
#' @export
#' @examples gene_immunescore_triangle("KLF7")
gene_immunescore_triangle=function(gene,method="pearson"){
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
  clist=t(clist)
  plist=t(plist)

  ### Plot
  ## Code from: YuLabSMU https://baijiahao.baidu.com/s?id=1655042135206685129

  UpColor <- circlize::colorRamp2(breaks = c(-1,0,1), colors = c("#164382","#FFFFFF","#eb1717"))
  DnColor <- circlize::colorRamp2(breaks = c(0,1), colors = c("#FFFFFF","#FFA500"))

  row_an <- HeatmapAnnotation(type = c(rep("Stromal", 1), rep("Immune", 1), rep("ESTIMATE", 1)),
                              show_annotation_name = F,
                              col = list(type = c("ESTIMATE" = "#5AC9FA","Immune" = "#FAC67A","Stromal" = "#51B743")),
                              show_legend = T,
                              annotation_legend_param = list(title = "Immunescore",nrow = 1),
                              which = "row")

  DiagFunc <- function(up, down){
    function(j, i, x, y, width, height, fill){
      grid.polygon(grid::unit.c(x - 0.5*width, x - 0.5*width, x + 0.5*width),
                   grid::unit.c(y - 0.5*height, y + 0.5*height, y + 0.5*height),
                   gp = grid::gpar(fill = DnColor(down[i, j]), col = "grey"))
      grid.polygon(grid::unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width),
                   grid::unit.c(y + 0.5*height, y - 0.5*height, y - 0.5*height),
                   gp = grid::gpar(fill = UpColor(up[i, j]), col = "grey"))
    }
  }

  p1 <- Heatmap(clist, column_title = paste0("Correlation between ",gene, " and immunescores"),
                rect_gp = gpar(type = "none"),
                show_heatmap_legend = F,
                cluster_rows = F,
                show_row_names = F,
                cluster_columns = F,
                left_annotation = row_an,
                cell_fun = DiagFunc(up = clist, down = plist))

  lgd <- list(Legend(title = "Correlation",
                     col_fun = UpColor,
                     at = c(-1,0,1),
                     direction = "horizontal"),
              Legend(title = "Pvalue",
                     col_fun = DnColor,
                     at = c(0,0.5,1),
                     direction = "horizontal"))

  draw(p1, annotation_legend_list = lgd,
       annotation_legend_side = "bottom",
       heatmap_legend_side = "bottom",
       merge_legend = TRUE)

}
