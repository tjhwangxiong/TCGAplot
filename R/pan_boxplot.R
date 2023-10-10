#' pan_boxplot
#'
#' @description Create a pan-cancer box plot for a single gene with symbols indicating statistical significance.
#' @param palette the color palette to be used for filling by groups. Allowed values include scientific journal palettes from ggsci R package, e.g.: "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty".
#' @param gene gene names likes "KLF7".
#' @param legend legend position. Allowed values include "top","bottom","left","right" and "none".
#' @import ggpubr
#' @importFrom dplyr select
#' @export
#' @examples pan_boxplot("KLF7")
pan_boxplot=function(gene,palette="jco",legend="right"){
  p <- ggboxplot(dplyr::select(tpm,Cancer,Group,all_of(gene)), x = "Cancer", y = gene,fill = "Group",
                 xlab='',
                 color = "black", palette = palette)+
    rotate_x_text(angle = 90)+
    grids(linetype="dashed")+
    theme(legend.title=element_text(size=14),
          legend.text=element_text(size=12),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size=14),
          axis.text.y = element_text(size = 12))+
    border("black")+
    theme(legend.position = legend)
  p + stat_compare_means(aes(group = Group),
                         label = "p.signif",
                         label.y.npc = "top",
                         hide.ns = T)
}
