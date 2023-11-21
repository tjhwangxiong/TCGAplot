#' pan_tumor_boxplot
#'
#' @description Create a pan-cancer box plot for a single gene in tumor samples.
#' @param gene gene names likes "KLF7".
#' @import ggpubr
#' @import RColorBrewer
#' @importFrom dplyr select
#' @export
#' @examples pan_tumor_boxplot("KLF7")
pan_tumor_boxplot=function(gene){
  data("tpm")
  ggboxplot(dplyr::select(subset(tpm,Group=="Tumor"),Cancer,all_of(gene)),
            x = "Cancer", y = gene,fill = "Cancer",
                 xlab='',
                 color = "black")+
    rotate_x_text(angle = 90)+
    grids(linetype="dashed")+
    theme(legend.title=element_text(size=14),
          legend.text=element_text(size=12),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size=14),
          axis.text.y = element_text(size = 12))+
    border("black")+
    theme(legend.position = "none")+
    scale_color_brewer(palette = "Set1")

}
