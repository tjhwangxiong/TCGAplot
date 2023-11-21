#' gene_methylation_scatter
#' @description Scatter plot of gene expression and gene promoter methylation correlation in a specific type of cancer. A pdf file will be generated in the working directory.
#' @param cancer cancer name likes "BLCA".
#' @param gene gene name likes "KLF7".
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom ggpubr ggscatter
#' @importFrom ggpubr ggscatterhist
#' @export
#' @examples gene_methylation_scatter("BLCA","KLF7")
gene_methylation_scatter=function(cancer,gene){
  data("tpm")
  if(gene %in% promoter_probe$gene){
    print('A pdf file will be generated in the working directory.')

    df1=subset(tpm,Group=="Tumor" & Cancer==cancer)%>%
      dplyr::select(all_of(gene))%>%
      dplyr::filter(rownames(.) %in% rownames(pro_methy))

    cg_list=promoter_probe$probe[promoter_probe$gene==gene]
    cg_list = cg_list[cg_list %in% colnames(pro_methy)]

    pdf(paste0(gene, "_methylation_correlation_",cancer,".pdf"),width = 5,height = 4)
    for (cg in cg_list){
      df2=pro_methy[rownames(df1),]%>%
        dplyr::select(all_of(cg))

      identical(rownames(df1),rownames(df2))
      df=cbind(df1,df2)
      colnames(df)=c("expression","methylation")

      p=ggpubr::ggscatter(df, x = "expression", y = "methylation",
                          title=cancer,
                          font.title = c(12, "bold", "black"),
                          add = "reg.line", conf.int = TRUE,
                          add.params = list(color = "blue", fill = "lightgray"),
                          cor.coef = TRUE, cor.method = "pearson",
                          color = "black", size = 3, alpha = 0.6,
                          xlab=paste0(gene," expression"),ylab=cg)
      print(p)
    }
    dev.off()
  } else{print('This gene is not included in the methylation matrix!')}

}
