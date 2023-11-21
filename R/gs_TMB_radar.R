#' gs_TMB_radar
#'
#' @description Create a pan-cancer radar chart for gene set and TMB correlation.
#' @param geneset gene set name likes "KEGG_APOPTOSIS".
#' @param method method="pearson" is the default value. The alternatives to be passed to correlation are "spearman" and "kendall".
#' @import dplyr
#' @import tibble
#' @importFrom magrittr %>%
#' @importFrom GSVA gsva
#' @importFrom fmsb radarchart
#' @export
#' @examples gs_TMB_radar("KEGG_APOPTOSIS")
gs_TMB_radar=function(geneset,method = "pearson"){

  df=dplyr::filter(tpm,Group=="Tumor")%>%
    dplyr::select(-Group)%>%
    dplyr::filter(rownames(.) %in% rownames(TMB))%>%
    dplyr::arrange(rownames(.))

  dat=TMB%>%
    dplyr::filter(rownames(.)%in% rownames(df))%>%
    dplyr::arrange(rownames(.))

  identical(rownames(df),rownames(dat))

  res=list()
  genelist=msigb_genesets[geneset]
  for (cancer in cancers) {
    dat1=subset(dat,rownames(dat) %in% rownames(subset(df,Cancer==cancer)))
    df1=df[rownames(dat1),]%>%
      dplyr::select(-Cancer)%>%
      t()

    #gsva
    df1 <- GSVA::gsva(expr=df1,
                           gset.idx.list=genelist,
                           kcdf="Gaussian",
                           verbose=T,
                           parallel.sz = parallel::detectCores())

    identical(colnames(df1),rownames(dat1))

    tmp <- cor.test(df1[names(genelist),], dat1$TMB,method = "pearson")
    res[[cancer]]=data.frame(rho=tmp[["estimate"]],pvalue=tmp[["p.value"]])
  }

  res=do.call(rbind,res)
  res$label=ifelse(res$pvalue<0.01,'**',
                                    ifelse(res$pvalue<0.05,'*',' '))
  res$group=paste0(rownames(res),res$label)
  res=res[,c(4,1)]%>%
    tibble::remove_rownames(.)%>%
    tibble::column_to_rownames("group")%>%
    tibble::add_column(Max=max(res$rho),.before="rho")%>%
    tibble::add_column(Min=min(res$rho),.before="rho")%>%
    t()%>%
    as.data.frame()%>%
    round(2)
  res=res[,c(1,33:2)]

  fmsb::radarchart(res,
             axistype = 1,
             # Customize the polygon
             pcol = "#00AFBB", pfcol = scales::alpha("#00AFBB", 0.1), plwd = 2, plty = 1,
             # Customize the grid
             cglcol = "grey", cglty = 1, cglwd = 0.8,
             # Customize the axis
             axislabcol = "black",
             title=paste0("Correlation between ",geneset, " and TMB"),
             # Variable labels
             vlcex = 1, vlabels = colnames(res),
             caxislabels = c(round(res["Min",1],1),
                             round(res["Min",1]+0.25*(res["Max",1]-res["Min",1]),1),
                             round(res["Min",1]+0.5*(res["Max",1]-res["Min",1]),1),
                             round(res["Min",1]+0.75*(res["Max",1]-res["Min",1]),1),
                             round(res["Max",1],1)))

}
