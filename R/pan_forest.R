#' pan_forest
#'
#' @description Create a pan-cancer Cox regression forest plot for a specific gene.
#' @param gene gene name likes "KLF7".
#' @param adjust whether the Cox regression analysis was adjusted by age and stage. adjust=F is the default value.
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom tibble add_column
#' @importFrom tibble remove_rownames
#' @importFrom tibble column_to_rownames
#' @importFrom survival Surv
#' @importFrom survival coxph
#' @importFrom forestplot forestplot
#' @importFrom forestplot fpColors
#' @importFrom forestplot fpTxtGp
#' @importFrom grid gpar
#' @export
#' @examples pan_forest("KLF7",adjust=F)

pan_forest=function(gene,adjust=F){

  if(adjust==F){
  cox_results <- list()
  for (cancer in cancers) {
    exprSet=subset(tpm,Group=="Tumor" & Cancer==cancer)%>%
      tibble::add_column(ID = stringr::str_sub(rownames(.),1,12),
                         .before="Cancer") %>%
      dplyr::filter(!duplicated(ID)) %>%
      tibble::remove_rownames(.) %>%
      tibble::column_to_rownames("ID")%>%
      dplyr::filter(rownames(.) %in% rownames(subset(meta,Cancer==cancer)))
    exprSet=exprSet[,-(1:2)]
    exprSet=as.matrix(t(exprSet))

    cl=meta[colnames(exprSet),]

    cl$symbol = exprSet[gene,]

    # code from Xiaojie Sun from Biotrainee.
    m = survival::coxph(survival::Surv(time, event) ~symbol, data = cl)
    beta <- coef(m)
    se <- sqrt(diag(vcov(m)))
    HR <- exp(beta)
    HRse <- HR * se

    #summary(m)
    tmp <- round(cbind(coef = beta,
                       se = se, z = beta/se,
                       p = 1 - pchisq((beta/se)^2, 1),
                       HR = HR, HRse = HRse,
                       HRz = (HR - 1) / HRse,
                       HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                       HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                       HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)

    cox_results[[cancer]]=(tmp['symbol',])
  }

  cox_results=do.call(rbind,cox_results)
  cox_results=as.data.frame(cox_results[,c(5,9:10,4)])

  np = paste0(cox_results$HR,' (', cox_results$HRCILL,'-',cox_results$HRCIUL,')')

  ## The rest of the columns in the table.
  tabletext <- cbind(c("Cancer",rownames(cox_results)),
                     c("HR (95%CI)",np),
                     c("P Value",cox_results$p))
  ##plot
  forestplot::forestplot(labeltext=tabletext, graph.pos=3,
             mean=c(NA,cox_results$HR),
             lower=c(NA,cox_results$HRCILL), upper=c(NA,cox_results$HRCIUL),
             title=paste0("Hazard Ratio Plot of ",gene),
             hrzl_lines=list("1" = grid::gpar(lwd=2, col="black"),"2" = grid::gpar(lwd=2, col="black"),"35" = grid::gpar(lwd=2, col="black")),
             is.summary=c(TRUE,rep(FALSE,33)),
             col=forestplot::fpColors(box="#1c61b6", lines="#1c61b6", zero = "gray50"),
             zero=1,
             cex=0.9, lineheight = "auto",
             colgap=unit(8,"mm"),
             txt_gp = forestplot::fpTxtGp(ticks=grid::gpar(cex=1)),
             boxsize=0.5,
             ci.vertices=TRUE, ci.vertices.height = 0.3)
  }

  if(adjust==T){
    cox_results <- list()
    for (cancer in cancers) {
      exprSet=subset(tpm,Group=="Tumor" & Cancer==cancer)%>%
        tibble::add_column(ID = stringr::str_sub(rownames(.),1,12),
                           .before="Cancer") %>%
        dplyr::filter(!duplicated(ID)) %>%
        tibble::remove_rownames(.) %>%
        tibble::column_to_rownames("ID")%>%
        dplyr::filter(rownames(.) %in% rownames(subset(meta,Cancer==cancer)))
      exprSet=exprSet[,-(1:2)]
      exprSet=as.matrix(t(exprSet))

      cl=meta[colnames(exprSet),]

      cl$symbol = exprSet[gene,]

      # code from Xiaojie Sun from Biotrainee.
      m = survival::coxph(survival::Surv(time, event) ~symbol+age, data = cl)
      beta <- coef(m)
      se <- sqrt(diag(vcov(m)))
      HR <- exp(beta)
      HRse <- HR * se

      #summary(m)
      tmp <- round(cbind(coef = beta,
                         se = se, z = beta/se,
                         p = 1 - pchisq((beta/se)^2, 1),
                         HR = HR, HRse = HRse,
                         HRz = (HR - 1) / HRse,
                         HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                         HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                         HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)

      cox_results[[cancer]]=(tmp['symbol',])
    }

    cox_results=do.call(rbind,cox_results)
    cox_results=as.data.frame(cox_results[,c(5,9:10,4)])

    np = paste0(cox_results$HR,' (', cox_results$HRCILL,'-',cox_results$HRCIUL,')')

    ## The rest of the columns in the table.
    tabletext <- cbind(c("Cancer",rownames(cox_results)),
                       c("HR (95%CI)",np),
                       c("P Value",cox_results$p))
    ##plot
    forestplot::forestplot(labeltext=tabletext, graph.pos=3,
                           mean=c(NA,cox_results$HR),
                           lower=c(NA,cox_results$HRCILL), upper=c(NA,cox_results$HRCIUL),
                           title=paste0("Hazard Ratio Plot of ",gene, " adjusted by age"),
                           hrzl_lines=list("1" = grid::gpar(lwd=2, col="black"),"2" = grid::gpar(lwd=2, col="black"),"35" = grid::gpar(lwd=2, col="black")),
                           is.summary=c(TRUE,rep(FALSE,33)),
                           col=forestplot::fpColors(box="#1c61b6", lines="#1c61b6", zero = "gray50"),
                           zero=1,
                           cex=0.9, lineheight = "auto",
                           colgap=unit(8,"mm"),
                           txt_gp = forestplot::fpTxtGp(ticks=grid::gpar(cex=1)),
                           boxsize=0.5,
                           ci.vertices=TRUE, ci.vertices.height = 0.3)
  }
}
