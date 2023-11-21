rm(list=ls())
setwd("./rawdata")
load("./tpm.Rdata")
library(cBioPortalData)
library(dplyr)
library(stringr)
cbio <- cBioPortal()
studies = getStudies(cbio)
df=studies[str_detect(studies$studyId,'tcga'),]

clist=unique(df$studyId)
MSI=list()
for (id in clist){
  clinical = clinicalData(cbio,id)
  if ('MSI_SCORE_MANTIS' %in% colnames(clinical)) {
    MSI[[id]] = clinical[,c('patientId','MSI_SCORE_MANTIS')]
    MSI[[id]]$MSI_SCORE_MANTIS=round(as.numeric(MSI[[id]]$MSI_SCORE_MANTIS),3)
  }
}

MSI=as.data.frame(do.call(rbind,MSI))%>%
  na.omit()%>%
  dplyr::filter(!duplicated(patientId))%>%
  dplyr::filter(patientId %in% stringr::str_sub(rownames(subset(tpm,Group=="Tumor")),1,12))%>%
  tibble::remove_rownames(.)%>%
  dplyr::rename_with(~c("ID","MSI"),1:2)%>%
  tibble::remove_rownames(.)%>%
  tibble::column_to_rownames('ID')

save(MSI,file="./MSI.Rdata")
setwd("../")
rm(list=ls())
