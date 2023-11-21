rm(list=ls())
setwd("./rawdata")
load("./tpm.Rdata")
library(dplyr)
dat=tpm%>%
  dplyr::mutate(ID=stringr::str_sub(rownames(.),1,12))

Normal=dplyr::filter(dat,Group=="Normal")
Tumor=dplyr::filter(dat,Group=="Tumor")

cancers=unique(dat$Cancer)

T_only=list()
for (cancer in cancers){
  df=dplyr::filter(Tumor,Cancer==cancer)%>%
    dplyr::filter(!duplicated(ID)) # preserve only one tissue for one tumor sample

  T_only[[cancer]]=df
}

T_only = as.data.frame(do.call(rbind,T_only))
rownames(T_only)=str_split(rownames(T_only),'[.]',simplify = T)[,2]

index <- intersect(Normal$ID,T_only$ID) # samples with both normal and tumor tissues
T1=dplyr::filter(T_only, ID %in% index)
N1=dplyr::filter(Normal, ID %in% index)
paired_tpm=rbind(T1,N1)
paired_tpm=paired_tpm[,c(1,2,ncol(paired_tpm),3:(ncol(paired_tpm)-1))]

## filter cancers with sample number <20

selected= paired_tpm %>%
  group_by(Cancer) %>%
  summarise(num = n()) %>%
  subset(num>20)
pcancers=selected$Cancer
paired_tpm = paired_tpm[paired_tpm$Cancer %in% pcancers,]
paired_tpm=paired_tpm[order(paired_tpm$Cancer,paired_tpm$ID),]
save(paired_tpm,file="paired_tpm.Rdata")

setwd("../")
rm(list=ls())
