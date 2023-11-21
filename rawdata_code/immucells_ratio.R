rm(list=ls())
setwd("./rawdata/")
load("./tpm.Rdata")
# https://gdc.cancer.gov/about-data/publications/panimmune contained TCGA.Kallisto.fullIDs.cibersort.relative.tsv file.
# download immune cell ratio fro tcga based on
# url="https://api.gdc.cancer.gov/data/b3df502e-3594-46ef-9f94-d041a20a0b9a"
immucells=read.table("TCGA.Kallisto.fullIDs.cibersort.relative.tsv",header = T,check.names = F)
colnames(immucells)=stringr::str_replace_all(colnames(immucells),'[.]',' ')
immucells$SampleID=stringr::str_replace_all(immucells$SampleID,'[.]','-')
immucells$SampleID=stringr::str_sub(immucells$SampleID,1,16)
immucells=immucells[!duplicated(immucells$SampleID),]
immucells=immucells[,-c(25:27)]
rownames(immucells)=immucells[,1]
immucells=immucells[,-(1:2)]
immucells=immucells[match(rownames(subset(tpm,Group=="Tumor")),rownames(immucells)),]
immucells=immucells[,order(colnames(immucells))]
immucells=round(immucells,4)
save(immucells,file="immucells.Rdata")
setwd("../")
rm(list=ls())
