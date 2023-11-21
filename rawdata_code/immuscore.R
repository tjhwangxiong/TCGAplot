rm(list=ls())
setwd("./rawdata/")
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos = rforge, dependencies = TRUE)
library(estimate)
library(stringr)
load("./tpm.Rdata")
if (!dir.exists("immuscore"))dir.create("immuscore")
setwd("./immuscore/")
estimate <- function(data,pro){
  dat=data[data$Cancer==pro,]
  dat=dat[,-(1:2)]
  dat=as.data.frame(t(dat))
  input.f=paste0(pro,'_estimate_input.txt')
  output.f=paste0(pro,'_estimate_gene.gct')
  output.ds=paste0(pro,'_estimate_score.gct')
  write.table(dat,file = input.f,sep = '\t',quote = F)
  filterCommonGenes(input.f=input.f,
                    output.f=output.f ,
                    id="GeneSymbol")
  estimateScore(input.ds = output.f,
                output.ds=output.ds,
                platform="illumina")   ## platform
  scores=read.table(output.ds,skip = 2,header = T,check.names = F)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  library(stringr)
  rownames(scores)=str_replace_all(rownames(scores),'[.]','-')
  return(scores)
}

projects = unique(tpm$Cancer)

immuscore=list()
for (project in projects) {
  immuscore[[project]]=round(estimate(tpm,project),2)
}

immuscore=do.call(rbind,lapply(immuscore, as.data.frame))
rownames(immuscore)=stringr::str_split(rownames(immuscore),'[.]',simplify = T)[,2]

immuscore=immuscore[rownames(immuscore) %in% rownames(subset(tpm,Group=="Tumor")),]

save(immuscore,file="../immuscore.Rdata")
setwd("../../")
rm(list=ls())
