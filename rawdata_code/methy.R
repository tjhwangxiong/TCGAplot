setwd("./rawdata/methy/")
load("../tpm.Rdata")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(tibble)
library(ChAMP)
projects <- getGDCprojects()$project_id
projects <- projects[grepl('^TCGA', projects, perl=TRUE)]
projects <- projects[order(projects)]

data(hm450.manifest.hg19)
data(probe.features)

# get probe with gene name
all_probe = probe.features[probe.features$gene!='',] %>%
  tibble::rownames_to_column("probe")%>%
  dplyr::select(probe,feature,gene)%>%
  dplyr::arrange(gene)

# get methylation beta value
methy_before_rbind=list()
for (project in projects) {
  query.met <- GDCquery(project = project,
                        data.category = "DNA Methylation",
                        data.type = "Methylation Beta Value",
                        platform = c("Illumina Human Methylation 450"))

  #GDCdownload(query.met)
  expdat=GDCprepare(query = query.met)
  count <- assay(expdat)%>%
    t()%>%
    as.data.frame()%>%
    tibble::add_column(ID = stringr::str_sub(rownames(.),1,16),.before=1)%>%
    filter(!duplicated(ID))%>%
    tibble::remove_rownames(.)%>%
    tibble::column_to_rownames("ID")%>%
    tibble::add_column(Cancer=stringr::str_split(project,'-',simplify=T)[,2],.before = 1)

  saveRDS(count,file=paste0(project,"_methy.rds"))

  methy_before_rbind[[project]]=count
}

save(methy_before_rbind,file="methy_before_rbind.Rdata")

methy=do.call(rbind,methy_before_rbind)
rownames(methy)=stringr::str_split(rownames(methy),'[.]',simplify = T)[,2]
methy=methy%>%
  tibble::add_column(Group = factor(ifelse(as.numeric(stringr::str_sub(rownames(.),14,15)) < 10,
                                           'Tumor','Normal'),
                                    levels=c("Normal","Tumor")),.before=2)

met=methy[,1:2]
df=methy[,3:ncol(methy)]%>%
  t()%>%
  na.omit()%>%
  t()%>%
  as.data.frame()%>%
  round(3)

df=df[,colnames(df) %in% all_probe$probe]

identical(rownames(met),rownames(df))
#[1] TRUE

methy=cbind(met,df);rm(df,met)

methy=methy[rownames(methy) %in% rownames(tpm),]

all_probe=all_probe[all_probe$probe %in% colnames(methy),]

promoter_probe = subset(all_probe,feature=="TSS1500")%>%
  dplyr::select(probe,gene)

save(methy,all_probe,file="../TCGA_Methylation.Rdata") # all probe methylation

pro_methy=subset(methy,Group=="Tumor")
pro_methy1=pro_methy[,1:2]
pro_methy2=pro_methy[,colnames(pro_methy) %in% promoter_probe$probe]
pro_methy=cbind(pro_methy1[,1],pro_methy2)
colnames(pro_methy)[1]="Cancer"

promoter_probe=promoter_probe[promoter_probe$probe %in% colnames(pro_methy),]

save(pro_methy,promoter_probe,file="../TCGA_Pomoter_Methylation.Rdata")

setwd("../../")
rm(list=ls())
