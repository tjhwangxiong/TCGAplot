rm(list=ls())
setwd("./rawdata/")
load("../tpm.Rdata")
if (!dir.exists("clinical"))dir.create("clinical")
setwd("./clinical")
library(TCGAbiolinks)
library(XML)
library(dplyr)
library(stringr)
library(tibble)
projects <- getGDCprojects()$project_id
projects <- projects[grepl('^TCGA', projects, perl=TRUE)]
projects <- projects[order(projects)]

projects1=projects[-which(projects %in% c("TCGA-GBM", "TCGA-LAML", "TCGA-LGG", "TCGA-PCPG", "TCGA-SARC"))]

# No stage information was found in project2.
projects2=c("TCGA-GBM", "TCGA-LAML", "TCGA-LGG", "TCGA-PCPG", "TCGA-SARC")

clinical1=list()
for (project in projects1) {
  query <- GDCquery(project = project,
                    data.category = "Clinical",
                    file.type = "xml")
  GDCdownload(query)
  dir1="/harmonized/Clinical/Clinical_Supplement/"
  xmls=dir(paste0('GDCdata/',project,dir1),pattern = "*.xml$",recursive = T)

  cl=list()

  use_col = c('bcr_patient_barcode',
              'vital_status',
              'days_to_death',
              'days_to_last_followup',
              'days_to_birth',
              'gender' ,
              'stage_event')

  for(i in 1:length(xmls)){
    result = xmlParse(paste0('GDCdata/',project,dir1,xmls[[i]]))
    rootnode = xmlRoot(result)
    dat=xmlToDataFrame(rootnode[2])
    if(sum(use_col %in% colnames(dat))==7)
      cl[[i]]=dat[,use_col]
  }

  clinical1[[project]] = do.call(rbind,cl)
}

clinical2=list()
for (project in projects2) {
  query <- GDCquery(project = project,
                    data.category = "Clinical",
                    file.type = "xml")
  GDCdownload(query)
  dir1="/harmonized/Clinical/Clinical_Supplement/"
  xmls=dir(paste0('GDCdata/',project,dir1),pattern = "*.xml$",recursive = T)

  cl=list()

  use_col= c('bcr_patient_barcode',
             'vital_status',
             'days_to_death',
             'days_to_last_followup',
             'days_to_birth',
             'gender')

  for(i in 1:length(xmls)){
    result = xmlParse(paste0('GDCdata/',project,dir1,xmls[[i]]))
    rootnode = xmlRoot(result)
    dat=xmlToDataFrame(rootnode[2])
    if(sum(use_col %in% colnames(dat))==6)
      cl[[i]]=dat[,use_col]
  }

  clinical2[[project]] = do.call(rbind,cl)
}

clinical1=as.data.frame(do.call(rbind,clinical1))
clinical2=as.data.frame(do.call(rbind,clinical2))
clinical2$stage_event=NA
identical(colnames(clinical1),colnames(clinical2))
meta=rbind(clinical1,clinical2)%>%
  tibble::add_column(Cancer = stringr::str_split(rownames(.),'-',simplify = T)[,2], .before="bcr_patient_barcode")
meta$Cancer=stringr::str_split(meta$Cancer,'[.]',simplify = T)[,1]

meta=meta[!duplicated(meta$bcr_patient_barcode),]
meta=meta[order(meta$Cancer,meta$bcr_patient_barcode),]
rownames(meta) <- meta$bcr_patient_barcode

# rename colname
colnames(meta)=c('Cancer','ID','event','death','last_followup','age','gender','stage')
# define missing value as NA
meta[meta==""]=NA

table(meta$Cancer)

# calculate survival time
table(meta$event)
# Alive  Dead
# 8405  2748
meta$time = ifelse(meta$event=="Alive",
                   meta$last_followup,
                   meta$death)
meta$time = round(as.numeric(meta$time)/30,2)
meta$age = round(-(as.numeric(meta$age)/365),0)

meta$gender=ifelse(meta$gender=="MALE","M","F")

# define event，Alive=0，Dead=1
meta$event=ifelse(meta$event=='Alive',
                  0,
                  1)
table(meta$event)
# 0     1
# 8405  2748

k1 = meta$time>=0.1;table(k1)
# FALSE  TRUE
# 598 10477
k2 = !(is.na(meta$time)|is.na(meta$event));table(k2)
# FALSE  TRUE
# 92 11075
meta = meta[k1&k2,]
nrow(meta)
# [1] 10477

# stage
meta$stage=str_split(meta$stage,"Stage ",simplify = T)[,2]
a = str_extract_all(str_sub(meta$stage,1,4),"I|V");head(a)
b = sapply(a,paste,collapse = "")
table(b)
#      I   II  III   IV
# 2562 2487 2206 2259  963

meta$stage = b
# set missing value as NA
table(meta$stage)
#      I   II  III   IV
# 2562 2487 2206 2259  963

meta[meta=="" | meta=="not reported"]=as.character(NA)
table(meta$stage,useNA = "always")
# I    II   III   IV  <NA>
# 2487 2206 2259  963 2562

meta=meta[,c('Cancer','event','time','age','gender','stage')]

# match or merge meta and expression
exprSet=dplyr::filter(tpm,Group=="Tumor")%>%
  tibble::add_column(ID = stringr::str_sub(rownames(.),1,12),
                     .before="Cancer") %>%
  dplyr::filter(!duplicated(ID)) %>%
  tibble::remove_rownames(.) %>%
  tibble::column_to_rownames("ID")%>%
  dplyr::select(-(1:2))

s = intersect(rownames(meta),rownames(exprSet));length(s)
# [1] 9678

meta = meta[s,]
dim(meta)
#[1] 9678    6

save(meta,file="../meta.Rdata")
setwd("../../")
