rm(list=ls())
setwd("./rawdata/TMB/")
load("../tpm.Rdata")
library(TCGAbiolinks)
library(dplyr)
library(stringr)
projects <- getGDCprojects()$project_id
projects <- projects[grepl('^TCGA', projects, perl=TRUE)]
projects <- projects[order(projects)]

TMB=list()
for (project in projects){
  query <- GDCquery(
    project = project,
    data.category = "Simple Nucleotide Variation",
    access = "open",
    data.type = "Masked Somatic Mutation",
    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")

  GDCdownload(query)

  mafFilePath2 = dir(path = paste0("GDCdata/",project),
                     pattern = "masked.maf.gz$",full.names = T,recursive=T)
  dat = lapply(mafFilePath2, data.table::fread, skip = "Hugo_Symbol")
  dat = data.table::rbindlist(l = dat, use.names = TRUE, fill = TRUE)

  ## TMB calculation
  # code reference https://zhuanlan.zhihu.com/p/394609586
  get_TMB <- function(file) {
    use_cols <- c("Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode",
                  "HGVSc", "t_depth", "t_alt_count")

    # read file
    df <- select(file, use_cols)
    data <- df %>%
      # calculate VAF
      mutate(vaf = t_alt_count / t_depth) %>%
      group_by(Tumor_Sample_Barcode) %>%
      summarise(mut_num = n(), TMB = mut_num / 30, MaxVAF = max(vaf))
    return(data)
  }

  TMB[[project]]=get_TMB(dat)
}

save(TMB,file="TMB_orig.Rdata")

TMB=as.data.frame(do.call(rbind,TMB))%>%
  dplyr::filter(!duplicated(Tumor_Sample_Barcode))%>%
  tibble::remove_rownames(.)%>%
  tibble::column_to_rownames('Tumor_Sample_Barcode')%>%
  dplyr::mutate(ID=stringr::str_sub(rownames(.),1,16))%>%
  dplyr::filter(!duplicated(ID))%>%
  dplyr::select(ID,TMB)%>%
  dplyr::filter(ID %in% rownames(subset(tpm,Group=="Tumor")))%>%
  tibble::remove_rownames(.)%>%
  tibble::column_to_rownames('ID')%>%
  round(2)

save(TMB,file="../TMB.Rdata")
setwd("../../")
rm(list=ls())
