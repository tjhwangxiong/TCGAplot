#' code to prepare `tpm_log2` dataset

setwd("./data-raw/count")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(tibble)

ltpm=function(project){
  # query
  query <- GDCquery(project = project,
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "STAR - Counts"
  )
  # download
  GDCdownload(query, method = "api", files.per.chunk = 100) # 100 files once

  # data prepare
  data=GDCprepare(query)
  # Only log2(tpm+1) of protein coding genes was used
  tpmdata=round(log2(SummarizedExperiment::assay(data,4)+1),2) # log2(tpm+1)
  gene_id=data.frame(id=rowData(data)@listData[["gene_id"]],
                     gene_name= rowData(data)@listData[["gene_name"]],
                     gene_type=rowData(data)@listData[["gene_type"]])
  identical(rownames(tpmdata),gene_id$id)
  dat=cbind(gene_id,tpmdata) %>%
    dplyr::filter(gene_type=="protein_coding") %>%     # only mRNA was selected
    dplyr::select(-c(1,3)) %>%                         # remove id and gene_type columns
    dplyr::filter(!duplicated(gene_name)) %>%          # remove duplicated gene
    tibble::remove_rownames(.) %>%
    tibble::column_to_rownames("gene_name") %>% # replace rowname from Ensembl id to gene name
    t() %>% as.data.frame() %>%
    tibble::add_column(Group = factor(ifelse(as.numeric(stringr::str_sub(rownames(.),14,15)) < 10,
                                             'Tumor','Normal'),
                                      levels=c("Normal","Tumor")),
                       .before="TSPAN6") %>% # add Group column before 'TSPAN6' column
    tibble::add_column(Cancer = stringr::str_split(project,'-',simplify = T)[,2],
                       .before="Group") %>%
    dplyr::mutate(ID=stringr::str_sub(rownames(.),1,16)) %>%
    dplyr::filter(!duplicated(ID)) %>% # remove duplicated samples randomly
    tibble::remove_rownames(.) %>%
    tibble::column_to_rownames("ID")
  return(dat)
}

projects <- getGDCprojects()$project_id
projects <- projects[grepl('^TCGA', projects, perl=TRUE)]
projects <- projects[order(projects)]

tpm=list()
for (project in projects) {
  tpm[[project]]=ltpm(project)
}

tpm=as.data.frame(do.call(rbind,tpm))
rownames(tpm)=stringr::str_split(rownames(tpm),'[.]',simplify = T)[,2]

tpm=tpm[,!colSums(tpm==0)==nrow(tpm)] # remove column of 0 cross all samples
save(tpm,file="../tpm.Rdata")

setwd("../../")

