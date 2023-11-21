rm(list=ls())
setwd("./rawdata/")
library(msigdbr)
library(dplyr)
all_gene_sets = msigdbr(species = "Homo sapiens")%>%
  select(gs_name,gs_exact_source,gene_symbol)

msigb_genesets=split(all_gene_sets$gene_symbol,all_gene_sets$gs_name)
save(msigb_genesets,file="msigb_genesets.Rdata")
setwd("../")
rm(list=ls())
