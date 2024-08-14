library(Seurat)
library(Signac)

library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
annotations=readRDS('/home/database/annotation/hg38/hg38_signac_ucsc_annotations.rds')


source('../../../frontal_cortex/scATAC/InferC.R')

d1=readRDS('../rds/f1_d53/geneAct_data.rds')
d2=readRDS('../rds/f2_d59/geneAct_data.rds')
d3=readRDS('../rds/f3_d78/geneAct_data.rds')
d4=readRDS('../rds/f4_d78_2/geneAct_data.rds')
d5=readRDS('../rds/f5_d89C/geneAct_data.rds')
d6=readRDS('../rds/f6_d89P/geneAct_data.rds')

colnames(d1)=paste0('d53_',colnames(d1))
colnames(d2)=paste0('d59_',colnames(d2))
colnames(d3)=paste0('d78_',colnames(d3))
colnames(d4)=paste0('d78.2_',colnames(d4))
colnames(d5)=paste0('d89C_',colnames(d5))
colnames(d6)=paste0('d89P_',colnames(d6))

D12=.simple_combine(d1, d2,FILL=TRUE)$combine
D34=.simple_combine(d3, d4,FILL=TRUE)$combine
D56=.simple_combine(d5, d6,FILL=TRUE)$combine

DATA=.simple_combine(D12,D34,FILL=TRUE)$combine
DATA=.simple_combine(DATA,D56,FILL=TRUE)$combine

saveRDS(DATA,'../rds/GeneActivity.rds')



















