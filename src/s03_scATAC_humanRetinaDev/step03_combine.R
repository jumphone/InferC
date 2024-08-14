library(Seurat)
library(Signac)

library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
annotations=readRDS('/home/database/annotation/hg38/hg38_signac_ucsc_annotations.rds')


source('../../../frontal_cortex/scATAC/InferC.R')


s1=readRDS('../rds/f1_d53/pbmc_merge.rds')
s2=readRDS('../rds/f2_d59/pbmc_merge.rds')
s3=readRDS('../rds/f3_d78/pbmc_merge.rds')
s4=readRDS('../rds/f4_d78_2/pbmc_merge.rds')
s5=readRDS('../rds/f5_d89C/pbmc_merge.rds')
s6=readRDS('../rds/f6_d89P/pbmc_merge.rds')


s1 <- RunTFIDF(s1)
s2 <- RunTFIDF(s2)
s3 <- RunTFIDF(s3)
s4 <- RunTFIDF(s4)
s5 <- RunTFIDF(s5)
s6 <- RunTFIDF(s6)

BATCH=c(rep('d53',ncol(s1)),
       rep('d59',ncol(s2)),
       rep('d78',ncol(s3)),
       rep('d78.2',ncol(s4)),
       rep('d89C',ncol(s5)),
       rep('d89P',ncol(s6))
        )

TIME=c(rep(1,ncol(s1)),
       rep(2,ncol(s2)),
       rep(3,ncol(s3)),
       rep(3,ncol(s4)),
       rep(4,ncol(s5)),
       rep(4,ncol(s6))
        )


d1=s1[['MERGE']]@data
d2=s2[['MERGE']]@data
d3=s3[['MERGE']]@data
d4=s4[['MERGE']]@data
d5=s5[['MERGE']]@data
d6=s6[['MERGE']]@data

rm(s1)
rm(s2)
rm(s3)
rm(s4)
rm(s5)
rm(s6)
gc()

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

saveRDS(DATA,'../rds/DATA.rds')
saveRDS(BATCH,'../rds/BATCH.rds')
saveRDS(TIME,'../rds/TIME.rds')



















