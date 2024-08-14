


library(Seurat)
library(Signac)
source('/home/toolkit/src/fitdevo.R')
source('../../frontal_cortex/scATAC/InferC.R')

pbmc=readRDS('./rds/scmul_pbmc.rds')

TYPE=pbmc$celltype
USED_TYPE=c('Astrocytes','EN','EN-fetal-early','EN-fetal-late',
            'IN-CGE','IN-fetal','IN-MGE','IPC',
            'Oligodendrocytes','OPC','RG')

pbmc_sub=subset(pbmc, cells=colnames(pbmc)[which(TYPE %in% USED_TYPE)])


###################################



AGE=readRDS('./rds/scmul_pbmc_usedType_age.rds')

GA=pbmc_sub[['RNA']]@data

ES=read.table('./GS/BHATTACHARYA_EMBRYONIC_STEM_CELL.gmt',sep='\t')[,1]
GA_ES=GA[which(rownames(GA) %in% ES),]

SS4=apply(GA_ES,2,mean)
cor(SS4,AGE, method='spearman')

saveRDS(SS4,file='./rds/geneAct_ES_score.rds')

saveRDS(pbmc_sub, file='./rds/scmul_pbmc_usedType.rds')
