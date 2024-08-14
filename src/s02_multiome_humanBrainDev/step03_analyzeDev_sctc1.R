
###############################
library(Seurat)
library(Signac)
source('/home/toolkit/src/fitdevo.R')
source('../../frontal_cortex/scATAC/InferC.R')

META=readRDS('./rds/scmul_pbmc_usedType_meta.rds')
UMAP=readRDS('./rds/scmul_pbmc_usedType_umap.rds')
PCA=readRDS('./rds/scmul_pbmc_usedType_pca.rds')
LSI=readRDS('./rds/scmul_pbmc_usedType_lsi.rds')
DP=readRDS('./rds/scmul_pbmc_usedType_dp.rds')
AGE=readRDS('./rds/scmul_pbmc_usedType_age.rds')


pbmc=readRDS('./rds/scmul_pbmc.rds')
TYPE=pbmc$celltype
USED_TYPE=c('Astrocytes','EN','EN-fetal-early','EN-fetal-late',
            'IN-CGE','IN-fetal','IN-MGE','IPC',
            'Oligodendrocytes','OPC','RG')
pbmc_sub=subset(pbmc, cells=colnames(pbmc)[which(TYPE %in% USED_TYPE)])
#########################################

set.seed(123)
random_index=sort(sample((1:ncol(pbmc_sub)),10000))
saveRDS(random_index,'./rds/sctc_random_index_10k.rds')

###########################################################
DATA=pbmc_sub[['RNA']]@data[,random_index]
rna=CreateSeuratObject(counts = DATA, project = "rna", min.cells = 1, min.features = 0)

###################################################################
library(SeuratData)
library(SeuratDisk)
sce_obj <- as.SingleCellExperiment(rna, assay = c("RNA"))
zellkonverter::writeH5AD(sce_obj, "./rds/sctc_10k.h5ad",'counts')
##################
# Run Python

###################################

library(SeuratData)
library(SeuratDisk)

DATA=pbmc_sub[['RNA']]@data
rna=CreateSeuratObject(counts = DATA, project = "rna", min.cells = 1, min.features = 0)

###################################################################
library(SeuratData)
library(SeuratDisk)
sce_obj <- as.SingleCellExperiment(rna, assay = c("RNA"))
zellkonverter::writeH5AD(sce_obj, "./rds/sctc_all.h5ad",'counts')






