library(Seurat)
library(Signac)

library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
annotations=readRDS('/home/database/annotation/hg38/hg38_signac_ucsc_annotations.rds')

source('../../../frontal_cortex/scATAC/InferC.R')

DATA=readRDS('../rds/DATA.rds')
UMAP=readRDS('../rds/pbmc_umap.rds')
TIME=readRDS('../rds/TIME.rds')

####################
genome.df=inferloop.getGenomeDF.hg38()
window=500000
sample_num=100

#####################
set.seed(123)
SI=sort(sample(1:ncol(DATA),5000))
indata=DATA[,SI]
used_coords=UMAP[SI,]

cicero_cds=inferc.ciceroFrame_step01_prepareInput(indata, used_coords=used_coords,genome.df=genome.df,window=window,sample_num=sample_num)

saveRDS(cicero_cds,'../rds/ciceroFrame_cicero_cds.rds')

#####################
used_function=inferc.rowCalD
conns <- inferc.ciceroFrame_step02_runUsedFuntion(cicero_cds, genome.df, window, sample_num, used_function)

saveRDS(conns,'../rds/ciceroFrame_conns_D.rds')

BATCH=readRDS('../rds/BATCH.rds')
AGG=inferc.aggMat(DATA, UMAP, clstNum=1000,seed=123, BATCH=BATCH)
saveRDS(AGG, file='../rds/AGG.rds')







