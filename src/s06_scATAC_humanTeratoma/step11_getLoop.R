library(Seurat)
library(Signac)
source('/home/toolkit/src/fitdevo.R')
source('../../frontal_cortex/scATAC/InferC.R')

pbmc=readRDS('../rds/pbmc.rds')

UMAP=pbmc@reductions$umap@cell.embeddings
BATCH=as.character(pbmc$batch)
AGG=inferc.aggMat(pbmc[['MERGE']]@data, UMAP, clstNum=1000,seed=123, BATCH=BATCH)
saveRDS(AGG, file=paste0('../rds/AGG.rds'))

atac_mat=pbmc[['MERGE']]@data
    set.seed(123)
    random_index=sample(1:ncol(atac_mat),5000)
    atac_mat_sub=atac_mat[,random_index]
    UMAP_sub=UMAP[random_index,]
    used_coords=UMAP_sub
    indata=atac_mat_sub
    genome.df=inferloop.getGenomeDF.hg19()

    window=500000
    sample_num=100
    set.seed(123)
    cicero_cds=inferc.ciceroFrame_step01_prepareInput(indata, used_coords=used_coords,genome.df=genome.df,window=window,sample_num=sample_num)
    saveRDS(cicero_cds,paste0('../rds/ciceroFrame_cicero_cds.rds'))
    used_function=inferc.rowCalD
    conns <- inferc.ciceroFrame_step02_runUsedFuntion(cicero_cds, genome.df, window, sample_num, used_function)
    saveRDS(conns,paste0('../rds/ciceroFrame_conns_D.rds'))






