library(Seurat)
library(Signac)

library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)


source('../../../frontal_cortex/scATAC/InferC.R')

pbmc=readRDS('../rds/pbmc.rds')

DATA=pbmc[['peaks']]@data
UMAP=pbmc@reductions$umap@cell.embeddings

####################
genome.df=inferloop.getGenomeDF.hg19()
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




