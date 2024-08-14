library(Seurat)
library(Signac)

library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
annotations=readRDS('/home/database/annotation/hg38/hg38_signac_ucsc_annotations.rds')

source('../../../frontal_cortex/scATAC/InferC.R')

DATA=readRDS('../rds/DATA.rds')
BATCH=readRDS('../rds/BATCH.rds')
TIME=readRDS('../rds/TIME.rds')

pbmc=CreateSeuratObject(counts = DATA)
gc()

pbmc$time=TIME
pbmc$batch=BATCH

pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q5')
pbmc <- RunSVD(pbmc,n = 30)
pbmc <- RunUMAP(object = pbmc,reduction = 'lsi', dims = 2:30 )

DimPlot(object = pbmc, label = T, group.by='batch' ) + NoLegend()

UMAP=pbmc@reductions$umap@cell.embeddings

saveRDS(UMAP,'../rds/pbmc_umap.rds')
saveRDS(pbmc,'../rds/pbmc.rds')
saveRDS(pbmc@reductions$lsi@cell.embeddings,'../rds/LSI.rds')
saveRDS(pbmc@meta.data,'../rds/META.rds')





