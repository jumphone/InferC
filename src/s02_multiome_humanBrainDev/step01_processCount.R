


library(Seurat)
library(Signac)
atac_count=readRDS('./rds/GSE204682_count_matrix_snATAC.RDS')
rna_count=readRDS('./rds/GSE204683_count_matrix_snRNA.RDS')
META=read.table('./rds/scp_metadata.csv',sep=',',row.names=1,header=T)
META=META[2:nrow(META),]


pbmc <- CreateSeuratObject(counts = rna_count,meta.data=META)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

pbmc[['ATAC']]=CreateAssayObject(counts=atac_count)

DefaultAssay(pbmc)='RNA'
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


DefaultAssay(pbmc)='ATAC'
pbmc <- FindTopFeatures(pbmc, min.cutoff = 5)
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)


pbmc <- FindMultiModalNeighbors(
  object = pbmc,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

pbmc <- RunUMAP(
  object = pbmc,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)


DimPlot(pbmc,group.by='celltype',label=T )+NoLegend()

saveRDS(pbmc,file='./rds/scmul_pbmc.rds')



set.seed(123)
random_cells=sample(colnames(pbmc),5000)
pbmc_sub=subset(pbmc, cells=random_cells)
DimPlot(pbmc_sub,group.by='celltype',label=T )+NoLegend()
saveRDS(pbmc_sub,file='./rds/scmul_pbmc_sub.rds')





















