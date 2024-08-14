



library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
set.seed(123)


library(Signac)
library(GenomicRanges)

opath='/home/database/data/InferX/data/frontal_cortex/scATAC/GSM5589415_UMB4540_snATAC_frontal_cortex_rep2_fragments.hg19.bed.gz.format'

fpath=paste0(opath, '/sorted.bed.gz')
#bpath=paste0('/home/database/reference/hg19/hg19.fa.size.100k_coolBin.bed')
bpath=paste0('/home/database/reference/hg19/hg19.fa.size.10k_coolBin.bed')
rpath=paste0(opath,'/mat_10k.rds')



MAT=readRDS(rpath)

chrom_assay <- CreateChromatinAssay(
  counts = MAT,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = fpath,
  min.cells = 20,
  min.features = 500
  )


pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
  )


################################
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
annotations=readRDS('/home/database/annotation/hg19/hg19_signac_ucsc_annotations.rds')

Annotation(pbmc) <- annotations
##########################

pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)

DimPlot(pbmc,label=TRUE)+NoLegend()

###########################################
gene.activities <- GeneActivity(pbmc)

pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)
#######################################
saveRDS(gene.activities,file=paste0(opath, '/gene.activities.rds'))


library(Seurat)
library(Signac)

META=read.table('/home/database/data/InferX/data/frontal_cortex/scATAC/GSE184462_metadata.tsv',header=TRUE,sep='\t')
TAG=t(matrix(unlist(stringr::str_split(META[,1],'\\+')),nrow=2))[,2]
META$tag=TAG


UM=META[which(TAG %in% colnames(pbmc)),]
head(UM)
TAB=sort(table(UM$sample),decreasing=TRUE)
SAMPLE=names(TAB)

####################################################
UM=UM[which(UM$sample == SAMPLE[1]),]
rownames(UM)=UM$tag

used.pbmc=subset(pbmc, cells=UM$tag)
used.pbmc=AddMetaData(used.pbmc,UM)

DimPlot(used.pbmc,group.by='cell.type',label=TRUE)+NoLegend()

#####################################

pbmc=used.pbmc

####################################

pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3,resolution=0.5)



pbmc$tmp=pbmc$cell.type
TAB=table(pbmc$cell.type)
CUT=100

pbmc$tmp[which(pbmc$cell.type %in% names(which(TAB<CUT)))]=NA
plot1=DimPlot(object = pbmc, group.by='tmp', label = TRUE) + NoLegend()
plot2=DimPlot(object = pbmc, group.by='seurat_clusters', label = TRUE) + NoLegend()
plot1+plot2


pbmc$type=pbmc$cell.type
pbmc$type[which(pbmc$seurat_clusters %in% c(7))]='Astro'
pbmc$type[which(pbmc$seurat_clusters %in% c(9))]='NA'
pbmc$type[which(pbmc$seurat_clusters %in% c(4))]='Microglia'
pbmc$type[which(pbmc$seurat_clusters %in% c(1,0))]='Neuron.Ex.Glu'
pbmc$type[which(pbmc$seurat_clusters %in% c(3,5))]='Neuron.In.GABA'
pbmc$type[which(pbmc$seurat_clusters %in% c(6,8))]='ODC'
pbmc$type[which(pbmc$seurat_clusters %in% c(2))]='OPC'


Idents(pbmc)=pbmc$type

DimPlot(pbmc,group.by='type',label=TRUE)+NoLegend()


used.pbmc=subset(pbmc, cells= colnames(pbmc)[which(pbmc$type!='NA')] )
DimPlot(used.pbmc,group.by='type',label=TRUE)+NoLegend()

pbmc=used.pbmc

saveRDS(pbmc, 'scATAC_pbmc.rds')



pbmc=readRDS('scATAC_pbmc.rds')

peaks <- CallPeaks(pbmc, macs2.path = "/home/toolkit/local/bin/macs2")

saveRDS(peaks, file='scATAC_pbmc_peaks_macs2.rds')

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg19, invert = TRUE)


macs2_counts <- FeatureMatrix(
  fragments = Fragments(pbmc),
  features = peaks,
  cells = colnames(pbmc)
  )

pbmc[["macs2"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fpath,
  annotation = readRDS('/home/database/annotation/hg19/hg19_signac_ucsc_annotations.rds')
  )

saveRDS(pbmc, 'scATAC_pbmc.rds')


























































saveRDS(used.pbmc, file='seurat_object_annotated.rds')



















































































