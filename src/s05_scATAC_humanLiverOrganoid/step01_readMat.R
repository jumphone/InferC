library(Seurat)
library(Signac)

library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)


ReadPeak<-function(PATH){
    counts <- Matrix::readMM(paste0(PATH, "_matrix.mtx"))
    barcodes <- readLines(paste0(PATH, "_barcodes.tsv"))
    peaks <- read.table(paste0(PATH,"_peaks.bed.restoreAll.bed"), sep="\t")
    peaknames <- paste(peaks$V1, peaks$V2, peaks$V3, sep="-")
    colnames(counts) <- barcodes
    rownames(counts) <- peaknames
    return(counts)
    }

hm1.data=ReadPeak('../data/GSE159557_RAW/GSM4832487_HM-1_hg19')
hm2.data=ReadPeak('../data/GSE159557_RAW/GSM4832488_HM-2_hg19')
dm1.data=ReadPeak('../data/GSE159557_RAW/GSM4832485_DM-1_hg19')
dm2.data=ReadPeak('../data/GSE159557_RAW/GSM4832486_DM-2_hg19')


.getUniqMat<-function(mat){
    #mat=hm1.data
    uniq_rname=names(table(rownames(mat)))
    umat=matrix(0,nrow=length(uniq_rname),ncol=ncol(mat))
    rownames(umat)=uniq_rname
    colnames(umat)=colnames(mat)
    umat=t(umat)
    mat=t(mat)
    mat=as.matrix(mat)
    i=1
    while(i<=ncol(umat)){
       this_index=which(colnames(mat)==uniq_rname[i])
       if(length(this_index)==1){
           umat[,i]=mat[,this_index]
           }else{
           umat[,i]=rowSums(mat[,this_index])
           }
       if(i%%10000==1){print(paste0(i,' / ',ncol(umat)))}
       i=i+1}
    out=t(umat)
    return(out)
    }


hm1.udata=.getUniqMat(hm1.data)
hm2.udata=.getUniqMat(hm2.data)
dm1.udata=.getUniqMat(dm1.data)
dm2.udata=.getUniqMat(dm2.data)

saveRDS(hm1.udata,'../rds/hm1.udata.rds')
saveRDS(hm2.udata,'../rds/hm2.udata.rds')
saveRDS(dm1.udata,'../rds/dm1.udata.rds')
saveRDS(dm2.udata,'../rds/dm2.udata.rds')

count2seurat<-function(counts){
    chrom_assay <- CreateChromatinAssay(
        counts = counts,
        sep = c(":", "-"),
        min.cells = 10,
        min.features = 200
        )
    pbmc <- CreateSeuratObject(
        counts = chrom_assay,
        assay = "peaks"
       )
    return(pbmc)
    }


hm1=count2seurat(hm1.udata)
hm2=count2seurat(hm2.udata)
dm1=count2seurat(dm1.udata)
dm2=count2seurat(dm2.udata)


hm1 <- RunTFIDF(hm1)
hm2 <- RunTFIDF(hm2)
dm1 <- RunTFIDF(dm1)
dm2 <- RunTFIDF(dm2)

saveRDS(hm1,'../rds/hm1.rds')
saveRDS(hm2,'../rds/hm2.rds')
saveRDS(dm1,'../rds/dm1.rds')
saveRDS(dm2,'../rds/dm2.rds')


BATCH=c(rep('hm1',ncol(hm1)),
       rep('hm2',ncol(hm2)),
       rep('dm1',ncol(dm1)),
       rep('dm2',ncol(dm2))
        )

pbmc=merge(hm1,c(hm2,dm1,dm2))
pbmc$batch=BATCH


pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q5')
pbmc <- RunSVD(pbmc,n = 30)
pbmc <- RunUMAP(object = pbmc,reduction = 'lsi', dims = 2:30 )

UMAP=pbmc@reductions$umap@cell.embeddings

DimPlot(pbmc,group.by='batch',label=TRUE)+NoLegend()

saveRDS(UMAP,'../rds/pbmc_umap.rds')
saveRDS(pbmc,'../rds/pbmc.rds')
saveRDS(pbmc@reductions$lsi@cell.embeddings,'../rds/LSI.rds')
saveRDS(pbmc@meta.data,'../rds/META.rds')


source('../../../frontal_cortex/scATAC/InferC.R')
BATCH=pbmc$batch
AGG=inferc.aggMat(pbmc[['peaks']]@data, UMAP, clstNum=1000,seed=123, BATCH=BATCH)
saveRDS(AGG, file='../rds/AGG.rds')



