
library(Signac)
library(GenomicRanges)
library(Seurat)
library(SeuratObject)


LST=read.table('LST.txt')
META=read.table('../data/GSE216323/GSM6668798_scATAC_metadata_new.txt',header=T,sep='\t',row.names=1)

pbmc_lst=c()
i=1
while(i<=nrow(LST)){
    this_pbmc=readRDS(paste0(LST[i,1],'.format/seuratWithMergePeak.rds'))
    this_pbmc@project.name=LST[i,2]
    OLD_NAME=colnames(this_pbmc)
    DefaultAssay(this_pbmc)='MERGE'
    this_pbmc[['peaks']]=NULL
    this_pbmc = RenameCells(this_pbmc, new.names = paste0(this_pbmc@project.name,'_',OLD_NAME) )
    this_pbmc=AddMetaData(this_pbmc, META)       
    if(this_pbmc@project.name=='h9esc'){this_pbmc$Clusters2=rep('H9_ESC',ncol(this_pbmc))}
    #this_pbmc=subset(this_pbmc, cells=colnames(this_pbmc)[which(!is.na(this_pbmc$Clusters2))])
    #this_pbmc <- RunTFIDF(this_pbmc)
    pbmc_lst=c(pbmc_lst,this_pbmc)
    print(i)
    i=i+1}

pbmc=merge(pbmc_lst[[1]],y=pbmc_lst[2:length(pbmc_lst)])
#######################################
tmp=t(matrix(unlist(stringr::str_split(colnames(pbmc),'_')),nrow=2))
pbmc$batch=tmp[,1]
#########################################

pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc,n=30)

pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = c(2:7))
#pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = c(2:10))

DimPlot(pbmc, group.by='Clusters2',label=T,raster=T)+NoLegend()

saveRDS(pbmc, file='../rds/pbmc.rds')


