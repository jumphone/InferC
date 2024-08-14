

library(Seurat)
library(Signac)

pbmc=readRDS('scATAC_pbmc.rds')
head(pbmc@meta.data)
DimPlot(pbmc,group.by='type')
MAT=as.matrix(pbmc[['macs2']]@data)
#MAT=as.matrix(pbmc[['peaks']]@data)
MAT[1:5,1:5]
VEC=pbmc@reductions$umap@cell.embeddings
TYPE=pbmc$type

source('InferC.R')
PEAKS=inferloop.splitLoop(rownames(MAT),'-',n=3)
write.table(PEAKS, file='./rds/scatac_macs2_peaks.bed',sep='\t',quote=F,row.names=F,col.names=F)



##################################################

source('InferC.R')
indata=MAT
used_coords=VEC
genome.df=inferloop.getGenomeDF.hg19()
window=500000
sample_num=100


####################################################

source('InferC.R')
cicero_cds=inferc.ciceroFrame_step01_prepareInput(indata, used_coords=used_coords,genome.df=genome.df,window=window,sample_num=sample_num)

saveRDS(cicero_cds,'./rds/ciceroFrame_cicero_cds.rds')


####################################################


source('InferC.R')
cicero_cds=readRDS('./rds/ciceroFrame_cicero_cds.rds')


used_function=used_function_cor
conns <- inferc.ciceroFrame_step02_runUsedFuntion(cicero_cds, genome.df, window, sample_num, used_function)
saveRDS(conns,'./rds/ciceroFrame_conns_cor.rds')


used_function=used_function_spearman
conns <- inferc.ciceroFrame_step02_runUsedFuntion(cicero_cds, genome.df, window, sample_num, used_function)
saveRDS(conns,'./rds/ciceroFrame_conns_spearman.rds')


conns <- inferc.ciceroFrame_step02_runCicero(cicero_cds, genome.df, window, sample_num)
saveRDS(conns,'./rds/ciceroFrame_conns_cicero.rds')


used_function=inferc.rowCalD
conns <- inferc.ciceroFrame_step02_runUsedFuntion(cicero_cds, genome.df, window, sample_num, used_function)
saveRDS(conns,'./rds/ciceroFrame_conns_D.rds')


#######################



conns_D=readRDS('./rds/ciceroFrame_conns_D.rds')
conns_cor=readRDS('./rds/ciceroFrame_conns_cor.rds')
conns_cicero=readRDS('./rds/ciceroFrame_conns_cicero.rds')
conns_spearman=readRDS('./rds/ciceroFrame_conns_spearman.rds')


conns_D_uniq=inferloop.getUniqLoop(conns_D)

.addRowName<-function(mat){
    rowName=apply(mat[,c(1:2)],1,paste0,collapse='.And.')
    rownames(mat)=rowName
    return(mat)
    }

conns_D_uniq=.addRowName(conns_D_uniq)
conns_cor=.addRowName(conns_cor)
conns_cicero=.addRowName(conns_cicero)
conns_spearman=.addRowName(conns_spearman)

conns_cor_uniq=conns_cor[rownames(conns_D_uniq),]
conns_cicero_uniq=conns_cicero[rownames(conns_D_uniq),]
conns_spearman_uniq=conns_spearman[rownames(conns_D_uniq),]


saveRDS(conns_D_uniq,'./rds/ciceroFrame_conns_D_uniq.rds')
saveRDS(conns_cor_uniq,'./rds/ciceroFrame_conns_cor_uniq.rds')
saveRDS(conns_cicero_uniq,'./rds/ciceroFrame_conns_cicero_uniq.rds')
saveRDS(conns_spearman_uniq,'./rds/ciceroFrame_conns_spearman_uniq.rds')


#####################################




BED1=inferloop.splitLoop(conns_D_uniq[,1],'_',3)
BED2=inferloop.splitLoop(conns_D_uniq[,2],'_',3)

BEDPE=cbind(BED1,BED2)

write.table(BEDPE, './rds/ciceroFrame_conns_uniq.bedpe',quote=F,row.names=F,col.names=F,sep='\t')




show_index=1:1000
plot(conns_cicero[show_index,3],conns_cor[show_index,3])

CONN=cbind(conns_cor[,3],conns_spearman[,3],conns_cicero[,3],conns_D[,3])





OUT_COR=inferloop.ciceroFrame(indata, used_function=used_function_cor, used_coords=used_coords,genome.df=genome.df)











