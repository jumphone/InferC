
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


###################################
P1=predict(lm(AGE~PCA[,1:50]+LSI[,2:40]))

NP1=P1
i=1
while(i<=max(AGE)){
   this_index=which(AGE==i)
   this_value=rank(P1[this_index])
   this_new=.norm1(this_value)
   NP1[this_index]=this_new
   i=i+1}



X=AGE+NP1*0.9
################################################


COL=rep('grey',nrow(META))
COL[which(META$celltype=='Astrocytes')]='pink'
COL[which(META$celltype=='EN-fetal-early')]='springgreen1'
COL[which(META$celltype=='EN-fetal-late')]='springgreen3'
COL[which(META$celltype=='EN')]='springgreen4'
COL[which(META$celltype=='Endothelial')]='grey60'
COL[which(META$celltype=='IN-fetal')]='lightgoldenrod1'
COL[which(META$celltype=='IN-CGE')]='lightgoldenrod3'
COL[which(META$celltype=='IN-MGE')]='lightgoldenrod4'
COL[which(META$celltype=='IPC')]='indianred1'
COL[which(META$celltype=='Microglia')]='violet'
COL[which(META$celltype=='OPC')]='skyblue1'
COL[which(META$celltype=='Oligodendrocytes')]='royalblue3'
COL[which(META$celltype=='Pericytes')]='sienna'
COL[which(META$celltype=='RG')]='red3'
COL[which(META$celltype=='VSMC')]='ivory3'
#####################################################


atac_mat=readRDS('./rds/scmul_pbmc_usedType_atacMat.rds')

set.seed(123)
random_index=sample(1:ncol(atac_mat),5000)
atac_mat_sub=atac_mat[,random_index]
UMAP_sub=UMAP[random_index,]

used_coords=UMAP_sub
indata=atac_mat_sub

genome.df=inferloop.getGenomeDF.hg38()
window=500000
sample_num=100

cicero_cds=inferc.ciceroFrame_step01_prepareInput(indata, used_coords=used_coords,genome.df=genome.df,window=window,sample_num=sample_num)

saveRDS(cicero_cds,'./rds/ciceroFrame_cicero_cds.rds')



source('../../frontal_cortex/scATAC/InferC.R')
cicero_cds=readRDS('./rds/ciceroFrame_cicero_cds.rds')


used_function=inferc.rowCalD
conns <- inferc.ciceroFrame_step02_runUsedFuntion(cicero_cds, genome.df, window, sample_num, used_function)
saveRDS(conns,'./rds/ciceroFrame_conns_D.rds')


#####################################################




source('../../frontal_cortex/scATAC/InferC.R')

conns=readRDS('./rds/ciceroFrame_conns_D.rds')

used_index=which(conns[,3]>0.5)
#used_index=which(conns[,3]>0.3)
length(used_index)

source('../../frontal_cortex/scATAC/InferC.R')
conns_uniq=inferloop.getUniqLoop(conns[used_index,])

saveRDS(conns_uniq,'./rds/ciceroFrame_conns_D_uniq.rds')
conns_uniq=readRDS('./rds/ciceroFrame_conns_D_uniq.rds')
#saveRDS(conns_uniq,'./rds/ciceroFrame_conns_D_uniq_0.3.rds')
#conns_uniq=readRDS('./rds/ciceroFrame_conns_D_uniq_0.3.rds')


AGG=inferc.aggMat(atac_mat, UMAP, clstNum=1000,seed=123)

saveRDS(AGG, file='./rds/scmul_pbmc_usedType_atacMat_agg.rds')


conns_uniq=readRDS('./rds/ciceroFrame_conns_D_uniq.rds')
AGG=readRDS('./rds/scmul_pbmc_usedType_atacMat_agg.rds')
NET=apply(conns_uniq[,c(1,2)],2,stringr::str_replace_all,'_','-')

MAT=AGG$agg

infercOut = inferc.calLoopScore(MAT, NET)

saveRDS(infercOut, './rds/scmul_pbmc_usedType_atacMat_infercOut.rds')


###############################
infercOut=readRDS('./rds/scmul_pbmc_usedType_atacMat_infercOut.rds')
AGG=readRDS('./rds/scmul_pbmc_usedType_atacMat_agg.rds')

LLS=infercOut$lls

LOOP=inferloop.splitLoop(rownames(LLS))
BED1=inferloop.splitLoop(LOOP[,1],'-',3)
BED2=inferloop.splitLoop(LOOP[,2],'-',3)
LLL=abs(as.numeric(BED1[,2])-as.numeric(BED2[,2]))




BED_ALL=inferloop.splitLoop(rownames(atac_mat),'-',3)
write.table(BED_ALL,'./rds/scmul_pbmc_usedType_atacMat.rds.bedAll.bed',quote=F,sep='\t',row.names=F,col.names=F)




























