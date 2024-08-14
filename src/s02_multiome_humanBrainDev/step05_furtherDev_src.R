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
PTMP=predict(lm(AGE~PCA[,1:50]+LSI[,2:40]))

NPTMP=PTMP
i=1
while(i<=max(AGE)){
   this_index=which(AGE==i)
   this_value=rank(PTMP[this_index])
   this_new=.norm1(this_value)
   NPTMP[this_index]=this_new
   i=i+1}
AGE_X=AGE+NPTMP*0.9
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

print('load data...')
atac_mat=readRDS('./rds/scmul_pbmc_usedType_atacMat.rds')
infercOut=readRDS('./rds/scmul_pbmc_usedType_atacMat_infercOut.rds')
AGG=readRDS('./rds/scmul_pbmc_usedType_atacMat_agg.rds')
print('finished!')


LLS=infercOut$lls
GLS=infercOut$gls

###################################
LOOP=inferloop.splitLoop(rownames(LLS))
BED1=inferloop.splitLoop(LOOP[,1],'-',3)
BED2=inferloop.splitLoop(LOOP[,2],'-',3)
LEN= abs((as.numeric(BED1[,2])+as.numeric(BED1[,3]))/2 - (as.numeric(BED2[,2])+as.numeric(BED2[,3]))/2)
summary(LEN)
###################################
###################################

PHYLOP=read.table('rds/scmul_pbmc_usedType_atacMat.rds.bedAll.bed.hg38.phyloP100way.txt',sep='\t')
rownames(PHYLOP)=PHYLOP[,1]

####################################################
S1=rep(0,ncol(LLS))
i=1
while(i<=length(S1)){
    this_s=cor(AGG$agg[,i], PHYLOP[,6])
    #this_s=inferloop.calABCD(AGG$agg[,i], PHYLOP[,6], 0, 0)$D
    S1[i]=this_s
    if(i%%100==1){print(paste0(i,' / ',length(S1)))}
    i=i+1
    }

####################################################
P1=PHYLOP[LOOP[,1],6]
P2=PHYLOP[LOOP[,2],6]






