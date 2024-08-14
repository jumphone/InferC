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
infercOut=readRDS('./rds/scmul_pbmc_usedType_atacMat_infercOut.rds')
AGG=readRDS('./rds/scmul_pbmc_usedType_atacMat_agg.rds')

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

############################
S2=rep(0,ncol(LLS))
i=1
while(i<=length(S2)){
    this_s=cor(LLS[,i]*P1, LLS[,i]*P2)
    #this_s=inferloop.calABCD(LLS[,i]*P1, LLS[,i]*P2, 0, 0)$D
    S2[i]=this_s
    if(i%%100==1){print(paste0(i,' / ',length(S2)))}
    i=i+1
    }

################
SS1=S1[AGG$clst]
SS2=S2[AGG$clst]

#############
cor(SS1,AGE, method='spearman')
cor(SS2,AGE, method='spearman')



plot(S1,S2)



P1=PHYLOP[LOOP[,1],6]
P2=PHYLOP[LOOP[,2],6]
W=apply(LLS,1,sd)
W=W/max(W)
SP1=P1*W
SP2=P2*W
############################
S3=rep(0,ncol(LLS))
i=1
while(i<=length(S3)){
    this_s=cor(LLS[,i]*SP1, LLS[,i]*SP2)
    #this_s=inferloop.calABCD(LLS[,i]*P1, LLS[,i]*P2, 0, 0)$D
    S3[i]=this_s
    if(i%%100==1){print(paste0(i,' / ',length(S3)))}
    i=i+1
    }

################
SS3=S3[AGG$clst]
cor(SS3,AGE, method='spearman')
















RP1=inferc.sigmoid(P1,a=0,b=1)
RP2=inferc.sigmoid(P2,a=0,b=1)
S3=cor(LLS, abs(RP1-RP2))
SS3=S3[AGG$clst]
cor(SS3,AGE, method='spearman')



R_PHYLOP=PHYLOP
R_PHYLOP[,6]=rank(PHYLOP[,6])
RP1=R_PHYLOP[LOOP[,1],6]
RP2=R_PHYLOP[LOOP[,2],6]

S3=cor(LLS, abs(RP1-RP2))




NS1=scale(S1)[,1]
NS2=scale(S2)[,1]

DS1=NS1-predict(loess(NS1~NS2))
DS2=NS2-predict(loess(NS2~NS1))

DD=sqrt(DS1**2+DS2**2)

CCC=inferc.colMap(DD,c(min(DD),median(DD),max(DD)),c('blue','white','red'))




S3=apply(cbind(scale(S1),scale(S2))
SS3=S3[AGG$clst]
cor(SS3,AGE, method='spearman')





AGE_B=AGE
AGE_B[which(AGE<=2)]=0
AGE_B[which(AGE>2)]=1

cor(SS1,AGE_B, method='spearman')
cor(SS2,AGE_B, method='spearman')




LIN1=which(META$celltype %in% c('RG','IPC','Astrocytes'))
LIN2=which(META$celltype %in% c('RG','IPC','OPC','Oligodendrocytes'))
LIN3=which(META$celltype %in% c('RG','IPC','EN-fetal-early','EN-fetal-late','EN'))
LIN4=which(META$celltype %in% c('RG','IPC','IN-fetal','IN-CGE','IN-MGE'))


SSS=SS1
cor(SSS, AGE, method='spearman')
cor(SSS[LIN1], AGE[LIN1], method='spearman')
cor(SSS[LIN2], AGE[LIN2], method='spearman')
cor(SSS[LIN3], AGE[LIN3], method='spearman')
cor(SSS[LIN4], AGE[LIN4], method='spearman')

SSS=SS2
cor(SSS, AGE, method='spearman')
cor(SSS[LIN1], AGE[LIN1], method='spearman')
cor(SSS[LIN2], AGE[LIN2], method='spearman')
cor(SSS[LIN3], AGE[LIN3], method='spearman')
cor(SSS[LIN4], AGE[LIN4], method='spearman')


SSS=DP
cor(SSS, AGE, method='spearman')
cor(SSS[LIN1], AGE[LIN1], method='spearman')
cor(SSS[LIN2], AGE[LIN2], method='spearman')
cor(SSS[LIN3], AGE[LIN3], method='spearman')
cor(SSS[LIN4], AGE[LIN4], method='spearman')



SSS=SS1
cor(SSS, AGE_B, method='spearman')
cor(SSS[LIN1], AGE_B[LIN1], method='spearman')
cor(SSS[LIN2], AGE_B[LIN2], method='spearman')
cor(SSS[LIN3], AGE_B[LIN3], method='spearman')
cor(SSS[LIN4], AGE_B[LIN4], method='spearman')

SSS=SS2
cor(SSS, AGE_B, method='spearman')
cor(SSS[LIN1], AGE_B[LIN1], method='spearman')
cor(SSS[LIN2], AGE_B[LIN2], method='spearman')
cor(SSS[LIN3], AGE_B[LIN3], method='spearman')
cor(SSS[LIN4], AGE_B[LIN4], method='spearman')


SSS=DP
cor(SSS, AGE_B, method='spearman')
cor(SSS[LIN1], AGE_B[LIN1], method='spearman')
cor(SSS[LIN2], AGE_B[LIN2], method='spearman')
cor(SSS[LIN3], AGE_B[LIN3], method='spearman')
cor(SSS[LIN4], AGE_B[LIN4], method='spearman')








