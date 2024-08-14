
###############################
library(Seurat)
library(Signac)
source('/home/toolkit/src/fitdevo.R')
source('../../frontal_cortex/scATAC/InferC.R')

CCI=read.table('./rds/sctc_all.cci',sep='\t',header=F,row.names=NULL)[,1]

META=readRDS('./rds/scmul_pbmc_usedType_meta.rds')
UMAP=readRDS('./rds/scmul_pbmc_usedType_umap.rds')
PCA=readRDS('./rds/scmul_pbmc_usedType_pca.rds')
LSI=readRDS('./rds/scmul_pbmc_usedType_lsi.rds')
DP=readRDS('./rds/scmul_pbmc_usedType_dp.rds')
AGE=readRDS('./rds/scmul_pbmc_usedType_age.rds')


cor(AGE,CCI,method='spearman')
#-0.3459126
cor(AGE,DP,method='spearman')
#-0.4158901

##########################
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



Y=CCI

set.seed(123)
O=sample(c(1:length(X)),length(X))

SS=smooth.spline(x=X,y=Y,df=10)

pdf('./plot/s2p02_age_sctc_cci.pdf',width=15,height=5)
plot(X[O],Y[O],col=COL[O],cex=0.5,pch=16)
abline(h=mean(Y),lwd=1,lty=1,col='black')
points(SS$x,SS$y,pch=16,cex=1,col='black',type='l',lty=2,lwd=4)
dev.off()


Y=DP

set.seed(123)
O=sample(c(1:length(X)),length(X))

SS=smooth.spline(x=X,y=Y,df=10)

pdf('./plot/s2p02_age_fitdevo_dp.pdf',width=15,height=5)
plot(X[O],Y[O],col=COL[O],cex=0.5,pch=16)
abline(h=mean(Y),lwd=1,lty=1,col='black')
points(SS$x,SS$y,pch=16,cex=1,col='black',type='l',lty=2,lwd=4)
dev.off()






