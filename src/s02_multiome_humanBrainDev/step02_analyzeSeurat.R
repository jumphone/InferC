

library(Seurat)
library(Signac)

pbmc=readRDS('./rds/scmul_pbmc.rds')


TYPE=pbmc$celltype
USED_TYPE=c('Astrocytes','EN','EN-fetal-early','EN-fetal-late',
            'IN-CGE','IN-fetal','IN-MGE','IPC',
            'Oligodendrocytes','OPC','RG')

pbmc_sub=subset(pbmc, cells=colnames(pbmc)[which(TYPE %in% USED_TYPE)])

META_ALL=pbmc@meta.data
UMAP_ALL=pbmc@reductions$umap@cell.embeddings

saveRDS(META_ALL, './rds/scmul_pbmc_meta.rds')
saveRDS(UMAP_ALL, './rds/scmul_pbmc_umap.rds')

###################################################
META=pbmc_sub@meta.data
PCA=pbmc_sub@reductions$pca@cell.embeddings
LSI=pbmc_sub@reductions$lsi@cell.embeddings
UMAP=pbmc_sub@reductions$umap@cell.embeddings


saveRDS(META, './rds/scmul_pbmc_usedType_meta.rds')
saveRDS(UMAP, './rds/scmul_pbmc_usedType_umap.rds')
saveRDS(PCA, './rds/scmul_pbmc_usedType_pca.rds')
saveRDS(LSI, './rds/scmul_pbmc_usedType_lsi.rds')

##################################

library(Seurat)
library(Signac)

pbmc=readRDS('./rds/scmul_pbmc.rds')
TYPE=pbmc$celltype
USED_TYPE=c('Astrocytes','EN','EN-fetal-early','EN-fetal-late',
            'IN-CGE','IN-fetal','IN-MGE','IPC',
            'Oligodendrocytes','OPC','RG')

pbmc_sub=subset(pbmc, cells=colnames(pbmc)[which(TYPE %in% USED_TYPE)])
META=pbmc_sub@meta.data
PCA=pbmc_sub@reductions$pca@cell.embeddings
LSI=pbmc_sub@reductions$lsi@cell.embeddings
UMAP=pbmc_sub@reductions$umap@cell.embeddings

#data <- readRDS("./rdsPancreas_10x_downsampled.rds")
# extract expression data
#expression_data <- data$expression_data
# running CytoTRACE 2 main function - cytotrace2 - with default parameters
#cytotrace2_result <- CytoTRACE2::cytotrace2(expression_data)


CT2=CytoTRACE2::cytotrace2(pbmc_sub[['RNA']]@data,species = 'human')
saveRDS(CT2,'./rds/scmul_pbmc_usedType_ct2.rds')

CT2_count=CytoTRACE2::cytotrace2(pbmc_sub[['RNA']]@counts,species = 'human')
saveRDS(CT2_count,'./rds/scmul_pbmc_usedType_ct2_count.rds')

source('/home/toolkit/src/fitdevo.R')
BGW=readRDS('/home/toolkit/src/BGW.rds')
MAT=pbmc_sub[['RNA']]@counts
DP=fitdevo(MAT=MAT, BGW=BGW, NORM=TRUE)

saveRDS(DP,'./rds/scmul_pbmc_usedType_dp.rds')

AGE=rep(0,length(META$donor_id))
AGE[which(META$donor_id %in% c('EaFet1','EaFet2'))]=1
AGE[which(META$donor_id %in% c('LaFet1','LaFet2'))]=2
AGE[which(META$donor_id %in% c('Inf1','Inf2'))]=3
AGE[which(META$donor_id %in% c('Child1','Child2'))]=4
AGE[which(META$donor_id %in% c('Adol1','Adol2'))]=5
AGE[which(META$donor_id %in% c('Adult1','Adult2'))]=6

saveRDS(AGE,'./rds/scmul_pbmc_usedType_age.rds')

#######################################################################
CT2=readRDS('./rds/scmul_pbmc_usedType_ct2.rds')
CT2_count=readRDS('./rds/scmul_pbmc_usedType_ct2_count.rds')
DP=readRDS('./rds/scmul_pbmc_usedType_dp.rds')
AGE=readRDS('./rds/scmul_pbmc_usedType_age.rds')

boxplot(CT2$CytoTRACE2_Score~AGE)

cor(CT2$CytoTRACE2_Score,AGE,method='spearman')
cor(CT2_count$CytoTRACE2_Score,AGE,method='spearman')
cor(DP,AGE,method='spearman')



#############################################################################


library(Seurat)
library(Signac)
source('/home/toolkit/src/fitdevo.R')
source('../../frontal_cortex/scATAC/InferC.R')

META_ALL=readRDS('./rds/scmul_pbmc_meta.rds')
UMAP_ALL=readRDS('./rds/scmul_pbmc_umap.rds')

COL_ALL=rep('grey',nrow(META_ALL))
COL_ALL[which(META_ALL$celltype=='Astrocytes')]='pink'
COL_ALL[which(META_ALL$celltype=='EN-fetal-early')]='springgreen1'
COL_ALL[which(META_ALL$celltype=='EN-fetal-late')]='springgreen3'
COL_ALL[which(META_ALL$celltype=='EN')]='springgreen4'
COL_ALL[which(META_ALL$celltype=='Endothelial')]='grey60'
COL_ALL[which(META_ALL$celltype=='IN-fetal')]='lightgoldenrod1'
COL_ALL[which(META_ALL$celltype=='IN-CGE')]='lightgoldenrod3'
COL_ALL[which(META_ALL$celltype=='IN-MGE')]='lightgoldenrod4'
COL_ALL[which(META_ALL$celltype=='IPC')]='indianred1'
COL_ALL[which(META_ALL$celltype=='Microglia')]='violet'
COL_ALL[which(META_ALL$celltype=='OPC')]='skyblue1'
COL_ALL[which(META_ALL$celltype=='Oligodendrocytes')]='royalblue3'
COL_ALL[which(META_ALL$celltype=='Pericytes')]='sienna'
COL_ALL[which(META_ALL$celltype=='RG')]='red3'
COL_ALL[which(META_ALL$celltype=='VSMC')]='ivory3'

plot(UMAP_ALL, col=COL_ALL, pch='+',cex=0.3)

TYPE_ALL=META_ALL$celltype
UMAP_ALL_TYPE=t(.generate_mean(t(UMAP_ALL), TYPE_ALL))


pdf('./plot/s2p01_umap.pdf',width=5,height=5)
plot(UMAP_ALL, col=COL_ALL, pch=16,cex=0.3)
text(x=UMAP_ALL_TYPE[,1],y=UMAP_ALL_TYPE[,2],labels=rownames(UMAP_ALL_TYPE))

plot(UMAP_ALL, col=COL_ALL, pch=16,cex=0.3)
dev.off()








META=readRDS('./rds/scmul_pbmc_usedType_meta.rds')
UMAP=readRDS('./rds/scmul_pbmc_usedType_umap.rds')
PCA=readRDS('./rds/scmul_pbmc_usedType_pca.rds')
LSI=readRDS('./rds/scmul_pbmc_usedType_lsi.rds')
DP=readRDS('./rds/scmul_pbmc_usedType_dp.rds')
AGE=readRDS('./rds/scmul_pbmc_usedType_age.rds')

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
Y=DP

set.seed(123)
O=sample(c(1:length(X)),length(X))



pdf('./plot/s2p02_age_fitdevo_dp.pdf',width=15,height=5)
fit=lm(Y~X)
plot(X[O],Y[O],col=COL[O],cex=0.5,pch=16)
abline(fit,lwd=3,lty=2,col='black')
abline(h=mean(Y),lwd=1,lty=1,col='black')
dev.off()




USED=which(META$celltype %in% c('RG','IPC','OPC','Oligodendrocytes'))

plot(X[O],Y[O],col='grey70',cex=0.3,pch=16)
points(X[USED],Y[USED],col=COL[USED],cex=0.5,pch=16)



source('../../frontal_cortex/scATAC/InferC.R')

atac_mat=pbmc_sub[['ATAC']]@data
saveRDS(atac_mat, './rds/scmul_pbmc_usedType_atacMat.rds')



###############################












































