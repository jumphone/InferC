


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

LLS=readRDS('./rds/LLS.rds')
P1=readRDS('./rds/P1.rds')
P2=readRDS('./rds/P2.rds')

AGG=readRDS('./rds/scmul_pbmc_usedType_atacMat_agg.rds')



MLLS=.generate_mean(LLS[,AGG$clst],AGE)
MLLS=MLLS[,order(as.numeric(colnames(MLLS)))]

pdf('./plot/s8p0x_demoDetail_agg.pdf',width=5,height=5)

plot(MLLS[,1]*P1,MLLS[,1]*P2,xlim=c(-6,6),ylim=c(-6,6),pch=16,cex=1)
abline(h=0,lty=2,lwd=2)
abline(v=0,lty=2,lwd=2)
abline(a=0,b=1,lty=2,lwd=2)


plot(MLLS[,6]*P1,MLLS[,6]*P2,xlim=c(-1.5,1.5),ylim=c(-1.5,1.5),pch=16,cex=1)
abline(h=0,lty=2,lwd=2)
abline(v=0,lty=2,lwd=2)
abline(a=0,b=1,lty=2,lwd=2)
dev.off()




absLLS=abs(MLLS)
A=absLLS*P1*P2
B=absLLS*(abs(P1)+abs(P2))
C=colSums(A)/colSums(B)
print(C)
#0.6022621 0.5354956 0.5031634 0.4859833 0.4843726 0.4813491







LIN1=which(META$celltype %in% c('RG','IPC','Astrocytes'))
LIN2=which(META$celltype %in% c('RG','IPC','OPC','Oligodendrocytes'))
LIN3=which(META$celltype %in% c('RG','IPC','EN-fetal-early','EN-fetal-late','EN'))
LIN4=which(META$celltype %in% c('RG','IPC','IN-fetal','IN-CGE','IN-MGE'))


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


SS1=readRDS('./rds/SS1.rds')


plot(UMAP, col=COL, pch=16,cex=0.3)

.norm1=function(x){
    if(min(x)!=max(x)){
        y=(x-min(x))/(max(x)-min(x))
       }else{
        y=rep(0,length(x))
       }
    return(y)
    }


alpha=0.2
beta=0.35
theta=0.05

NSS1=.norm1(SS1)
X0=.norm1(UMAP[,1])
Y0=.norm1(UMAP[,2])

X1=X0+Y0
Y1=Y0*alpha
Y2=NSS1+beta+theta
XLIM=c(0,2)
YLIM=c(0,1+beta+theta)


pdf('./plot/s8p0x_demoDetail_umap.pdf',width=5,height=6.5)
plot(X1,Y1,xlim=XLIM,ylim=YLIM,col=COL,pch=16,cex=0.2)
segments(x0=0,x1=1,y0=0,y1=alpha)
segments(x0=1,x1=2,y0=alpha,y1=alpha)
segments(x0=0,x1=1,y0=0,y1=0)
segments(x0=1,x1=2,y0=0,y1=alpha)

points(X1,Y2,pch=16,col=COL,cex=0.5)
abline(h=beta,lty=2,lwd=1,col='grey60')

######################################################
START_INDEX=which(AGE==1 & SS1==max(SS1[which(AGE==1)]))[1]
END_INDEX=which(AGE==6 & SS1<quantile(SS1[which(AGE==6)],0.1))[1]
saveRDS(START_INDEX,file='START_INDEX.rds')
saveRDS(END_INDEX,file='END_INDEX.rds')

points(x=X1[START_INDEX],y=Y2[START_INDEX],cex=2,lwd=2,col='black')
points(x=X1[START_INDEX],y=Y2[START_INDEX],cex=1.5,pch='x',col='black')
points(x=X1[END_INDEX],y=Y2[END_INDEX],cex=2,lwd=2,col='black')
points(x=X1[END_INDEX],y=Y2[END_INDEX],cex=1.5,pch='x',col='black')
dev.off()

rownames(META)[START_INDEX]
#11_ACCGCAATCTTGCAGG-1
rownames(META)[END_INDEX]
#150656_AAACCGCGTTTACTTG-1

META$celltype[START_INDEX]
#RG

META$celltype[END_INDEX]
#EN

START_INDEX=readRDS(file='START_INDEX.rds')
END_INDEX=readRDS(file='END_INDEX.rds')


pdf('./plot/s8p0x_demoDetail_start.pdf',width=5,height=5)
#####################################################################
print(SS1[START_INDEX])
#0.6638889
this_index=as.numeric(names(START_INDEX))
X=LLS[,this_index]*P1
Y=LLS[,this_index]*P2
plot(X,Y,pch=16,cex=1, xlim=c(-50,50),ylim=c(-50,50))
abline(h=0,lty=2,lwd=2)
abline(v=0,lty=2,lwd=2)
abline(a=0,b=1,lty=2,lwd=2)
show_index=which(X>30 & Y>30)
rownames(LLS)[show_index]
#"chr3-27980390-27980726.And.chr3-27992697-27993868"
points(X[show_index],Y[show_index],cex=3,lwd=2,col='red')
###############
this_den=density(LLS[show_index,])
plot(this_den,lwd=3,col='grey80',type='h')
points(this_den,lwd=5,col='grey40',type='l')
abline(v=LLS[show_index,this_index],lwd=5,col='red')
LLS[show_index,this_index]
#19.73779
ecdf(LLS[show_index,])(LLS[show_index,this_index])
#0.999
###############
this_den=density(c(P1,P2))
plot(this_den,lwd=3,col='grey80',type='h')
points(this_den,lwd=5,col='grey40',type='l')
abline(v=P1[show_index],lwd=2,col='red')
abline(v=P2[show_index],lwd=2,col='red')
P1[show_index]
#2.35228
P2[show_index]
#2.28381
ecdf(c(P1,P2))(P1[show_index])
#0.93054
ecdf(c(P1,P2))(P2[show_index])
#0.923385
#################
dev.off()



pdf('./plot/s8p0x_demoDetail_end.pdf',width=5,height=5)
#####################################################################
print(SS1[END_INDEX])
#0.4649822
this_index=as.numeric(names(END_INDEX))
X=LLS[,this_index]*P1
Y=LLS[,this_index]*P2
plot(X,Y,pch=16,cex=1, xlim=c(-20,20),ylim=c(-20,20))
abline(h=0,lty=2,lwd=2)
abline(v=0,lty=2,lwd=2)
abline(a=0,b=1,lty=2,lwd=2)
show_index=which(X>15 & Y<10)
rownames(LLS)[show_index]
#"chr15-63612246-63612605.And.chr15-63869807-63870969"
points(X[show_index],Y[show_index],cex=3,lwd=2,col='blue')
###############
this_den=density(LLS[show_index,])
plot(this_den,lwd=3,col='grey80',type='h')
points(this_den,lwd=5,col='grey40',type='l')
abline(v=LLS[show_index,this_index],lwd=5,col='blue')
LLS[show_index,this_index]
#3.931666
ecdf(LLS[show_index,])(LLS[show_index,this_index])
#0.995
###############
this_den=density(c(P1,P2))
plot(this_den,lwd=3,col='grey80',type='h')
points(this_den,lwd=5,col='grey40',type='l')
abline(v=P1[show_index],lwd=5,col='blue')
abline(v=P2[show_index],lwd=5,col='blue')
P1[show_index]
#4.66371
P2[show_index]
#0.202695
ecdf(c(P1,P2))(P1[show_index])
#0.995135
ecdf(c(P1,P2))(P2[show_index])
#0.00497
#################
dev.off()


















#####################################################################
print(SS1[START_INDEX])




































