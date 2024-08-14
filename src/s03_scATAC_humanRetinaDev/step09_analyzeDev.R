library(Seurat)
library(Signac)

library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
annotations=readRDS('/home/database/annotation/hg38/hg38_signac_ucsc_annotations.rds')

source('../../../frontal_cortex/scATAC/InferC.R')

UMAP=readRDS('../rds/pbmc_umap.rds')
BATCH=readRDS('../rds/BATCH.rds')
TIME=readRDS('../rds/TIME.rds')
AGG=readRDS('../rds/AGG.rds')
conns=readRDS('../rds/ciceroFrame_conns_D.rds')
conns[,1]=as.character(conns[,1])
conns[,2]=as.character(conns[,2])

AGE=TIME


COL=rep('NA',nrow(UMAP))
COL[which(BATCH=='d53')]='red3'
COL[which(BATCH=='d59')]='pink'
COL[which(BATCH=='d78')]='lightgoldenrod1'
COL[which(BATCH=='d78.2')]='lightgoldenrod1'
COL[which(BATCH=='d89C')]='springgreen3'
COL[which(BATCH=='d89P')]='skyblue1'


plot(UMAP,col=COL, pch=16,cex=0.3)

pdf('../plot/s9p01_umap.pdf',width=5,height=5)
plot(UMAP,col=COL, pch=16,cex=0.3)
set.seed(123)
O=sample(c(1:nrow(UMAP)),nrow(UMAP))
plot(UMAP[O,],col=COL[O], pch=16,cex=0.3)
dev.off()


################################################
PHYLOP=read.table('../data/GSE184386/ALL_PEAK.bed.hg38.phyloP100way.txt',sep='\t')
rownames(PHYLOP)=PHYLOP[,1]
PHYLOP_=PHYLOP
rownames(PHYLOP_)=stringr::str_replace_all(PHYLOP[,1],'-','_')

PP_=rank(PHYLOP_[,6],ties.method='average')
names(PP_)=rownames(PHYLOP_)

####################
conns_used=conns
R1=rank(conns_used[,3],ties.method='average')
R2=rank(PP_[conns_used[,1]]+PP_[conns_used[,2]],ties.method='average')
RR=R1+R2
this_cut=sort(RR,decreasing=T)[100000*2]

#################
used_index=which(RR>=this_cut)
conns_used=conns[used_index,]
print(dim(conns_used))
conns_uniq=inferloop.getUniqLoop(conns_used)

#################
NET=apply(conns_uniq[,c(1,2)],2,stringr::str_replace_all,'_','-')
PEAK_USED=c(NET[,1],NET[,2])
MAT=AGG$agg[which(rownames(AGG$agg) %in% PEAK_USED),]

###################
infercOut = inferc.calLoopScore(MAT, NET)
LLS=infercOut$lls
GLS=infercOut$gls

####################
LOOP=inferloop.splitLoop(rownames(LLS))
P1=PHYLOP[LOOP[,1],6]
P2=PHYLOP[LOOP[,2],6]

absLLS=abs(LLS)
A=absLLS*P1*P2
B=absLLS*(abs(P1)+abs(P2))
C=colSums(A)/colSums(B)

S1=C
SS1=S1[AGG$clst]
cor(SS1, TIME, method='spearman')
# -0.3455438

boxplot(SS1~TIME)
########################################
START_INDEX=which(TIME==1 & SS1==max(SS1[which(TIME==1)]))[1]
END_INDEX=which(TIME==4 & SS1<quantile(SS1[which(TIME==4)],0.1))[1]


rownames(UMAP)[START_INDEX]
#d53_AACCAACGTCGCTAGC-1
rownames(UMAP)[END_INDEX]
#d89C_AAACGAAAGCAGAAAG-1



BATCH[START_INDEX]
BATCH[END_INDEX]

#####################


pdf('../plot/s9p0x_demoDetail.pdf',width=5,height=5)
set.seed(123)
O=sample(c(1:nrow(UMAP)),nrow(UMAP))
plot(UMAP[O,],col=COL[O], pch=16,cex=0.3)
points(x=UMAP[START_INDEX,1],y=UMAP[START_INDEX,2],cex=2,lwd=2,col='black')
points(x=UMAP[START_INDEX,1],y=UMAP[START_INDEX,2],cex=1.5,pch='x',col='black')
points(x=UMAP[END_INDEX,1],y=UMAP[END_INDEX,2],cex=2,lwd=2,col='black')
points(x=UMAP[END_INDEX,1],y=UMAP[END_INDEX,2],cex=1.5,pch='x',col='black')
#####################################################################
print(SS1[START_INDEX])
#0.5868679
this_index=as.numeric(names(START_INDEX))
X=LLS[,this_index]*P1
Y=LLS[,this_index]*P2
plot(X,Y,pch=16,cex=1, xlim=c(-250,250),ylim=c(-250,250))
abline(h=0,lty=2,lwd=2)
abline(v=0,lty=2,lwd=2)
abline(a=0,b=1,lty=2,lwd=2)
show_index=which(X>100 & Y>100)
rownames(LLS)[show_index]
#"chr1-216705622-216706145.And.chr1-216907028-216907504"
points(X[show_index],Y[show_index],cex=3,lwd=2,col='red')
###############
this_den=density(LLS[show_index,])
plot(this_den,lwd=3,col='grey80',type='h')
points(this_den,lwd=5,col='grey40',type='l')
abline(v=LLS[show_index,this_index],lwd=5,col='red')
LLS[show_index,this_index]
#64.41517
ecdf(LLS[show_index,])(LLS[show_index,this_index])
#1
###############
this_den=density(c(P1,P2))
plot(this_den,lwd=3,col='grey80',type='h')
points(this_den,lwd=5,col='grey40',type='l')
abline(v=P1[show_index],lwd=5,col='red')
abline(v=P2[show_index],lwd=5,col='red')
P1[show_index]
#1.74911
P2[show_index]
#2.84873
ecdf(c(P1,P2))(P1[show_index])
#0.850495
ecdf(c(P1,P2))(P2[show_index])
#0.959545
#################
dev.off()




pdf('../plot/s9p0x_demoDetail_end.pdf',width=5,height=5)
this_col='blue'
print(SS1[END_INDEX])
#0.4038064
this_index=as.numeric(names(END_INDEX))
X=LLS[,this_index]*P1
Y=LLS[,this_index]*P2
plot(X,Y,pch=16,cex=1, xlim=c(-10,10),ylim=c(-10,10))
abline(h=0,lty=2,lwd=2)
abline(v=0,lty=2,lwd=2)
abline(a=0,b=1,lty=2,lwd=2)
show_index=which(X<5 & Y>5)
rownames(LLS)[show_index]
#"chr3-101255222-101255984.And.chr3-101325359-101326192"
points(X[show_index],Y[show_index],cex=3,lwd=2,col=this_col)
###############
this_den=density(LLS[show_index,])
plot(this_den,lwd=3,col='grey80',type='h')
points(this_den,lwd=5,col='grey40',type='l')
abline(v=LLS[show_index,this_index],lwd=5,col=this_col)
LLS[show_index,this_index]
#4.81307
ecdf(LLS[show_index,])(LLS[show_index,this_index])
#0.988
###############
###############
this_den=density(c(P1,P2))
plot(this_den,lwd=3,col='grey80',type='h')
points(this_den,lwd=5,col='grey40',type='l')
abline(v=P1[show_index],lwd=5,col='blue')
abline(v=P2[show_index],lwd=5,col='blue')
P1[show_index]
#0.232545
P2[show_index]
#1.64412
ecdf(c(P1,P2))(P1[show_index])
#0.008825
ecdf(c(P1,P2))(P2[show_index])
#0.82903
#################

#################
dev.off()







##################################################




X=LLS[,as.numeric(names(END_INDEX))]*P1
Y=LLS[,as.numeric(names(END_INDEX))]*P2
plot(X,Y,pch=16,cex=1, xlim=c(-20,20),ylim=c(-20,20))
abline(h=0,lty=2,lwd=2)
abline(v=0,lty=2,lwd=2)













#################



########################################

CUMAP=t(.generate_mean(t(UMAP),AGG$clst))
CUMAP=CUMAP[order(as.numeric(rownames(CUMAP))),]
DIST=as.matrix(dist(CUMAP))
RMAT=apply(DIST,2,rank)

N=nrow(CUMAP)
Q=0.99
S2=( S1 %*% ( Q**RMAT ) ) * (1-Q)/Q
SS2=S2[AGG$clst]
cor(SS2, TIME, method='spearman')
# -0.5022852


###########################################
CUMAP=t(.generate_mean(t(UMAP),AGG$clst))
CUMAP=CUMAP[order(as.numeric(rownames(CUMAP))),]
VVV=S1
CCC=inferc.colMap(VVV,c(min(VVV),median(VVV),max(VVV)),c('blue','grey90','red'))



pdf('../plot/s9p00_trajectory_thisStudy.pdf',width=4,height=4.5)

plot(UMAP,col=CCC[AGG$clst],pch=16,cex=0.3)

plot(x=VVV,y=rep(1,length(VVV)),col=CCC,type='h',lwd=2)

source('/home/toolkit/src/fitdevo.R')
FIELD=fitdevo.field(DP=SS1, VEC=UMAP, SHOW=TRUE,N=11,P=0.9,AL=0.5,CUT=0,CEX=0.3,COL=CCC[AGG$clst])

dev.off()



META=readRDS('../rds/META.rds')

cor(META$nCount_RNA,TIME, method='spearman')
#-0.1115383

META$nFeature_RNA
cor(META$nFeature_RNA,TIME, method='spearman')
#-0.1610665


pdf('../plot/s9p0x_time_N.pdf',width=4,height=4)
boxplot(META$nFeature_RNA~TIME,pch='+',cex=0.5,col='grey60',main='nFeature')
boxplot(META$nCount_RNA~TIME,pch='+',cex=0.5,col='grey60',main='nCount')
dev.off()


tapply(META$nFeature_RNA,TIME,mean)
#1        2        3        4
#7224.699 6201.750 3713.469 4031.564

tapply(META$nCount_RNA,TIME,mean)
#       1        2        3        4
#19880.60 17162.45 13125.73 14533.35



GA=readRDS('../rds/GeneActivity.rds')
ES=read.table('./GS/BHATTACHARYA_EMBRYONIC_STEM_CELL.gmt',sep='\t')[,1]
GA_ES=GA[which(rownames(GA) %in% ES),]

SS5=apply(GA_ES,2,mean)
cor(SS5,TIME, method='spearman')
#-0.3156221

USED=which(TIME %in% c(3,4))
cor(SS4[USED],TIME[USED],method='spearman')
cor(SS3[USED],TIME[USED],method='spearman')


saveRDS(SS1,'../rds/SS1.rds')
saveRDS(SS5,'../rds/SS5.rds')
###################################


#########

############################

############################
pdf('../plot/s9p0x_marker.pdf',width=3,height=3.5)
CEX=0.5
PCH=16

GENE='PAX6'
TMP=GA[GENE,]
TMP_COL=inferc.colMap(TMP,c(0,1.5,3),c('blue','grey90','red'))
plot(UMAP[order(TMP),], col=TMP_COL[order(TMP)], cex=CEX, pch=PCH,main=GENE)
plot(x=TMP,y=rep(1,length(TMP)),type='h',col=TMP_COL,lwd=2)

GENE='NES'
TMP=GA[GENE,]
TMP_COL=inferc.colMap(TMP,c(0,1.4,3),c('blue','grey90','red'))
plot(UMAP[order(TMP),], col=TMP_COL[order(TMP)], cex=CEX, pch=PCH,main=GENE)
plot(x=TMP,y=rep(1,length(TMP)),type='h',col=TMP_COL,lwd=2)

GENE='SOX2'
TMP=GA[GENE,]
TMP_COL=inferc.colMap(TMP,c(0,1,3),c('blue','grey90','red'))
plot(UMAP[order(TMP),], col=TMP_COL[order(TMP)], cex=CEX, pch=PCH,main=GENE)
plot(x=TMP,y=rep(1,length(TMP)),type='h',col=TMP_COL,lwd=2)

GENE='NOTCH1'
TMP=GA[GENE,]
TMP_COL=inferc.colMap(TMP,c(0,1.9,3),c('blue','grey90','red'))
plot(UMAP[order(TMP),], col=TMP_COL[order(TMP)], cex=CEX, pch=PCH,main=GENE)
plot(x=TMP,y=rep(1,length(TMP)),type='h',col=TMP_COL,lwd=2)


GENE='ABCG2'
TMP=GA[GENE,]
TMP_COL=inferc.colMap(TMP,c(0,1.9,3),c('blue','grey90','red'))
plot(UMAP[order(TMP),], col=TMP_COL[order(TMP)], cex=CEX, pch=PCH,main=GENE)
plot(x=TMP,y=rep(1,length(TMP)),type='h',col=TMP_COL,lwd=2)

GENE='POU5F1'
TMP=GA[GENE,]
TMP_COL=inferc.colMap(TMP,c(0,1,3),c('blue','grey90','red'))
plot(UMAP[order(TMP),], col=TMP_COL[order(TMP)], cex=CEX, pch=PCH,main=GENE)
plot(x=TMP,y=rep(1,length(TMP)),type='h',col=TMP_COL,lwd=2)


dev.off()







#plot(UMAP, col=TMP_COL, cex=0.3, pch=16)




#####################################

X=TIME

Y=SS1

pdf('../plot/s9p02_time_thisStudy.pdf',width=4,height=4)
boxplot(Y~X,pch='+',cex=0.5,col='grey60')
dev.off()

tapply(Y,X,mean)
#  1         2         3         4
#0.5069504 0.4998529 0.4861146 0.4760211


X=TIME

Y=SS5

pdf('../plot/s9p03_time_ESC.pdf',width=4,height=4)
boxplot(Y~X,pch='+',cex=0.5,col='grey60')
dev.off()

tapply(Y,X,mean)
#  1         2         3         4
#0.18745847 0.14431457 0.07116677 0.10433890



#####################################
















