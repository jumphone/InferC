library(Seurat)
library(Signac)

library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)


source('../../../frontal_cortex/scATAC/InferC.R')

pbmc=readRDS('../rds/pbmc.rds')

BATCH=pbmc$batch
TIME=rep(0,length(BATCH))
TIME[which(BATCH %in% c('hm1','hm2'))]=1

UMAP=readRDS('../rds/pbmc_umap.rds')
AGG=readRDS('../rds/AGG.rds')
conns=readRDS('../rds/ciceroFrame_conns_D.rds')
conns[,1]=as.character(conns[,1])
conns[,2]=as.character(conns[,2])


################################################
PHYLOP=read.table('../data/GSE159557_RAW/ALL_PEAK.bed.hg19.phyloP100way.txt',sep='\t')
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
#0.05802713

CUMAP=t(.generate_mean(t(UMAP),AGG$clst))
CUMAP=CUMAP[order(as.numeric(rownames(CUMAP))),]
DIST=as.matrix(dist(CUMAP))
RMAT=apply(DIST,2,rank)

N=nrow(CUMAP)
Q=0.99
S2=( S1 %*% ( Q**RMAT ) ) * (1-Q)/Q
SS2=S2[AGG$clst]
cor(SS2, TIME, method='spearman')
#0.08440593

cor(pbmc$nCount_peaks,TIME,method='spearman')
#-0.02332737
cor(pbmc$nFeature_peaks,TIME,method='spearman')
#-0.03062504

pdf('../plot/liver_s2p01_boxplot.pdf',width=2,height=5.5)
boxplot(SS1[which(TIME==1)],SS1[which(TIME==0)],col=c('indianred1','grey60'),pch='+',cex=0.5)
dev.off()

t.test(SS1[which(TIME==1)],SS1[which(TIME==0)])
#0.01497



