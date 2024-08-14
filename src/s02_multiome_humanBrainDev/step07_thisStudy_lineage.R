
setwd('/home/database/data/InferX/data/humanBrainDev/scMulti')

source('step05_furtherDev_src.R')

PHYLOP=read.table('rds/scmul_pbmc_usedType_atacMat.rds.bedAll.bed.hg38.phyloP100way.txt',sep='\t')
rownames(PHYLOP)=PHYLOP[,1]

conns=readRDS('./rds/ciceroFrame_conns_D.rds')
conns[,1]=as.character(conns[,1])
conns[,2]=as.character(conns[,2])

rownames(PHYLOP)=PHYLOP[,1]
PHYLOP_=PHYLOP
rownames(PHYLOP_)=stringr::str_replace_all(PHYLOP[,1],'-','_')

PP_=rank(PHYLOP_[,6],ties.method='average')
names(PP_)=rownames(PHYLOP_)

####################
R1=rank(conns[,3],ties.method='average')
R2=rank(PP_[conns[,1]]+PP_[conns[,2]],ties.method='average')
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

saveRDS(LLS,'./rds/LLS.rds')

####################
LOOP=inferloop.splitLoop(rownames(LLS))
P1=PHYLOP[LOOP[,1],6]
P2=PHYLOP[LOOP[,2],6]


saveRDS(P1,'./rds/P1.rds')
saveRDS(P2,'./rds/P2.rds')
#####################


absLLS=abs(LLS)
A=absLLS*P1*P2
B=absLLS*(abs(P1)+abs(P2))
C=colSums(A)/colSums(B)

S1=C
SS1=S1[AGG$clst]
cor(SS1, AGE, method='spearman')
########################



CUMAP=t(.generate_mean(t(UMAP),AGG$clst))
CUMAP=CUMAP[order(as.numeric(rownames(CUMAP))),]
DIST=as.matrix(dist(CUMAP))
RMAT=apply(DIST,2,rank)

N=nrow(CUMAP)
Q=0.99
S2=( S1 %*% ( Q**RMAT ) ) * (1-Q)/Q
SS2=S2[AGG$clst]
cor(SS2, AGE, method='spearman')

############################################

saveRDS(SS1,'./rds/SS1.rds')
saveRDS(SS2,'./rds/SS2.rds')



CUMAP=t(.generate_mean(t(UMAP),AGG$clst))
CUMAP=CUMAP[order(as.numeric(rownames(CUMAP))),]
V=S1
CCC=inferc.colMap(V,c(min(V),median(V),max(V)),c('blue','grey90','red'))
plot(CUMAP,col=CCC,pch=16)

CUMAP=t(.generate_mean(t(UMAP),AGG$clst))
CUMAP=CUMAP[order(as.numeric(rownames(CUMAP))),]
V=S2
CCC=inferc.colMap(V,c(min(V),median(V),max(V)),c('blue','grey90','red'))
plot(CUMAP,col=CCC,pch=16)




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




Y=SS1

set.seed(123)
O=sample(c(1:length(X)),length(X))

SS=smooth.spline(x=X,y=Y,df=3)

pdf('./plot/s2p02_age_thisStudy.pdf',width=15,height=5)
plot(X[O],Y[O],col=COL[O],cex=0.5,pch=16)
abline(h=mean(Y),lwd=1,lty=1,col='black')
points(SS$x,SS$y,pch=16,cex=1,col='black',type='l',lty=2,lwd=4)
dev.off()




####################################

LIN1=which(META$celltype %in% c('RG','IPC','Astrocytes'))
LIN2=which(META$celltype %in% c('RG','IPC','OPC','Oligodendrocytes'))
LIN3=which(META$celltype %in% c('RG','IPC','EN-fetal-early','EN-fetal-late','EN'))
LIN4=which(META$celltype %in% c('RG','IPC','IN-fetal','IN-CGE','IN-MGE'))


pdf('./plot/s2p03_age_thisStudy_lineages.pdf',width=5,height=5)
###########
this_x=X[LIN1]
this_y=Y[LIN1]
this_col=COL[LIN1]
this_ss=smooth.spline(x=this_x,y=this_y,df=10)
plot(this_x,this_y,col=this_col,cex=0.3,pch=16)
#abline(h=mean(this_y),lwd=1,lty=1,col='black')
points(this_ss$x,this_ss$y,pch=16,cex=1,col='black',type='l',lty=2,lwd=3)
###########
this_x=X[LIN2]
this_y=Y[LIN2]
this_col=COL[LIN2]
this_ss=smooth.spline(x=this_x,y=this_y,df=10)
plot(this_x,this_y,col=this_col,cex=0.3,pch=16)
#abline(h=mean(this_y),lwd=1,lty=1,col='black')
points(this_ss$x,this_ss$y,pch=16,cex=1,col='black',type='l',lty=2,lwd=3)
###########
this_x=X[LIN3]
this_y=Y[LIN3]
this_col=COL[LIN3]
this_ss=smooth.spline(x=this_x,y=this_y,df=10)
plot(this_x,this_y,col=this_col,cex=0.3,pch=16)
#abline(h=mean(this_y),lwd=1,lty=1,col='black')
points(this_ss$x,this_ss$y,pch=16,cex=1,col='black',type='l',lty=2,lwd=3)
###########
this_x=X[LIN4]
this_y=Y[LIN4]
this_col=COL[LIN4]
this_ss=smooth.spline(x=this_x,y=this_y,df=10)
plot(this_x,this_y,col=this_col,cex=0.3,pch=16)
#abline(h=mean(this_y),lwd=1,lty=1,col='black')
points(this_ss$x,this_ss$y,pch=16,cex=1,col='black',type='l',lty=2,lwd=3)

dev.off()




SSS=SS1
c1=cor(SSS, AGE, method='spearman')
c2=cor(SSS[LIN1], AGE[LIN1], method='spearman')
c3=cor(SSS[LIN2], AGE[LIN2], method='spearman')
c4=cor(SSS[LIN3], AGE[LIN3], method='spearman')
c5=cor(SSS[LIN4], AGE[LIN4], method='spearman')
C1=c(c1,c2,c3,c4,c5)


SSS=DP
c1=cor(SSS, AGE, method='spearman')
c2=cor(SSS[LIN1], AGE[LIN1], method='spearman')
c3=cor(SSS[LIN2], AGE[LIN2], method='spearman')
c4=cor(SSS[LIN3], AGE[LIN3], method='spearman')
c5=cor(SSS[LIN4], AGE[LIN4], method='spearman')
C2=c(c1,c2,c3,c4,c5)


CCI=read.table('./rds/sctc_all.cci',sep='\t',header=F,row.names=NULL)[,1]

SSS=CCI
c1=cor(SSS, AGE, method='spearman')
c2=cor(SSS[LIN1], AGE[LIN1], method='spearman')
c3=cor(SSS[LIN2], AGE[LIN2], method='spearman')
c4=cor(SSS[LIN3], AGE[LIN3], method='spearman')
c5=cor(SSS[LIN4], AGE[LIN4], method='spearman')
C3=c(c1,c2,c3,c4,c5)


SSS=META$nCount_RNA
c1=cor(SSS, AGE, method='spearman')
c2=cor(SSS[LIN1], AGE[LIN1], method='spearman')
c3=cor(SSS[LIN2], AGE[LIN2], method='spearman')
c4=cor(SSS[LIN3], AGE[LIN3], method='spearman')
c5=cor(SSS[LIN4], AGE[LIN4], method='spearman')
C4=c(c1,c2,c3,c4,c5)


SSS=META$nFeature_RNA
c1=cor(SSS, AGE, method='spearman')
c2=cor(SSS[LIN1], AGE[LIN1], method='spearman')
c3=cor(SSS[LIN2], AGE[LIN2], method='spearman')
c4=cor(SSS[LIN3], AGE[LIN3], method='spearman')
c5=cor(SSS[LIN4], AGE[LIN4], method='spearman')
C5=c(c1,c2,c3,c4,c5)


SSS=META$nCount_ATAC
c1=cor(SSS, AGE, method='spearman')
c2=cor(SSS[LIN1], AGE[LIN1], method='spearman')
c3=cor(SSS[LIN2], AGE[LIN2], method='spearman')
c4=cor(SSS[LIN3], AGE[LIN3], method='spearman')
c5=cor(SSS[LIN4], AGE[LIN4], method='spearman')
C6=c(c1,c2,c3,c4,c5)

SSS=META$nFeature_ATAC
c1=cor(SSS, AGE, method='spearman')
c2=cor(SSS[LIN1], AGE[LIN1], method='spearman')
c3=cor(SSS[LIN2], AGE[LIN2], method='spearman')
c4=cor(SSS[LIN3], AGE[LIN3], method='spearman')
c5=cor(SSS[LIN4], AGE[LIN4], method='spearman')
C7=c(c1,c2,c3,c4,c5)




ES=readRDS(file='./rds/geneAct_ES_score.rds')
SSS=ES
c1=cor(SSS, AGE, method='spearman')
c2=cor(SSS[LIN1], AGE[LIN1], method='spearman')
c3=cor(SSS[LIN2], AGE[LIN2], method='spearman')
c4=cor(SSS[LIN3], AGE[LIN3], method='spearman')
c5=cor(SSS[LIN4], AGE[LIN4], method='spearman')
C8=c(c1,c2,c3,c4,c5)


MAT=cbind(C1,C2,C8,C3,C4,C5,C6,C7)
MAT=-MAT

pdf('./plot/s2p04_age_thisStudy_lineages_boxplot.pdf',width=5,height=3.5)
boxplot(MAT, pch='+',col=c('indianred1',rep('grey60',7)),lwd=1)
abline(h=0,lty=2,lwd=1)
dev.off()


t.test(MAT[,1],MAT[,2],paired=T)
#0.01487
t.test(MAT[,1],MAT[,3],paired=T)
#0.02296

apply(MAT,2,mean)

#       C1           C2           C8           C3           C4           C5
# 0.730277537  0.527093020  0.402411415  0.031507863  0.058137013 -0.004104048
#          C6           C7
#-0.195835538 -0.217098871




print(MAT)
#           C1        C2        C8         C3          C4          C5
#[1,] 0.7910820 0.4158901 0.4134077  0.3459126  0.09574268  0.01353653
#[2,] 0.5504066 0.4204766 0.3338671  0.5221800  0.11337556  0.11458036
#[3,] 0.5919859 0.3604371 0.5537996  0.3219227  0.44632593  0.39972971
#[4,] 0.9005172 0.7092931 0.3406960 -0.4937173 -0.12859045 -0.26212821
#[5,] 0.8173960 0.7293683 0.3702866 -0.5387587 -0.23616866 -0.28623863
#             C6         C7
#[1,] -0.1990504 -0.2159282
#[2,] -0.1652922 -0.1842995
#[3,]  0.1167110  0.1036840
#[4,] -0.4163768 -0.4428904
#[5,] -0.3151693 -0.3460602






