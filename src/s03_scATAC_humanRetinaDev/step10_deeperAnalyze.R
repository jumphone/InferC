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


SS1=readRDS('../rds/SS1.rds')
SS5=readRDS('../rds/SS5.rds')


COL=rep('NA',nrow(UMAP))
COL[which(BATCH=='d53')]='red3'
COL[which(BATCH=='d59')]='pink'
COL[which(BATCH=='d78')]='lightgoldenrod1'
COL[which(BATCH=='d78.2')]='lightgoldenrod1'
COL[which(BATCH=='d89C')]='springgreen3'
COL[which(BATCH=='d89P')]='skyblue1'


GA=readRDS('../rds/GeneActivity.rds')
ES=read.table('./GS/BHATTACHARYA_EMBRYONIC_STEM_CELL.gmt',sep='\t')[,1]
GA_ES=GA[which(rownames(GA) %in% ES),]

####################
GA_MAT=as.matrix(GA)
COR=cor(t(GA_MAT),TIME,method='spearman')
saveRDS(COR,'../rds/GA_MAT_corTime.rds')
###################
COR=readRDS('../rds/GA_MAT_corTime.rds')
UROW=which(!is.na(COR[,1]))
UCOR=COR[UROW,1]
names(UCOR)=rownames(COR)[UROW]
UCOR=-UCOR
#####################

length(UCOR)
#19477

this_cor=cor(SS1,-TIME,method='spearman')

length(which(UCOR>this_cor))
#[1] 20

print(names(which(UCOR>this_cor)))
# [1] "AGMO"          "COL14A1"       "COL24A1"       "EMCN"
# [5] "ESR1"          "GABRB1"        "GRM1"          "HBG2"
# [9] "KCND2"         "KCNS3"         "KLHL1"         "MSR1"
#[13] "PIK3C2G"       "RP11-385J1.3"  "RP11-545J16.1" "RXFP1"
#[17] "SLC25A21"      "TRPC6"         "TTC6"          "ZPLD1"



plot(density(UCOR))
DEN=density(UCOR)

#########################

pdf('../plot/s10p01_dentisy.pdf',width=8,height=3)
###############
plot(DEN)
points(DEN$x,DEN$y,type='h',col='grey80',lwd=2)
points(DEN,type='l',col='grey60',lwd=4)
######
this_cor=cor(SS1,-TIME,method='spearman')
abline(v=this_cor,lwd=4,col='indianred1')
ecdf(UCOR)(this_cor)
#0.9989731
#########
this_cor=cor(SS5,-TIME,method='spearman')
abline(v=this_cor,lwd=4,col='royalblue1')
ecdf(UCOR)(this_cor)
#0.9976896
################
used_index=which(names(UCOR) %in% ES)
points(x=UCOR[used_index],y=rep(3,length(used_index)),lwd=1,col='royalblue1',type='h')
summary(UCOR[used_index])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.02201 0.06075 0.07904 0.08701 0.10845 0.23515

abline(v=0,lty=2,lwd=2)
####################
dev.off()

###################################################

GA_MAT=as.matrix(GA)

COR=readRDS('../rds/GA_MAT_corTime.rds')
UROW=which(!is.na(COR[,1]))


CMAT=c()
this_used=which(TIME%in%c(1,2))
this_cor=cor(t(GA_MAT[,this_used]),TIME[this_used],method='spearman')
CMAT=cbind(CMAT, this_cor)

this_used=which(TIME%in%c(2,3))
this_cor=cor(t(GA_MAT[,this_used]),TIME[this_used],method='spearman')
CMAT=cbind(CMAT, this_cor)

this_used=which(TIME%in%c(3,4))
this_cor=cor(t(GA_MAT[,this_used]),TIME[this_used],method='spearman')
CMAT=cbind(CMAT, this_cor)

saveRDS(CMAT, file='../rds/GA_MAT_CMAT.rds')


UCMAT=-CMAT
UCMAT=UCMAT[UROW,]
UCMAT[which(is.na(UCMAT))]=0


C1=c()
C5=c()
this_used=which(TIME%in%c(1,2))
this_cor1=cor(SS1[this_used],TIME[this_used],method='spearman')
this_cor5=cor(SS5[this_used],TIME[this_used],method='spearman')
C1=c(C1,this_cor1)
C5=c(C5,this_cor5)

this_used=which(TIME%in%c(2,3))
this_cor1=cor(SS1[this_used],TIME[this_used],method='spearman')
this_cor5=cor(SS5[this_used],TIME[this_used],method='spearman')
C1=c(C1,this_cor1)
C5=c(C5,this_cor5)

this_used=which(TIME%in%c(3,4))
this_cor1=cor(SS1[this_used],TIME[this_used],method='spearman')
this_cor5=cor(SS5[this_used],TIME[this_used],method='spearman')
C1=c(C1,this_cor1)
C5=c(C5,this_cor5)

C1=-C1
C5=-C5



this_index=which(UCOR>cor(SS1,-TIME,method='spearman'))
this_mat=as.matrix(GA[this_index,])

t(apply(this_mat,1,tapply,TIME,mean))
























































