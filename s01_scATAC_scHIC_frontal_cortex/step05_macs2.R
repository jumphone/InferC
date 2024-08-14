

library(Seurat)
library(Signac)

pbmc=readRDS('scATAC_pbmc.rds')

head(pbmc@meta.data)

DimPlot(pbmc,group.by='type')

MAT=as.matrix(pbmc[['macs2']]@data)
MAT[1:5,1:5]
VEC=pbmc@reductions$umap@cell.embeddings

#############################
set.seed(123)
N=100
KM=kmeans(VEC, centers=N)
CLST=KM$cluster
saveRDS(KM, file='./rds/KM.rds')
############################

source('/home/database/data/InferX/data/frontal_cortex/scATAC/source.R')
KM=readRDS('./rds/KM.rds')
CLST=KM$cluster
pbmc$clst=CLST
DimPlot(pbmc,group.by='clst',label=TRUE)+NoLegend()

source('source.R')
AGG=.generate_mean(MAT, CLST)
AGG=AGG[,order(as.numeric(colnames(AGG)))]
saveRDS(AGG, file='./rds/macs2_AGG.rds')

############################
AGG=readRDS('AGG.rds')
############################

source('source.R')
OUT=.getChrLoc(rownames(AGG),SPLIT1='-',SPLIT2='-')
CHR=OUT$chr
START=OUT$start
END=OUT$end
LOC=(START+END)/2


# GJA1 CHR6 120M-124M
# PDGFRA CHR4 55100000   53100000 57100000
#USED_INDEX1=which(CHR=='chr4' & LOC > 50100000 & LOC < 60100000)
#USED_INDEX2=which(CHR=='chr4' & LOC > 50100000 & LOC < 60100000)

USED_INDEX1=which(CHR=='chr6' & LOC > 120000000 & LOC < 124000000)
USED_INDEX2=which(CHR=='chr6' & LOC > 120000000 & LOC < 124000000)

length(USED_INDEX1)
length(USED_INDEX2)

UMAT1=AGG[USED_INDEX1,]
UMAT2=AGG[USED_INDEX2,]
############################
COR=cor(t(UMAT1),t(UMAT2))
################################
COR[which(is.na(COR))]=0
COR_POS=COR
COR_POS[which(COR<0)]=0

####################################
UCHR1=CHR[USED_INDEX1]
ULOC1=LOC[USED_INDEX1]
UNAME1=rownames(UMAT1)
UCHR2=CHR[USED_INDEX2]
ULOC2=LOC[USED_INDEX2]
UNAME2=rownames(UMAT2)

source('/home/toolkit/src/InferLoop.R')

CUT=0.1
mat=AGG[unique(c(USED_INDEX1,USED_INDEX2)),]
vvv=as.numeric(COR)
p1=rep(rownames(COR),times=ncol(COR),each=1)
p2=rep(colnames(COR),times=1,each=nrow(COR))
over_cut=which(vvv>CUT)

################################################
net=cbind(p1[over_cut],p2[over_cut])
net_uniq=inferloop.getUniqLoop(net)
ILS=inferloop.inferLoopSignal(mat, net_uniq, r=0)

TAB=table(CLST,pbmc$type)

TAB_MAT=matrix(TAB, nrow=nrow(TAB),ncol=ncol(TAB))
rownames(TAB_MAT)=rownames(TAB)
colnames(TAB_MAT)=colnames(TAB)
TAB_MAT=as.data.frame(TAB_MAT)

USED_BIN=which(TAB_MAT$Astro / rowSums(TAB_MAT)> 0.8)
#USED_BIN=which(TAB_MAT$Microglia / rowSums(TAB_MAT)> 0.8)
#USED_BIN=which(TAB_MAT$Neuron.Ex.Glu / rowSums(TAB_MAT)> 0.8)
#USED_BIN=which(TAB_MAT$Neuron.In.GABA / rowSums(TAB_MAT)> 0.8)
#USED_BIN=which(TAB_MAT$OPC / rowSums(TAB_MAT)> 0.8)
#USED_BIN=which(TAB_MAT$ODC / rowSums(TAB_MAT)> 0.8)

ASTRO_ILS=rowMeans(ILS[,USED_BIN])
LOOP=inferloop.splitLoop(names(ASTRO_ILS))
OUT1=.getChrLoc(LOOP[,1])
OUT2=.getChrLoc(LOOP[,2])
CHR1=OUT1$chr
CHR2=OUT2$chr
LOC1=(OUT1$start+OUT1$end )/2
LOC2=(OUT2$start+OUT2$end )/2

Z=ASTRO_ILS
X=(LOC1+LOC2)/2
Y=abs(LOC2-LOC1)/2
Xall=c(LOC1,LOC2)

source('source.R')
#############
lineCol='grey50'
baseCol='grey90'
maxCol='red1'
maxCol.Z=2
#############
Zmax=2
Hh=1
Zh=1
Ch=1
#############
Hbase=Ch+Zh
Y=Y/max(Y)*Hh
#############
plot(c(min(Xall),max(Xall)),c(0,0),col='white',
     ylim=c(0,Hh+Zh+Ch ),xlim=c(min(Xall),max(Xall)),
     yaxt='n',xlab=CHR1[1],ylab='')
#############
COL=.colMap(Z,c(0,maxCol.Z),c(baseCol,maxCol))
points(X[order(Z)],Y[order(Z)]+Hbase,col=COL[order(Z)],pch=18)
#############
Zpos=Z-0
Zpos[which(Zpos<0)]=0
Zpos[which(Zpos>Zmax)]=Zmax
ZZZ= Hbase - Zpos/Zmax * Zh
segments(x0=LOC1[order(Z)],y0=ZZZ[order(Z)],x1=LOC2[order(Z)],y1=ZZZ[order(Z)],col=COL[order(Z)])
##############
segments(x0=min(Xall),y0=Hbase, x1=max(Xall), y1=Hbase, col=lineCol)
segments(x0=min(Xall),y0=Hbase, x1=(max(Xall)+min(Xall))/2, y1= Hbase+Hh, col=lineCol)
segments(x0=max(Xall),y0=Hbase, x1=(max(Xall)+min(Xall))/2, y1= Hbase+Hh, col=lineCol)
abline(h=Hbase-Zh,lty=2,col=lineCol)
##############

POS=121756745

D1=abs(LOC1-POS)
D2=abs(LOC2-POS)
MIN=min(c(D1,D2))
MIN_INDEX=which(c(D1,D2)==MIN)
MIN_LOC=c(LOC1,LOC2)[MIN_INDEX][1]

XL=c(min(c(LOC1,LOC2)),max(c(LOC1,LOC2)))
ZL=c(0,0)
Xp=c()
Zp=c()
i=1
while(i<=length(LOC1)){
    this_loc1=LOC1[i]
    this_loc2=LOC2[i]
    this_z=Zpos[i]
    if(this_loc1==MIN_LOC){
        this_x=this_loc2
        XL=c(XL, this_x)
        ZL=c(ZL, this_z)
        #points(this_loc1,Hbase-this_z/Zmax * Zh,pch=16,col='black')
        points(this_loc2,Hbase-this_z/Zmax * Zh,pch=1,col='black')
        #points(this_loc2, this_z/Zmax * Ch,pch=1,col='black')
        Xp=c(Xp,this_loc2)
        Zp=c(Zp,this_z/Zmax * Ch)
        }
    if(this_loc2==MIN_LOC){
        this_x=this_loc1
        XL=c(XL, this_x)
        ZL=c(ZL, this_z)
        #points(this_loc2,Hbase-this_z/Zmax * Zh,pch=16,col='black')
        points(this_loc1,Hbase-this_z/Zmax * Zh,pch=1,col='black')
        #points(this_loc1, this_z/Zmax * Ch,pch=1,col='black')
        Xp=c(Xp,this_loc1)
        Zp=c(Zp,this_z/Zmax * Ch)
        }
    i=i+1}


XL=c(XL,POS)
ZL=c(ZL,Zmax)

ZL=ZL[order(XL)]
XL=XL[order(XL)]

ZL[which(ZL<0)]=0
ZL=ZL/Zmax*Ch

NNN=1000
Lapp=c(1:NNN)/NNN * (max(XL)-min(XL))+min(XL)  #c(min(XL):max(XL))
Lapp=c(Lapp,XL)
Lapp=Lapp[order(Lapp)]
Zapp=approxfun(x=XL,y=ZL,ties='mean')(Lapp)
points(Lapp,Zapp,type='h',lwd=2,col='grey80')#col=.colMap(Zapp,c(0,maxCol.Z ),c(baseCol,maxCol)))
points(XL,ZL,type='b',col='black',cex=0)
#points(XL,ZL,type='p',pch=1,col='black')
points(Xp, Zp,pch=1,col='black')

segments(x0=POS,x1=POS,y0=0,y1=Ch+Zh,lty=1,col='black')
#points(x=POS, y=Ch+Zh, pch=16, col='black')
#points(x=POS, y=Ch, pch=16, col='black')












