
library(Seurat)
library(Signac)

pbmc=readRDS('scATAC_pbmc.rds')

head(pbmc@meta.data)

DimPlot(pbmc,group.by='type')

MAT=as.matrix(pbmc[['peaks']]@data)
MAT[1:5,1:5]
VEC=pbmc@reductions$umap@cell.embeddings

#############################
set.seed(123)
N=100
KM=kmeans(VEC, centers=N)
CLST=KM$cluster
saveRDS(KM, file='KM.rds')
############################

source('/home/database/data/InferX/data/frontal_cortex/scATAC/source.R')
KM=readRDS('KM.rds')
CLST=KM$cluster
pbmc$clst=CLST
DimPlot(pbmc,group.by='clst',label=TRUE)+NoLegend()

AGG=.generate_mean(MAT, CLST)
AGG=AGG[,order(as.numeric(colnames(AGG)))]
saveRDS(AGG, file='AGG.rds')

#####################################
AGG=readRDS('AGG.rds')

TAGG=t(AGG)

TMP1=rownames(AGG)
TMP2=t(matrix(unlist(stringr::str_split(TMP1,'-')),nrow=3))
CHR=TMP2[,1]
LOC=(as.numeric(TMP2[,2])+as.numeric(TMP2[,3]))/2




# GJA1 CHR6 120M-124M
USED_INDEX=which(CHR=='chr6' & LOC > 120000000 & LOC < 124000000)
length(USED_INDEX)


UMAT=AGG[USED_INDEX,]
COR=cor(t(UMAT))
COR[which(is.na(COR))]=0
COR_POS=COR
COR_POS[which(COR<0)]=0

UCHR=CHR[USED_INDEX]
ULOC=LOC[USED_INDEX]
UNAME=rownames(UMAT)

source('/home/toolkit/src/InferLoop.R')

CUT=0.1
mat=UMAT
vvv=COR_POS[upper.tri(COR_POS,diag=F)]
net=matrix(0,nrow=length(which(vvv>CUT)),ncol=2)
net_index=rep(0,length(which(vvv>CUT)))
i=1
t=1
while(i<nrow(COR_POS)){
    j=i+1
    while(j<=ncol(COR_POS)){
        if(COR_POS[i,j]>CUT){
            net[t,]=c(UNAME[i],UNAME[j])#paste0(UNAME[i],'.And.',UNAME[j])
            net_index[t]=(j-1)*nrow(COR_POS)+i
            t=t+1
            }
        j=j+1}
    print(i)
    i=i+1}

net_uniq=net
ILS=inferloop.inferLoopSignal(mat, net_uniq, r=0)

TAB=table(CLST,pbmc$type)

TAB_MAT=matrix(TAB, nrow=nrow(TAB),ncol=ncol(TAB))
rownames(TAB_MAT)=rownames(TAB)
colnames(TAB_MAT)=colnames(TAB)
TAB_MAT=as.data.frame(TAB_MAT)

USED_BIN=which(TAB_MAT$Astro / rowSums(TAB_MAT)> 0.8)
#USED_BIN=which(TAB_MAT$Microglia / rowSums(TAB_MAT)> 0.8)
#USED_BIN=which(TAB_MAT$Neuron.Ex.Glu / rowSums(TAB_MAT)> 0.8)

ASTRO_ILS=rowMeans(ILS[,USED_BIN])
ASTRO_ILS_POS=ASTRO_ILS
ASTRO_ILS_POS[which(ASTRO_ILS<0)]=0
ASTRO_ILS_POS_MAT=matrix(0,nrow=nrow(COR_POS),ncol=ncol(COR_POS))
ASTRO_ILS_POS_MAT[net_index]=ASTRO_ILS_POS


#######################
FINAL=ASTRO_ILS_POS_MAT * COR_POS
diag(FINAL)=1
FINAL[lower.tri(FINAL)]=t(FINAL)[lower.tri(t(FINAL))]
#FINAL=COR_POS
#FINAL=ASTRO_ILS_POS_MAT

.drawHeatmap<-function(o.mat, ...){
    library('ComplexHeatmap')
    library('circlize')
    o.mat=as.matrix(o.mat)
    ###################
    vvv=as.numeric(o.mat)
    
    col.fun=col_fun =colorRamp2(c(0,0.1),
              c('grey90','red3'))
    ################
    ht=Heatmap(o.mat,row_title='',name="v",
        cluster_rows=F, cluster_columns=F,
        show_column_dend = F, show_row_dend = F,
        show_column_names= F, show_row_names= F,
        border = TRUE,
        col=col.fun,
        row_names_side='right',
        row_names_gp = gpar(fontsize = 10, lineheight=NULL),
        column_names_gp = gpar(fontsize = 10, lineheight=NULL),
        row_names_max_width = max_text_width(rownames(o.mat),gp = gpar(fontsize = 15)),
        column_names_max_height = max_text_width(colnames(o.mat),gp = gpar(fontsize = 15)),
        ...)
    pdf('tmp_initialize.pdf',width=7,height=7)
    ht=draw(ht)
    dev.off()
    return(ht)
    }

o.mat=FINAL
.drawHeatmap(o.mat)



















UMAT=AGG
TMP1=rownames(UMAT)
TMP2=t(matrix(unlist(stringr::str_split(TMP1,'-')),nrow=3))

CHR=TMP2[,1]
LOC=(as.numeric(TMP2[,2])+as.numeric(TMP2[,3]))/2


REGION


WIN=1000
STEP=100







library(reticulate)
use_python("/home/toolkit/local/bin/python3",required=TRUE)
np = import("numpy",convert=FALSE)








TTT=AGG[1:10,] %*% t(AGG[1:10,])






COV = as.matrix(np$cov(AGG[1:10000,]))

TMP=cov(t(AGG[1:1000,]))


UMAT=AGG
REGION=rownames(UMAT)
TMP=t(matrix(unlist(stringr::str_split(REGION,'-')),nrow=3))







