

library(Seurat)
library(Signac)

pbmc=readRDS('scATAC_pbmc.rds')
head(pbmc@meta.data)
DimPlot(pbmc,group.by='type')
MAT=as.matrix(pbmc[['macs2']]@data)
#MAT=as.matrix(pbmc[['peaks']]@data)
MAT[1:5,1:5]
VEC=pbmc@reductions$umap@cell.embeddings
TYPE=pbmc$type

source('old/InferC_20240125.R')
s1_agg = inferc.aggMat(MAT, VEC, clstNum=100, seed=123)


TSS=read.table('/home/database/annotation/hg19/Homo_sapiens.GRCh37.75.chr.gtf.tran.bed.TSS.sort.bed.geneUniq.txt',sep='\t',row.names=NULL,header=FALSE)
this_index=1


SLOP=2000000
#SLOP=500000

CHR='chr4'
VIEW_POINT=55095264
START=VIEW_POINT-SLOP
END=VIEW_POINT+SLOP



CHR='chr18'
VIEW_POINT=74729224
START=VIEW_POINT-SLOP
END=VIEW_POINT+SLOP



source('old/InferC_20240125.R')

s2_loop = inferc.calCorLoop(s1_agg$agg, CHR, START, END,
                            Dcut=0.5, SPLIT1='-', SPLIT2='-')


LOOP = inferc.agg2cell(AGG=s2_loop$loop, CLST=s1_agg$clst)
LOOP_TYPE = inferc.generate_mean(LOOP, TYPE)

colnames(LOOP_TYPE)

source('old/InferC_20240125.R')
i=5
LOOP_SIGNAL=LOOP_TYPE[,i]
s3_draw=inferc.preContact(LOOP_SIGNAL,exMin=0,MAIN=colnames(LOOP_TYPE)[i],Zmax=2.5)
s3_draw=inferc.preV4C(s3_draw, VIEW_POINT)

#inferc.drawContact(s3_draw)

inferc.drawContact(s3_draw)
inferc.addV4C(s3_draw)


inferc.drawV4C(s3_draw)


#############
pdf('tmp.pdf',width=18,height=5)
par(mfrow=c(1,ncol(LOOP_TYPE)))
par(mar = c(4, 1, 2, 1)) #par(mar = c(bottom, left, top, right))
i=1
while(i<=ncol(LOOP_TYPE)){
    LOOP_SIGNAL=LOOP_TYPE[,i]
    s3_draw=inferc.preContact(LOOP_SIGNAL,exMin=0,MAIN=colnames(LOOP_TYPE)[i],Zmax=2.5)
    s3_draw=inferc.preV4C(s3_draw,VIEW_POINT=VIEW_POINT)
    inferc.drawContact(s3_draw)
    inferc.addV4C(s3_draw)
    print(paste0(i,' / ',ncol(LOOP_TYPE)))
    print(colnames(LOOP_TYPE)[i])
    i=i+1
    }
dev.off()













VP_LOC=TSS[,3]
VP_CHR=TSS[,1]
VP_STRAND=TSS[,4]
VP_TYPE=TSS[,7]

MACS=inferc.getChrLoc(rownames(MAT))

TMP_LOC=0
VP_MIN=c()
VP_USE=c()
i=1
while(i<=nrow(TSS)){
    this_chr=VP_CHR[i]
    this_loc=VP_LOC[i]
    this_index=which(MACS$chr==this_chr)
    this_macs_loc=(MACS$end[this_index]+MACS$start[this_index])/2
    this_dist=abs(this_macs_loc-this_loc)
    this_min=min(this_dist)
    VP_MIN=c(VP_MIN,this_min)
    if(abs(this_loc-TMP_LOC)>SLOP & this_min<1000){
        TMP_LOC=this_loc;this_use=1}else{this_use=0}
    VP_USE=c(VP_USE,this_use)
    if(i%%1000==1){print(paste0(i,' / ',nrow(TSS)))}
    i=i+1}

#saveRDS(VP_MIN, file='./rds/VP_MIN.rds')
#saveRDS(VP_USE, file='./rds/VP_USE.rds')

VP_MIN=readRDS('./rds/VP_MIN.rds')
VP_USE=readRDS('./rds/VP_USE.rds')

VP_USED_INDEX=which(VP_USE==1)


USED_TSS=TSS[VP_USED_INDEX,]

OUT=list()
i=5
this_index=1
while(this_index<=nrow(USED_TSS)){    
    this_chr=USED_TSS[this_index,1]
    this_vp=USED_TSS[this_index,3]
    this_strand=USED_TSS[this_index,4]
    this_start=this_vp-SLOP
    this_end=this_vp+SLOP
    ##############################
    s2_loop = inferc.calCorLoop(s1_agg$agg, this_chr, this_start, this_end, corCut=0.1)
    LOOP = inferc.agg2cell(AGG=s2_loop$loop, CLST=s1_agg$clst)
    LOOP_TYPE = inferc.generate_mean(LOOP, TYPE)
    this_signal=LOOP_TYPE[,i]
    s3_draw=inferc.preContact(this_signal,exMin=0,MAIN=colnames(LOOP_TYPE)[i],Zmax=3)
    s4_4c=inferc.preV4C(s3_draw, VIEW_POINT=this_vp) 
    s4_4c$v4c_strand=this_strand
    OUT[[this_index]] <- s4_4c
    ###############################
    #if(this_index%%1000==1){print(paste0(this_index,' / ',nrow(USED_TSS)))}
    print(paste0(this_index,' / ',nrow(USED_TSS),' !!!'))
    this_index=this_index+1
    }


saveRDS(OUT, file='./rds/OUT5.rds')

inferc.drawV4C(OUT[[3]])

NNN=1000

MMM=c()
i=1
while(i<=300){
    this_out=OUT[[i]]
    this_x=this_out$v4c_x
    this_z=this_out$v4c_z
    this_strand=this_out$v4c_strand
    if(this_strand=='-'){max(this_x)-this_x}
    this_x=(this_x-min(this_x))/(max(this_x)-min(this_x))
    this_z_pred=approxfun(x=this_x, y=this_z,ties=mean)(c(1:NNN)/NNN)
    MMM=cbind(MMM, this_z_pred)
    i=i+1}


CCC=cor(MMM, c(1:NNN),method='pearson')
CCC[which(is.na(CCC))]=0


SMMM=apply(MMM,2,scale)
SMMM[which(is.na(SMMM))]=0
SMMM=SMMM
.smooth<-function(x,df){return(smooth.spline(x,df=df)$y)}
SMMM=apply(SMMM,2,.smooth,10)


source('/home/database/repli/stepClst.R')

CLST=as.integer(c(1:nrow(SMMM))/((nrow(SMMM)+1)/10))

SSS=inferc.generate_mean(t(SMMM),CLST)

SC=.stepClst(SSS)






library('ComplexHeatmap')
library('circlize')
o.mat=SMMM[,order(SC$split)]
###################
vvv=as.numeric(o.mat)
col_fun =colorRamp2( c(-1, 0, 1),
              c('royalblue3','white','red3'))
################
ht=Heatmap(o.mat,row_title='',name="v",
        cluster_rows=F, cluster_columns=F,
        show_column_dend = F, show_row_dend = F,
        show_column_names= F, show_row_names= F,
        col=col_fun, border = TRUE,
        row_names_side='right',
        row_names_gp = gpar(fontsize = 10, lineheight=NULL),
        column_names_gp = gpar(fontsize = 10, lineheight=NULL),
        row_names_max_width = max_text_width(rownames(o.mat),gp = gpar(fontsize = 15)),
        column_names_max_height = max_text_width(colnames(o.mat),gp = gpar(fontsize = 15)),
        )

print(ht)


inferc.drawV4C(OUT[[order(CCC)[10]]])



CHR='chr4'
VIEW_POINT=55095264
START=VIEW_POINT-2000000
END=VIEW_POINT+2000000


#############
source('InferC.R')

s2_loop = inferc.calCorLoop(s1_agg$agg, CHR, START, END,
                            corCut=0.1, SPLIT1='-', SPLIT2='-')

LOOP = inferc.agg2cell(AGG=s2_loop$loop, CLST=s1_agg$clst)
LOOP_TYPE = inferc.generate_mean(LOOP, TYPE)

colnames(LOOP_TYPE)

i=5
LOOP_SIGNAL=LOOP_TYPE[,i]
s3_draw=inferc.preContact(LOOP_SIGNAL,exMin=0,MAIN=colnames(LOOP_TYPE)[i],Zmax=3)
s3_draw=inferc.preV4C(s3_draw, VIEW_POINT)

#inferc.drawContact(s3_draw)

inferc.drawContact(s3_draw)
inferc.addV4C(s3_draw)


inferc.drawV4C(s3_draw)





LOOP=inferloop.splitLoop(names(LOOP_SIGNAL))
    


#############
pdf('tmp.pdf',width=18,height=5)
par(mfrow=c(1,ncol(LOOP_TYPE)))
par(mar = c(4, 1, 2, 1)) #par(mar = c(bottom, left, top, right))
i=1
while(i<=ncol(LOOP_TYPE)){
    LOOP_SIGNAL=LOOP_TYPE[,i]
    s3_draw=inferc.preContact(LOOP_SIGNAL,exMin=0,MAIN=colnames(LOOP_TYPE)[i],Zmax=3)
    s3_draw=inferc.preV4C(s3_draw,VIEW_POINT=VIEW_POINT)
    inferc.drawContact(s3_draw)
    inferc.addV4C(s3_draw)
    print(paste0(i,' / ',ncol(LOOP_TYPE)))
    print(colnames(LOOP_TYPE)[i])
    i=i+1
    }
dev.off()











LOOP_SIGNAL=LOOP_TYPE[,5]
s3_draw=inferc.preContact(LOOP_SIGNAL,exMin=0,MAIN=colnames(LOOP_TYPE)[i])
plot(s3_draw$X,s3_draw$Z)



points(s3_draw$X,predict(loess(s3_draw$Z~s3_draw$X,span=0.1)),col='red')












LOOP_SIGNAL=LOOP_TYPE[,'Astro']
source('InferC.R')
LOOP_SIGNAL=LOOP_TYPE[,'Astro']
s3_draw=inferc.preContact(LOOP_SIGNAL,exMin=0)
inferc.drawContact(s3_draw)
inferc.add4C(s3_draw, VIEW_POINT=VIEW_POINT)


saveRDS(s1_agg,  './rds/s1_agg.rds')
saveRDS(s2_loop, './rds/s2_loop.rds')
saveRDS(s3_draw, './rds/s3_draw.rds')











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












