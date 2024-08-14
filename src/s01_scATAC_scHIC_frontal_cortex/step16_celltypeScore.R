
source('InferC.R')


COUNT=readRDS('./rds/COUNT_ALL.rds')
CONNS=readRDS('./rds/CONNS_ALL.rds')


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


LOOP=rownames(CONNS)

NET=inferloop.splitLoop(LOOP)
CELLBIN=inferloop.generateBin(MAT,VEC, n=100)
saveRDS(CELLBIN, file='./rds/inferloop_CELLBIN.rds')
saveRDS(NET, file='./rds/inferloop_NET.rds')



ILS_MAT=inferloop.inferLoopSignal(CELLBIN$mat, NET, r=0)
ILS_MAT=ILS_MAT[,order(as.numeric(colnames(ILS_MAT)))]
saveRDS(ILS_MAT, file='./rds/inferloop_ILS_MAT.rds')
#####################################






source('InferC.R')
library(Seurat)
library(Signac)

COUNT=readRDS('./rds/COUNT_ALL.rds')
CONNS=readRDS('./rds/CONNS_ALL.rds')
META=readRDS('scATAC_pbmc_meta.rds')
CELLBIN=readRDS('./rds/inferloop_CELLBIN.rds')
ILS_MAT=readRDS('./rds/inferloop_ILS_MAT.rds')
PV_MAT=pnorm(ILS_MAT)

TAB=table(META$type, CELLBIN$clst)
PMAT=t(t(TAB)/colSums(TAB))
MAX_ROW_INDEX=apply(-PMAT,2,order)[1,]
TAG=rownames(PMAT)[MAX_ROW_INDEX]



TYPE_ILS_MAT=inferc.generate_mean(ILS_MAT, TAG)
TYPE_ILS_MAT=TYPE_ILS_MAT[,order(colnames(TYPE_ILS_MAT))]
saveRDS(TYPE_ILS_MAT, file='./rds/inferloop_TYPE_ILS_MAT.rds')





######################
# Reanalysis
######################


source('InferC.R')
library(Seurat)
library(Signac)

COUNT=readRDS('./rds/COUNT_ALL.rds')
CONNS=readRDS('./rds/CONNS_ALL.rds')
META=readRDS('scATAC_pbmc_meta.rds')
CELLBIN=readRDS('./rds/inferloop_CELLBIN.rds')
ILS_MAT=readRDS('./rds/inferloop_ILS_MAT.rds')
TYPE_ILS_MAT=readRDS('./rds/inferloop_TYPE_ILS_MAT.rds')



.normRPM <- function(x){
    x.sum=sum(x)
    if(x.sum==0){
        y=x
    }else{
        y=x / x.sum *1000000
        }
    return(y)
    }


COUNT_nonAll=COUNT[,c(1:6)]

COUNT_RPM=apply(COUNT_nonAll,2,.normRPM)
COUNT_RPM_S=Seurat::FastRowScale(COUNT_RPM,scale_max = 100)
COUNT_RPM_S[which(is.na(COUNT_RPM_S))]=0



USED_INDEX=which(rowSums(COUNT_nonAll)>0)

inferloop_mat=TYPE_ILS_MAT[USED_INDEX,]
COR1=cor(inferloop_mat,COUNT_RPM_S[USED_INDEX,],method='spearman')

#thisstudy_mat=TYPE_ILS_MAT[USED_INDEX,] * CONNS[USED_INDEX,1]


.ssfunc<-function(x,a=0.5,b=10){return(1/(1+exp(-(x-a)*b)))}

x=c(1:1000)/1000
y=.ssfunc(x,a=0.5,b=10)
plot(x,y)



P=CONNS[USED_INDEX,1]
#W=pnorm(logit(P))
W=.ssfunc(P,a=0.5,b=10)

thisstudy_mat = TYPE_ILS_MAT[USED_INDEX,] * W

COR2=cor(thisstudy_mat,COUNT_RPM_S[USED_INDEX,],method='spearman')













t.test(diag(COR1),c(COR1[upper.tri(COR1)],COR1[lower.tri(COR1)]))
t.test(diag(COR2),c(COR2[upper.tri(COR2)],COR2[lower.tri(COR2)]))

min(diag(COR1))-max(c(COR1[upper.tri(COR1)],COR1[lower.tri(COR1)]))
min(diag(COR2))-max(c(COR2[upper.tri(COR2)],COR2[lower.tri(COR2)]))

t.test( diag(COR1) , diag(COR2) ,paired=T)
#0.001447
boxplot(diag(COR1) , diag(COR2))


boxplot(diag(COR1) , diag(COR2), diag(COR3))


mean(diag(COR1))
#0.06621812
mean(diag(COR2))
#0.08056688

saveRDS(COR1,'./rds/schic_zvalue_vs_lls_inferloop_cor.rds')
saveRDS(COR2,'./rds/schic_zvalue_vs_lls_thisstudy_cor.rds')


COR1=readRDS('./rds/schic_zvalue_vs_lls_inferloop_cor.rds')
COR2=readRDS('./rds/schic_zvalue_vs_lls_thisstudy_cor.rds')


pdf('./plot/p06_lls_vs_zvalue_boxplot.pdf',width=3,height=3.5)
boxplot(diag(COR2),diag(COR1),col=c('indianred1','grey60'),ylim=c(0.02,0.22))
boxplot(diag(COR1),diag(COR2),col=c('grey60','indianred1'),ylim=c(0.02,0.22))
dev.off()



library('ComplexHeatmap')
library('circlize')

o.mat=COR1
###############################
vvv=as.numeric(o.mat)
col_fun =colorRamp2(
             c(0, 0.04, 0.05, 0.10),
             c('royalblue3','skyblue1','grey95','indianred1')
             )

COL_LIST=list(CType=c('Astro'='red','MG'='gold1','Neuron.Ex'='darkgreen','Neuron.In'='darkolivegreen1','ODC'='blue3','OPC'='skyblue'))
CType=c('Astro','MG','Neuron.Ex','Neuron.In','ODC','OPC')

ha_top = HeatmapAnnotation( CType=CType,
                            col=COL_LIST
                          )
ha_left = rowAnnotation( CType=CType,
                         col=COL_LIST
                          )



ht1=Heatmap(o.mat,row_title='',name="v",
        cluster_rows=F, cluster_columns=F,
        show_column_dend = F, show_row_dend = F,
        show_column_names= T, show_row_names= T,
        col=col_fun, border = TRUE,
        row_names_side='right',
        row_names_gp = gpar(fontsize = 10, lineheight=NULL),
        column_names_gp = gpar(fontsize = 10, lineheight=NULL),
        row_names_max_width = max_text_width(rownames(o.mat),gp = gpar(fontsize = 15)),
        column_names_max_height = max_text_width(colnames(o.mat),gp = gpar(fontsize = 15)),
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.2f", o.mat[i, j]), x, y, gp = gpar(fontsize = 10))
            }
          )



ht2=Heatmap(o.mat,row_title='',name="v",
        cluster_rows=F, cluster_columns=F,
        show_column_dend = F, show_row_dend = F,
        show_column_names= T, show_row_names= T,
        col=col_fun, border = TRUE,
        row_names_side='right',
        row_names_gp = gpar(fontsize = 10, lineheight=NULL),
        column_names_gp = gpar(fontsize = 10, lineheight=NULL),
        row_names_max_width = max_text_width(rownames(o.mat),gp = gpar(fontsize = 15)),
        column_names_max_height = max_text_width(colnames(o.mat),gp = gpar(fontsize = 15)),
        top_annotation = ha_top,
        left_annotation = ha_left,
        )


pdf('./plot/p06_lls_vs_zvalue_heatmap_inferloop.pdf',width=7,height=5)
print(ht1)
print(ht2)
dev.off()




library('ComplexHeatmap')
library('circlize')

o.mat=COR2
###############################
vvv=as.numeric(o.mat)
col_fun =colorRamp2(
             c(0, 0.04, 0.05, 0.10),
             c('royalblue3','skyblue1','grey95','indianred1')
             )

ht1=Heatmap(o.mat,row_title='',name="v",
        cluster_rows=F, cluster_columns=F,
        show_column_dend = F, show_row_dend = F,
        show_column_names= T, show_row_names= T,
        col=col_fun, border = TRUE,
        row_names_side='right',
        row_names_gp = gpar(fontsize = 10, lineheight=NULL),
        column_names_gp = gpar(fontsize = 10, lineheight=NULL),
        row_names_max_width = max_text_width(rownames(o.mat),gp = gpar(fontsize = 15)),
        column_names_max_height = max_text_width(colnames(o.mat),gp = gpar(fontsize = 15)),
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.2f", o.mat[i, j]), x, y, gp = gpar(fontsize = 10))
            }
          )



ht2=Heatmap(o.mat,row_title='',name="v",
        cluster_rows=F, cluster_columns=F,
        show_column_dend = F, show_row_dend = F,
        show_column_names= T, show_row_names= T,
        col=col_fun, border = TRUE,
        row_names_side='right',
        row_names_gp = gpar(fontsize = 10, lineheight=NULL),
        column_names_gp = gpar(fontsize = 10, lineheight=NULL),
        row_names_max_width = max_text_width(rownames(o.mat),gp = gpar(fontsize = 15)),
        column_names_max_height = max_text_width(colnames(o.mat),gp = gpar(fontsize = 15)),
        top_annotation = ha_top,
        left_annotation = ha_left,
        )


pdf('./plot/p06_lls_vs_zvalue_heatmap_thisstudy.pdf',width=7,height=5)
print(ht1)
print(ht2)
dev.off()









zvalue_mat=COUNT_RPM_S[USED_INDEX,]
inferloop_mat=TYPE_ILS_MAT[USED_INDEX,]
gls_v=CONNS[USED_INDEX,1]



X=c()
Y=c()
i=0
step=0.1
while(i<=0.9){
   this_index=which(ecdf(gls_v)(gls_v)>=i & ecdf(gls_v)(gls_v)<i+step)
   this_cor=mean(diag(cor(zvalue_mat[this_index,], inferloop_mat[this_index,], method='spearman')))
   X=c(X,i)
   Y=c(Y,this_cor)
   print(i)
   i=i+step}




plot(quantile(W,X),Y)








